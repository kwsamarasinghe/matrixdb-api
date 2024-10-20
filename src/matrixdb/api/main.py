import hashlib
import traceback

from flask import Flask, request, Response
from flask_caching import Cache
from flask_cors import CORS
import json
from itertools import groupby

from gevent.pywsgi import WSGIServer
import requests as requests

from src.matrixdb.services.biomolecules.protein_data_manager import ProteinDataManager
from src.matrixdb.services.interactome.interaction_manager import InteractionDataManager
from src.matrixdb.services.interactome.network_manager import NetworkManager
from src.matrixdb.utils.cache.cache_manager import CacheManager
from src.matrixdb.utils.database.database_manager import DatabaseManager
from src.matrixdb.utils.solr.solr_query_manager import SolrQueryManager

from dotenv import load_dotenv
import os

# Load environment variables from .env file
load_dotenv()

app_config = {
    "database_url": os.getenv('DATABASE_URL'),
    "primary_database_name": os.getenv("PRIMARY_DATABASE_NAME"),
    "secondary_database_name": os.getenv("SECONDARY_DATABASE_NAME"),
    "solr_url": os.getenv('SOLR_URL'),
    "pdb_api_url": os.getenv('PDB_API_URL'),
}

config = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "SimpleCache",  # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 300
}

app = Flask(__name__)
cors = CORS(app, resources={
    r"/api/*": {
        "origins": ["http://localhost:3000"],
        "methods": ["GET", "POST", "PUT", "DELETE"],
        "allow_headers": ["Content-Type", "Authorization"],
        "supports_credentials": True
    }
})
app.config.from_mapping(config)
cache = Cache(app)

database_manager = DatabaseManager(app_config)
cache_manager = CacheManager(app_config, database_manager)
protein_data_manager = ProteinDataManager(app_config, cache_manager, database_manager)
interaction_data_manager = InteractionDataManager(app_config, database_manager)
network_manager = NetworkManager(database_manager, cache_manager, protein_data_manager, interaction_data_manager)
solr_manager = SolrQueryManager(app_config)


@app.route('/api/biomolecules/<id>', methods=['GET'])
@cache.cached(timeout=50000, query_string=True)
def get_biomolecule_by_id(id):
    core_database_connection = database_manager.get_primary_connection()
    biomolecule = core_database_connection["biomolecules"].find_one(
        {
            "id": id
        },
        {
            "_id": False
        }
    )

    # Include GO, uniprot keyword and interpro definitions
    meta_data_cache = cache_manager.get_meta_data()
    if "annotations" in biomolecule and "go" in biomolecule["annotations"]:
        go_terms = list()
        for go in biomolecule["annotations"]["go"]:
            go = {
                "id": go,
                "term": meta_data_cache["go"][go]["term"],
                "definition": meta_data_cache["go"][go]["description"],
                "category": meta_data_cache["go"][go]["category"]
            }
            go_terms.append(go)
        biomolecule["annotations"]["go"] = go_terms

    if biomolecule["type"] == 'protein' or biomolecule["type"] == 'multimer':
        keywords = list()
        if "annotations" in biomolecule and "keywords" in biomolecule["annotations"]:
            for keyword in biomolecule["annotations"]["keywords"]:
                # Have to handle unipro/trembl separately for the moment
                if type(keyword) == dict:
                    keyword_id = keyword['id']
                    keyword = {
                        "id": keyword_id,
                        "definition": meta_data_cache["uniprotKeywords"][keyword_id]["definition"],
                        "term": meta_data_cache["uniprotKeywords"][keyword_id]["name"]
                    }
                else:
                    keyword = {
                        "id": keyword,
                        "definition": meta_data_cache["uniprotKeywords"][keyword]["definition"],
                        "term": meta_data_cache["uniprotKeywords"][keyword]["name"]
                    }
                    keywords.append(keyword)
            biomolecule["annotations"]["keywords"] = keywords

        if "xrefs" in biomolecule and "interpro" in biomolecule["xrefs"]:
            if biomolecule["xrefs"]["interpro"]:
                interpro_terms = list()
                for interpro in biomolecule["xrefs"]["interpro"]:
                    if 'IPR' in interpro:
                        if interpro in meta_data_cache["interpro"]:
                            interpro_terms.append({
                                'id': meta_data_cache["interpro"][interpro]["id"],
                                'value': meta_data_cache["interpro"][interpro]["name"]
                            })
                biomolecule["xrefs"]["interpro"] = interpro_terms

        if biomolecule["type"] == 'multimer' and "xrefs" in biomolecule and "reactome" in biomolecule["xrefs"]:
            if biomolecule["xrefs"]["reactome"]:
                reactome_terms = list()
                for reactome in biomolecule["xrefs"]["reactome"]:
                        if reactome["id"] in meta_data_cache["reactome"]:
                            reactome_data = meta_data_cache["reactome"][reactome["id"]]
                            reactome_terms.append({
                                'id': reactome_data["id"],
                                'value': reactome_data["term"].strip('"')
                            })
                biomolecule["xrefs"]["reactome"] = reactome_terms

    return Response(json.dumps(biomolecule), mimetype='application/json')


@app.route('/api/biomolecules/', methods=['POST'])
def get_biomolecules_by_id():
    core_database_connection = database_manager.get_primary_connection()
    biomoelcule_ids = json.loads(request.data)["ids"]
    biomolecules = list(core_database_connection["biomolecules"].find({
        "id": {
            "$in": biomoelcule_ids
        }
    },
    {
        "_id": False
    }))

    for biomolecule in biomolecules:
        if biomolecule["type"] == 'protein':

            expressions_by_protein = core_database_connection["geneExpression"].find_one({
                "uniprot": biomolecule["id"]
            })

            proteomics_expression_by_protein = core_database_connection["proteomicsExpressions"].find_one({
                "uniprot": biomolecule["id"]
            })

            gene_expressions = list()
            if expressions_by_protein is not None:

                expression_values = list()
                for expression in expressions_by_protein["expressions"]:
                    for e in list(expression.values()):
                        expression_values.append(e[0])

                grouped_expressions = groupby(sorted(list(expression_values), key=lambda x: x["uberonName"]), key=lambda x: x["uberonName"])

                for tissue_uberon_name, group in grouped_expressions:
                    max_group = max(group, key=lambda x: x["tpm"])
                    gene_expressions.append({
                        "tissueUberonName": tissue_uberon_name.strip('"'),
                        "tpm": max_group["tpm"]
                    })

            prot_expressions = list()
            if proteomics_expression_by_protein is not None:
                for e in proteomics_expression_by_protein["expressions"]:
                    prot_expressions.append(
                        {
                            "tissueId": e["tissueId"],
                            "sampleName": e["sampleName"] if "sampleName" in e else "",
                            "sample": e["sample"] if "sample" in e else "",
                            "score": e["confidenceScore"] if "confidenceScore" in e else 0,
                        })

            biomolecule["geneExpression"] = gene_expressions
            biomolecule["proteomicsExpression"] = prot_expressions

    return Response(json.dumps(biomolecules), mimetype='application/json')


@app.route('/api/biomolecules/<id>/interactors/', methods=['GET'])
def get_biomolecule_interactors_by_id(id):
    core_database_connection = database_manager.get_primary_connection()
    interactions = list(core_database_connection["interactions"].find(
        {
            "participants": id
        },
        {
            "_id": 0
        }
    ))

    neighborhood = {}
    direct = set()
    predictions = 0
    inferred = 0
    publications = set()
    if len(interactions) > 0:
        for interaction in interactions:
            participant_ids = interaction["participants"]
            participant_ids.remove(id)
            if len(participant_ids) == 0:
                continue
            partner = participant_ids[0]
            if partner not in neighborhood:
                neighborhood[partner] = {}
                neighborhood[partner]["association"] = interaction["id"]

            if "experiments" in interaction:
                if "spoke_expanded_from" in interaction["experiments"]["direct"]:
                    if "spoke_expanded_from" not in neighborhood[partner]:
                        neighborhood[partner]["spokeExpandedFrom"] = list()
                    neighborhood[partner]["spokeExpandedFrom"].extend(interaction["experiments"]["direct"]["spoke_expanded_from"])
                    for se in interaction["experiments"]["direct"]["spoke_expanded_from"]:
                        direct.add(se)

                if "binary" in interaction["experiments"]["direct"]:
                    if "directly_supported_by" not in neighborhood[partner]:
                        neighborhood[partner]["directlySupportedBy"] = list()
                    neighborhood[partner]["directlySupportedBy"].extend(interaction["experiments"]["direct"]["binary"])
                    for b in interaction["experiments"]["direct"]["binary"]:
                        direct.add(b)

            if "prediction" in interaction:
                if "predictions" not in neighborhood[partner]:
                    neighborhood[partner]["predictions"] = list()
                neighborhood[partner]["predictions"].append(partner)
                predictions += 1

    print()
    return {
        "count": len(neighborhood.keys()),
        "direct": len(direct),
        "predictions": predictions,
        "inferred": 0,
        "details": neighborhood,
        "publications": list(publications)
    }


@app.route('/api/biomolecules/proteins/expressions/', methods=['POST'])
def get_protein_expression():
    expression_data = dict()
    for protein in json.loads(request.data):
        gene = protein_data_manager.get_gene_name(protein)
        expression_data[protein] = {
            'gene': gene,
            'geneExpression': protein_data_manager.get_gene_expressions(protein),
            'proteomicsExpression': protein_data_manager.get_proteomics_expressions(protein)
        }
    return expression_data


@app.route('/api/biomolecules/<id>/binding-regions', methods=['GET'])
def get_binding_by_id(id):
    core_database_connection = database_manager.get_primary_connection()

    # Get pdb to uniprot mapping
    biomolecule = core_database_connection["biomolecules"].find_one({
        'id': id
    })
    pdb_entities = dict()
    for pdb in biomolecule['molecular_details']['pdb']:
        if "properties" in pdb:
            for prop in pdb["properties"]:
                if prop["type"] == "chains":
                    chain_names = prop["value"].split('=')[0]

                    if '/' in chain_names:
                        chain_names = chain_names.split('/')
                    else:
                        chain_names = [chain_names]

                    for chain_name in chain_names:
                        chain_start = prop["value"].split('=')[1].split('-')[0]
                        chain_end = prop["value"].split('=')[1].split('-')[1]
                        if pdb['id'] not in pdb_entities:
                            pdb_entities[pdb['id']] = list()

                        pdb_entities[pdb['id']].append({
                            'chain': chain_name,
                            'start': chain_start,
                            'end': chain_end
                        })

    # Get the pdb to uniprot mapping
    try:
        response = requests.get(f'{app_config["pdb_api_url"]}/{id}', timeout=20)
        pdb_response = response.json()
        if "mappings" in pdb_response[id]:
            pdb_mappings = pdb_response[id]["mappings"]
            for pdb_mapping in pdb_mappings:
                pdb_id = pdb_mapping['entry_id']
                for segment in pdb_mapping['segments']:
                    for chain in segment["chains"]:
                        if pdb_id.upper() in pdb_entities:
                            for pdb_chain in pdb_entities[pdb_id.upper()]:
                                if chain == pdb_chain['chain']:
                                    pdb_chain['pdb_start'] = segment['pdb_start']
                                    pdb_chain['pdb_end'] = segment['pdb_end']
                                    pdb_chain['unp_start'] = segment['unp_start']
                                    pdb_chain['unp_end'] = segment['unp_end']

    except requests.exceptions.Timeout:
        print("No pdb to uniprot mapping retrieved")

    interactions = core_database_connection["interactions"].find({"participants": id})

    # Experiment ids
    experiment_ids = set()
    for interaction in interactions:
        if "experiments" not in interaction:
            continue

        if len(interaction["experiments"]["direct"]["binary"]) > 0:
            for e in interaction["experiments"]["direct"]["binary"]:
                experiment_ids.add(e)
        if len(interaction["experiments"]["direct"]["spoke_expanded_from"]) > 0:
            for e in interaction["experiments"]["direct"]["spoke_expanded_from"]:
                experiment_ids.add(e)

    experiments = core_database_connection["experiments"].find({"id": {
        "$in": list(experiment_ids)
    }})

    mapping_regions = list()
    for experiment in experiments:
        for participant in experiment["participants"]:
            if participant['id'] == id:
                for feature in participant["features"]:
                    if (feature["feature_name"] == "binding-associated region" or
                            feature["feature_name"] == "sufficient binding region" or
                            feature["feature_name"] == "direct binding region"):
                        if "(" in feature["featur_value"]:
                            feature_value = feature["featur_value"].split("(")[0]
                        else:
                            feature_value = feature["featur_value"]
                        if '..' in feature_value:
                            feature_value = feature_value.split('..')[1]
                        feature["featur_value"] = feature_value
                        mapping_regions.append(feature)

    for pdb_id in pdb_entities:
        pdb_chains = pdb_entities[pdb_id]
        for chain in pdb_chains:
            for mapping_region in mapping_regions:
                feature_value = mapping_region['featur_value']
                if '?' in feature_value: continue
                if '-' not in feature_value: continue
                region_start = feature_value.split('-')[0]
                region_end = feature_value.split('-')[1]

                if region_start == '?' or region_end == '?':
                    continue

                if int(chain['unp_start']) <= int(region_start) and int(region_end) <= int(chain['unp_end']):
                    if 'mapping_regions' not in chain:
                        chain['binding_regions'] = list()

                    chain['binding_regions'].append({
                        'region_start': region_start,
                        'region_end': region_end
                    })
    return pdb_entities


@app.route('/api/associations/', methods=['POST'])
def get_associations_by_biomolecules():
    try:
        core_database_connection = database_manager.get_primary_connection()
        biomolecule_ids = json.loads(request.data)["biomolecules"]
        interactions = list(core_database_connection["interactions"].find(
            {
                "participants": {
                    "$in": biomolecule_ids
                }
            },
            {
                "_id": 0
            }
        ))

        participants = set()
        for interaction in interactions:
            for participant in interaction["participants"]:
                participants.add(participant)

        participants = list(core_database_connection["biomolecules"].find(
            {
                "id": {
                    "$in": list(participants)
                }
            },
            {
                "_id": 0
            }))
        unique_participants = dict()
        for participant in participants:
            unique_participants[participant["id"]] = participant

        for expression in list(core_database_connection["proteomicsExpression"].find(
            {
                "uniprot": {
                    "$in": list(key for key in unique_participants)
                }
            },
            {
                "_id": 0
            }
        )):
            selected_expressions = list(filter(lambda e: "60-year-old stage (human)" in e["sex"], expression["expressions"]))
            if len(selected_expressions) > 0:
                unique_participants[expression["uniprot"]]["expressions"] = list()
                for selected_expression in selected_expressions:
                    unique_participants[expression["uniprot"]]["expressions"].append({
                        "tissueId": selected_expression["uberonId"],
                        "tpm": selected_expression["tpm"]
                    })

        return json.dumps({
            "interactions": interactions,
            "participants": list(unique_participants.values())
        })
    except:
        traceback.print_tb()
        return json.dumps({
            "interactions": []
        })


@app.route('/api/associations/<id>', methods=['GET'])
def get_association_by_id(id):
    core_database_connection = database_manager.get_primary_connection()
    association = core_database_connection["interactions"].find_one(
        {
            "id": id
        },
        {
            "_id": 0
        }
    )

    if association is not None:
        return {
            "message": "Association found",
            "association": association
        }
    else:
        return {
            "message": "Association not found",
        }


@app.route('/api/experiments/<id>', methods=['GET'])
def get_experiments_by_id(id):
    core_database_connection = database_manager.get_primary_connection()
    meta_data_cache =cache_manager.get_meta_data()
    experiments = core_database_connection["experiments"].find_one(
        {
            "id": id
        },
        {
            "_id": 0
        }
    )

    if experiments is not None:

        for participant in experiments['participants']:
            if 'biological_role' in participant:
                participant['biological_role'] = {
                    'id': meta_data_cache['psimi'][participant['biological_role']]['id'],
                    'name': meta_data_cache['psimi'][participant['biological_role']]['name']
                }
            if 'experimental_role' in participant:
                participant['experimental_role'] = {
                    'id': meta_data_cache['psimi'][participant['experimental_role']]['id'],
                    'name': meta_data_cache['psimi'][participant['experimental_role']]['name']
                }
            if 'identification_method' in participant:
                participant['participant_detection_method'] = {
                    'id': meta_data_cache['psimi'][participant['identification_method']]['id'],
                    'name': meta_data_cache['psimi'][participant['identification_method']]['name']
                }

        experiments['source'] = {
            'id': meta_data_cache['psimi'][experiments['source']]['id'],
            'name': meta_data_cache['psimi'][experiments['source']]['name']
        }
        experiments['interaction_type'] = {
            'id': meta_data_cache['psimi'][experiments['interaction_type']]['id'],
            'name': meta_data_cache['psimi'][experiments['interaction_type']]['name']
        }
        experiments['interaction_detection_method'] = {
            'id': meta_data_cache['psimi'][experiments['interaction_detection_method']]['id'],
            'name': meta_data_cache['psimi'][experiments['interaction_detection_method']]['name']
        }

        return {
            "message": "Association found",
            "experiment": experiments
        }
    else:
        return {
            "message" : "No assoication found for " + id
        }


@app.route('/api/experiments/', methods=['POST'])
def get_experiments_by_ids():
    experiment_ids = json.loads(request.data)["ids"]
    core_database_connection = database_manager.get_primary_connection()
    experiments = core_database_connection["experiments_new"].find(
        {
            "id": {
                '$in': experiment_ids
            }
        },
        {
            "_id": 0
        }
    )

    if experiments is not None:
        return {
            "message": "Association found",
            "experiments": list(e for e in experiments)
        }
    else:
        return {
            "message": "No assoication found for " + id
        }


@app.route('/api/xrefs', methods=['GET'])
def get_xrefs_by_ids():
    ids_to_search = request.args.getlist('id')
    core_database_connection = database_manager.get_primary_connection()
    xrefs = core_database_connection["Keywrds"].find(
        {
            "id": {
                "$in": ids_to_search
            }
        },
        {
            "_id": 0
        }
    )

    return json.dumps(list(xrefs))


def search_query_key():
    query = request.args.get("query")
    hash_obj = hashlib.new('sha256')
    hash_obj.update(query.encode('utf-8'))
    return hash_obj.hexdigest()


@app.route('/api/search', methods=['GET'])
@cache.cached(timeout=120, make_cache_key=search_query_key)
def search_with_text_solr():
    args = request.args
    search_text = args['query']
    search_mode = "0"
    if "mode" in args:
        search_mode = str(args["mode"])

    biomolecules = solr_manager.search_biomolecules(search_text, search_mode)

    publications = list()
    if search_mode == "0":
        publications = solr_manager.search_publications(search_text)

    return json.dumps({
        "biomolecules": biomolecules,
        "publications": publications
    })


def convert_name(name):
    if len(name.split('_')) > 0:
        return name.split('_')[0]
    else:
        return name


@app.route('/api/statistics/', methods=['GET'])
def get_stats():
    core_database_connection = database_manager.get_primary_connection()
    statistics = list(core_database_connection["statistics"].find())
    statistics_reply = dict()
    for statistic in statistics:
        if statistic["category"] == "biomolecule_counts":
            statistics_reply["biomolecules"] = statistic["statistics"]

        if statistic["category"] == "interaction_counts":
            statistics_reply["interactions"] = statistic["statistics"]

        if statistic["category"] == "experiment_counts":
            statistics_reply["experiments"] = statistic["statistics"]
    return json.dumps(statistics_reply)


@app.route('/api/biomolecules/suggestions/<search_query>', methods=['GET'])
def get_biomolcules_suggestions(search_query):
    return {
        "biomolecules": solr_manager.query_biomolecules(search_query)
    }


def network_request_key():
    """A function which is called to derive the key for a computed value.
      The key in this case is the concat value of all the json request
      parameters. Other strategy could to use any hashing function.
    :returns: unique string for which the value should be cached.
    """
    user_data = request.get_json()
    option = True
    if "onlyDirectPartners" in user_data:
        option = user_data["onlyDirectPartners"]
    key = ",".join([f"{value}" for value in user_data['biomolecules']])
    return f'{key}_{option}'


@app.route('/api/network', methods=['POST'])
@cache.cached(timeout=120, make_cache_key=network_request_key)
def generate_network():
    request_body = json.loads(request.data)
    biomolecules = request_body["biomolecules"]
    second_neighborhood = False
    if "onlyDirectPartners" in request_body:
        second_neighborhood = not request_body["onlyDirectPartners"]

    return network_manager.generate_network(biomolecules, second_neighborhood)


@app.route('/api/publications/<pubmed_id>', methods=['GET'])
@cache.cached(timeout=60)
def get_publication_details_by_id(pubmed_id):
    core_database_connection = database_manager.get_primary_connection()
    experiments_by_pubmed = core_database_connection['experiments'].find({
        'pmid': pubmed_id
    })

    if experiments_by_pubmed is None:
        return json.dumps({
            'pubmed': pubmed_id,
            'interactions': []
        })

    experiment_ids = list(experiment['id'] for experiment in experiments_by_pubmed)
    interactions_by_pubmed = core_database_connection['interactions'].find({
        '$or' : [
            {
                'experiments.direct.spoke_expanded_from': {
                    '$in' : experiment_ids
                }
            },
            {
                'experiments.direct.binary': {
                    '$in': experiment_ids
                }
            }
        ]
    })

    interactions_to_return = list()
    participants_to_return = set()
    for interaction in interactions_by_pubmed:
        # Filterout non-relevant evidences from the list of interactions
        relevant_binary_evidences = list()
        relavant_spoke_evidence = list()
        if 'experiments' in interaction:
            if 'direct' in interaction['experiments']:
                if 'binary' in interaction['experiments']['direct']:
                    for binary_evidence in interaction['experiments']['direct']['binary']:
                        if pubmed_id in binary_evidence:
                            relevant_binary_evidences.append(binary_evidence)
                if 'spoke_expanded_from' in interaction['experiments']['direct']:
                    for spoke_evidence in interaction['experiments']['direct']['spoke_expanded_from']:
                        if pubmed_id in spoke_evidence:
                            relavant_spoke_evidence.append(spoke_evidence)

        interaction['experiments']['direct']['binary'] = relevant_binary_evidences
        interaction['experiments']['direct']['spoke_expanded_from'] = relavant_spoke_evidence

        interactions_to_return.append({
            'id': interaction['id'],
            'participants': interaction['participants'],
            'experiments': interaction['experiments'],
            'score': interaction['score']
        })

        for participant in  interaction['participants']:
            participants_to_return.add(participant)

    publication_to_return = {
        'publication': pubmed_id,
        'interactions': interactions_to_return,
        'participants': list(participants_to_return)
    }

    # publication details
    publication = core_database_connection['publications'].find_one({
        'id': pubmed_id
    })

    if publication:
        publication_to_return['title'] = publication['title']
        publication_to_return['authors'] =  publication['authors']
        if 'journal' in publication:
            publication_to_return['journal'] =  publication['journal']
        if 'epubdate' in publication:
            publication_to_return['epubdate'] = publication['epubdate']
        if 'pubdate' in publication:
            publication_to_return['pubdate'] = publication['pubdate']
        publication_to_return['source'] = publication['source']
        publication_to_return['volume'] = publication['volume']
        publication_to_return['issue'] = publication['issue']
        publication_to_return['pages'] = publication['pages']

    return json.dumps(publication_to_return)


@app.route('/api/metadata/<metadata_type>', methods=['GET'])
@cache.cached(timeout=60)
def get_meta_data(metadata_type):
    meta_data_response = dict()
    meta_data = cache_manager.get_meta_data()
    if metadata_type == 'ncbi':
        if 'ncbiTaxonomy' in meta_data:
            meta_data_response['ncbi'] = meta_data['ncbiTaxonomy']

    return meta_data_response


if __name__ == '__main__':
    http_server = WSGIServer(('127.0.0.1', 8000), app)
    print("Server started at 8000")
    http_server.serve_forever()
