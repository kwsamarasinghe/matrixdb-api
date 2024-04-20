import re
import traceback

from flask import Flask, request, Response
from flask_caching import Cache

from gevent.pywsgi import WSGIServer
from pymongo import MongoClient
import pandas as pd
import json
from functools import reduce
from itertools import groupby

from src.matrixdb.biomolecule_services.protein_data_manager import ProteinDataManager
from src.matrixdb.interactome.network_manager import NetworkManager
from src.matrixdb.utils.solr.solr_query_controller import query_solr

config = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "SimpleCache",  # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 300
}

app = Flask(__name__)
app.config.from_mapping(config)
cache = Cache(app)

database = None
biomolecule_registry = list()
meta_data_cache = {
    "psimi": dict(),
    "go": dict(),
    "interpro": dict(),
    "uniprotKeywords": dict(),
    "uberon": dict(),
    "bto": dict()
}


# Build biomolecule registry
def build_biomolecule_registry():
    biomolecules = list()
    for biomolecule in database["biomolecules"].find():
        biomolecules.append(biomolecule["id"])

    return sorted(biomolecules)


def build_meta_data_cache():
    # psimi
    for psimi in database["psimi"].find():
        meta_data_cache["psimi"][psimi["id"]] = psimi

    # go
    for go in database["go"].find():
        meta_data_cache["go"][go["id"]] = go

    # Interpro
    for interpro in database["interpro"].find():
        meta_data_cache["interpro"][interpro["id"]] = interpro

    # Uniprot keywords
    for uniprot_keyword in database["uniprotKeywords"].find():
        meta_data_cache["uniprotKeywords"][uniprot_keyword["id"]] = uniprot_keyword

    for uberon in database["uberon"].find():
        meta_data_cache["uberon"][uberon["id"]] = uberon

    for bto in database["brenda"].find():
        meta_data_cache["bto"][bto["id"]] = bto


@app.route('/api/biomolecules/<id>', methods=['GET'])
@cache.cached(timeout=50000, query_string=True)
def get_biomolecule_by_id(id):
    biomolecule = database["biomolecules"].find_one({
        "id": id
    },
    {
        "_id": False
    })

    # Include GO, uniprot keyword and interpro definitions
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

    if biomolecule["type"] == 'protein':
        keywords = list()
        if "annotations" in biomolecule and "keywords" in biomolecule["annotations"]:
            for keyword in biomolecule["annotations"]["keywords"]:
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
                                'name': meta_data_cache["interpro"][interpro]["name"]
                            })
                biomolecule["xrefs"]["interpro"] = interpro_terms

    return Response(json.dumps(biomolecule), mimetype='application/json')


@app.route('/api/biomolecules/', methods=['POST'])
def get_biomolecules_by_id():
    biomoelcule_ids = json.loads(request.data)["ids"]
    biomolecules = list(database["biomolecules"].find({
        "id": {
            "$in": biomoelcule_ids
        }
    },
    {
        "_id": False
    }))

    for biomolecule in biomolecules:
        if biomolecule["type"] == 'protein':

            expressions_by_protein = database["geneExpression"].find_one({
                "uniprot": biomolecule["id"]
            })

            proteomics_expression_by_protein = database["proteomicsExpressions"].find_one({
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
    interactions = list(database["interactions"].find(
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
    '''
    database_url = "mongodb://localhost:27018/"
    try:
        database_client = MongoClient(database_url)
        database = database_client["matrixdb-4_0-pre-prod"]
    except Exception:
        print("Problem connecting to db " + database_url)

    protein_ids = json.loads(request.data)
    proteins = database["biomolecules"].find({
        "id": {
            '$in': protein_ids
        }
    })

    expression_data = dict()
    for protein in proteins:
        expressions_by_protein = database["geneExpression"].find_one({
            "uniprot": protein['id']
        })

        proteomics_expression_by_protein = database["proteomicsExpression"].find_one({
            "uniprot": protein['id']
        })

        gene_expressions = list()
        if expressions_by_protein is not None:

            expression_values = list()
            for expression in expressions_by_protein["expressions"]:
                for e in list(expression.values()):
                    expression_values.append(e[0])

            # filtered_expression = filter(lambda e: "life cycle" in e["sex"], expressions_by_protein["expressions"])
            grouped_expressions = groupby(sorted(list(expression_values), key=lambda x: x["uberonName"]),
                                          key=lambda x: x["uberonName"])

            for tissue_uberon_name, group in grouped_expressions:
                max_group = max(group, key=lambda x: x["tpm"])
                gene_expressions.append({
                    "tissueUberonName": tissue_uberon_name.strip('"'),
                    "tpm": max_group["tpm"]
                })

        prot_expressions = list()
        if proteomics_expression_by_protein is not None:
            for e in proteomics_expression_by_protein["expressions"]:
                prot_expressions.append({
                    "tissueId": e["tissueId"],
                    "name": e["sampleName"] if "sampleName" in e else "",
                    "sample": e["sample"] if "sample" in e else "",
                    "score": e["confidenceScore"] if "confidenceScore" in e else 0,
                })

        expression_data[protein['id']] ={
            "gene": protein["relations"]["gene_name"] if protein is not None and "gene_name" in protein["relations"] else None,
            "geneExpression": gene_expressions,
            "proteomicsExpression": prot_expressions
        }

    '''
    protein_data_manager = ProteinDataManager(meta_data_cache)
    expression_data = dict()
    for protein in json.loads(request.data):
        expression_data[protein] = {
            'geneExpression': protein_data_manager.get_gene_expressions(protein),
            'proteomicsExpression': protein_data_manager.get_proteomics_expressions(protein)
        }
    return expression_data


@app.route('/api/associations/', methods=['POST'])
def get_associations_by_biomolecules():
    try:
        biomolecule_ids = json.loads(request.data)["biomolecules"]
        interactions = list(database["interactions"].find(
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

        participants = list(database["biomolecules"].find(
            {
                "id": {
                    "$in": list(participants)
                }
            },
            {
                "_id": 0
            }))
        unique_participants = dict()
        for p in participants:
            unique_participants[p["id"]] = p

        for expression in list(database["proteomicsExpression"].find(
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
    association = database["interactions"].find_one(
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
    experiments = database["experiments"].find_one(
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
                participant['biological_role'] = meta_data_cache['psimi'][participant['biological_role']]['name']
            if 'experimental_role' in participant:
                participant['experimental_role'] = meta_data_cache['psimi'][participant['experimental_role']]['name']
            if 'participant_detection_method' in participant:
                participant['participant_detection_method'] = meta_data_cache['psimi'][participant['participant_detection_method']]['name']

        experiments['source'] = meta_data_cache['psimi'][experiments['source']]['name']
        experiments['interaction_type'] = meta_data_cache['psimi'][experiments['interaction_type']]['name']
        experiments['interaction_detection_method'] = meta_data_cache['psimi'][experiments['interaction_detection_method']]['name']


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

    experiments = database["experiments_new"].find(
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
    xrefs = database["Keywrds"].find(
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

'''
@app.route('/api/search', methods=['GET'])
@cache.cached(timeout=50000, query_string=True)
def search_with_text():
    args = request.args
    search_text = args['text']
    # Search biomolecules for the moment
    biomolecules_by_id = list(database["biomolecules"].find(
        {
            "$or": [
                {
                    "id":  search_text
                },
                {
                    "id": search_text.lower()
                },
                {
                    "names.name": {"$regex": search_text, '$options': 'i'}
                },
                {
                    "names.common_name": {"$regex": search_text, '$options': 'i'}
                },
            ]
        },
        {
            "_id" : 0
        }
    ))

    if len(biomolecules_by_id) < 10:
        biomolecules_by_other = database["biomolecules"].find(
            {

                    "$or":[
                        {
                            "id": { "$regex" : search_text, '$options': 'i' }
                        },
                        {
                            "names.name": {"$regex": search_text, '$options': 'i'}
                        },
                        {
                            "names.common_name": {"$regex": search_text, '$options': 'i'}
                        },
                        {
                            "names.other_name": {"$regex": search_text, '$options': 'i'}
                        },
                        {
                            "names.scientific_name": {"$regex": search_text, '$options': 'i'}
                        },
                        {
                            "relations.belongs_to": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "description": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "annotations.keywords": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "moleculer_details": {"$regex": search_text, '$options': 'i'}
                        },
                        {
                            "xrefs.chebi": search_text
                        },
                        {
                            "xrefs.uniprot": search_text
                        },
                        {
                            "xrefs.kegg": search_text
                        }
                    ]

            },
            {
                "_id": 0
            }
        )
        biomolecules_by_id.extend(biomolecules_by_other)

    results = []
    for biomolecule in list(biomolecules_by_id):
        result = dict()
        result["id"] = biomolecule["id"]
        result["names"] = biomolecule["names"]
        if "description" in biomolecule:
            result["description"] = biomolecule["description"]
        if "xrefs" in biomolecule:
            result["xrefs"] = biomolecule["xrefs"]
        results.append(result)

    return json.dumps(results)
'''

@app.route('/api/search', methods=['GET'])
def search_with_text_solr():
    args = request.args
    search_text = args['text']

    biomolecules_core_url = 'http://localhost:8983/solr/biomolecules'
    publications_core_url = 'http://localhost:8983/solr/publications'

    biomolecule_query_params = {
        'q': f'{search_text}',
        'defType': 'dismax',
        'qf': 'biomolecule_id^10.0 name^5 common_name^4 recommended_name^3 other_name^2 species^2 description^2 chebi complex_portal go_names go_ids keyword_ids keyword_names',
        'rows': 100
    }
    biomolecule_solr_docs = query_solr(biomolecules_core_url, biomolecule_query_params)
    biomolecule_solr_docs = sorted(biomolecule_solr_docs, key=lambda doc: doc['interaction_count'])

    publication_query_params = {
        'q': f'{search_text}',
        'defType': 'dismax',
        'qf': 'publication_id title authors journal',
        'rows': 100
    }
    publication_solr_docs = query_solr(publications_core_url, publication_query_params)
    publication_solr_docs = sorted(publication_solr_docs, key=lambda doc: doc['interaction_count'])

    return json.dumps({
        "biomolecules": biomolecule_solr_docs,
        "publications": publication_solr_docs
    })

def convert_name(name):
    if len(name.split('_')) > 0:
        return name.split('_')[0]
    else:
        return name

@app.route('/api/statistics/', methods=['GET'])
def get_stats():
    statistics = list(database["statistics"].find())
    statistics_reply = dict()
    for statistic in statistics:
        if statistic["category"] == "biomolecule_counts":
            statistics_reply["biomolecules"] = statistic["statistics"]

        if statistic["category"] == "interaction_counts":
            statistics_reply["interactions"] = statistic["statistics"]

        if statistic["category"] == "experiment_counts":
            statistics_reply["experiments"] = statistic["statistics"]
    return json.dumps(statistics_reply)


@app.route('/api/statistics/biomolecules', methods=['GET'])
@cache.cached(timeout=50000)
def get_biomol_stats():
    biomols = list(database["biomolecules"].find({}))

    # Group by type
    df = pd.DataFrame(biomols)
    biomol_groups = df.groupby('type').size().to_dict()

    biomol_types = {b["id"] : b["type"] for b in biomols}

    # For proteins check which associations have at least one ECM
    matrisome_keywords = database["matrisome"].find({})
    matrisome_dict = {matrisome_keyword["protein"]: matrisome_keyword["keyword"] for matrisome_keyword in
                      list(matrisome_keywords)}

    associations = database["cleaned_associations"].find({})
    biomols_in_associations = dict()

    ecm_stats = {
        "both": 0,
        "one": 0,
        "neither":  0
    }

    for association in list(associations):
        if len(association["biomolecules"]) != 2:
            print("non binary")
            b0 = association["biomolecules"][0]
            type_0 = biomol_types[b0]
            if b0 in matrisome_dict:
                ecm_stats["one"] += 1

            if type_0 not in biomols_in_associations:
                biomols_in_associations[type_0] = 0
            continue

        b0 = association["biomolecules"][0]
        b1 = association["biomolecules"][1]
        if b0 in biomol_types:
            type_0 = biomol_types[b0]

        if b1 in biomol_types:
            type_1 = biomol_types[b1]


        if type_0 is not None:
            if type_0 not in biomols_in_associations:
                biomols_in_associations[type_0] = set()
            biomols_in_associations[type_0].add(b0)

            if type_0 == 'protein':
                b0_obj = list(filter(lambda b: b["id"] == b0, biomols))[0]
                if "annotations" in b0_obj:
                    if "subcell_location" in b0_obj["annotations"]:
                        print(b0_obj["annotations"]["subcell_location"])
                    if "keywords" in b0_obj["annotations"]:
                        print(b0_obj["annotations"]["keywords"])
                else:
                    print(b0 + " No annotations")

        if type_1 is not None:
            if type_1 not in biomols_in_associations:
                biomols_in_associations[type_1] = set()
            biomols_in_associations[type_1].add(b1)

            if type_1 == 'protein':
                b1_obj = list(filter(lambda b: b["id"] == b1, biomols))[0]
                if "annotations" in b1_obj:
                    if "subcell_location" in b1_obj["annotations"]:
                        print(b1_obj["annotations"]["subcell_location"])
                    if "keywords" in b1_obj["annotations"]:
                        print(b1_obj["annotations"]["keywords"])
                else:
                    print(b1_obj)

        if b0 in matrisome_dict and b1 in matrisome_dict:
            ecm_stats["both"] += 1

        elif b0 not in matrisome_dict and b1 not in matrisome_dict:
            ecm_stats["neither"] += 1
        else:
            ecm_stats["one"] += 1

    for biomol_type in biomols_in_associations:
        biomols_in_associations[biomol_type] = len(biomols_in_associations[biomol_type])

    return json.dumps({
        "biomolecules": biomol_groups,
        "in_assocs" : biomols_in_associations,
        "ecm_stats" : ecm_stats
    })


@app.route('/api/statistics/experiments', methods=['GET'])
@cache.cached(timeout=50000)
def get_experiment_stats():
    # Count experiments
    experiments = list(database["experiments_new"].find({}))

    all_exps_supported = 0
    for experiment in experiments:
        if "directly_supports" in experiment:
            all_exps_supported += len(experiment["directly_supports"])

        if "spoke_expanded_into" in experiment:
            all_exps_supported += len(experiment["spoke_expanded_into"])

    # Exps and assocs with source matrixdb
    matrixdb_experiments = [me for me in experiments if "source" in me and me["source"] == "MI:0917"]
    matrixdb_experiments_supporting_associations = set()
    for me in matrixdb_experiments:
        if "directly_supports" in me:
            for exp in me["directly_supports"]:
                matrixdb_experiments_supporting_associations.add(exp)

        if "spoke_expanded_into" in me:
            for exp in me["spoke_expanded_into"]:
                matrixdb_experiments_supporting_associations.add(exp)

    experiment_counts = {
        "total": {
            "experiments": len(experiments),
            "experimentally_supported_assocs": all_exps_supported
        },
        "matrixdb_core": {
            "experiments": len(matrixdb_experiments),
            "experimentally_supported_assocs": len(list(matrixdb_experiments_supporting_associations))
        }
    }

    return json.dumps({
        "experiments": experiment_counts
    })


@app.route('/api/biomolecules/suggestions/<search_query>', methods=['GET'])
def get_biomolcules_suggestions(search_query):
    biomolecules_core_url = 'http://localhost:8983/solr/biomolecules'

    biomolecule_query_params = {
        'q': '*:*',
        'qf': 'biomolecule_id^10.0 name^5 common_name^4 recommended_name^3 description^2 keywords',
        'fq': f'biomolecule_id:*{search_query}* OR name:*{search_query}* OR common_name:*{search_query}* OR '
              f'recommended_name:*{search_query}* OR description:*{search_query}* OR keywords:*{search_query}*',
        'rows': 1000
    }
    biomolecule_solr_docs = query_solr(biomolecules_core_url, biomolecule_query_params)
    biomolecule_solr_docs = sorted(biomolecule_solr_docs, key=lambda doc: doc['interaction_count'], reverse=True)

    if biomolecule_solr_docs is not None:
        suggestions = []
        for bd in biomolecule_solr_docs:
            if 'name' not in bd:
                continue
            '''
            suggestion = {
                'id': bd['biomolecule_id'][0],
                'name': bd['name'][0],
            }
            if 'xrefs' in bd:
                suggestion['xref'] = bd['xrefs'][0]
            '''


            suggestions.append(bd['name'])
        return {
            "suggestions": suggestions
        }
    else:
        return {
            "suggestions": []
        }


def network_request_key():
   """A function which is called to derive the key for a computed value.
      The key in this case is the concat value of all the json request
      parameters. Other strategy could to use any hashing function.
   :returns: unique string for which the value should be cached.
   """
   user_data = request.get_json()
   return ",".join([f"{value}" for value in user_data['biomolecules']])

@app.route('/api/network', methods=['POST'])
@cache.cached(timeout=60, make_cache_key=network_request_key)
def generate_network():
    biomolecules = json.loads(request.data)["biomolecules"]
    return network_manager.generate_network(biomolecules)

@app.route('/api/publications/<pubmed_id>', methods=['GET'])
@cache.cached(timeout=60)
def get_publication_details_by_id(pubmed_id):
    experiments_by_pubmed = database['experiments'].find({
        'pmid': pubmed_id
    })

    if experiments_by_pubmed is None:
        return json.dumps({
            'pubmed': pubmed_id,
            'interactions': []
        })

    experiment_ids = list(experiment['id'] for experiment in experiments_by_pubmed)
    interactions_by_pubmed = database['interactions'].find({
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

    publication_to_return = {
        'publication': pubmed_id,
        'interactions': interactions_to_return
    }

    # publication details
    publication = database['publications'].find_one({
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

    return json.dumps(publication_to_return)


if __name__ == '__main__':

    # Connects to the db
    database_url = "mongodb://localhost:27018/"
    datasets = {}
    try:
        database_client = MongoClient(database_url)
        database = database_client["matrixdb_4_0"]
    except Exception:
        print("Problem connecting to db " + database_url)

    print("Building biomolecule registry")
    biomolecule_registry = build_biomolecule_registry()

    print("Building meta data cache")
    build_meta_data_cache()

    network_manager = NetworkManager(database_connection=database, meta_data_cache=meta_data_cache)

    # Serve the src with gevent
    http_server = WSGIServer(('127.0.0.1', 8000), app)
    print("Server started at 8000")
    http_server.serve_forever()
