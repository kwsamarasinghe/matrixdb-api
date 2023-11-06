import traceback

from flask import Flask, request, Response
from flask import Flask
from flask_caching import Cache
from flask_cors import CORS, cross_origin
from gevent.pywsgi import WSGIServer
from pymongo import MongoClient
import pandas as pd
import json
from functools import reduce
from itertools import groupby

config = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "SimpleCache",  # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 300
}

app = Flask(__name__)
app.config.from_mapping(config)
cache = Cache(app)
#CORS(app, support_credentials=True)

database = None
biomolecule_registry = list()


# Build biomolecule registry
def build_biomolecule_registry():
    biomolecules = list()
    for biomolecule in database["biomolecules"].find():
        biomolecules.append(biomolecule["id"])

    for biomolecule in database["biomolecules_missing"].find():
        biomolecules.append(biomolecule["id"])
    return sorted(biomolecules)


@app.route('/api/biomolecules/<id>', methods=['GET'])
@cache.cached(timeout=50000, query_string=True)
def get_biomolecule_by_id(id):
    biomolecule = database["biomolecules"].find_one({
        "id": id
    },
    {
        "_id": False
    })
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

    biomolecules_from_missing = list(database["biomolecules_new"].find({
        "id": {
            "$in": biomoelcule_ids
        }
    },
        {
            "_id": False
        }))
    if len(biomolecules_from_missing) > 0:
        biomolecules.append(biomolecules_from_missing)

    for biomolecule in biomolecules:
        if biomolecule["type"] == 'protein':

            expressions_by_protein = database["proteinExpression"].find_one({
                "uniprot": biomolecule["id"]
            })

            proteomics_expression_by_protein = database["proteomicsExpressions"].find_one({
                "uniprot": biomolecule["id"]
            })

            gene_expressions = list()
            if expressions_by_protein is not None:
                filtered_expression = filter(lambda e: "life cycle" in e["sex"], expressions_by_protein["expressions"])
                grouped_expressions = groupby(sorted(list(filtered_expression), key=lambda x: x["uberonName"]), key=lambda x: x["uberonName"])

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
    interactions = list(database["interactions_new"].find(
        {
            "participants": id
        },
        {
            "_id": 0
        }
    ))

    neighborhood = {}
    direct = 0
    inferred = 0
    publications = set()
    if len(interactions) > 0:
        for interaction in interactions:
            publications.add(interaction["pubmed"])
            participant_ids = interaction["participants"]
            participant_ids.remove(id)
            partner = participant_ids[0]
            if partner not in neighborhood:
                neighborhood[partner] = {}
                neighborhood[partner]["association"] = interaction["id"]

            if "spoke_expanded_from" in interaction:
                if "spoke_expanded_from" not in neighborhood[partner]:
                    neighborhood[partner]["spokeExpandedFrom"] = list()
                direct += interaction["spoke_expanded_from"]
                neighborhood[partner]["spokeExpandedFrom"].extend(interaction["spoke_expanded_from"])
            else:
                if "directly_supported_by" not in neighborhood[partner]:
                    neighborhood[partner]["directlySupportedBy"] = list()
                direct += len(interaction["experiments"])
                neighborhood[partner]["directlySupportedBy"].extend(interaction["experiments"])

    print()
    return {
        "count": len(neighborhood.keys()),
        "direct": direct,
        "inferred": 0,
        "details": neighborhood,
        "publications": list(publications)
    }


@app.route('/api/biomolecules/proteins/expressions/<id>', methods=['GET'])
#@cache.cached(timeout=50000, query_string=True)
def get_protein_expression(id):

    protein = database["biomolecules"].find_one({
        "id": id
    })

    expressions_by_protein = database["proteinExpression"].find_one({
        "uniprot": id
    })

    proteomics_expression_by_protein = database["proteomicsExpressions"].find_one({
        "uniprot": id
    })

    gene_expressions = list()
    if expressions_by_protein is not None:
        filtered_expression = filter(lambda e: "life cycle" in e["sex"], expressions_by_protein["expressions"])
        grouped_expressions = groupby(sorted(list(filtered_expression), key=lambda x: x["uberonName"]),
                                      key=lambda x: x["uberonName"])

        for tissue_uberon_name, group in grouped_expressions:
            max_group = max(group, key=lambda x: x["tpm"])
            gene_expressions.append({
                "tissueUberonName": tissue_uberon_name.strip('"'),
                "tpm": max_group["tpm"]
            })

    prot_expressions = {}
    if proteomics_expression_by_protein is not None:
        for e in proteomics_expression_by_protein["expressions"]:
            if "sampleName" in e:
                sample_name = e["sampleName"]
                sample = e["sample"]
                tissue_id = e["tissueId"]
                score = e["confidenceScore"]

                if tissue_id not in prot_expressions:
                    prot_expressions[tissue_id] = {
                        "expressionValues": []
                    }

                prot_expressions[tissue_id]["expressionValues"].append({
                    "tissueId": tissue_id,
                    "score": score,
                    "name": sample_name,
                    "sample": sample,
                })

    return {
        "protein": id,
        "gene": protein["relations"]["gene_name"] if "gene_name" in protein["relations"] else None,
        "geneExpression": gene_expressions,
        "proteomicsExpression": prot_expressions
    }


@app.route('/api/associations/', methods=['POST'])
#@cache.cached(timeout=50000, query_string=True)
def get_associations_by_biomolecules():
    try:
        biomolecule_ids = json.loads(request.data)["biomolecules"]
        interactions = list(database["interactions_new"].find(
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

        for expression in list(database["proteinExpression"].find(
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
    association = database["interactions_new"].find_one(
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
    experiments = database["experiments_new"].find_one(
        {
            "id": id
        },
        {
            "_id": 0
        }
    )

    if experiments is not None:
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


def convert_name(name):
    if len(name.split('_')) > 0:
        return name.split('_')[0]
    else:
        return name

@app.route('/api/statistics/associations', methods=['GET'])
@cache.cached(timeout=50000)
def get_interaction_stats():
    associations = database["interactions_new"].find({})

    associations_by_biomolecule_type = {}
    associations_by_biomolecule_type["protein"] = {}
    associations_by_biomolecule_type["protein"]["protein"] = 0
    associations_by_biomolecule_type["protein"]["multimer"] = 0
    associations_by_biomolecule_type["protein"]["gag"] = 0
    associations_by_biomolecule_type["protein"]["pfrag"] = 0
    associations_by_biomolecule_type["multimer"]  = {}
    associations_by_biomolecule_type["multimer"]["multimer"] = 0
    associations_by_biomolecule_type["multimer"]["pfrag"] = 0
    associations_by_biomolecule_type["multimer"]["gag"] = 0
    associations_by_biomolecule_type["pfrag"] = {}
    associations_by_biomolecule_type["pfrag"]["pfrag"] = 0
    associations_by_biomolecule_type["pfrag"]["gag"] = 0
    associations_by_biomolecule_type["gag"] = {}
    associations_by_biomolecule_type["gag"]["gag"] = 0

    # Same association is duplicated as to distinguish the experiments
    # Furthermore, there can be symmetric form of the association
    # Need to get the unique biomolecule pairs removing the symmetry and then count

    sorted_associations = set()
    all_assoc = 0
    for association in list(associations):
        all_assoc+= 1
        if "inferred_from" in association and len(association["inferred_from"]) > 0:
            continue

        biomolecules = association["biomolecules"]
        if len(biomolecules) > 1:
            if biomolecules[0] > biomolecules[1]:
                sorted_associations.add((biomolecules[1], biomolecules[0], association["id"]))
            else:
                sorted_associations.add((biomolecules[0], biomolecules[1], association["id"]))
        else:
            sorted_associations.add((biomolecules[0],association["id"]))

    print("All assoc " + str(all_assoc))
    print("sorted " + str(len(sorted_associations)))
    # Count associations by biomolecule pairs, i.e protein - protein, protein - gag , etc
    biomol_cache = {}
    all_exp_supported = 0
    for biomolecules in list(sorted_associations):

        if len(biomolecules) > 2:
            interactant1 = biomolecules[0]
            interactant2 = biomolecules[1]
            association_id = biomolecules[2]
        else:
            interactant1 = biomolecules[0]
            interactant2 = biomolecules[0]
            association_id = biomolecules[1]

        if interactant1 not in biomol_cache:
            biomolecule = database["cleaned_biomolecules"].find_one({
                "id": interactant1
            })
            if biomolecule is None:
                #print("No biomol " + " " + interactant1)
                continue
            interactant1_type = biomolecule["type"].lower()
            if interactant1_type is None:
                print("No type? " + interactant1)
            biomol_cache[interactant1] = interactant1_type
        else:
            interactant1_type = biomol_cache[interactant1].lower()
            if interactant1_type is None:
                print("No type?" + interactant1)

        if interactant2 not in biomol_cache:
            biomolecule = database["cleaned_biomolecules"].find_one({
                "id": interactant2
            })
            if biomolecule is None:
                #print("No biomol "  + interactant2)
                continue
            interactant2_type = biomolecule["type"].lower()
            if interactant2_type is None:
                print("No type " + interactant2)
            biomol_cache[interactant2] = interactant2_type
        else:
            interactant2_type = biomol_cache[interactant2].lower()
            if interactant2_type is None:
                print("No type " + interactant2)

        all_exp_supported += 1
        if "protein" in interactant1_type and "protein" in interactant2_type:
            print(interactant1 + "," + interactant2+"," + association_id +",protein-protien")
            associations_by_biomolecule_type["protein"]["protein"] += 1

        elif "protein" in interactant1_type and "gag" in interactant2_type:
            print(interactant1 + "," + interactant2+",protein-gag")
            associations_by_biomolecule_type["protein"]["gag"] += 1

        elif "gag" in interactant1_type and "protein" in interactant2_type:
            print(interactant1 + "," + interactant2+",gag-protien")
            associations_by_biomolecule_type["protein"]["gag"] += 1

        elif "protein" in interactant1_type and "mult" in interactant2_type:
            print(interactant1 + "," + interactant2+",protein-mult")
            associations_by_biomolecule_type["protein"]["multimer"] += 1

        elif 'mult' in interactant1_type and "protein" in interactant2_type:
            print(interactant1 + "," + interactant2+",mult-protine")
            associations_by_biomolecule_type["protein"]["multimer"] += 1

        elif "protein" in interactant1_type and "pfrag" in interactant2_type:
            print(interactant1 + "," + interactant2+",protein-pfrag")
            associations_by_biomolecule_type["protein"]["pfrag"] += 1

        elif "pfrag" in interactant1_type and "protein" in interactant2_type:
            print(interactant1 + "," + interactant2+",pfrag-protine")
            associations_by_biomolecule_type["protein"]["pfrag"] += 1

        elif "mult" in interactant1_type and "mult" in interactant2_type:
            print(interactant1 + "," + interactant2+",mult-mult")
            associations_by_biomolecule_type["multimer"]["multimer"] += 1

        elif "mult" in interactant1_type and "pfrag" in interactant2_type:
            print(interactant1 + "," + interactant2+",mult-pfrag")
            associations_by_biomolecule_type["multimer"]["pfrag"] += 1

        elif "pfrag" in interactant1_type and "mult" in interactant2_type:
            print(interactant1 + "," + interactant2+",pfrag-mult")
            associations_by_biomolecule_type["multimer"]["pfrag"] += 1

        elif "mult" in interactant1_type and "gag" in interactant2_type:
            print(interactant1 + "," + interactant2+",mult-gag")
            associations_by_biomolecule_type["multimer"]["gag"] += 1

        elif "gag" in interactant1_type and "mult" in interactant2_type:
            print(interactant1 + "," + interactant2+",gag-mult")
            associations_by_biomolecule_type["multimer"]["gag"] += 1

        elif "pfrag" in interactant1_type and "pfrag" in interactant2_type:
            print(interactant1 + "," + interactant2+",pfrag-pfrag")
            associations_by_biomolecule_type["pfrag"]["pfrag"] += 1

        elif "pfrag" in interactant1_type and "gag" in interactant2_type:
            print(interactant1 + "," + interactant2+",pfrag-gag")
            associations_by_biomolecule_type["pfrag"]["gag"] += 1

        elif "gag" in interactant1_type and "pfrag" in interactant2_type:
            print(interactant1 + "," + interactant2+",gag-pfrag")
            associations_by_biomolecule_type["pfrag"]["gag"] += 1

        elif "gag" in interactant1_type and "gag" in interactant2_type:
            print(interactant1 + "," + interactant2+",gag-gag")
            associations_by_biomolecule_type["gag"]["gag"] += 1

    return json.dumps({
        "associations": associations_by_biomolecule_type
    })


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
@cache.cached(timeout=50000)
def get_biomolcules_suggestions(search_query):
    if search_query in biomolecule_registry:
        return {
            "suggestions": [search_query]
        }

    index = next((i for i, item in enumerate(biomolecule_registry) if item.lower().startswith(search_query.lower())),
                 None)

    if index is not None:
        biomolecules_found = biomolecule_registry[index + 1:index + 11]
        return {
            "suggestions": biomolecules_found
        }
    else:
        return {
            "suggestions": []
        }


@app.route('/api/networks', methods=['POST'])
def generate_network():
    biomolecules = json.loads(request.data)["biomolecules"]
    # Get the associations of biomolecuels
    neighbors_by_biomolecule = dict()
    associations = list()
    participants = set()

    for biomolecule in biomolecules:
        participants.add(biomolecule)

        interactions = database["interactions_new"].find({
            "participants": biomolecule
        })
        for interaction in interactions:
            if interaction["id"] == "GAG_2__P08581":
                print()
            associations.append({
                "id": interaction["id"],
                "participants": interaction["participants"],
                "experiments": interaction["experiments"]
            })
            if biomolecule not in neighbors_by_biomolecule:
                neighbors_by_biomolecule[biomolecule] = set()

            # Check
            for participant in interaction["participants"]:
                if participant is not biomolecule:
                    neighbors_by_biomolecule[biomolecule].add(participant)
                    participants.add(participant)

    # Compute the interaction
    if len(biomolecules) > 1 and len(neighbors_by_biomolecule.keys()) > 1:
        common_neighbors = reduce(lambda x, y: x.intersection(y), list(neighbors_by_biomolecule.values()),
                                  list(neighbors_by_biomolecule.values())[0])

        interactions = database["interactions_new"].find({
            "participants": {
                "$in": list(common_neighbors)
            }
        })
        for interaction in interactions:
            associations.append({
                "id": interaction["id"],
                "participants": interaction["participants"],
                "experiments": interaction["experiments"]
            })
            for participant in interaction["participants"]:
                participants.add(participant)

    # Remove duplicates, must be fixed in data
    unique_associations = dict()
    for association in associations:
        unique_associations[association["id"]] = association

    return json.dumps({
        "associations": list(unique_associations.values()),
        "participants": list(participants)
    })


if __name__ == '__main__':

    # Connects to the db
    database_url = "mongodb://localhost:27018/"
    datasets = {}
    try:
        database_client = MongoClient(database_url)
        database = database_client["matrixdb-pre-prod"]
    except Exception:
        print("Problem connecting to db " + database_url)

    print("Building biomolecule registry")
    biomolecule_registry = build_biomolecule_registry()

    # Serve the src with gevent
    http_server = WSGIServer(('127.0.0.1', 8000), app)
    print("Server started at 8000")
    http_server.serve_forever()
