from flask import Flask, request, Response
from flask import Flask
from flask_caching import Cache
from flask_cors import CORS, cross_origin
from gevent.pywsgi import WSGIServer
from pymongo import MongoClient
import pandas as pd
import json

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

@app.route('/api/biomolecules/<id>', methods=['GET'])
def get_biomolecule_by_id(id):
    biomolecule = database["biomolecules"].find_one({
        "id": id
    },
    {
        "_id": False
    })
    return Response(json.dumps(biomolecule), mimetype='application/json')

@app.route('/api/biomolecules/<id>/interactors/', methods=['GET'])
def get_biomolecule_interactors_by_id(id):
    '''
    assocations = database["associations"].find(
        {
            "biomolecule": id
        },
        {
            "_id" :0
        }
    )
    neighborhood_map = {}
    for assocation in list(assocations):
        #print(id + " " + assocation["id"])
        assocation['biomolecule'].remove(id)
        neighbor_id = assocation['biomolecule'][0]
        if neighbor_id not in neighborhood_map:
            neighborhood_map[neighbor_id] = {}
        neighbor = neighborhood_map[neighbor_id]


        if "directlysupportedby" in assocation:
            if assocation["directlysupportedby"] is not None:
                if "directlysupportedby" not in neighbor:
                    neighbor["directlysupportedby"] = []
                if type(assocation["directlysupportedby"]).__name__ != 'list':
                    neighbor["directlysupportedby"].append(assocation["directlysupportedby"])
                else:
                    neighbor["directlysupportedby"] += assocation["directlysupportedby"]

        if "spokeexpandedfrom" in assocation:
            if assocation["spokeexpandedfrom"] is not None:
                if "spokeexpandedfrom" not in neighbor:
                    neighbor["spokeexpandedfrom"] = []
                if type(assocation["spokeexpandedfrom"]).__name__ != 'list':
                    neighbor["spokeexpandedfrom"].append(assocation["spokeexpandedfrom"])
                else:
                    neighbor["spokeexpandedfrom"] += assocation["spokeexpandedfrom"]

        if "inferredfrom" in assocation:
            if assocation["inferredfrom"] is not None:
                if "inferredfrom" not in neighbor:
                    neighbor["inferredfrom"] = []
                if type(assocation["inferredfrom"]).__name__ != 'list':
                    neighbor["inferredfrom"].append(assocation["inferredfrom"])
                else:
                    neighbor["inferredfrom"] += assocation["inferredfrom"]

        # Extract inferred from experiments

    direct = 0
    inferred = 0
    for k in neighborhood_map.keys():
        if "directlysupportedby" in neighborhood_map[k] or "spokeexpandedfrom" in neighborhood_map[k]:
            direct += 1
        else:
            inferred += 1


    return json.dumps({
        "count": len(neighborhood_map.keys()),
        "direct": direct,
        "inferred": inferred,
        "interactors": neighborhood_map
    })
    '''
    associations = database["associationsByBiomolecule"].find_one(
        {
            "biomolecule": id
        },
        {
            "_id": 0
        }
    )

    neighborhood = {}
    direct = 0
    inferred = 0
    if associations is not None:
        for association_id in associations["associations"]:
            association = associations["associations"][association_id]
            if len(association["directlysupportedby"]) > 0 or len(association["spokeexpandedfrom"]) > 0:
                direct = direct + 1
            #if len(association["inferredfrom"]) > 0:
            else:
                inferred = inferred + 1

            neighborhood[association["partner"]] = {
                "association" : association_id,
                "directlysupportedby" : association["directlysupportedby"],
                "spokeexpandedfrom" : association["spokeexpandedfrom"],
                "inferredfrom" : association["inferredfrom"]
            }

    parterns_from_matrixdb_3_5 = len(neighborhood.keys())

    # Build a neighborhood map from intact associations
    intact_associations = database["intactAssociationsByBiomolecule"].find_one(
        {
            "biomolecule": id
        },
        {
            "_id": 0
        }
    )
    if intact_associations is not None:
        for association_id in intact_associations["associations"]:
            intact_association = intact_associations["associations"][association_id]
            if intact_association["partner"] not in neighborhood:
                neighborhood[intact_association["partner"]] = {}
                neighborhood[intact_association["partner"]]["association"] = association_id
                neighborhood[intact_association["partner"]]["direct_from_intact"] = []
                neighborhood[intact_association["partner"]]["direct_from_intact"] = intact_association["directlysupportedby"]
            else:
                neighborhood[intact_association["partner"]]["direct_from_intact"] = []
                new_from_intact = list(set(intact_association["directlysupportedby"]).difference(
                    set(neighborhood[intact_association["partner"]]["directlysupportedby"])))
                if len(new_from_intact) > 0:
                    neighborhood[intact_association["partner"]]["direct_from_intact"] = new_from_intact

    all_partners = len(neighborhood.keys())
    return {
        "count" : len(neighborhood.keys()),
        "countMatrixdb": parterns_from_matrixdb_3_5,
        "countIntact": all_partners - parterns_from_matrixdb_3_5,
        "direct": direct,
        "inferred": inferred,
        "interactors": neighborhood
    }

@app.route('/api/associations/<id>', methods=['GET'])
def get_association_by_id(id):
    association = database["associations"].find_one(
        {
            "id": {
                "$regex": id,
                "$options": "i"
            }
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
        association = database["intactAssociations"].find_one(
            {
                "id": {
                    "$regex": id,
                    "$options": "i"
                }
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
                "message" : "No assoication found for " + id
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
        return {
            "message": "Association found",
            "experiment": experiments
        }
    else:
        experiments = database["intactExperiments"].find_one(
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
    biomolecules_by_id = database["cleaned_biomolecules"].find(
        {
            "or" : [
                {
                    "id":  search_text
                },
                {
                    "id": search_text.lower()
                }
            ]
        },
        {
            "_id" : 0
        }
    )

    if biomolecules_by_id.retrieved == 0:
        biomolecules_by_id = database["cleaned_biomolecules"].find(
            {

                    "$or":[
                        {
                            "id": { "$regex" : search_text, '$options' : 'i' }
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
                            "annotation.keywords": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "moleculer_details": {"$regex": search_text, '$options': 'i'}
                        },
                        {
                            "xref.chebi": search_text
                        },
                        {
                            "xref.uniprot": search_text
                        },
                        {
                            "xref.kegg": search_text
                        }
                    ]

            },
            {
                "_id": 0
            }
        )

    response = list({"id": b["id"] } for b in list(biomolecules_by_id))

    return json.dumps(response)


def convert_name(name):
    if len(name.split('_')) > 0:
        return name.split('_')[0]
    else:
        return name

@app.route('/api/statistics', methods=['GET'])
@cache.cached(timeout=50000)
def get_statistics():
    associations_by_biomolecules = database["associationsByBiomolecule"].find({})

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

    # Count associations by biomolecule pairs, i.e protein - protein, protein - gag , etc
    biomol_cache = {}
    for association in list(associations_by_biomolecules):
        biomolecule_id = association["biomolecule"]

        if biomolecule_id not in biomol_cache:
            biomolecule = database["biomolecules"].find_one({
                "id": biomolecule_id
            })
            if biomolecule is None:
                print()
            biomolecule_type = biomolecule["type"].lower()
            if biomolecule_type is None:
                print()
            biomol_cache[biomolecule_id] = biomolecule_type
        else:
            biomolecule_type = biomol_cache[biomolecule_id].lower()
            if biomolecule_type is None:
                print()

        for association_id in association["associations"]:
            if len(association["associations"][association_id]["directlysupportedby"]) == 0  and len(association["associations"][association_id]["spokeexpandedfrom"]) == 0:
                continue

            partner_id = association["associations"][association_id]["partner"]
            if partner_id is None:
                print()

            if partner_id not in biomol_cache:
                partner = database["biomolecules"].find_one({
                    "id": partner_id
                })

                if partner is None:
                    print()

                if "type" in partner:
                    partner_type = partner["type"].lower()

                if partner_type is None:
                    print()

                biomol_cache[partner_id] = partner_type
            else:
                partner_type = biomol_cache[partner_id].lower()
                if partner_type is None:
                    print()

            if "protein" in biomolecule_type and "protein" in partner_type:
                associations_by_biomolecule_type["protein"]["protein"] += 1

            elif "protein" in biomolecule_type and "gag" in partner_type:
                associations_by_biomolecule_type["protein"]["gag"] += 1

            elif "gag" in biomolecule_type and "protein" in partner_type:
                associations_by_biomolecule_type["protein"]["gag"] += 1

            elif "protein" in biomolecule_type and "mult" in partner_type:
                associations_by_biomolecule_type["protein"]["multimer"] += 1

            elif 'mult' in biomolecule_type and "protein" in partner_type:
                associations_by_biomolecule_type["protein"]["multimer"] += 1

            elif "protein" in biomolecule_type and "pfrag" in partner_type:
                associations_by_biomolecule_type["protein"]["pfrag"] += 1

            elif "pfrag" in biomolecule_type and "protein" in partner_type:
                associations_by_biomolecule_type["protein"]["pfrag"] += 1

            elif "mult" in biomolecule_type and "mult" in partner_type:
                associations_by_biomolecule_type["multimer"]["multimer"] += 1

            elif "mult" in biomolecule_type and "pfrag" in partner_type:
                associations_by_biomolecule_type["multimer"]["pfrag"] += 1

            elif "pfrag" in biomolecule_type and "mult" in partner_type:
                associations_by_biomolecule_type["multimer"]["pfrag"] += 1

            elif "mult" in biomolecule_type and "gag" in partner_type:
                associations_by_biomolecule_type["multimer"]["gag"] += 1

            elif "gag" in biomolecule_type and "mult" in partner_type:
                associations_by_biomolecule_type["multimer"]["gag"] += 1

            elif "pfrag" in biomolecule_type and "pfrag" in partner_type:
                associations_by_biomolecule_type["pfrag"]["pfrag"] += 1

            elif "pfrag" in biomolecule_type and "gag" in partner_type:
                associations_by_biomolecule_type["pfrag"]["gag"] += 1

            elif "gag" in biomolecule_type and "pfrag" in partner_type:
                associations_by_biomolecule_type["pfrag"]["gag"] += 1

            elif "gag" in biomolecule_type and "gag" in partner_type:
                associations_by_biomolecule_type["gag"]["gag"] += 1

            #else:
                #print("biomolecule type " + biomolecule_type)
                #print("partner type " + partner_type)
                #print()

    # Since all our interactions are symetric, number needs to be divided by 2
    for k in associations_by_biomolecule_type.keys():
        l = associations_by_biomolecule_type[k]
        for m in l.keys():
            associations_by_biomolecule_type[k][m] = associations_by_biomolecule_type[k][m] / 2

    # Count the biomolecules
    biomol_counts = pd.Series(biomol_cache).value_counts().to_dict()

    aggregated_biomol_counts = {}
    for i , j in list(map(lambda x: (convert_name(x), biomol_counts[x]), list(biomol_counts))):
        if i in aggregated_biomol_counts:
            aggregated_biomol_counts[i] = aggregated_biomol_counts[i] + j
        else:
            aggregated_biomol_counts[i] = j

    return json.dumps({
        'associations' : associations_by_biomolecule_type,
        'biomolecules' : biomol_counts
    })


if __name__ == '__main__':

    # Connects to the db
    database_url = "mongodb://localhost:27018/"
    datasets = {}
    try:
        database_client = MongoClient(database_url)
        database = database_client["matrixdb"]
    except Exception:
        print("Problem connecting to db " + database_url)

    # Serve the src with gevent
    http_server = WSGIServer(('127.0.0.1', 8000), app)
    print("Server started at 8000")
    http_server.serve_forever()
