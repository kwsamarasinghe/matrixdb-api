from flask import Flask, request, Response
from flask_cors import CORS, cross_origin
from gevent.pywsgi import WSGIServer
from pymongo import MongoClient
import json

app = Flask(__name__)
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
    for biomolecule in associations["associations"]:
        association = associations["associations"][biomolecule]
        if len(association["directlysupportedby"]) > 0 or len(association["spokeexpandedfrom"]) > 0:
            direct = direct + 1
        #if len(association["inferredfrom"]) > 0:
        else:
            inferred = inferred + 1

        neighborhood[association["partner"]] = {
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
        for biomolecule in intact_associations["associations"]:
            intact_association = intact_associations["associations"][biomolecule]
            if intact_association["partner"] not in neighborhood:
                neighborhood[intact_association["partner"]] = {}
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
        return {
            "message" : "No assoication found for " + id
        }


@app.route('/api/search', methods=['GET'])
def search_with_text():
    args = request.args
    search_text = args['text']
    # Search biomolecules for the moment
    biomolecules_by_id = database["biomolecules"].find(
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
        biomolecules_by_id = database["biomolecules"].find(
            {

                    "$or":[
                        {
                            "id": { "$regex" : search_text, '$options' : 'i' }
                        },
                        {
                            "Belongs_to": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "GAG_Comments": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "SmallMol_Definition": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "FragmentName": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "GAG_Name": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "Personal_Keyword": {"$regex": search_text, '$options' : 'i'}
                        },
                        {
                            "GlycoCT_String": {"$regex": search_text, '$options': 'i'}
                        },
                        {
                            "CheBI_identifier": search_text
                        }
                    ]

            },
            {
                "_id": 0
            }
        )

    return json.dumps(list(biomolecules_by_id))

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
