from src.matrixdb.biomolecule_services.protein_data_manager import ProteinDataManager


class NetworkManager:

    def __init__(self, database_connection, meta_data_cache):
        # initialize database connection
        self.database_connection = database_connection
        self.biomolecule_cache = meta_data_cache
        self.protein_data_manager = ProteinDataManager(meta_data_cache)

    def merge_associations(self, *interactions):
        result = interactions[0].copy()

        for d in interactions[1:]:
            for key, value in d.items():
                if key in result and isinstance(result[key], list) and isinstance(value, list):
                    result[key] += value
                elif key in result and isinstance(result[key], dict) and isinstance(value, dict):
                    result[key] = self.merge_associations(result[key], value)
                else:
                    result[key] = value

        return result

    def transform_interactor(self, interactor):
        transformed_interactor = {
            'id': interactor['id'],
            'type': interactor['type']
        }

        if 'molecular_details' in interactor and 'pdb' in interactor['molecular_details']:
            transformed_interactor['pdb'] = len(interactor['molecular_details']) > 0

        if 'ecm' in interactor:
            transformed_interactor['ecm'] = True

        if 'relations' in interactor and 'gene_name' in interactor['relations']:
            transformed_interactor['geneName'] = interactor['relations']['gene_name']

        # Get expressions and add to interactor
        if interactor['type'] == 'protein':
            transformed_interactor['geneExpression'] = self.protein_data_manager.get_gene_expressions(interactor['id'])
            transformed_interactor['proteomicsExpression'] = self.protein_data_manager.get_proteomics_expression(interactor['id'])

        return transformed_interactor

    def transform_interaction(self, interaction):
        transformed_interaction = {
            'id': interaction['id'],
            'participants': interaction['participants'],
        }

        if "experiments" in interaction:
            transformed_interaction['experiments'] = interaction['experiments']

        if 'prediction' in interaction:
            transformed_interaction['prediction'] = True

        if 'score' in interaction:
            transformed_interaction['score'] = interaction['score']

        return transformed_interaction

    def generate_network(self, biomolecules):
        '''
        Generates the interaction network of one or more biomolecules
        Network consists of all biomolecule partners of the give biomolecule list, found in interactions
        :param biomolecules:
        :return: {
                    biomolecules: list of all biomoleules
                    interactions: list of all interactions
                 }
        '''
        # Get the associations with participants in given biomolecules
        interactions = list(self.database_connection["interactions"].find({
            'participants': {
                '$in': biomolecules
            }
        }))

        # Get the participants from interactions
        network_interactions = dict()
        network_interactors_map = dict()
        if len(interactions) > 0:
            for interaction in interactions:

                if interaction['id'] in interactions:
                    # Multiple associations with same id
                    network_interactions[interaction['id']] = \
                        self.merge_associations(network_interactions[interaction['id']], interaction)
                else:
                    network_interactions[interaction['id']] = interaction

                network_interactors_map[interaction['participants'][0]] = interaction['participants'][0]
                if len(interaction['participants']) > 1:
                    network_interactors_map[interaction['participants'][1]] = interaction['participants'][1]

            # Populate interactors
            network_interactors = list(self.database_connection["biomolecules"].find({
                'id': {
                    '$in': list(network_interactors_map.keys())
                }
            }))

            return {
                'interactions': list(map(lambda interactions_id: self.transform_interaction(network_interactions[interactions_id]),
                                   network_interactions)),
                'interactors': list(map(lambda network_interactor: self.transform_interactor(network_interactor),
                                   network_interactors))
            }
        else:
            return {
                'interactions': [],
                'interactors': []
            }