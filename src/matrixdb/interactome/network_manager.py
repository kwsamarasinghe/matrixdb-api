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
            transformed_interactor['proteomicsExpression'] = self.protein_data_manager.get_proteomics_expressions(interactor['id'])

        return transformed_interactor

    def transform_interactors(self, interactors):
        # Extract the expresssions for interactors
        gene_expressions = self.protein_data_manager.get_gene_expression_for_proteins(list(i['id'] for i in interactors))
        proteomics_expressions = self.protein_data_manager.get_proteomics_expressions_for_proteins(list(i['id'] for i in interactors))

        transformed_interactors = list()
        for interactor in interactors:
            interactor_id = interactor['id']
            transformed_interactor = {
                'id': interactor_id,
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
                if interactor_id in gene_expressions:
                    transformed_interactor['geneExpression'] = gene_expressions[interactor_id]
                if interactor_id in proteomics_expressions:
                    transformed_interactor['proteomicsExpression'] = proteomics_expressions[interactor_id]

            transformed_interactors.append(transformed_interactor)

        return transformed_interactors

    def transform_interaction(self, interaction):
        transformed_interaction = {
            'id': interaction['id'],
            'participants': interaction['participants'],
            'type': ''
        }

        if "experiments" in interaction:
            transformed_interaction['experiments'] = interaction['experiments']

        if 'prediction' in interaction:
            transformed_interaction['prediction'] = True

        if 'score' in interaction:
            transformed_interaction['score'] = interaction['score']

        if 'prediction' in interaction and 'experiments' in interaction:
            transformed_interaction['type'] = 'Experimental and Predicted'

        if 'prediction' in interaction and 'experiments' not in interaction:
            transformed_interaction['type'] = 'Only Predicted'

        if 'prediction' not in interaction and 'experiments' in interaction:
            transformed_interaction['type'] = 'Only Experimental'



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

            unique_interactors = dict()
            for i in network_interactors:
                if i['id'] not in unique_interactors:
                    unique_interactors[i['id']] = i
                else:
                    print()
            network_interactors = unique_interactors.values()

            # Populate expressions for interactors
            ev = list()
            for i in network_interactions:
                interaction = network_interactions[i]
                if 'experiments' in interaction:
                    if 'direct' in interaction['experiments']:
                        if 'binary' in interaction['experiments']['direct']:
                            ev.extend(interaction['experiments']['direct']['binary'])
                        if 'spoke_expanded_from' in interaction['experiments']['direct']:
                            ev.extend(interaction['experiments']['direct']['spoke_expanded_from'])
            print()
            return {
                'interactions': list(map(lambda interactions_id: self.transform_interaction(network_interactions[interactions_id]),
                                   network_interactions)),
                'interactors': self.transform_interactors(network_interactors)
            }
        else:
            return {
                'interactions': [],
                'interactors': []
            }