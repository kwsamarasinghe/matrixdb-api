
class NetworkManager:

    def __init__(self, database_connection,
                 meta_data_cache,
                 protein_data_manager,
                 interaction_data_manager=None):
        # initialize database connection
        self.database_connection = database_connection
        self.meta_data_cache = meta_data_cache
        self.protein_data_manager = protein_data_manager
        self.interaction_data_manager = interaction_data_manager

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

    def transform_interaction(self, interaction, experiment_map):
        transformed_interaction = {
            'id': interaction['id'],
            'participants': interaction['participants'],
            'type': ''
        }

        if "experiments" in interaction:
            transformed_interaction['experiments'] = {
                "direct" : {
                    "binary": list(),
                    "spoke_expanded_from": list()
                }
            }

            transformed_interaction["detection_method"] = list()
            if "binary" in interaction['experiments']["direct"]:
                for binary in interaction['experiments']['direct']['binary']:
                    experiment = experiment_map[binary]
                    transformed_interaction['experiments']['direct']['binary'].append(experiment["id"])
                    method = self.meta_data_cache["psimi"][experiment["interaction_detection_method"]]["name"]
                    transformed_interaction["detection_method"].append(method)

            if "spoke_expanded_from" in interaction['experiments']['direct']:
                for se in interaction['experiments']['direct']['spoke_expanded_from']:
                    experiment = experiment_map[se]
                    transformed_interaction["experiments"]["direct"]["spoke_expanded_from"].append(experiment["id"])
                    method = self.meta_data_cache["psimi"][experiment["interaction_detection_method"]]["name"]
                    transformed_interaction["detection_method"].append(method)

        if 'prediction' in interaction:
            transformed_interaction['prediction'] = True

        if 'score' in interaction:
            if interaction['score'] != '-':
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
        # Get the 1-neighborhood of each biomolecule and derive the interactions and partners
        interactor_ids = set()
        interaction_ids = set()

        for biomolecule in biomolecules:
            # Add current biomolecule
            interactor_ids.add(biomolecule)

            # Partners
            interactor_ids.update(partner for partner in self.interaction_data_manager.get_neighborhood(biomolecule))

            # Interactions
            for partner in interactor_ids:
                sorted_partners = sorted([biomolecule, partner])
                interaction_id = f"{sorted_partners[0]}__{sorted_partners[1]}"
                interaction_ids.add(interaction_id)

        if len(interaction_ids) > 0:

            # Extract the interactors
            interactors = list(self.database_connection["biomolecules"].find({
                "id": {"$in": list(interactor_ids)}
            }))
            transformed_interactors = self.transform_interactors(interactors)

            # Extract the interactions
            interactions = list(self.database_connection["interactions"].find({
                "id": {"$in": list(interaction_ids)}
            }))

            # Extracts a map of experiments
            experiment_ids = set()
            for interaction in interactions:
                if "experiments" in interaction:
                    if "binary" in interaction["experiments"]["direct"]:
                        for binary_experiments in interaction["experiments"]["direct"]["binary"]:
                            experiment_ids.add(binary_experiments)
                    if "spoke_expanded_from" in interaction["experiments"]["direct"]:
                        for spoke_exapnded in interaction["experiments"]["direct"]["spoke_expanded_from"]:
                            experiment_ids.add(spoke_exapnded)

            experiments = dict()
            for experiment in list(self.database_connection["experiments"].find({
                "id": {"$in": list(experiment_ids)}
            })):
                del experiment["_id"]
                experiments[experiment["id"]] = experiment

            transformed_interactions = list()
            for interaction in interactions:
                transformed_interactions.append(self.transform_interaction(interaction, experiments))

            return {
                'interactions': transformed_interactions,
                'interactors': transformed_interactors
            }
        else:
            return {
                'interactions': [],
                'interactors': []
            }

