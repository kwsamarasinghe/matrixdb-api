
class InteractionDataManager:

    def __init__(self, app_config, database_manager):
        self.app_config = app_config
        self.database_manager = database_manager

        # Generates interaction cache
        self.neighborhood_cache = dict()
        print("Building interaction cache")
        core_database_connection = self.database_manager.get_primary_connection()
        for interaction in core_database_connection["interactions"].find({}):
            if len(interaction["participants"]) == 2:
                participant_0 = interaction["participants"][0]
                participant_1 = interaction["participants"][1]

            # Self interactions
            if len(interaction["participants"]) == 1:
                participant_0 = interaction["participants"][0]
                participant_1 = interaction["participants"][0]

            if participant_0 not in self.neighborhood_cache:
                self.neighborhood_cache[participant_0] = set()
            self.neighborhood_cache[participant_0].add(participant_1)

            if participant_1 not in self.neighborhood_cache:
                self.neighborhood_cache[participant_1] = set()
            self.neighborhood_cache[participant_1].add(participant_0)

        print(f"Neighborhood cache built : entries {len(self.neighborhood_cache.keys())}")

    def get_neighborhood(self, biomolecule):
        return self.neighborhood_cache[biomolecule]