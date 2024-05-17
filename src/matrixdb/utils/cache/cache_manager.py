class CacheManager:
    """
    Cache manager which maintains a in-memory cache of meta-data
    """

    def __init__(self, app_config, database_manager):
        self.meta_data_cache = {
            "psimi": dict(),
            "go": dict(),
            "interpro": dict(),
            "uniprotKeywords": dict(),
            "uberon": dict(),
            "bto": dict()
        }
        self.app_config = app_config
        self.database_manager = database_manager

        core_database_connection = self.database_manager.get_primary_connection()
        # psimi
        for psimi in core_database_connection["psimi"].find():
            self.meta_data_cache["psimi"][psimi["id"]] = psimi

        # go
        for go in core_database_connection["go"].find():
            self.meta_data_cache["go"][go["id"]] = go

        # Interpro
        for interpro in core_database_connection["interpro"].find():
            self.meta_data_cache["interpro"][interpro["id"]] = interpro

        # Uniprot keywords
        for uniprot_keyword in core_database_connection["uniprotKeywords"].find():
            self.meta_data_cache["uniprotKeywords"][uniprot_keyword["id"]] = uniprot_keyword

        for uberon in core_database_connection["uberon"].find():
            self.meta_data_cache["uberon"][uberon["id"]] = uberon

        for bto in core_database_connection["brenda"].find():
            self.meta_data_cache["bto"][bto["id"]] = bto

        print("Cache manager initialized")

    def get_meta_data(self):
        return self.meta_data_cache
