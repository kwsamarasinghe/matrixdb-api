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
            "bto": dict(),
            "ncbiTaxonomy": dict(),
            "reactome": dict()
        }
        self.app_config = app_config
        self.database_manager = database_manager

        core_database_connection = self.database_manager.get_primary_connection()
        # psimi
        for psimi in core_database_connection["psimi"].find():
            del psimi["_id"]
            self.meta_data_cache["psimi"][psimi["id"]] = psimi

        # go
        for go in core_database_connection["go"].find():
            del go["_id"]
            self.meta_data_cache["go"][go["id"]] = go

        # Interpro
        for interpro in core_database_connection["interpro"].find():
            del interpro["_id"]
            self.meta_data_cache["interpro"][interpro["id"]] = interpro

        # Uniprot keywords
        for uniprot_keyword in core_database_connection["uniprotKeywords"].find():
            del uniprot_keyword["_id"]
            self.meta_data_cache["uniprotKeywords"][uniprot_keyword["id"]] = uniprot_keyword

        # Uberon
        for uberon in core_database_connection["uberon"].find():
            del uberon["_id"]
            self.meta_data_cache["uberon"][uberon["id"]] = uberon

        #BTO
        for bto in core_database_connection["brenda"].find():
            del bto["_id"]
            self.meta_data_cache["bto"][bto["id"]] = bto

        # NCBI
        for ncbi in core_database_connection["ncbiTaxonomy"].find():
            del ncbi["_id"]
            self.meta_data_cache["ncbiTaxonomy"][ncbi["id"]] = ncbi

        # Reactome
        for reactome in core_database_connection["reactome"].find():
            del reactome["_id"]
            self.meta_data_cache["reactome"][reactome["id"]] = reactome

        print("Cache manager initialized")

    def get_meta_data(self):
        return self.meta_data_cache
