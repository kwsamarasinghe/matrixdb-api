from pymongo import MongoClient


class DatabaseManager:

    def __init__(self, app_config):
        self.app_config = app_config
        self.database_url = app_config['database_url']

    def get_primary_connection(self):
        """
        Gets a database connection to the primary database
        :return: database connection
        """
        if self.app_config["primary_database_name"] is None:
            raise Exception("Configration has no primary_database_name")

        return MongoClient(self.database_url)[self.app_config["primary_database_name"]]

    def get_secondary_connection(self):
        """
        Gets a database connection to the secondary database
        :return: database connection
        """
        if self.app_config["secondary_database_name"] is None:
            raise Exception("Configration has no secondary_database_name")
        return MongoClient(self.database_url)[self.app_config["secondary_database_name"]]
