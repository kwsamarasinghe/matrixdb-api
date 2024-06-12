import urllib
from urllib.parse import urlencode

from numpy import percentile
import requests as requests


class SolrQueryManager:

    def __init__(self, app_config):
        self.solr_url = app_config["solr_url"]
        self.biomolecules_core_url = f'{self.solr_url}/solr/biomolecules'
        self.publications_core_url = f'{self.solr_url}/solr/publications'

    def query_biomolecules(self, search_query):
        biomolecule_query_params = {
            'q': f'{search_query}',
            "q.op": "OR",
            'defType': 'dismax',
            'qf': 'biomolecule_id^5 name^5 common_name^5 recommended_name^9 other_name^2 species^2 description^2 gene^6 chebi complex_portal go_names go_ids keyword_ids keyword_names',
            'pf': 'name^2 recommended_name^2 go_names^1.2 common_name^2',
            'fl': "*,score",
            'rows': 1000
        }

        biomolecule_solr_docs = self.query_solr(self.biomolecules_core_url, biomolecule_query_params)

        # Check if the query exactly matched to id, name, common_name
        #exact_matches = self.get_exact_matches(search_query, biomolecule_solr_docs)
        #if len(exact_matches) > 0:
        #    biomolecule_solr_docs = exact_matches
        # Get the 75% percentile
        if len(biomolecule_solr_docs) >= 5:
            scores = list(b["score"] for b in biomolecule_solr_docs)
            percentile_cutoff = min(scores) + (max(scores) - min(scores)) * 0.60
            outlier_cutoff = self.get_outlier_cutoff(scores)

            for document in biomolecule_solr_docs:

                if outlier_cutoff != -1:
                    if document["score"] >= outlier_cutoff:
                        document["outlier"] = True

                if document["score"] >= percentile_cutoff:
                    document["most_relevant"] = True

        '''
        # Sort by interactions
        biomolecule_solr_docs_with_interactions = list(
            filter(lambda doc: doc['interaction_count'] > 0, biomolecule_solr_docs))
        
                if len(biomolecule_solr_docs_with_interactions) >= 100:
            biomolecule_solr_docs = sorted(biomolecule_solr_docs_with_interactions,
                                           key=lambda doc: doc['interaction_count'])[0:99]
        else:
            biomolecule_solr_docs = sorted(biomolecule_solr_docs, key=lambda doc: doc['interaction_count'],
                                           reverse=True)[0:99]
        '''
        return biomolecule_solr_docs

    def query_publications(self, search_query):

        publication_query_params = {
            'q': f'{search_query}',
            'defType': 'dismax',
            'qf': 'publication_id title authors journal',
            'rows': 50
        }
        publication_solr_docs = self.query_solr(self.publications_core_url, publication_query_params)
        return sorted(publication_solr_docs, key=lambda doc: doc['interaction_count'])

    @staticmethod
    def query_solr(core_url, params):
        """
        Perform a Solr query on the specified core with the given parameters.

        Parameters:
        - core_url: The URL of the Solr core (e.g., http://localhost:8983/solr/core1)
        - params: A dictionary of query parameters

        Returns:
        - The Solr query response in JSON format
        """
        full_url = f"{core_url}/select?{urllib.parse.urlencode(params)}"
        response = requests.get(full_url)
        return response.json()["response"]["docs"]

    def get_exact_matches(self, search_text, result_documents):
        return list(filter(lambda doc: doc['biomolecule_id'].lower() == search_text.lower()
                    or 'gene' in doc and search_text.lower() in doc['gene'].lower()
                    or 'name' in doc and doc['name'].lower() == search_text.lower()
                    or 'recommended_name' in doc and doc['recommended_name'].lower() == search_text.lower()
                    or 'chebi' in doc and doc['chebi'].lower() == search_text.lower()
                    or 'complex_portal' in doc and doc['complex_portal'].lower() == search_text.lower()
                    or 'species' in doc and doc['species'].lower() == search_text.lower(), result_documents))

    def get_outlier_cutoff(self, scores):
        q_95 = percentile(scores, 95)
        q_5 = percentile(scores, 5)
        iqr = q_95 - q_5
        lower_bound = q_5 - 1.5 * iqr
        upper_bound = q_95 + 1.5 * iqr

        outliers = [x for x in scores if x > upper_bound]
        if len(outliers) > 0:
            return sorted(outliers)[0]
        else:
            return -1