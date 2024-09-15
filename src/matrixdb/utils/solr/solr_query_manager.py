import urllib
from urllib.parse import urlencode

from numpy import percentile
import requests as requests


class SolrQueryManager:

    def __init__(self, app_config):
        self.solr_url = app_config["solr_url"]
        self.biomolecules_core_url = f'{self.solr_url}/solr/biomolecules'
        self.publications_core_url = f'{self.solr_url}/solr/publications'

        self.advanced_query_field_mappings = {
            'id': ['biomolecule_id_exact'],
            'uniprot': ['biomolecule_id_exact'],
            'name': ['recommended_name_exact', 'name', 'common_name'],
            'gene': ['exact_gene'],
            'go': [
                'go_names'
            ],
            'keywords': [
                'keyword_names'
            ],
            'matrixdb_ecm': ['matrixdb_ecm'],
            'matrisome_subcategory': ['matrisome_subcategory'],
            'matrisome_category': ['matrisome_category'],
            'species': ['species'],
            'reactome': ['reactome']
        }

    def search_biomolecules(self, search_query, mode="0"):
        if mode == "0":
            return self.do_simple_search(search_query)

        if mode == "1":
            return self.do_advanced_search(search_query)

    def do_simple_search(self, search_query):
        biomolecule_query_params = {
            'q': f'{search_query}',
            "q.op": "OR",
            'defType': 'dismax',
            'qf': 'biomolecule_id_exact^20 biomolecule_id^8 name^8 common_name^5 recommended_name_exact^15 recommended_name^8 other_name^2 species^2 description^2 gene^2 chebi^10 complex_portal^10 go_names go_ids keyword',
            'pf': 'name^4 recommended_name^4 go_names^4',
            'bq': '-dataset:TrEMBL^10',
            'fl': "*,score",
            'rows': 100
        }

        biomolecule_solr_docs = self.query_solr(self.biomolecules_core_url, biomolecule_query_params)

        # Check if the query exactly matched to id, name, common_name
        #exact_matches = self.get_exact_matches(search_query, biomolecule_solr_docs)
        #if len(exact_matches) > 0:
        #    biomolecule_solr_docs = exact_matches
        #    return biomolecule_solr_docs

        # Get the 75% percentile
        if len(biomolecule_solr_docs) >= 5:
            scores = list(b["score"] for b in biomolecule_solr_docs)
            percentile_cutoff = min(scores) + (max(scores) - min(scores)) * 0.60
            outlier_cutoff = self.get_outlier_cutoff(scores)

            for document in biomolecule_solr_docs:

                if 'dataset' in document and document['dataset'] == 'TrEMBL':
                    continue

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

    def parse_advanced_query(self, query_part):
        parsed_query_part = {
            "field": query_part[:query_part.index(':')].strip(),
        }

        query_values = query_part[query_part.index(":") + 1:].strip()
        if "," in query_values:
            query_values = query_values.split(",")
        else:
            query_values = [query_values]

        parsed_query_part["values"] = query_values
        return parsed_query_part

    def do_advanced_search(self, search_query):
        if ":" not in search_query:
            return list()

        # Check if it has conjunction AND
        parsed_query = []
        if "and" in search_query.lower():
            query_parts = search_query.split("and")
            if len(query_parts) > 2:
                return []

            for query_part in query_parts:
                parsed_query.append(self.parse_advanced_query(query_part))
        else:
            parsed_query.append(self.parse_advanced_query(search_query))

        #for query_field in query_fields:
        #    if query_field not in self.advanced_query_field_mappings:
        #        return []

        biomolecule_query_param_string = "q=*:*&rows=1000"

        fq = []
        for parsed_query_part in parsed_query:
            query_field = parsed_query_part["field"]
            query_parameters = parsed_query_part["values"]

            if query_field not in self.advanced_query_field_mappings:
                return []

            query_field_mappings = self.advanced_query_field_mappings[query_field]
            sub_fqs = []
            for query_field_mapping in query_field_mappings:
                field_level_fq = [f'{query_field_mapping}:"{query_parameter}"' for query_parameter in query_parameters]
                if len(field_level_fq) > 1:
                    sub_fq = ' OR '.join(field_level_fq)
                else:
                    sub_fq = field_level_fq[0]
                sub_fqs.append(sub_fq)
            fq.append(f'({" OR ".join(sub_fqs)})')

        biomolecule_query_param_string += f'&fq=({" AND ".join(fq)})'
        biomolecule_solr_docs = self.query_solr_with_url(f'{self.biomolecules_core_url}/select?{biomolecule_query_param_string}')
        return biomolecule_solr_docs

    def search_publications(self, search_query):

        publication_query_params = {
            'q': f'{search_query}',
            'defType': 'dismax',
            'qf': 'publication_id title authors journal',
            'rows': 50
        }
        publication_solr_docs = self.query_solr(self.publications_core_url, publication_query_params)
        return sorted(publication_solr_docs, key=lambda doc: doc['interaction_count'])

    def query_solr(self, core_url, params):
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

    def query_solr_with_url(self, full_url):
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
        upper_bound = q_95 + 1.5 * iqr

        outliers = [x for x in scores if x > upper_bound]
        if len(outliers) > 0:
            return sorted(outliers)[0]
        else:
            return -1
