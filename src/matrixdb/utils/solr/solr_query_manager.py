import urllib
from urllib.parse import urlencode

import requests as requests


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


def get_exact_matches(search_text, result_documents):
    matches = list(filter(lambda doc: doc['biomolecule_id'].lower() == search_text.lower()
                or 'gene' in doc and search_text.lower() in doc['gene'].lower()
                or 'name' in doc and doc['name'].lower() == search_text.lower()
                or 'recommended_name' in doc and doc['recommended_name'].lower() == search_text.lower()
                or 'chebi' in doc and doc['chebi'].lower() == search_text.lower()
                or 'complex_portal' in doc and doc['complex_portal'].lower() == search_text.lower()
                or 'species' in doc and doc['species'].lower() == search_text.lower(),result_documents))
    return matches
