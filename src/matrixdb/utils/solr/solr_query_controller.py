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
    full_url = f"{core_url}/select"
    response = requests.get(full_url, params=params)
    return response.json()["response"]["docs"]