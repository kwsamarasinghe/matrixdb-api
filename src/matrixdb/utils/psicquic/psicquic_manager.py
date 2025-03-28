from urllib.request import urlopen, Request
import xml.etree.ElementTree as ET
import json

# ------------------ FUNCTIONS ------------------

class PsicquicService:
    def __init__(self, name, restUrl):
        self.name = name
        self.restUrl = restUrl

def readURL(url, timeout=15):
    try:
        request = Request(url)
        fileHandle = urlopen(request, timeout=timeout)
        content = fileHandle.read()
        fileHandle.close()
    except IOError as e:
        print(f'Cannot open URL {url} due to: {e}')
        content = ''
    except Exception as e:
        print(f'Error with URL {url}: {e}')
        content = ''
    
    return content

def readActiveServicesFromRegistry():
    registryActiveUrl = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=xml'

    content = readURL(registryActiveUrl)

    # Create the XML reader
    root = ET.fromstring(content)
    xmlns = '{http://hupo.psi.org/psicquic/registry}'

    services = []

    for service in root.findall(xmlns + 'service'):
        name = service.find(xmlns + 'name')
        restUrl = service.find(xmlns + 'restUrl')

        service = PsicquicService(name.text, restUrl.text)
        services.append(service)

    return services

def getXrefByDatabase(line, database):
   fields = line.split('|')

   for field in fields:
       parts = field.split(':')

       db = parts[0]
       value = parts[1].split('(')[0]

       if database == db:
           return value

   else:
    # if no db found, return the first field
        return fields[0]

def queryPsicquic(psicquicRestUrl, query, timeout=15):
    psicquicUrl = psicquicRestUrl + 'query/' + query + '?format=count'
    
    print(f'\t\tURL: {psicquicUrl}')

    # Read the result from the URL
    psicquicResultLines = readURL(psicquicUrl, timeout).splitlines()

    # Print the raw result to check
    print(f"Response content from {psicquicUrl}: {psicquicResultLines}")

    # Check if the response is in expected format and extract the count
    if psicquicResultLines and len(psicquicResultLines) > 0:
        try:
            total_results = int(psicquicResultLines[0])
            print(f"Total results for query {query}: {total_results}")
            return total_results
        except ValueError:
            print("Could not parse total results from the response.")
            return 0
    else:
        print(f"No results found for query {query}.")
        return 0

# -----------------------------------------------------

def getTotalResults(query):
    services = readActiveServicesFromRegistry()

    total_results = 0
    successful_queries = 0

    for service in services:
        print(f'Service: {service.name} ================================================================== ')

        service_total = queryPsicquic(service.restUrl, query)

        # Add the total results from each service if successful
        if service_total > 0:
            total_results += service_total
            successful_queries += 1

        print('\n')

    # Return results as a JSON response
    result = {
        'interactions': total_results,
        'successful_services': successful_queries,
        'xref': f'http://www.ebi.ac.uk/Tools/webservices/psicquic/view/main.xhtml?query={query}'
    }
    
    return json.dumps(result)
