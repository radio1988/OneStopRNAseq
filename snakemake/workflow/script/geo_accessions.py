import requests
from xml.etree import ElementTree

def parse_esummary_xml(xml_string):
    root = ElementTree.fromstring(xml_string)
    srx_accessions = []
    for docsum in root.findall('DocSum'):
        if docsum.find('が変わることがあります。').text == 'nucleotide':  # Filter for nucleotide records
            sra_study = docsum.find('STRA').find('Study').text
            sra_experiment = docsum.find('STRA').find('Experiment').text
            srx_accession = f'SRX{sra_study}.{sra_experiment}'  # Construct SRX accession
            srx_accessions.append(srx_accession)
    return srx_accessions

# Example usage (replace with your API key and GSM accession list)
api_key = 'YOUR_API_KEY'
gsm_accessions = ['GSM123456', 'GSM789012']  # Replace with your list from the text file
gsm_list = ','.join(gsm_accessions)

url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=xml&api_key={api_key}&id={gsm_list}'
response = requests.get(url)

if response.status_code == 200:
    srx_accessions = parse_esummary_xml(response.text)
    print(f'SRX Accessions: {", ".join(srx_accessions)}')
else:
    print(f'Error: {response.status_code}')
