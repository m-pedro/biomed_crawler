import os
import sys
import time
from Bio import Entrez
from Bio import Medline
import requests

# Set your email for Entrez
Entrez.email = "your_email@example.com"  # Replace with your email

def search_pubmed(query, max_results=10):
    """
    Search PubMed for articles based on the given query.
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        results = Entrez.read(handle)
        handle.close()
        return results["IdList"]
    except Exception as e:
        print(f"Error searching PubMed: {e}")
        return []

def fetch_pubmed_details(id_list):
    """
    Fetch details of PubMed articles given their IDs.
    """
    try:
        ids = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        return list(records)
    except Exception as e:
        print(f"Error fetching PubMed details: {e}")
        return []

def check_pmc_availability(pmid):
    """
    Check if an article is available in PubMed Central.
    """
    try:
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", linkname="pubmed_pmc", id=pmid)
        result = Entrez.read(handle)
        handle.close()
        
        if result[0]["LinkSetDb"]:
            return result[0]["LinkSetDb"][0]["Link"][0]["Id"]
        else:
            return None
    except Exception as e:
        print(f"Error checking PMC availability: {e}")
        return None

def download_pmc_article(pmcid, output_dir):
    """
    Download an article from PubMed Central.
    """
    base_url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC"
    pdf_url = f"{base_url}{pmcid}/pdf/main.pdf"

    
    try:
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
        }
        response = requests.get(pdf_url, headers=headers)
        response.raise_for_status()
        
        filename = f"PMC{pmcid}.pdf"
        filepath = os.path.join(output_dir, filename)
        
        with open(filepath, 'wb') as f:
            f.write(response.content)
        
        print(f"Article downloaded: {filepath}")
        return filepath
    except requests.RequestException as e:
        print(f"Error downloading article: {e}")
        return None

def main():
    query = input("Enter your PubMed search query: ")
    max_results = int(input("Enter the maximum number of results (default 10): ") or 10)
    
    print(f"Searching PubMed for: {query}")
    id_list = search_pubmed(query, max_results)
    
    if not id_list:
        print("No results found.")
        return
    
    print(f"Found {len(id_list)} results.")
    records = fetch_pubmed_details(id_list)
    
    available_pmc = []
    for record in records:
        pmid = record.get("PMID", "N/A")
        title = record.get("TI", "No title")
        print(f"PMID: {pmid}")
        print(f"Title: {title}")
        
        pmcid = check_pmc_availability(pmid)
        if pmcid:
            print("Available in PMC")
            available_pmc.append((pmid, pmcid, title))
        else:
            print("Not available in PMC")
        print("-" * 50)
    
    if not available_pmc:
        print("No articles available in PubMed Central.")
        return
    
    print("\nArticles available in PubMed Central:")
    for i, (pmid, pmcid, title) in enumerate(available_pmc, 1):
        print(f"{i}. PMID: {pmid}, PMCID: {pmcid}")
        print(f"   Title: {title}")
    
    for pmid, pmcid, title in available_pmc:
        print(f"Downloading article: {title}")
        output_dir = "downloaded_articles"
        os.makedirs(output_dir, exist_ok=True)
        
        downloaded_file = download_pmc_article(pmcid, output_dir)
        if downloaded_file:
            print(f"Article successfully downloaded to {downloaded_file}")
        else:
            print(f"Failed to download the article: {title}")
    
    print("All articles downloaded.")

if __name__ == "__main__":
    main()
