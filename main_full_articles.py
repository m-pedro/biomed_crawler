from Bio import Entrez
import pandas as pd
import time
import requests
import os

# --------------------- Configuration --------------------- #

# Replace with your actual email
Entrez.email = "mediratta@gmail.com"

# (Optional) Replace with your actual NCBI API key
# Entrez.api_key = "YOUR_NCBI_API_KEY"

# Define the search term
cancer_type = "uveal melanoma"

# Define the search query
search_query = f'"{cancer_type}"[Title/Abstract]'

# Define batch size
batch_size = 100  # Number of records to fetch per request

# Output filenames
metadata_output = "uveal_melanoma_pubmed_articles.csv"
pmc_output = "uveal_melanoma_pmc_fulltext.csv"

# Directory to save PMC full-texts
fulltext_dir = "PMC_Fulltexts"
os.makedirs(fulltext_dir, exist_ok=True)

# ---------------------------------------------------------- #

def search_pubmed(query):
    """Search PubMed and return the total count, WebEnv, and QueryKey."""
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=0,  # Get only metadata
        usehistory="y"
    )
    record = Entrez.read(handle)
    handle.close()
    count = int(record["Count"])
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    return count, webenv, query_key

def fetch_pubmed_articles(webenv, query_key, count, batch_size):
    """Fetch articles from PubMed in batches."""
    num_batches = (count // batch_size) + 1
    pmids = []
    titles = []
    abstracts = []
    authors_list = []
    journals = []
    publication_dates = []
    pmcids = []
    fulltext_links = []
    
    for batch in range(num_batches):
        start = batch * batch_size
        end = min(start + batch_size, count)
        
        print(f"Fetching records {start + 1} to {end}...")
        
        fetch_handle = Entrez.efetch(
            db="pubmed",
            rettype="xml",
            retmode="xml",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key
        )
        
        records = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        for article in records["PubmedArticle"]:
            # PMID
            pmid = article["MedlineCitation"]["PMID"]
            pmids.append(pmid)
            
            # Title
            title = article["MedlineCitation"]["Article"].get("ArticleTitle", "")
            titles.append(title)
            
            # Abstract
            abstract_sections = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [])
            if isinstance(abstract_sections, list):
                abstract = " ".join(abstract_sections)
            elif isinstance(abstract_sections, str):
                abstract = abstract_sections
            else:
                abstract = ""
            abstracts.append(abstract)
            
            # Authors
            authors = article["MedlineCitation"]["Article"].get("AuthorList", [])
            author_names = []
            for author in authors:
                if "LastName" in author and "Initials" in author:
                    author_names.append(f"{author['LastName']} {author['Initials']}")
                elif "CollectiveName" in author:
                    author_names.append(author["CollectiveName"])
            authors_list.append(", ".join(author_names))
            
            # Journal
            journal = article["MedlineCitation"]["Article"]["Journal"].get("Title", "")
            journals.append(journal)
            
            # Publication Date
            pub_date_info = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].get("PubDate", {})
            year = pub_date_info.get("Year", "")
            month = pub_date_info.get("Month", "")
            day = pub_date_info.get("Day", "")
            publication_date = f"{year} {month} {day}".strip()
            publication_dates.append(publication_date)
            
            # PMC ID (if available)
            pmcid = ""
            article_ids = article["MedlineCitation"].get("ArticleIds", [])
            for id_info in article_ids:
                if id_info.attributes["IdType"] == "pmc":
                    pmcid = id_info
                    break
            pmcid = pmcid if isinstance(pmcid, str) else pmcid.get("_") if pmcid else ""
            pmcids.append(pmcid)
            
            # Full-text link (PMC or external)
            if pmcid:
                fulltext_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/"
            else:
                # Attempt to get URL from PubMed data
                # This may not always be available
                # Alternatively, use Europe PMC for external links
                article_links = article["PubmedData"].get("ArticleIdList", [])
                ext_url = ""
                for link in article_links:
                    if link.attributes["IdType"] == "doi":
                        ext_url = f"https://doi.org/{link}"
                        break
                fulltext_url = ext_url
            fulltext_links.append(fulltext_url)
        
        # Respect NCBI's rate limits
        time.sleep(0.3)  # Adjust as needed based on your API key usage
    
    return (pmids, titles, abstracts, authors_list, journals, publication_dates, pmcids, fulltext_links)

def fetch_pmc_fulltext(pmcid):
    """Fetch full-text PDF from PMC given a PMC ID."""
    if not pmcid.startswith("PMC"):
        pmcid = f"PMC{pmcid}"
    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/"
    response = requests.get(pdf_url)
    if response.status_code == 200:
        return response.content  # Binary content of PDF
    else:
        return None

def main():
    print(f"Searching PubMed for '{cancer_type}'...")
    count, webenv, query_key = search_pubmed(search_query)
    print(f"Total records found: {count}")
    
    if count == 0:
        print("No records found. Exiting.")
        return
    
    print("Starting to fetch records...")
    (pmids, titles, abstracts, authors_list, journals, publication_dates, pmcids, fulltext_links) = fetch_pubmed_articles(
        webenv, query_key, count, batch_size
    )
    
    # Create a DataFrame for metadata
    metadata = {
        "PMID": pmids,
        "Title": titles,
        "Abstract": abstracts,
        "Authors": authors_list,
        "Journal": journals,
        "Publication_Date": publication_dates,
        "PMCID": pmcids,
        "Fulltext_Link": fulltext_links
    }
    
    df_metadata = pd.DataFrame(metadata)
    
    # Save metadata to CSV
    df_metadata.to_csv(metadata_output, index=False)
    print(f"Metadata saved to {metadata_output}")
    
    # Fetch and save full-text PDFs from PMC
    # Note: Only for articles with PMCIDs
    df_pmc = df_metadata[df_metadata["PMCID"] != ""].copy()
    print(f"Total PMC articles with full text available: {len(df_pmc)}")
    
    fulltext_records = []
    
    for idx, row in df_pmc.iterrows():
        pmcid = row["PMCID"]
        title = row["Title"]
        print(f"Fetching full text for {pmcid} - {title}")
        pdf_content = fetch_pmc_fulltext(pmcid)
        if pdf_content:
            # Save PDF to file
            pdf_filename = f"{pmcid}.pdf"
            pdf_path = os.path.join(fulltext_dir, pdf_filename)
            with open(pdf_path, "wb") as f:
                f.write(pdf_content)
            fulltext_records.append({
                "PMCID": pmcid,
                "Title": title,
                "PDF_Path": pdf_path
            })
            print(f"Saved PDF to {pdf_path}")
        else:
            print(f"Failed to fetch PDF for {pmcid}")
        # To avoid overloading PMC servers
        time.sleep(0.3)
    
    # Create DataFrame for PMC full-texts
    df_pmc_fulltext = pd.DataFrame(fulltext_records)
    df_pmc_fulltext.to_csv(pmc_output, index=False)
    print(f"PMC full-text records saved to {pmc_output}")

if __name__ == "__main__":
    main()
