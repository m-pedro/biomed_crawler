from Bio import Entrez
import pandas as pd
import time

# Replace with your actual email
Entrez.email = "mediratta@gmail.com"

# (Optional) Replace with your actual NCBI API key
# Entrez.api_key = "YOUR_NCBI_API_KEY"

# Define the search term
cancer_type = "uveal melanoma"

# Define the search query
# Here, we're searching for articles that mention "uveal melanoma" in the title or abstract
search_query = f'"{cancer_type}"[Title/Abstract]'

# Perform the search on PubMed
handle = Entrez.esearch(
    db="pubmed",
    term=search_query,
    retmax=100000,  # Maximum number of records to retrieve (adjust as needed)
    usehistory="y"   # Use history to enable fetching results in batches
)

# Parse the search results
record = Entrez.read(handle)
handle.close()

# Extract relevant information
count = int(record["Count"])
print(f"Total records found: {count}")

# Extract WebEnv and QueryKey for fetching records in batches
webenv = record["WebEnv"]
query_key = record["QueryKey"]


# Define batch size
batch_size = 100  # Number of records to fetch per request

# Calculate the number of batches
num_batches = (count // batch_size) + 1

# Initialize lists to store the fetched data
pmids = []
titles = []
abstracts = []
authors_list = []
journals = []
publication_dates = []

for batch in range(num_batches):
    # Calculate the start and end indices for the current batch
    start = batch * batch_size
    end = min(start + batch_size, count)
    
    print(f"Fetching records {start + 1} to {end}...")
    
    # Fetch the records
    # fetch_handle = Entrez.efetch(
    #     db="pubmed",
    #     rettype="medline",
    #     retmode="text",
    #     retstart=start,
    #     retmax=batch_size,
    #     webenv=webenv,
    #     query_key=query_key
    # )
    
    # # Parse the fetched data
    # records = Entrez.read(fetch_handle)
    # fetch_handle.close()
    
    # Note: The above efetch with 'medline' format returns data in a different structure.
    # To simplify, we'll use the 'xml' format instead.
    # Let's adjust the efetch call to use 'xml'.
    
    # Adjusted fetch
    fetch_handle = Entrez.efetch(
        db="pubmed",
        rettype="xml",
        retmode="xml",
        retstart=start,
        retmax=batch_size,
        webenv=webenv,
        query_key=query_key
    )
    
    # Parse the XML data
    records = Entrez.read(fetch_handle)
    fetch_handle.close()
    
    # Iterate through each article
    for article in records["PubmedArticle"]:
        # Extract PMID
        pmid = article["MedlineCitation"]["PMID"]
        pmids.append(pmid)
        
        # Extract Article Title
        try:
            title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        except KeyError:
            title = ""
        titles.append(title)
        
        # Extract Abstract
        try:
            abstract_sections = article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            abstract = " ".join(abstract_sections)
        except KeyError:
            abstract = ""
        abstracts.append(abstract)
        
        # Extract Authors
        try:
            authors = article["MedlineCitation"]["Article"]["AuthorList"]
            author_names = []
            for author in authors:
                # Some authors might not have a LastName or Initials
                if "LastName" in author and "Initials" in author:
                    author_names.append(f"{author['LastName']} {author['Initials']}")
                elif "CollectiveName" in author:
                    author_names.append(author["CollectiveName"])
            authors_list.append(", ".join(author_names))
        except KeyError:
            authors_list.append("")
        
        # Extract Journal Name
        try:
            journal = article["MedlineCitation"]["Article"]["Journal"]["Title"]
        except KeyError:
            journal = ""
        journals.append(journal)
        
        # Extract Publication Date
        try:
            pub_date_info = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
            year = pub_date_info.get("Year", "")
            month = pub_date_info.get("Month", "")
            day = pub_date_info.get("Day", "")
            publication_date = f"{year} {month} {day}".strip()
        except KeyError:
            publication_date = ""
        publication_dates.append(publication_date)
    
    # To comply with NCBI rate limits, especially if not using an API key
    # it's good practice to wait between requests
    time.sleep(0.3)  # Sleep for 300 milliseconds

# Create a DataFrame
data = {
    "PMID": pmids,
    "Title": titles,
    "Abstract": abstracts,
    "Authors": authors_list,
    "Journal": journals,
    "Publication_Date": publication_dates
}

df = pd.DataFrame(data)

# Save to CSV
output_filename = "uveal_melanoma_pubmed_articles.csv"
df.to_csv(output_filename, index=False)

print(f"Data saved to {output_filename}")
