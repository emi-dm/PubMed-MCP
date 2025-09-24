from Bio import Entrez
from fastmcp import FastMCP

# Configure Entrez email (required by NCBI)
Entrez.email = 'user@example.com'  # Replace with your actual email

mcp = FastMCP(name="PubMed-MCP")

@mcp.tool
def search_pubmed(query, max_results=10):
    #TODO: Add error handling and logging
    """
    Searches PubMed for a given query and returns a list of articles as JSON objects.
    """
    try:
        # Search PubMed for article IDs using Entrez.esearch
        print(f"Searching for: {query}")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=str(max_results))
        search_record = Entrez.read(handle)
        handle.close()
        
        # Get the list of PMIDs
        pmid_list = search_record["IdList"]
        total_count = search_record["Count"]
        
        print(f"Se encontraron {total_count} artículos. Los primeros {len(pmid_list)} PMIDs son: {pmid_list}")
        
        if not pmid_list:
            print("No articles found for your query.")
            return []
        
        # Fetch detailed information for each PMID using Entrez.efetch
        handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="xml", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        papers_list = []
        
        # Process each article
        for record in records['PubmedArticle']:
            # Get PMID
            pmid = str(record['MedlineCitation']['PMID'])
            
            # Get article details
            article = record['MedlineCitation']['Article']
            
            # Get title
            title = str(article.get('ArticleTitle', 'No title found'))
            
            # Get authors
            authors = []
            if 'AuthorList' in article:
                for author in article['AuthorList']:
                    if 'ForeName' in author and 'LastName' in author:
                        authors.append(f"{author['ForeName']} {author['LastName']}")
                    elif 'LastName' in author:
                        authors.append(str(author['LastName']))
            
            # Get abstract
            abstract = "No abstract found"
            if 'Abstract' in article and 'AbstractText' in article['Abstract']:
                abstract_parts = article['Abstract']['AbstractText']
                if isinstance(abstract_parts, list):
                    # Handle multiple abstract sections
                    abstract_texts = []
                    for part in abstract_parts:
                        if hasattr(part, 'get') and part.get('@Label'):
                            # Structured abstract with labels
                            abstract_texts.append(f"{part['@Label']}: {part}")
                        else:
                            abstract_texts.append(str(part))
                    abstract = " ".join(abstract_texts)
                else:
                    abstract = str(abstract_parts)
            
            # Get publication date
            pub_date_info = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
            year = str(pub_date_info.get('Year', 'Unknown'))
            month = str(pub_date_info.get('Month', 'Unknown'))
            
            # Get journal name
            journal = str(article.get('Journal', {}).get('Title', 'No journal found'))
            
            # Create paper dictionary
            paper = {
                "pmid": pmid,
                "title": title,
                "authors": authors,
                "abstract": abstract,
                "journal": journal,
                "publication_year": year,
                "publication_month": month,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            }
            
            papers_list.append(paper)
        
        print(f"Successfully processed {len(papers_list)} articles.")
        return papers_list
        
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()
        return []

if __name__ == '__main__':
    #parser = argparse.ArgumentParser(description="Search PubMed for articles.")
    #parser.add_argument("query", type=str, help="The search query for PubMed.")
    #parser.add_argument("--max_results", type=int, default=10, help="Maximum number of results to return.")
    
    #args = parser.parse_args()
    
    #search_pubmed(args.query, args.max_results)
    mcp.run(transport="stdio")