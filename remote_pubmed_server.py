from Bio import Entrez
from fastmcp import FastMCP
from typing import List
from urllib.error import HTTPError

# Configure Entrez email (required by NCBI)
Entrez.email = 'user@example.com'  # Replace with your actual email

mcp = FastMCP(name="PubMed-MCP")

async def _build_field_query(user_query: str, in_title: bool, in_abstract: bool, in_keywords: bool) -> str:
    """Build a PubMed (Entrez) query applying field restrictions.

    Fields mapping used:
    - Title: [ti]
    - Abstract: [ab]
    - Title/Abstract convenience: [tiab]
    - Keywords (Other Term): [ot] (Author provided keywords)
    - MeSH Headings: [mh]

    Strategy:
    - If only one of title or abstract is selected, use that specific field tag.
    - If both selected, use tiab (lets PubMed optimize) plus ot/mh if requested.
    - Keywords option expands with OR clauses for ot and mh.
    - Parentheses ensure proper boolean grouping.
    """
    core_clauses: List[str] = []

    if in_title and in_abstract:
        # tiab covers both Title and Abstract text
        core_clauses.append(f"({user_query})[tiab]")
    elif in_title:
        core_clauses.append(f"({user_query})[ti]")
    elif in_abstract:
        core_clauses.append(f"({user_query})[ab]")
    else:
        # No field restriction for title/abstract selected, let user_query as-is (PubMed default: all fields)
        core_clauses.append(f"({user_query})")

    if in_keywords:
        # Include author keywords (ot) and MeSH terms (mh) as expansion
        keywords_clause = f"({user_query})[ot] OR ({user_query})[mh]"
        # Combine with previous core clauses using OR to broaden search
        core_group = " OR ".join(core_clauses)
        combined = f"({core_group}) OR ({keywords_clause})"
        return combined

    return " OR ".join(core_clauses)


@mcp.tool
async def search_pubmed(query: str,
                  max_results: int = 10,
                  title: bool = True,
                  abstract: bool = True,
                  keywords: bool = True):
    """Search PubMed and return a list of article JSON objects.

    Parameters:
        query: Free-text user query; boolean operators (AND/OR/NOT) supported by PubMed.
        max_results: Maximum number of records to retrieve (retmax).
        title: If True, include Title field in search restriction (ti / tiab).
        abstract: If True, include Abstract field in search restriction (ab / tiab).
        keywords: If True, expand search to Author Keywords (ot) and MeSH Headings (mh).

    Field logic:
        - title and abstract both True => core search uses [tiab]
        - only title True => uses [ti]
        - only abstract True => uses [ab]
        - neither title nor abstract True => no restriction (all fields)
        - keywords True => additionally OR with [ot] and [mh] versions of the query

    Returns:
        List[dict]: Each dict contains pmid, title, authors, abstract, journal, publication_year,
                    publication_month, url.
    """
    try:
        if not isinstance(query, str) or not query.strip():
            print("Empty query provided; returning empty result list.")
            return []
        if max_results <= 0:
            max_results = 10

        # Build refined query with field tags
        refined_query = _build_field_query(query.strip(), title, abstract, keywords)

        # Search PubMed for article IDs using Entrez.esearch
        print(f"Searching for: {refined_query}")
        handle = Entrez.esearch(db="pubmed", term=refined_query, retmax=str(max_results))
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
        
    except HTTPError as http_err:
        print(f"HTTP error during Entrez request: {http_err}")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()
        return []

@mcp.prompt(name="precise_pubmed_query", description="Craft a precise PubMed Boolean query from a natural language information need.")
async def precise_pubmed_query(information_need: str) -> str:
    """Returns a very brief prompt (in English) so a model generates a SINGLE concise query
    to search the topic in PubMed. Advanced options are ignored; signature kept for compatibility.

    Main parameter:
        information_need: Topic or information need.
    """
    tema = information_need.strip() or "(specify a topic)"
    prompt = (
        "Generate one brief, precise PubMed query for the given topic. "
        "Then run the search (tool: search_pubmed) and return ONLY a list of articles in the following format, "
        "one per line, with no extra text and no JSON:\n\n"
        "1. Title: <article title> Author: <first author> Journal: <journal> Year: <year> URL: <url>\n"
        f"Topic: {tema}\n\n"
    )
    return prompt

if __name__ == '__main__':
    #parser = argparse.ArgumentParser(description="Search PubMed for articles.")
    #parser.add_argument("query", type=str, help="The search query for PubMed.")
    #parser.add_argument("--max_results", type=int, default=10, help="Maximum number of results to return.")
    
    #args = parser.parse_args()
    
    #search_pubmed(args.query, args.max_results)
    mcp.run(transport="streamable-http", host="0.0.0.0", port=8000)