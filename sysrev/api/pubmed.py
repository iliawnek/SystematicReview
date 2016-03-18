import Bio.Entrez

def pubmed_query(query):
    return Bio.Entrez.esearch(db="pubmed",
                            sort='relevance',
                            retmax=1000,
                            term=query)