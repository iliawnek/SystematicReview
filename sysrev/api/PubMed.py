from Bio import Entrez

# Abuse contact for if we query too much/cause a problem
Entrez.email = 'systematicreview@nallar.me'


def query(query):
    return Entrez.read(Entrez.esearch(db="pubmed",
                                      sort='relevance',
                                      retmax=1000,
                                      term=query,
                                      field="title",
                                      rettype='xml'))
