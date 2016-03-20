from Bio import Entrez

# Abuse contact for if we query too much/cause a problem
Entrez.email = 'systematicreview@nallar.me'

db = "pubmed"


def query(query):
    return Entrez.read(Entrez.esearch(db=db,
                                      sort='relevance',
                                      retmax=1000,
                                      term=query,
                                      field="title",
                                      rettype='json'))


def paper(ids):
    # Convert list of ids to comma-delimited list
    idList = ",".join(map(str, ids))
    return Entrez.read(Entrez.efetch(db=db, id=idList, retmode='xml'))
