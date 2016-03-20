from Bio import Entrez

# Abuse contact for if we query too much/cause a problem
Entrez.email = 'systematicreview@nallar.me'

db = "pubmed"

def get_data_from_query(query):
    """Returns raw data from pubmed query from API"""
    return Entrez.read(Entrez.esearch(db=db,
                                      sort='relevance',
                                      retmax=1000,
                                      term=query,
                                      field="All Fields",
                                      rettype='json'))


def get_ids_from_query(query):
    """Returns list of ids matched by a query from API"""
    return get_data_from_query(query)[u'IdList']


def read_papers_from_ids(ids):
    """Reads data for all ids in list from API"""
    # Convert list of ids to comma-delimited list
    idList = ",".join(map(str, ids))
    return Entrez.read(Entrez.efetch(db=db, id=idList, retmode='xml'))


def _get_authors(article):
    """Gets comma delimited list of authors from article"""
    if not u'AuthorList' in article:
        return ""

    authorList = article[u'AuthorList']
    return ", ".join(filter(lambda it: len(it) > 0, map(_author_to_string, authorList)))


def _author_to_string(author):
    try:
        return author[u'ForeName'] + " " + author[u'LastName']
    except KeyError:
        return ""


def url_from_id(id):
    """Returns URL of article with given ID"""
    return "https://www.ncbi.nlm.nih.gov/pubmed/" + id


def _get_date(medlineCitation):
    """Returns date of the given medline citation"""
    date = medlineCitation[u'DateCreated']
    return date[u'Year'] + '-' + date[u'Month'] + '-' + date[u'Day']
