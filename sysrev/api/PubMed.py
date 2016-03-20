from Bio import Entrez

from sysrev.models import Paper

# Abuse contact for if we query too much/cause a problem
Entrez.email = 'systematicreview@nallar.me'

db = "pubmed"


def get_data_from_query(query):
    """Returns raw data from pubmed query from API"""
    return Entrez.read(Entrez.esearch(db=db,
                                      sort='relevance',
                                      retmax=1000,
                                      term=query,
                                      field="title",
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
    authorList = article[u'AuthorList']
    return ", ".join(map(lambda author: author[u'ForeName'] + " " + author[u'LastName'], authorList))


def url_from_id(id):
    """Returns URL of article with given ID"""
    return "https://www.ncbi.nlm.nih.gov/pubmed/" + id


def _get_date(medlineCitation):
    """Returns date of the given medline citation"""
    date = medlineCitation[u'DateCreated']
    return date[u'Year'] + '-' + date[u'Month'] + '-' + date[u'Day']


def create_paper_from_data(data, review, pool):
    """Creates Paper model from given data, review and pool"""
    medlineCitation = data[u'MedlineCitation']
    article = medlineCitation[u'Article']
    paper = Paper.objects.get_or_create(review=review, title=article[u'ArticleTitle'])[0]
    paper.review = review
    paper.authors = _get_authors(article)

    # TODO: label for section headings is lost
    # eg. StringElement('some text here', attributes={u'NlmCategory': u'METHODS', u'Label': u'METHODS'})
    abstractText = ""
    for stringElement in article[u'Abstract'][u'AbstractText']:
        abstractText += stringElement

    paper.abstract = abstractText
    paper.publish_date = _get_date(medlineCitation)
    paper.url = url_from_id(medlineCitation[u'PMID'])
    paper.notes = ""
    paper.pool = pool
    paper.save()
    return paper


def create_papers_from_ids(ids, review, pool='A'):
    """Creates papers from all of the given ids, in the given review and pool"""
    papers = read_papers_from_ids(ids)
    return map(lambda data: create_paper_from_data(data, review, pool), papers)
