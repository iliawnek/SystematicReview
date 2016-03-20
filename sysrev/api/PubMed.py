from Bio import Entrez

# Abuse contact for if we query too much/cause a problem
from sysrev.models import Paper

Entrez.email = 'systematicreview@nallar.me'

db = "pubmed"

def query(query):
    return Entrez.read(Entrez.esearch(db=db,
                                      sort='relevance',
                                      retmax=1000,
                                      term=query,
                                      field="title",
                                      rettype='json'))


def read_papers_from_ids(ids):
    # Convert list of ids to comma-delimited list
    idList = ",".join(map(str, ids))
    return Entrez.read(Entrez.efetch(db=db, id=idList, retmode='xml'))


def get_authors(article):
    authorList = article[u'AuthorList']
    return ", ".join(map(lambda author: author[u'ForeName'] + " " + author[u'LastName'], authorList))


def url_from_id(id):
    return "https://www.ncbi.nlm.nih.gov/pubmed/" + id


def get_date(medlineCitation):
    date = medlineCitation[u'DateCreated']
    return date[u'Year'] + '-' + date[u'Month'] + '-' + date[u'Day']


def create_paper_from_data(data, review, pool):
    medlineCitation = data[u'MedlineCitation']
    article = medlineCitation[u'Article']
    paper = Paper.objects.get_or_create(review=review, title=article[u'ArticleTitle'])[0]
    paper.review = review
    paper.authors = get_authors(article)
    paper.abstract = article[u'Abstract'][u'AbstractText']
    paper.publish_date = get_date(medlineCitation)
    paper.url = url_from_id(medlineCitation[u'PMID'])
    paper.notes = ""
    paper.pool = pool
    paper.save()
    return paper


def create_papers_from_ids(ids, review, pool='A'):
    papers = read_papers_from_ids(ids)
    return map(lambda data: create_paper_from_data(data, review, pool), papers)
