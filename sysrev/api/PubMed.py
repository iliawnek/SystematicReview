from Bio import Entrez

from sysrev.models import Paper

# Abuse contact for if we query too much/cause a problem
Entrez.email = 'systematicreview@nallar.me'

db = "pubmed"


def pubmed_query(query):
    return Entrez.read(Entrez.esearch(db=db,
                                      sort='relevance',
                                      retmax=1000,
                                      term=query,
                                      field="title",
                                      rettype='json'))


def get_data_from_query(query):
    response = pubmed_query(query)
    data = {"count": response[u'Count'],
            "ids": response[u'IdList']}
    return data


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

    # TODO: label for section headings is lost
    # eg. StringElement('some text here', attributes={u'NlmCategory': u'METHODS', u'Label': u'METHODS'})
    abstractText = ""
    for stringElement in article[u'Abstract'][u'AbstractText']:
        abstractText += stringElement

    paper.abstract = abstractText
    paper.publish_date = get_date(medlineCitation)
    paper.url = url_from_id(medlineCitation[u'PMID'])
    paper.notes = ""
    paper.pool = pool
    paper.save()
    return paper


def create_papers_from_ids(ids, review, pool='A'):
    papers = read_papers_from_ids(ids)
    return map(lambda data: create_paper_from_data(data, review, pool), papers)

