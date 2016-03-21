from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse
from django.db                      import models
from django.template.defaultfilters import slugify

from sysrev.api.PubMed import _get_authors, _get_date, url_from_id, read_papers_from_ids

from sysrev.api import PubMed


class Review(models.Model):
    participants       = models.ManyToManyField(User)

    title              = models.CharField(max_length=128, unique=False)
    slug               = models.SlugField()
    description        = models.TextField(default="")
    date_created       = models.DateTimeField(auto_now_add=True)
    last_modified      = models.DateTimeField(auto_now=True)
    completed          = models.BooleanField(default=False)
    date_completed     = models.DateTimeField(default=None, null=True)
    query              = models.TextField(default="")


    def perform_query(self):
        # TODO: discard existing papers if there are any
        data = PubMed.get_data_from_query(self.query)
        Paper.create_papers_from_pubmed_ids(data[u'IdList'], self)


    def paper_pool_percentages(self):
        counts = self.paper_pool_counts()
        total = counts["total"]

        if total is not 0:
            abstract = (float(counts["abstract"]) / float(total)) * 100.0
            document = (float(counts["document"]) / float(total)) * 100.0
            final    = (float(counts["final"])    / float(total)) * 100.0
            rejected = (float(counts["rejected"]) / float(total)) * 100.0
            return {"abstract": abstract,
                    "document": document,
                    "final": final,
                    "rejected": rejected,
                    "progress": final + rejected}
        else:
            return



    def paper_pool_counts(self):
        relevant_papers = Paper.objects.filter(review=self)
        abstract_count = relevant_papers.filter(pool="A").count()
        document_count = relevant_papers.filter(pool="D").count()
        final_count = relevant_papers.filter(pool="F").count()
        rejected_count = relevant_papers.filter(pool="R").count()
        return {"abstract": abstract_count,
                "document": document_count,
                "final": final_count,
                "rejected": rejected_count,
                "remaining": abstract_count + document_count,
                "total": abstract_count + document_count + final_count + rejected_count}

    def invite(self, invitees):
        for invitee in invitees:
            user = None
            if invitee.find("@") == -1:
                user = User.objects.get(username=invitee)
            else:
                user = User.objects.get(email=invitee)
            self.participants.add(user)

    def clean(self):
        if (not self.participants) or self.participants.count() < 1:
            raise ValidationError('Need at least one participant')

    def save(self, *args, **kwargs):
        self.slug = slugify(self.title)
        super(Review, self).save()

    def get_absolute_url(self):
        return reverse('review_detail', args=[str(self.pk)])[:-1] + "-" + self.slug

    def __unicode__(self):
        return str(self.pk) + ": " + self.title


class Paper(models.Model):
    ABSTRACT_POOL = 'A'
    DOCUMENT_POOL = 'D'
    FINAL_POOL    = 'F'
    REJECTED      = 'R'
    POOLS         = (
        (ABSTRACT_POOL, 'Abstract pool'),
        (DOCUMENT_POOL, 'Document pool'),
        (FINAL_POOL,    'Final pool'),
        (REJECTED,      'Rejected')
    )

    review       = models.ForeignKey(Review)
    title        = models.CharField(max_length=128)
    authors      = models.CharField(max_length=128)
    abstract     = models.TextField(default="")
    publish_date = models.DateField(null=True)
    url          = models.URLField(default="")
    notes        = models.TextField(default="")
    pool         = models.CharField(max_length=1, choices=POOLS, default=ABSTRACT_POOL)

    @staticmethod
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
        try:
            for stringElement in article[u'Abstract'][u'AbstractText']:
                abstractText += stringElement
        except KeyError:
            pass

        paper.abstract = abstractText
        paper.publish_date = _get_date(medlineCitation)
        paper.url = url_from_id(medlineCitation[u'PMID'])
        paper.notes = ""
        paper.pool = pool
        paper.save()
        return paper

    @staticmethod
    def create_papers_from_pubmed_ids(ids, review, pool='A'):
        """Creates papers from all of the given ids, in the given review and pool"""
        papers = read_papers_from_ids(ids)

        # Commit all papers in single transaction
        # Improves performance, as django won't automatically commit after every save call when creating lots of papers
        from django.db import transaction

        with transaction.atomic():
            return map(lambda data: Paper.create_paper_from_data(data, review, pool), papers)

    def get_absolute_url(self):
        return self.review.get_absolute_url() + "/" + str(self.pk)

    def __unicode__(self):
        return str(self.review) + " - " + self.title
