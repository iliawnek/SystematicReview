from django.test import TestCase

from api import PubMed
from sysrev.models import Review


class PubmedQueryTestCase(TestCase):
    def test_query(self):
        result = PubMed.query("smoking")
        self.assertGreater(result[u'Count'], 25000, "Expected >25000 results for smoking")

    def test_paper(self):
        result = PubMed.read_papers_from_ids([25929677])

        self.assertEquals(len(result[0][u'MedlineCitation'][u'Article'][u'AuthorList']),
                          7,
                          "25929677 should have 7 authors")

    def test_create_papers_from_ids(self):
        review = Review.objects.get_or_create(title="Investigating the effects of acupuncture on children with ADHD")[0]
        result = PubMed.create_papers_from_ids([26502548], review)[0]
        print result.title
        self.assertEquals("[A Meta-analysis on Acupuncture Treatment of Attention Deficit/Hyperactivity Disorder].",
                          result.title)

    def test_adhd_query(self):
        query = """(adhd OR adhs OR addh) AND (child OR adolescent) AND acupuncture"""
        result = PubMed.get_ids_from_query(query)
        self.assertGreater(len(result), 0, "Expected some results for ADHD query")
