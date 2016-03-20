from django.test import TestCase

from api import PubMed
from sysrev.models import Paper, Review


class PubmedQueryTestCase(TestCase):
    test_ids = [
        26502548,
        25929677,
    ]

    def test_query(self):
        result = int(PubMed.get_data_from_query("smoking")[u'Count'])
        min_ = 25000
        self.assertGreater(result, min_, "Expected >%d results for smoking, got %d" % (min_, result))

    def test_paper(self):
        result = PubMed.read_papers_from_ids([25929677])

        self.assertEquals(len(result[0][u'MedlineCitation'][u'Article'][u'AuthorList']),
                          7,
                          "25929677 should have 7 authors")

    def test_create_papers_from_ids(self):
        review = Review.objects.get_or_create(title="Investigating the effects of acupuncture on children with ADHD")[0]
        result = Paper.create_papers_from_pubmed_ids(self.test_ids, review)[0]
        self.assertEquals("[A Meta-analysis on Acupuncture Treatment of Attention Deficit/Hyperactivity Disorder].",
                          result.title)

    def test_adhd_query(self):
        query = """(adhd OR adhs OR addh) AND (child OR adolescent) AND acupuncture"""
        self._query_min_results(query)

    def test_lung_cancer_query(self):
        query = """(solder OR soldering) AND (lung AND cancer)"""
        self._query_min_results(query)

    def _query_min_results(self, query, min=0):
        """Tests given query, ensures at least min results are found"""
        result = PubMed.get_ids_from_query(query)
        self.assertGreater(len(result), 0, "Expected at lest " + str(min) + " results for query '" + query + "'")
