from django.test import TestCase
from api import PubMed


class PubmedQueryTestCase(TestCase):
    def test_query(self):
        result = PubMed.query("smoking")
        self.assertGreater(result[u'Count'], 25000, "Expected >25000 results for smoking")

    def test_paper(self):
        result = PubMed.paper([25929677])

        self.assertEquals(len(result[0][u'MedlineCitation'][u'Article'][u'AuthorList']),
                          7,
                          "25929677 should have 7 authors")
