from django.test import TestCase
from api import PubMed

class PubmedQueryTestCase(TestCase):
    def test_query(self):
        result = PubMed.query("smoking")
        self.assertGreater(result[u'Count'], 25000, "Expected >25000 results for smoking")