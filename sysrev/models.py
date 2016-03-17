from django.db import models
from django.contrib.auth.models import User


class Review(models.Model):
    user = models.ForeignKey(User)
    title = models.CharField(max_length=128, unique=True)
    description = models.TextField(default="")
    date_created = models.DateTimeField(auto_now_add=True)
    last_modified = models.DateTimeField(auto_now=True)
    completed = models.BooleanField(default=False)

    query = models.TextField(default="")
    abstract_pool_size = models.IntegerField(default=0)
    document_pool_size = models.IntegerField(default=0)
    final_pool_size = models.IntegerField(default=0)
    rejected_pool_size = models.IntegerField(default=0)

    def __unicode__(self):
        return str(self.user) + " - " + self.title


class Paper(models.Model):
    review = models.ForeignKey(Review)
    title = models.CharField(max_length=128)
    authors = models.CharField(max_length=128)
    abstract = models.TextField(default="")
    publish_date = models.DateField(null=True)
    url = models.URLField(default="")
    notes = models.TextField(default="")

    ABSTRACT_POOL = 'A'
    DOCUMENT_POOL = 'D'
    FINAL_POOL = 'F'
    REJECTED = 'R'
    POOLS = ((ABSTRACT_POOL, 'Abstract pool'),
             (DOCUMENT_POOL, 'Document pool'),
             (FINAL_POOL, 'Final pool'),
             (REJECTED, 'Rejected'))

    pool = models.CharField(max_length=1, choices=POOLS, default=ABSTRACT_POOL)

    def __unicode__(self):
        return str(self.review) + " - " + self.title
