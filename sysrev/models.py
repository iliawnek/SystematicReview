from django.db import models


class Researcher(models.Model):
    # TODO: username
    # TODO: password
    # TODO: email
    pass


class Review(models.Model):
    user = models.ForeignKey(Researcher)
    title = models.CharField(max_length=128)
    description = models.TextField()
    date_created = models.DateTimeField(auto_now_add=True)
    last_modified = models.DateTimeField(auto_now=True)

    query = models.TextField()
    abstract_pool_size = models.IntegerField()
    document_pool_size = models.IntegerField()
    final_pool_size = models.IntegerField()


class Paper(models.Model):
    review = models.ForeignKey(Review)
    title = models.CharField(max_length=128)
    authors = models.CharField(max_length=128)
    abstract = models.TextField()
    publish_date = models.DateField()
    url = models.URLField()
    notes = models.TextField()

    ABSTRACT_POOL = 'A'
    DOCUMENT_POOL = 'D'
    FINAL_POOL = 'F'
    REJECTED = 'R'
    POOLS = ((ABSTRACT_POOL, 'Abstract pool'),
             (DOCUMENT_POOL, 'Document pool'),
             (FINAL_POOL, 'Final pool'),
             (REJECTED, 'Rejected'))

    pool = models.CharField(max_length=1, choices=POOLS, default=ABSTRACT_POOL)
