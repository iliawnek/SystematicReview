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
    query = models.CharField()
    date_created = models.DateTimeField(auto_now_add=True)
    last_modified = models.DateTimeField(auto_now=True)


class Paper(models.Model):
    title = models.CharField()
    abstract = models.TextField()
    author = models.CharField()
    publish_date = models.DateField()
    url = models.URLField()

    ABSTRACT_STAGE = 'A'
    DOCUMENT_STAGE = 'D'
    FINAL_STAGE = 'F'
    STAGES = ((ABSTRACT_STAGE, 'Abstract stage'),
              (DOCUMENT_STAGE, 'Document stage'),
              (FINAL_STAGE, 'Final stage'))

    stage = models.CharField(max_length=1, choices=STAGES, default=ABSTRACT_STAGE)
