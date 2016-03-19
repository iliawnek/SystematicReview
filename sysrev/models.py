from django.db import models
from django.contrib.auth.models import User
from django.template.defaultfilters import slugify
from django.core.exceptions import ValidationError


class Review(models.Model):
    participants       = models.ManyToManyField(User)

    title              = models.CharField(max_length=128, unique=False)
    slug               = models.SlugField()
    description        = models.TextField(default="")
    date_created       = models.DateTimeField(auto_now_add=True)
    last_modified      = models.DateTimeField(auto_now=True)
    completed          = models.BooleanField(default=False)

    query              = models.TextField(default="")
    abstract_pool_size = models.IntegerField(default=0)
    document_pool_size = models.IntegerField(default=0)
    final_pool_size    = models.IntegerField(default=0)
    rejected_pool_size = models.IntegerField(default=0)

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

    def __unicode__(self):
        return str(self.review) + " - " + self.title
