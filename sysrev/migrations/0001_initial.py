# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(max_length=128)),
                ('authors', models.CharField(max_length=128)),
                ('abstract', models.TextField(default=b'')),
                ('publish_date', models.DateField(null=True)),
                ('url', models.URLField(default=b'')),
                ('pubmed_id', models.CharField(max_length=16)),
                ('notes', models.TextField(default=b'')),
                ('pool', models.CharField(default=b'A', max_length=1, choices=[(b'A', b'Abstract pool'), (b'D', b'Document pool'), (b'F', b'Final pool'), (b'R', b'Rejected')])),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Review',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(max_length=128)),
                ('slug', models.SlugField()),
                ('description', models.TextField(default=b'')),
                ('date_created', models.DateTimeField(auto_now_add=True)),
                ('last_modified', models.DateTimeField(auto_now=True)),
                ('completed', models.BooleanField(default=False)),
                ('date_completed', models.DateTimeField(default=None, null=True)),
                ('query', models.TextField(default=b'')),
                ('participants', models.ManyToManyField(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='paper',
            name='review',
            field=models.ForeignKey(to='sysrev.Review'),
            preserve_default=True,
        ),
    ]
