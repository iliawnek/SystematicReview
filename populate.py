import os
import django
import datetime
from random import randint
from sysrev.models import User, Review, Paper

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tango-with-django.settings')
django.setup()


def populate():
    jill = add_user(username="jill", email="jill@jill.com", password="jill")
    bob = add_user(username="bob", email="bob@bob.com", password="bob")
    jen = add_user(username="jen", email="jen@jen.com", password="jen")

    adhd = add_review(user=jill,
                      title="Investigating the effects of acupuncture on children with ADHD.",
                      description="Something about ADHD.",
                      date_created=generate_random_date(),
                      last_modified=generate_random_date(recent=True),
                      query="(adhd OR adhs OR addh) AND (child OR adolescent) AND (acupuncture)")

    stress = add_review(user=jill,
                        title="Stress experienced by students during examinations",
                        description="Something about stress.",
                        date_created=generate_random_date(),
                        last_modified=generate_random_date(recent=True),
                        query="(stress) AND (student) AND (exam OR examination OR test)")

    lung_cancer = add_review(user=jill,
                             title="Development of lung cancer from inhalation of soldering fumes",
                             description="Something about lung cancer.",
                             date_created=generate_random_date(),
                             last_modified=generate_random_date(recent=True),
                             query="(solder OR soldering) AND (lung AND cancer)")

    rsi = add_review(user=jen,
                     title="RSI in drummers and guitarists",
                     description="Something about RSI.",
                     date_created=generate_random_date(),
                     last_modified=generate_random_date(recent=True),
                     query="(RSI OR repetitive OR strain OR injury) AND (drums OR drumming OR drummer OR guitar OR guitarist)")

    # TODO: get papers from PubMed using the above reviews/queries

    for user in User.objects.all():
        for review in Review.objects.filter(user=user):
            for paper in Paper.objects.filter(review=review):
                print "{0} - {1} - {2}".format(user, review, paper)


def add_user(username, email, password):
    user = User.objects.get_or_create(username=username, email=email, password=password)[0]
    user.save()
    return user


def add_review(user, title, description, date_created, last_modified, query):
    review = Review.objects.get_or_create(title=title)[0]
    review.user = user
    review.description = description
    review.date_created = date_created
    review.last_modified = last_modified
    review.query = query
    review.save()
    return review


def add_paper(review, title, authors, abstract, publish_date, url, notes, pool):
    paper = Paper.objects.get_or_create(title=title)[0]
    paper.review = review
    paper.authors = authors
    paper.abstract = abstract
    paper.publish_date = publish_date
    paper.url = url
    paper.notes = notes
    paper.pool = pool
    paper.save()
    return paper


def generate_random_date(recent=False):
    if recent:
        year = 2016
        month = randint(1, 2)
    else:
        year = randint(2010, 2015)
        month = randint(1, 12)
    day = randint(1, 28)
    hour = randint(0, 23)
    minute = randint(0, 59)

    return datetime.datetime(year, month, day, hour, minute)


# Start execution here!
if __name__ == '__main__':
    print "Starting population script..."
    populate()
