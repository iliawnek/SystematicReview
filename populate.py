import datetime
import os
import random

import django

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'systematic_review.settings')
django.setup()

from sysrev.api import PubMed
from sysrev.models import User, Review


def populate():
    jill = add_user("jill")
    bob = add_user("bob")
    jen = add_user("jen")

    adhd = add_review(users=[jill, bob],
                      title="Investigating the effects of acupuncture on children with ADHD",
                      description="""This is a sentence about why this particular review was created. It involves
                      something to do with ADHD and acupuncture. These are a few more test sentences to make the description
                      longer. This is another sentence.""",
                      date_created=generate_random_date(),
                      last_modified=generate_random_date(recent=True),
                      query="""(adhd OR adhs OR addh) AND (child OR adolescent) AND acupuncture""")

    stress = add_review(users=[jill],
                        title="Stress experienced by students during examinations",
                        description="Something about stress.",
                        date_created=generate_random_date(),
                        last_modified=generate_random_date(recent=True),
                        query="(stress) AND (student) AND (exam OR examination OR test)")

    rsi = add_review(users=[jen, bob],
                     title="RSI in drummers and guitarists",
                     description="Something about RSI.",
                     date_created=generate_random_date(),
                     last_modified=generate_random_date(recent=True),
                     query="(RSI OR repetitive OR strain OR injury) AND (drums OR drumming OR drummer OR guitar OR guitarist)")

    # TODO: find out why lung_cancer and adhd queries return 0 results
    reviews = [stress, rsi]
    for review in reviews:
        PubMed.create_papers_from_ids(PubMed.get_ids_from_query(review.query)["ids"], review)

    add_paper_by_id(review=adhd, id=26502548)
    add_paper_by_id(review=adhd, id=26990084)


def add_review(users, title, description, date_created, last_modified, query, completed=False):
    review = Review.objects.get_or_create(title=title)[0]
    review.description = description
    review.date_created = date_created
    review.last_modified = last_modified
    review.query = query
    review.completed = completed
    review.save()
    for user in users:
        review.participants.add(user)
    return review


def add_paper_by_id(review, id):
    PubMed.create_papers_from_ids([id], review)


def add_user(name, email=None, password=None):
    if email is None:
        email = name + "@example.com"

    if password is None:
        password = name

    try:
        user = User.objects.get(username=name)
    except User.DoesNotExist:
        user = User.objects.create_user(username=name, email=email, password=password)
        user.save()

    return user


def generate_random_date(recent=False):
    if recent:
        year = 2016
        month = random.randint(1, 2)
    else:
        year = random.randint(2010, 2015)
        month = random.randint(1, 12)
    day = random.randint(1, 28)
    hour = random.randint(0, 23)
    minute = random.randint(0, 59)

    return datetime.datetime(year, month, day, hour, minute)


# Start execution here!
if __name__ == '__main__':
    print "Starting population script..."
    populate()
