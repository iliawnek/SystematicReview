import datetime
import os
import random

import django

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'systematic_review.settings')
django.setup()

from sysrev.api import PubMed
from sysrev.models import Paper, Review, User


def populate():
    jill = add_user("jill")
    bob = add_user("bob")
    jen = add_user("jen")

    # TODO: add a few more reviews
    adhd = add_review(users=[jill, bob],
                      title="Investigating the effects of acupuncture on children with ADHD",
                      description="The purpose of this review is to identify any potential positive or negative "
                                  "effects when applying the practice of acupuncture to young children suffering from "
                                  "attention deficit hyperactivity disorder, which is more commonly referred to by its abbreviation, ADHD.",
                      date_created=generate_random_date(),
                      last_modified=generate_random_date(recent=True),
                      query="""(adhd OR adhs OR addh) AND (child OR adolescent) AND acupuncture""")

    lung_cancer = add_review(users=[jen],
                             title="Development of lung cancer from inhalation of soldering fumes",
                             description="This review will retrieve all papers discussing links between soldering fume "
                                         "inhalation and lung cancer.",
                             date_created=generate_random_date(),
                             last_modified=generate_random_date(recent=True),
                             query="(solder OR soldering) AND (lung AND cancer)")

    rsi = add_review(users=[jen, bob],
                     title="RSI in musicians such as drummers and guitarists",
                     description="This review explores any connection between the development of repetitive strain "
                                 "injury, which is more commonly known as RSI, in musicians, or more specifically, "
                                 "people who play drums or play guitar. This issue is more frequently discussed in "
                                 "the world of computer use and office life, but rarely in the performance of music.",
                     date_created=generate_random_date(),
                     last_modified=generate_random_date(recent=True),
                     query="(RSI OR repetitive OR strain OR injury) AND (drums OR drumming OR drummer OR guitar OR guitarist)")

    stress = add_review(users=[jill],
                        title="Stress experienced by students during examinations",
                        description="This review will retrieve all medical papers which discuss the issues regarding "
                                    "stress, anxiety, and a general lack of wellbeing of students while their "
                                    "educational establishment undergo examinations.",
                        date_created=generate_random_date(),
                        last_modified=generate_random_date(recent=True),
                        query="stress AND anxiety AND student AND (exam OR examination OR test) AND (university OR college)")

    reviews = [
        adhd,
        lung_cancer,
        rsi,
        stress,
    ]

    import time

    for review in reviews:
        start = time.time()

        Paper.create_papers_from_pubmed_ids(PubMed.get_ids_from_query(review.query), review)

        end = time.time()
        duration = (end - start)
        print "Populating %s took %f seconds" % (review.title, duration)


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
    Paper.create_papers_from_pubmed_ids([id], review)


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
    from django.db import transaction

    with transaction.atomic():
        populate()
