import os
import django
import datetime
from random import randint

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'systematic_review.settings')
django.setup()

from sysrev.models import User, Review, Paper


def populate():
    jill = User.objects.create_user(username="jill", email="jill@example.com", password="jill")
    jill.save()
    bob = User.objects.create_user(username="bob", email="bob@example.com", password="bob")
    bob.save()
    jen = User.objects.create_user(username="jen", email="jen@example.com", password="jen")
    jen.save()

    adhd = add_review(users=[jill, bob],
                      title="Investigating the effects of acupuncture on children with ADHD",
                      description="""This is a sentence about why this particular review was created. It involves
                      something to do with ADHD and acupuncture. These are a few more test sentences to make the description
                      longer. This is another sentence.""",
                      date_created=generate_random_date(),
                      last_modified=generate_random_date(recent=True),
                      query="(adhd OR adhs OR addh) AND (child OR adolescent) AND (acupuncture)")

    stress = add_review(users=[jill],
                        title="Stress experienced by students during examinations",
                        description="Something about stress.",
                        date_created=generate_random_date(),
                        last_modified=generate_random_date(recent=True),
                        query="(stress) AND (student) AND (exam OR examination OR test)",
                        completed=True)

    lung_cancer = add_review(users=[jill, jen],
                             title="Development of lung cancer from inhalation of soldering fumes",
                             description="Something about lung cancer.",
                             date_created=generate_random_date(),
                             last_modified=generate_random_date(recent=True),
                             query="(solder OR soldering) AND (lung AND cancer)")

    rsi = add_review(users=[jen, bob],
                     title="RSI in drummers and guitarists",
                     description="Something about RSI.",
                     date_created=generate_random_date(),
                     last_modified=generate_random_date(recent=True),
                     query="(RSI OR repetitive OR strain OR injury) AND (drums OR drumming OR drummer OR guitar OR guitarist)")

    for pool in ["A", "A", "A", "D", "D", "F", "F", "F", "F", "R"]:
        add_paper(review=adhd,
                  title="A Meta-analysis on Acupuncture Treatment of Attention Deficit/Hyperactivity Disorder",
                  authors="Ni XQ, Zhang JY, Han XM, Yin DQ",
                  abstract="""OBJECTIVE:
    To assess the efficacy and safety of acupuncture in treating attention-deficit/hyperactivity disorder (ADHD) children.
    METHODS:
    A literature search was conducted to retrieve randomized cotrolled clinical trials of acupuncture in treating ADHD covering the period of the years of establishment of the databases to January 2014 from database of CBM, CNKI, PubMed, Cochrane Library by using key words "attention deficit hyperactivity disorder" "hyperactivity""minimal brain dysfunction" "acupuncture". Two independent researchers extracted data from located articles in a pre-defined structured way, and consulted the third researcher if necessary.
    RESULTS:
    Thirteen original trials including 1 304 cases of children with ADHD were obtained in this study according to our included criteria and excluded criteria. In these trials, acupuncture intervention alone, or acupuncture plus pharmacotherapy (methylphenidate, haloperidol) or acupuncture plus behavioral therapy were compared with simple pharmacotherapy or behavioral therapy alone. Results of Meta-analysis indicated that the total effective rate and Conners' index of hyperactivity (CIH) score-reduction rate in the acupuncture group were significantly superior to those of the other treatment groups [OR = 2.22, 95% CI (1.65, 3.00), Z = 5.22, P < 0.00001] [SMD = -0.94, 95% CI (-1.41, -0.47), Z = 3.89, P < 0.0001]. Acupuncture treatment is more effective than haloperidol in reducing the score of Conners' Rating Scale for ADHD [SMD = -7.28, 95% CI (-8.32, -6.23), Z = 13.62, P < 0.00001]. Acupuncture is similarly effective as Methylphenidate (Ritalin) in improving the Chinese medicine syndrome (liver-kidney yin hypoactivity) of children with ADHD [SMD = -1.14, 95% CI (-2.53, 0.25), Z = 1.60, P = 0.11]. Less severe adverse effects were reported with acupuncture therapy than the pharmacotherapy (poor appetite, dry mouth, nausea and constipation). These effects were not likely due to publication bias (approximately symmetry funnel plot, Egger's test P > 0.1).
    CONCLUSION:
    Acupuncture is an effective and safe therapy in treating ADHD, combined administration of acupuncture and pharmacotherapy or behavioral therapy is more effective than the pharmacotherapy or behavioral therapy alone. However, more rigorously designed and high-quality RCTs are needed to confirm the above conclusion.""",
                  publish_date=datetime.date(2015, 8, 1),
                  url="http://www.ncbi.nlm.nih.gov/pubmed/26502548",
                  notes="this is a note",
                  pool=pool)

    # TODO: get papers from PubMed using the above reviews/queries

    for user in User.objects.all():
        for review in Review.objects.filter(participants=user):
            for paper in Paper.objects.filter(review=review):
                print "{0} - {1} - {2}".format(user, review, paper)


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


def add_paper(review, title, authors, abstract, publish_date, url, notes, pool="A"):
    paper = Paper.objects.get_or_create(review=review, title=title)[0]
    paper.review = review
    paper.authors = authors
    paper.abstract = abstract
    paper.publish_date = publish_date
    paper.url = url
    paper.notes = notes
    paper.pool = pool
    paper.save()
    if pool is "A":
        review.abstract_pool_size += 1
    elif pool is "D":
        review.document_pool_size += 1
    elif pool is "F":
        review.final_pool_size += 1
    elif pool is "R":
        review.rejected_pool_size += 1

    review.save()
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
