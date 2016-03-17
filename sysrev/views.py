from django.shortcuts import render
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from itertools import chain

from sysrev.models import Review


@login_required
def index(request):
    context_dict = {}

    in_progress_reviews = Review.objects.order_by('-last_modified')\
        .filter(user=request.user, completed=False)
    completed_reviews = Review.objects.order_by('-last_modified')\
        .filter(user=request.user, completed=True)
    reviews = list(chain(in_progress_reviews, completed_reviews))

    for i in range(0, len(reviews)):
        abstract = reviews[i].abstract_pool_size
        document = reviews[i].document_pool_size
        final = reviews[i].final_pool_size
        rejected = reviews[i].rejected_pool_size
        total = abstract + document + final + rejected

        if total is not 0:
            abstract = int((float(abstract) / float(total)) * 100.0)
            document = int((float(document) / float(total)) * 100.0)
            final = int((float(final) / float(total)) * 100.0)
            rejected = int((float(rejected) / float(total)) * 100.0)

            total = abstract + document + final + rejected
            if total < 100:
                abstract += 100 - total

        reviews[i] = {"review": reviews[i], "abstract": abstract, "document": document, "final": final, "rejected": rejected}

    context_dict["reviews"] = reviews

    return render(request, 'sysrev/index.html', context_dict)


@login_required
def create(request):
    pass
