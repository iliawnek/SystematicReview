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
    context_dict["reviews"] = list(chain(in_progress_reviews, completed_reviews))
    return render(request, 'sysrev/index.html', context_dict)


@login_required
def create(request):
    pass
