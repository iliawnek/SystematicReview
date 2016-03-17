from django.shortcuts import render
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required

from sysrev.models import Review


@login_required
def index(request):
    context_dict = {"reviews": Review.objects.filter(user=request.user)}
    return render(request, 'sysrev/index.html', context_dict)


@login_required
def create(request):
    pass
