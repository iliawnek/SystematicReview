from django.shortcuts               import render
from django.http                    import HttpResponse
from django.contrib.auth.decorators import login_required
from django.utils.decorators        import method_decorator
from django.views.generic.edit      import CreateView, UpdateView
from django.views.generic           import ListView, DetailView

from itertools import chain

from sysrev.models import *
from sysrev.forms  import *


class ProfileView(UpdateView):
    template_name = "sysrev/profile_form.html"
    form_class    = ProfileForm
    model         = User
    success_url   = "#"

    def get_object(self, queryset=None):
        return self.request.user

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ProfileView, self).dispatch(*args, **kwargs)


class ReviewListView(ListView):
    model = Review

    def get_context_data(self, **kwargs):
        context = super(ReviewListView, self).get_context_data(**kwargs)

        in_progress_reviews = Review.objects.order_by('-last_modified').filter(user=self.request.user, completed=False)
        completed_reviews   = Review.objects.order_by('-last_modified').filter(user=self.request.user, completed=True)
        reviews             = list(chain(in_progress_reviews, completed_reviews))

        # progress bar
        for i in range(0, len(reviews)):
            abstract = reviews[i].abstract_pool_size
            document = reviews[i].document_pool_size
            final    = reviews[i].final_pool_size
            rejected = reviews[i].rejected_pool_size
            total    = abstract + document + final + rejected

            if total is not 0:
                abstract = int((float(abstract) / float(total)) * 100.0)
                document = int((float(document) / float(total)) * 100.0)
                final    = int((float(final)    / float(total)) * 100.0)
                rejected = int((float(rejected) / float(total)) * 100.0)

                total = abstract + document + final + rejected
                if total < 100:
                    abstract += 100 - total

            reviews[i] = {"review": reviews[i], "abstract": abstract, "document": document, "final": final, "rejected": rejected}

        context["reviews"] = reviews
        return context

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewListView, self).dispatch(*args, **kwargs)


class ReviewDetailView(DetailView):
    model = Review

    def get_context_data(self, object=None):
        context = {}
        try:
            if object.user == self.request.user:
                context["review"] = object
        except Review.DoesNotExist:
            pass
        return context

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewDetailView, self).dispatch(*args, **kwargs)


class ReviewCreateView(CreateView):
    model = Review

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewCreateView, self).dispatch(*args, **kwargs)
