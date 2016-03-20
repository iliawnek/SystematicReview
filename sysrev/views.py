from django.shortcuts               import render
from django.http                    import HttpResponse, HttpResponseRedirect, Http404
from django.contrib.auth.decorators import login_required
from django.contrib.formtools.wizard.views import SessionWizardView
from django.utils.decorators        import method_decorator
from django.views.generic.edit      import CreateView, UpdateView
from django.views.generic           import ListView, DetailView
from django.core.urlresolvers       import reverse


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

        in_progress_reviews = Review.objects.order_by('-last_modified').filter(participants=self.request.user, completed=False)
        completed_reviews   = Review.objects.order_by('-last_modified').filter(participants=self.request.user, completed=True)
        reviews             = list(chain(in_progress_reviews, completed_reviews))

        for i in range(0, len(reviews)):
            reviews[i] = {"review": reviews[i],
                          "count": paper_pool_counts(reviews[i]),
                          "percent": paper_pool_percentages(reviews[i])}

        context["reviews"] = reviews
        return context

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewListView, self).dispatch(*args, **kwargs)


def paper_pool_percentages(review):
    counts = paper_pool_counts(review)
    total = sum(counts.values()) - counts["remaining"]

    if total is not 0:
        abstract = (float(counts["abstract"]) / float(total)) * 100.0
        document = (float(counts["document"]) / float(total)) * 100.0
        final    = (float(counts["final"])    / float(total)) * 100.0
        rejected = (float(counts["rejected"]) / float(total)) * 100.0
        return {"abstract": abstract,
                "document": document,
                "final": final,
                "rejected": rejected,
                "progress": final + rejected}
    else:
        return



def paper_pool_counts(review):
    relevant_papers = Paper.objects.filter(review=review)
    abstract_count = relevant_papers.filter(pool="A").count()
    document_count = relevant_papers.filter(pool="D").count()
    final_count = relevant_papers.filter(pool="F").count()
    rejected_count = relevant_papers.filter(pool="R").count()
    return {"abstract": abstract_count,
            "document": document_count,
            "final": final_count,
            "rejected": rejected_count,
            "remaining": abstract_count + document_count}


class ReviewDetailView(DetailView):
    model = Review

    def get_context_data(self, object=None):
        context = {}
        try:
            if self.request.user in object.participants.all():
                context["review"] = object
                context["count"] = paper_pool_counts(object)
                context["percent"] = paper_pool_percentages(object)
            else:
                raise Http404("Review not found")
        except Review.DoesNotExist:
            raise Http404("Review not found")
        return context

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewDetailView, self).dispatch(*args, **kwargs)


class PaperDetailView(DetailView):
    model = Paper

    def get_context_data(self, object=None):
        context = {}
        try:
            if (object.review == Review.objects.get(pk=self.kwargs['pk'])) and (self.request.user in object.review.participants.all()):
                paper = Paper.objects.get(pk=self.kwargs['pk2'])
                context["paper"] = paper
                context["review"] = object.review

                titles = {'A': 'Abstract screening', 'D': 'Document screening', 'F': 'Final document', 'R': 'Rejected document'};
                context["title"] = titles[paper.pool]

                context["to_judge"] = ('A', 'D')
                context["to_embed_full"] = ('D', 'F')

                context["count"] = paper_pool_counts(object.review)
                context["percent"] = paper_pool_percentages(object.review)
            else:
                raise Http404("Paper not found")
        except Review.DoesNotExist:
            raise Http404("Paper not found")
        return context

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(PaperDetailView, self).dispatch(*args, **kwargs)


class ReviewCreateWizard(SessionWizardView):
    form_list     = [ReviewCreateStep1, ReviewCreateStep2]
    template_name = "sysrev/review_create_wizard.html"

    def done(self, form_list, **kwargs):
        s1 = form_list[0].cleaned_data
        s2 = form_list[1].cleaned_data

        review = Review()


        review.title       = s1["title"]
        review.description = s1["description"]
        review.query       = s2["query"]

        review.save()

        review.participants.add(self.request.user)
        invited = filter(lambda i: i, map(lambda l: str.strip(str(l)), s1["invited"].splitlines()))
        review.invite(invited)

        review.save()

        return HttpResponseRedirect(reverse("review", args=(review.pk,)))

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewCreateWizard, self).dispatch(*args, **kwargs)
