from datetime import datetime
from itertools import chain
from random import randint

from django.contrib.auth.decorators import login_required
from django.contrib.formtools.wizard.views import SessionWizardView
from django.http import HttpResponseRedirect, Http404
from django.utils.decorators        import method_decorator
from django.views.decorators.cache  import cache_control
from django.views.generic           import ListView, DetailView, RedirectView
from django.views.generic.edit import UpdateView, DeleteView
from registration.backends.simple.views import RegistrationView

from sysrev.forms  import *


class SRRegistrationView(RegistrationView):
    def get_success_url(self, user=None):
        return "/"


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
            reviews[i] = {"review":  reviews[i],
                          "count":   reviews[i].paper_pool_counts(),
                          "percent": reviews[i].paper_pool_percentages()}

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
            if self.request.user in object.participants.all():
                context["review"]  = object
                context["count"]   = object.paper_pool_counts()
                context["percent"] = object.paper_pool_percentages()
                context["final_papers"] = Paper.objects.filter(review=object, pool="F")
                context["rejected_papers"] = Paper.objects.filter(review=object, pool="R")
            else:
                raise Http404("Review not found")
        except Review.DoesNotExist:
            raise Http404("Review not found")
        return context

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewDetailView, self).dispatch(*args, **kwargs)


class ReviewUpdateView(UpdateView):
    model  = Review
    fields = ['title', 'description', 'participants', 'query']

    def get_success_url(self):
        return self.request.path[:-7]

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewUpdateView, self).dispatch(*args, **kwargs)


class ReviewDeleteView(DeleteView):
    model       = Review
    success_url = "/"

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewDeleteView, self).dispatch(*args, **kwargs)


class ReviewWorkView(RedirectView):
    permanent = False

    def get(self, request, *args, **kwargs):
        try:
            review = Review.objects.get(pk=self.kwargs['pk'])
            papers = Paper.objects.filter(review=review)
            counts = review.paper_pool_counts()

            if counts["abstract"] == 0 and counts["document"] == 0:
                review.completed = True
                review.date_completed = datetime.now()
                review.save()
                self.url = review.get_absolute_url()
                return super(ReviewWorkView, self).get(request, args, **kwargs)
            elif counts["abstract"] > 0:
                papers = papers.filter(pool="A")
            elif counts["document"] > 0:
                papers = papers.filter(pool="D")

            paper = papers.all()[randint(0, papers.count()-1)]
            self.url = paper.get_absolute_url()

            return super(ReviewWorkView, self).get(request, args, **kwargs)
        except Review.DoesNotExist:
            raise Http404("Paper not found")
        except Paper.DoesNotExist:
            raise Http404("Paper not found")

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewWorkView, self).dispatch(*args, **kwargs)


class PaperDetailView(DetailView):
    model = Paper

    def get_context_data(self, object=None):
        context = {}
        try:
            review = Review.objects.get(pk=self.kwargs['pk'])
            if self.request.user in review.participants.all():
                paper = Paper.objects.get(pk=self.kwargs['pk2'])
                context["paper"] = paper
                context["review"] = review

                titles = {'A': 'Abstract screening', 'D': 'Document screening', 'F': 'Final document', 'R': 'Rejected document'}
                context["title"] = titles[paper.pool]

                context["to_judge"] = ('A', 'D')
                context["to_embed_full"] = ('D', 'F')

                context["count"]   = review.paper_pool_counts()
                context["percent"] = review.paper_pool_percentages()
            else:
                raise Http404("Paper not found")
        except Review.DoesNotExist:
            raise Http404("Paper not found")
        except Paper.DoesNotExist:
            raise Http404("Paper not found")
        return context

    @method_decorator(login_required)
    @cache_control(no_cache=True, must_revalidate=True, no_store=True)
    def dispatch(self, *args, **kwargs):
        return super(PaperDetailView, self).dispatch(*args, **kwargs)


class PaperChoiceView(RedirectView):
    permanent = False

    def get(self, request, *args, **kwargs):
        try:
            review = Review.objects.get(pk=self.kwargs['pk'])
            paper = Paper.objects.get(pk=self.kwargs['pk2'], review=review)
            choice = self.kwargs['choice']
            if choice == "yes":
                if paper.pool == "A":
                    paper.pool = "D"
                elif paper.pool == "D":
                    paper.pool = "F"
            elif choice == "no":
                paper.pool = "R"
            else:
                raise Http404("Invalid choice")
            paper.save()
            self.url = review.get_absolute_url() + "/work/"
            return super(PaperChoiceView, self).get(request, args, **kwargs)
        except Review.DoesNotExist:
            raise Http404("Review not found")
        except Paper.DoesNotExist:
            raise Http404("Paper not found")

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(PaperChoiceView, self).dispatch(*args, **kwargs)


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

        ids = PubMed.get_data_from_query(review.query)["ids"]
        Paper.create_papers_from_pubmed_ids(ids, review)

        return HttpResponseRedirect(review.get_absolute_url())

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(ReviewCreateWizard, self).dispatch(*args, **kwargs)
