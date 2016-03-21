from django.conf.urls          import patterns, include, url
from django.contrib            import admin
from django.views.generic.base import RedirectView
from sysrev.views              import *

urlpatterns = patterns(
    '',

    # Index redirect
    url(r'^$',                                              RedirectView.as_view(url='/reviews'), name='index'),

    # About
    url(r'^about/$',                                        AboutView.as_view(),          name='about'),

    # Reviews
    url(r'^reviews$',                                       ReviewListView.as_view(),     name='review_list'),
    url(r'^reviews/(?P<pk>\d+)(-([\w\-]+))?/$',             ReviewDetailView.as_view(),   name='review_detail'),
    url(r'^reviews/(?P<pk>\d+)(-([\w\-]+))?/edit$',         ReviewUpdateView.as_view(),   name='review_update'),
    url(r'^reviews/(?P<pk>\d+)(-([\w\-]+))?/delete$',       ReviewDeleteView.as_view(),   name='review_delete'),
    url(r'^reviews/(?P<pk>\d+)(-([\w\-]+))?/dl$',           ReviewDownloadView.as_view(), name='review_download'),
    url(r'^reviews/create$',                                ReviewCreateWizard.as_view(), name='review_create'),
    url(r'^reviews/(?P<pk>\d+)(-([\w\-]+))?/work/$',        ReviewWorkView.as_view(),     name='review_work'),

    # Papers
    url(r'^reviews/(?P<pk>\d+)(-([\w\-]+))?/(?P<pk2>\d+)/$', PaperDetailView.as_view(),    name='paper_detail'),
    url(r'^reviews/(?P<pk>\d+)(-([\w\-]+))?/(?P<pk2>\d+)/(?P<choice>\w+)$', PaperChoiceView.as_view(), name='paper_choice'),

    # Authentication
    url(r'^profile/',                                       ProfileView.as_view(),        name='profile'),
    url(r'^accounts/register/$',                            SRRegistrationView.as_view(), name='registration_register'),
    url(r'^accounts/',                                      include('registration.backends.simple.urls')),
    url(r'^accounts/',                                      include('django.contrib.auth.urls')),

    # Administration
    url(r'^admin/',                                         include(admin.site.urls))
)
