from django.conf.urls import patterns, include, url
from django.contrib   import admin
from sysrev.views     import *

urlpatterns = patterns(
    '',
    url(r'^$',                                              ReviewListView.as_view(),     name='index'),
    url(r'^review/(?P<pk>\d+)(-([\w\-]+))?/$',              ReviewDetailView.as_view(),   name='review'),
    url(r'^review/(?P<pk>\d+)(-([\w\-]+))?/update$',        ReviewUpdateView.as_view(),   name='review_update'),
    url(r'^review/(?P<pk>\d+)(-([\w\-]+))?/delete$',        ReviewDeleteView.as_view(),   name='review_delete'),
    url(r'^create/',                                        ReviewCreateWizard.as_view(), name='create'),

    url(r'^review/(?P<pk>\d+)(-([\w\-]+))?/work/$',         WorkView.as_view(),           name='work'),
    url(r'^review/(?P<pk>\d+)(-([\w\-]+))?/(?P<pk2>\d+)/$', PaperDetailView.as_view(),    name='paper'),

    url(r'^profile/',                                       ProfileView.as_view(),        name='profile'),
    url(r'^accounts/register/$',                            SRRegistrationView.as_view(), name='registration_register'),
    url(r'^accounts/',                                      include('registration.backends.simple.urls')),
    url(r'^accounts/',                                      include('django.contrib.auth.urls')),

    url(r'^admin/',                                         include(admin.site.urls)),
)
