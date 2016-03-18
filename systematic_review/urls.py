from django.conf.urls import patterns, include, url
from django.contrib import admin
from sysrev import views
from sysrev.views import *
from registration.backends.simple.views import RegistrationView


class MyRegistrationView(RegistrationView):
    def get_success_url(self, user=None):
        return "/"

urlpatterns = patterns(
    '',
    url(r'^$',                          ReviewListView.as_view(),     name='index'),
    url(r'^review/(?P<pk>\d+)(-([\w\-]+))?/$', ReviewDetailView.as_view(),   name='review'),
    url(r'^create/',                    ReviewCreateWizard.as_view(), name='create'),
    url(r'^profile/',                   ProfileView.as_view(),        name='profile'),
    url(r'^admin/',                     include(admin.site.urls)),
    url(r'^accounts/register/$',        MyRegistrationView.as_view(), name='registration_register'),
    url(r'^accounts/',                  include('registration.backends.simple.urls')),
    url(r'^accounts/',                  include('django.contrib.auth.urls')),
)
