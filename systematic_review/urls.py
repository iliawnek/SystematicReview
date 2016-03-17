from django.conf.urls import patterns, include, url
from django.contrib import admin
from sysrev import views
from registration.backends.simple.views import RegistrationView


class MyRegistrationView(RegistrationView):
    def get_success_url(self, user=None):
        return "/"

urlpatterns = patterns('',
                       url(r'^$', views.index, name='index'),
                       url(r'^create/', views.create, name='create'),
                       url(r'^admin/', include(admin.site.urls)),
                       url(r'^accounts/register/$', MyRegistrationView.as_view(), name='registration_register'),
                       url(r'^accounts/', include('registration.backends.simple.urls')),
                       url(r'^accounts/', include('django.contrib.auth.urls')),
                       )
