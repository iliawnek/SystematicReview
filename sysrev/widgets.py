
from django.forms import widgets
from django.utils.safestring import *

class QueryWidget(widgets.Textarea):
    def render(self, name, value, attrs=None):
        textAreaHtml = super(QueryWidget, self).render(name, value, attrs)
        return mark_safe("I am the query editor widget, with some custom behaviour.") + textAreaHtml