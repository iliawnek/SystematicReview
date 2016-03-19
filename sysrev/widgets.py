from django.forms import widgets
from django.utils.safestring import *

# See static/js/querywidget.js
class QueryWidget(widgets.Textarea):
    def __init__(self, attrs=None):
        default_attrs = {'class': "queryWidget"}
        if attrs:
            default_attrs.update(attrs)
        super(QueryWidget, self).__init__(default_attrs)

    def script(self):
        return mark_safe("<noscript>Please enable javascript to use the query editor. Without it enabled, you can only "
                         "use the advanced editor</noscript>"
                         "<script src='/static/js/querywidget.js' type='text/javascript' defer></script>")

    def render(self, name, value, attrs=None):
        textAreaHtml = super(QueryWidget, self).render(name, value, attrs)
        return self.script() + textAreaHtml
