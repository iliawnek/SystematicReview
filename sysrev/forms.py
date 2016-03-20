from django import forms
from django.contrib.auth.models import User
from registration.forms import RegistrationForm
from sysrev.models import *
from sysrev.api import PubMed
from widgets import *

class ProfileForm(forms.ModelForm):
    class Meta:
        model  = User
        fields = ("email",)

    def clean_email(self):
        email = self.cleaned_data.get('email')
        if self.instance and self.instance.email == email:
            return email
        if User.objects.filter(email=email).count():
            raise forms.ValidationError('That email address is already taken.')
        return email




class ReviewCreateStep1(forms.Form):
    title       = forms.CharField(max_length=128, label="Review Title")
    description = forms.CharField(widget=forms.Textarea, required=False)
    invited     = forms.CharField(widget=forms.Textarea, required=False, label="Invite Participants", help_text="Enter the email address or username of each participant, one per line")

    def clean_invited(self):
        invited = self.cleaned_data.get('invited')
        invlist = filter(lambda i: i, map(lambda l: str.strip(str(l)), invited.splitlines()))

        for invitee in invlist:
            try:
                if invitee.find("@") == -1:
                    User.objects.get(username=invitee)
                else:
                    User.objects.get(email=invitee)
            except User.DoesNotExist:
                raise forms.ValidationError('User '+invitee+' not found.')

        return invited


class ReviewCreateStep2(forms.Form):
    query       = forms.CharField(widget=QueryWidget)

    def clean_query(self):
        query = self.cleaned_data.get('query')
        data = PubMed.get_data_from_query(query)
        count = data["Count"]
        if count >= 1000:
            raise forms.ValidationError("""Your query returned %s papers.\n
                                        It must return fewer than 1000 papers.\n
                                        Modify your query and try again.""" % str(count))
        elif count == 0:
            raise forms.ValidationError("Your query did not return any papers.")

        return query
