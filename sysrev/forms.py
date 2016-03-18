from django import forms
from django.contrib.auth.models import User
from registration.forms import RegistrationForm
from sysrev.models import *

class ProfileForm(forms.ModelForm):
    class Meta:
        model  = User
        fields = ('email',)

    def clean_email(self):
        email = self.cleaned_data.get('email')
        if self.instance and self.instance.email == email:
            return email
        if User.objects.filter(email=email).count():
            raise forms.ValidationError('That email address is already taken.')
        return email




class ManageReview(forms.ModelForm):
    class Meta:
        model  = Review
