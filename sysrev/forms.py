from django import forms
from django.contrib.auth.models import User
from registration.forms import RegistrationForm

class ManageAccount(forms.ModelForm):
    email = forms.EmailField(required=True)

    class Meta:
        model = User
        fields = ('email',)

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop('user', None)
        super(ManageAccount, self).__init__(*args, **kwargs)

    def clean_email(self):
        email = self.cleaned_data.get('email')
        if self.user and self.user.email == email:
            return email
        if User.objects.filter(email=email).count():
            raise forms.ValidationError('That email address is already taken.')
        return email

    def save(self, commit=True):
        self.user.email = self.cleaned_data['email']

        if commit:
            self.user.save()

        return self.user