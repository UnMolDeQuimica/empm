from django.shortcuts import render
from django.views.generic import TemplateView, View


class Home(TemplateView):
    template_name = "index.html"


class HowDoesItWork(TemplateView):
    template_name = "how_does_it_work.html"


class CookiesPolicy(TemplateView):
    template_name = "cookies_policy.html"

