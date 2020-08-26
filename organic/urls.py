from django.urls import path
from . import views

urlpatterns = [

    path("", views.index, name="index"),
    path("formula/<str:chemicalFormula>", views.organic_identifier, name="organic_identifier"),
    path("redirect", views.organic_identifier_form, name="organic_identifier_form")
]


