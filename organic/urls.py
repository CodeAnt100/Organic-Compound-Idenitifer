from django.urls import path
from . import views

urlpatterns = [

    path("", views.index, name="index"),
    path("<str:chemicalFormula>", views.organic_identifier, name="organic_identifier")
]


