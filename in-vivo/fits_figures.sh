#!/usr/bin/env bash

# Fit models
jupyter nbconvert --execute --inplace --allow-errors fits/CTCF*.ipynb

# Draw figures
jupyter nbconvert --execute --inplace --allow-errors fits/FigureS4.ipynb