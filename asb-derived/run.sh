#!/usr/bin/env bash

# Create tables
python3 table.py

# Fit models
jupyter nbconvert --execute --inplace --allow-errors CTCF_ASB.ipynb

# Draw figures
jupyter nbconvert --execute --inplace --allow-errors Figure*.ipynb