#! /bin/sh

# Cleans up existing installation, then creates fresh DB and starts server

# Installs pip requirements
pip install -r requirements.txt

# Remove any existing DB
rm db.sqlite3

# Make migrations and migrate
python manage.py makemigrations
python manage.py migrate

# Populate fresh DB with test data
python populate.py

# Start local server
python manage.py runserver --insecure
