pip install -r requirements.txt
@echo off
( type nul >> db.sqlite3 ) 2>nul || (echo The database is currently in use by another process. && exit /b)
@echo on
del db.sqlite3
python manage.py makemigrations
python manage.py migrate
python populate.py
python manage.py runserver --insecure