pip install -r requirements.txt
del db.sqlite3
python manage.py makemigrations
python manage.py migrate
python populate.py
python manage.py runserver