# Setup

Ensure a C++ compiler is available. On Windows, you need to grab https://www.microsoft.com/en-us/download/details.aspx?id=44266.

Set up by running `fresh-start.bat` or `fresh-start.sh`. The script performs the following actions:

1. Installs modules from requirements.txt with `pip install -r requirements.txt`.
2. Removes any existing database, if it exists.
3. Creates a new set of migrations and runs them.
4. Populates the new database by running the population script `populate.py`.
5. Starts up the server on http://127.0.0.1:8000/.

The population script `populate.py` adds the following users:

| username | password | email address    |
| -------- | -------- | ---------------- |
| jill     | jill     | jill@example.com |
| jen      | jen      | jen@example.com  |
| bob      | bob      | bob@example.com  |

# Debug Mode

By default, DEBUG is enabled in this master branch.

To start the server in release mode run:

`python manage.py runserver --settings=no_debug`

# Tests

Tests can be ran with `./manage.py test` or `python manage.py test`.

# Group members

- David Robertson - [DavidJRobertson](https://github.com/DavidJRobertson): 2133246
- Ewan Morton - [ewanmorton](https://github.com/ewanmorton): 2135980
- Ken Li - [iliawnek](https://github.com/iliawnek): 2131620
- Ross Allan - [nallar](https://github.com/nallar): 2133868
