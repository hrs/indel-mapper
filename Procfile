web: gunicorn -w 4 -b "0.0.0.0:$PORT" app:app
worker: celery worker -A app.celery --loglevel=info
