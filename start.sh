gunicorn -b 127.0.0.1:8001 --pythonpath dash_app/apps/dash-singlecell-vr app:server --timeout 300
gunicorn -b 127.0.0.1:8000 --pythonpath dash_app/apps/singlecell-vr-api app:server
