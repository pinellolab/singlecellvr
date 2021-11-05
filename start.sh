#!/usr/bin/env bash

url=$1
frontend_port=$2
backend_port=$3

sed -i "s%API = .*$%API = \'${url}:${backend_port}\'%g" /app/dash_app/apps/dash-singlecell-vr/app.py
sed -i "s%API = .*$%API = \'${url}:${backend_port}\'%g" /app/dash_app/apps/dash-singlecell-vr/static/index_v4.js

bind_url=$(echo $url | sed "s%https://%%")

gunicorn --bind $bind_url:$frontend_port --bind $bind_url:$backend_port --certfile=/localhost.pem --keyfile=/localhost-key.pem -w 2 --pythonpath dash_app/apps/dash-singlecell-vr app:server --timeout 300
