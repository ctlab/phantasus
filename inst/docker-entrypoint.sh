#!/bin/bash

_term() {
  echo "Caught signal!"
  nginx -s stop
  apachectl -k graceful-stop
}

trap _term SIGINT SIGTERM SIGWINCH


chown -R $OCPU_USER /var/log/apache2
chown -R $OCPU_USER /var/run/apache2
chown -R $OCPU_USER /var/log/opencpu
mkdir -p /run/apache2
chown $OCPU_USER /run/apache2

chown -R $OCPU_USER /var/log/nginx
chown -R $OCPU_USER /var/lib/nginx

touch /run/nginx.pid
chown $OCPU_USER /run/nginx.pid

gosu $OCPU_USER nginx

gosu $OCPU_USER apachectl start

tail -f /var/log/nginx/access.log
