#!/bin/bash

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

exec gosu $OCPU_USER apachectl -D FOREGROUND


