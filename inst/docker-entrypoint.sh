#!/bin/bash

_term() {
  echo "Caught signal!"
  nginx -s stop
  apachectl -k graceful-stop
}

_term_config(){
    echo "Failed to complete setup! Bad configuration file or user: $OCPU_USER doesn't have permissions to write cache or configuration file."
    exit 1
}

trap _term SIGINT SIGTERM SIGWINCH

trap _term_config ERR
set -eE

printf '%s\n' "#### ++++ ----"\
              "#### ++++ ----"\
              "#### ++++ ----"\
              "====      &&&&      _____    _                       _"\
              "====      &&&&     |  __ \  | |                     | |"\
              "====      &&&&     | |__) | | |__     __ _   _ __   | |_    __ _   ___   _   _   ___"\
              "&&&& #### ----     |  ___/  | '_ \   / _\` | | '_ \  | __|  / _\` | / __| | | | | / __|"\
              "&&&& #### ----     | |      | | | | | (_| | | | | | | |_  | (_| | \__ \ | |_| | \__ \ "\
              "&&&& #### ----     |_|      |_| |_|  \__,_| |_| |_|  \__|  \__,_| |___/  \__,_| |___/ "\
              "----"\
              "----"\
              "----"


# fixes permission error on old docker versions (workaround from https://github.com/moby/moby/issues/6047#issuecomment-68608697)
mv /var/log/apache2 /var/log/apache2.bak && mv /var/log/apache2.bak /var/log/apache2

chown -R $OCPU_USER /var/log/apache2
chown -R $OCPU_USER /var/run/apache2
chown -R $OCPU_USER /var/log/opencpu
chown -R $OCPU_USER /var/phantasus/ocpu-root

mkdir -p /run/apache2
chown $OCPU_USER /run/apache2

chown -R $OCPU_USER /var/log/nginx
chown -R $OCPU_USER /var/lib/nginx
chown -R $OCPU_USER /var/cache/nginx

touch /run/nginx.pid
chown $OCPU_USER /run/nginx.pid

# Let's not change user files ownership?
# chown -R $OCPU_USER $R_USER_CONFIG_DIR/R/phantasus

gosu $OCPU_USER R -e "phantasus:::createDockerConf(); phantasus::setupPhantasus()" #|| _term_config

gosu $OCPU_USER nginx

gosu $OCPU_USER apachectl start

printf '%s\n'  "■ ▦ ▥   __"\
               "▤   ▩  |__) |_   _   _  |_  _    _      _" \
               "▩ ■ ▥  |    | ) (_| | ) |_ (_| _) |_| _)     is up"\
               "▥ "

tail -f /var/log/nginx/access.log
