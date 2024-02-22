#!/bin/bash

_term() {
  trap ERR
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


color() {
    STARTCOLOR="\e[38;5;$2";
    ENDCOLOR="\e[0m";
    export "$1"="$STARTCOLOR%s"
}
color darkred 196m
color pinkish 211m
color blue 27m
color lightblue 68m
color coral 203m
color lessblue 26m
color normal 0m

printf $darkred "#### "
printf $pinkish "++++ "
printf $blue'\n' "----"
printf $darkred "#### "
printf $pinkish "++++ "
printf $blue'\n' "----"
printf $darkred "#### "
printf $pinkish "++++ "
printf $blue'\n' "----"

printf $lightblue "====      "
printf $coral "&&&&"
printf "\e[0m%b"'\n' "      _____    _                       _                                   _                         "
printf $lightblue "====      "
printf $coral "&&&&"
printf "\e[0m%s\n" "     |  __ \  | |                     | |                                 (_)                        "
printf $lightblue   "====      "
printf $coral"&&&&     "
printf "\e[0m%s\n" "| |__) | | |__     __ _   _ __   | |_    __ _   ___   _   _   ___     _   ___     _   _   _ __  "
printf $coral "&&&& "
printf $darkred "#### "
printf $lessblue "----     "
printf "\e[0m%s\n" "|  ___/  | '_ \   / _\` | | '_ \  | __|  / _\` | / __| | | | | / __|   | | / __|   | | | | | '_ \ "
printf $coral "&&&& "
printf $darkred "#### "
printf $lessblue "----     "
printf "\e[0m%s\n" "| |      | | | | | (_| | | | | | | |_  | (_| | \__ \ | |_| | \__ \   | | \__ \   | |_| | | |_) |"
printf $coral "&&&& "
printf $darkred "#### "
printf $lessblue "----     "
printf "\e[0m%s\n" "|_|      |_| |_|  \__,_| |_| |_|  \__|  \__,_| |___/  \__,_| |___/   |_| |___/    \__,_| | .__/ "
printf $blue "----"
printf "\e[0m%s\n" "                                                                                                        | |    "
printf $blue "----"
printf "\e[0m%s\n" "                                                                                                        |_|    "
printf $blue "----"
printf "\e[0m%s\n" ""


tail -f /var/log/nginx/access.log
