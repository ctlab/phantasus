#!/bin/bash
PHANTASUS_VERSION=$(grep -r 'Version:' ./DESCRIPTION | awk {'print $2'})
PHANTASUS_BUILD=$(git rev-parse HEAD)
printf "window.PHANTASUS_VERSION='$PHANTASUS_VERSION';\nwindow.PHANTASUS_BUILD='none-${PHANTASUS_BUILD::10}';\n" | tee inst/www/phantasus.js/RELEASE.js


