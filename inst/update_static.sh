#!/bin/bash

cd ../inst-raw

git submodule init && git submodule update

cp phantasus.js/js/phantasus-latest.min.js phantasus.js/js/phantasus-external-latest.min.js ../inst/www/phantasus.js/js/
cp -r phantasus.js/jasmine ../inst/www/phantasus.js/
cp phantasus.js/my.conf.js phantasus.js/package.json ../inst/www/phantasus.js/
cp -r phantasus.js/css ../inst/www/phantasus.js
