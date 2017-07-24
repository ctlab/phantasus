#!/bin/bash
cd inst/www/phantasus.js/
npm install

R -e "phantasus::servePhantasus('0.0.0.0', 8000, cacheDir = file.path(getwd(), 'jasmine', 'cache'))" &
karma start my.conf.js --single-run

cd

