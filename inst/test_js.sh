#!/bin/bash
cd inst/www/phantasus.js/
npm install

R -e "phantasus::servePhantasus('0.0.0.0', 8000, cacheDir = file.path(getwd(), 'jasmine', 'cache'))" &
karma start my.conf.js --single-run
RETVAL=$?

cd

if [ $RETVAL -eq 0 ]
then
    echo "Jasmine tests are successfully passed"
else
    echo "Tests failed" >&2
fi
