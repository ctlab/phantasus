#!/bin/bash
R -e "devtools::install('.')"

cd inst/www/phantasus.js/
npm install

R -e "phantasus::servePhantasus('0.0.0.0', 8000, cacheDir = file.path(getwd(), 'jasmine', 'cache'))" &
PH_PID=$!

karma start my.conf.js --single-run
RETVAL=$?

cd

function finish {
    kill $PH_PID
}
trap finish EXIT

if [ $RETVAL -eq 0 ]
then
    echo "Jasmine tests are successfully passed"
    exit 0
else
    echo "Tests failed" >&2
    exit 1
fi



