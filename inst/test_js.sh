#!/bin/bash
R -e "devtools::install('.')"

cd inst/www/phantasus.js/
npm install karma --save-dev
npm install

mkdir "jasmine/cache"
R -e "phantasus::getES('GSE53986', destdir='jasmine/cache')"

R -e "phantasus::servePhantasus('0.0.0.0', 8000, cacheDir = 'jasmine/cache', preloadedDir = 'jasmine/cache', openInBrowser=FALSE)" &
PH_PID=$!

sleep 2

./node_modules/karma/bin/karma start my.conf.js --single-run
RETVAL=$?

cd

function finish {
    kill $PH_PID
    echo "Killed OpenCPU-server"
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



