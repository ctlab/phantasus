#!/bin/bash
R --vanilla -e "devtools::install('.')"

mkdir -p inst/www/phantasus.js/jasmine/cache
cp "./inst/testdata/GSE27112-GPL6103.rda" inst/www/phantasus.js/jasmine/cache
cd inst/www/phantasus.js/
npm install karma --save-dev
npm install

mkdir -p jasmine/preloaded
R --vanilla -e "phantasus::getES('GSE53986', destdir = 'jasmine/preloaded')" 

touch server.log
R --vanilla -e "phantasus::servePhantasus('0.0.0.0', 8000, cacheDir = 'jasmine/cache', preloadedDir = 'jasmine/preloaded', openInBrowser=FALSE)"  > server.log 2>&1 &
PH_PID=$!


# Waiting for the server to start
while ! grep -q started server.log ; do sleep 0.1; done

sleep 0.1

curl http://localhost:8000/ocpu/library/phantasus/R/loadPreloaded -d "name='GSE53986'" | tee curl.log
curl http://localhost:8000`grep stdout curl.log`/text

./node_modules/karma/bin/karma start my.conf.js --single-run
RETVAL=$?

cd

function finish {
    kill $PH_PID
    echo Killed $PH_PID
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



