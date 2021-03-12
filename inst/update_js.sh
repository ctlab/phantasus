#!/bin/bash

PHANTASUS_JS_DIR="$1"

if [ -z "${PHANTASUS_JS_DIR}" ]; then
    echo "Please specify git dir with phantasus.js"
    exit 1
    # PHANTASUS_JS_DIR=inst-raw/phantasus.js
fi

PHANTASUS_JS_DIR=$(readlink -f "${PHANTASUS_JS_DIR}")

echo "Using phantasus.js from ${PHANTASUS_JS_DIR}"

grunt --base "${PHANTASUS_JS_DIR}" dist

git --work-tree="${PHANTASUS_JS_DIR}" log -1 > ./inst/www/phantasus.js/commit

cp ${PHANTASUS_JS_DIR}/index.html ./inst/www/phantasus.js/
cp ${PHANTASUS_JS_DIR}/js/phantasus.js ./inst/www/phantasus.js/js/
cp ${PHANTASUS_JS_DIR}/js/phantasus-external-*.min.js ./inst/www/phantasus.js/js/
cp -r ${PHANTASUS_JS_DIR}/jasmine ./inst/www/phantasus.js/
cp ${PHANTASUS_JS_DIR}/my.conf.js ./inst/www/phantasus.js/
cp ${PHANTASUS_JS_DIR}/package.json ./inst/www/phantasus.js/
cp -r ${PHANTASUS_JS_DIR}/css/images ./inst/www/phantasus.js/css/
cp -r ${PHANTASUS_JS_DIR}/css/phantasus-latest.min.css ./inst/www/phantasus.js/css/
