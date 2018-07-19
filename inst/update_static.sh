#!/bin/bash

R -e "rmarkdown::render('./vignettes/phantasus-tutorial.Rmd')"
cp ./vignettes/phantasus-tutorial.html inst/www/phantasus.js/

cd inst-raw/phantasus.js

grunt dist

git log -1 > ../../inst/www/phantasus.js/commit

cp ./index.html ../../inst/www/phantasus.js/
cp ./js/phantasus.js ../../inst/www/phantasus.js/js/
cp ./js/phantasus-external-*.min.js ../../inst/www/phantasus.js/js/
cp -r ./jasmine ../../inst/www/phantasus.js/
cp ./my.conf.js ./package.json ../../inst/www/phantasus.js/
cp -r ./css/images ../../inst/www/phantasus.js/css/
cp -r ./css/phantasus-latest.min.css ../../inst/www/phantasus.js/css/
