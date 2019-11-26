FROM assaron/phantasus-preimage

ARG TARGET_BRANCH=master
ARG PHANTASUS_BUILD
ARG GITHUB_PAT
ENV OCPU_MASTER_HOME=/var/phantasus/ocpu-root

#RUN git clone -b ${TARGET_BRANCH} --recursive https://github.com/ctlab/phantasus /root/phantasus
COPY . /root/phantasus

RUN R -e 'devtools::install("/root/phantasus", dependencies=TRUE, upgrade=FALSE, build_vignettes=TRUE)'
RUN printf "window.PHANTASUS_BUILD='$PHANTASUS_BUILD';" >> /root/phantasus/inst/www/phantasus.js/RELEASE.js
RUN cp -r /root/phantasus/inst/www/phantasus.js /var/www/html/phantasus
RUN cp /root/phantasus/inst/configs/default /etc/nginx/sites-available/default
RUN mkdir /etc/opencpu
RUN cp /root/phantasus/inst/configs/opencpu.conf /etc/opencpu/server.conf
RUN cp /root/phantasus/inst/configs/index.html /var/www/html/
RUN rm -rf /root/phantasus/inst

EXPOSE 80

RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN mkdir -p /var/phantasus/cache
RUN mkdir -p /var/phantasus/preloaded
RUN mkdir -p /var/phantasus/ocpu-root

CMD service nginx start && \
   R -e 'library(phantasus); servePhantasus("0.0.0.0", 8001, openInBrowser = F, cacheDir="/var/phantasus/cache", preloadedDir="/var/phantasus/preloaded")'


