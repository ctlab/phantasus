FROM assaron/phantasus-preimage

ARG TARGET_BRANCH=master
ARG PHANTASUS_BUILD
ARG GITHUB_PAT
ENV OCPU_MASTER_HOME=/var/phantasus/ocpu-root

# RUN apt install -y git && git clone -b ${TARGET_BRANCH} --recursive https://github.com/ctlab/phantasus /root/phantasus

COPY . /root/phantasus

RUN R -e 'devtools::install("/root/phantasus", dependencies=TRUE, upgrade=FALSE, build_vignettes=TRUE); remove.packages("BH")'


RUN printf "window.PHANTASUS_BUILD='$PHANTASUS_BUILD';" >> /root/phantasus/inst/www/phantasus.js/RELEASE.js
RUN cp -r /root/phantasus/inst/www/phantasus.js /var/www/html/phantasus

RUN cp /root/phantasus/inst/configs/default /etc/nginx/sites-available/default
RUN mkdir -p /etc/opencpu
RUN cp /root/phantasus/inst/configs/opencpu.conf /etc/opencpu/server.conf
RUN cp /root/phantasus/inst/configs/index.html /var/www/html/
RUN cp -f /root/phantasus/inst/configs/apache_conf/apache_ports.conf /etc/apache2/ports.conf 
RUN cp -f /root/phantasus/inst/configs/apache_conf/000-default.conf /etc/apache2/sites-available/000-default.conf
RUN cp -f /root/phantasus/inst/configs/apache_conf/apache_opencpu.conf /etc/apache2/sites-available/opencpu.conf
RUN cp -f /root/phantasus/inst/configs/apache_conf/apache_mpm_prefork.conf /etc/apache2/mods-available/mpm_prefork.conf
RUN cp -f /root/phantasus/inst/configs/apache_conf/Rprofile /etc/opencpu/Rprofile
RUN rm -rf /root/phantasus/inst

RUN a2dissite default-ssl.conf
RUN a2dissite 000-default.conf
RUN a2dissite rstudio.conf

EXPOSE 8000

RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN mkdir -p /var/phantasus/cache && chown www-data /var/phantasus/cache 
RUN mkdir -p /var/phantasus/preloaded && chown www-data /var/phantasus/preloaded 
RUN mkdir -p /var/phantasus/ocpu-root && chown www-data /var/phantasus/ocpu-root 


ENV OCPU_USER=www-data
CMD /usr/bin/docker-entrypoint.sh

