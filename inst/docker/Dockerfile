FROM ubuntu:16.04

ARG TARGET_BRANCH=master
ENV OCPU_MASTER_HOME=/var/phantasus/ocpu-root

RUN apt-get -y update && \
    apt-get -y dist-upgrade && \
    apt-get -y install \
        software-properties-common \
        git \
        libcairo2-dev \
        libxt-dev \
        libssl-dev \
        libssh2-1-dev \
        libcurl4-openssl-dev \
        apache2 \
        libapparmor-dev \
        libxml2-dev \
        locales

RUN apt-add-repository -y 'deb http://ppa.launchpad.net/maarten-fonville/protobuf/ubuntu zesty main' && \
    apt-get -y update

RUN apt-get -y --allow-unauthenticated install \
    libprotobuf-dev \
    protobuf-compiler


#RUN apt-add-repository -y ppa:marutter/rrutter3.5 && \
#    apt-get update

RUN apt-add-repository -y ppa:marutter/rrutter && \
    apt-get update

RUN apt-get install -y r-base

RUN R -e 'install.packages("protolite", repo = "https://cran.rstudio.com/")'

RUN apt-add-repository -y ppa:opencpu/opencpu-2.0 && \
    apt-get update && \
    apt-get install -y --allow-unauthenticated \
        opencpu-lib


#RUN R -e 'install.packages("opencpu", repo = "https://cran.rstudio.com/")'
#
#
RUN touch /etc/apache2/sites-available/opencpu2.conf
RUN printf "ProxyPass /ocpu/ http://localhost:8001/ocpu/\nProxyPassReverse /ocpu/ http://localhost:8001/ocpu\n" >> /etc/apache2/sites-available/opencpu2.conf
RUN a2ensite opencpu2

#RUN sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list'
# protobuf 3.5, with 2GB byte limit
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -
    
RUN R -e 'install.packages("devtools", repo = "https://cran.rstudio.com/")'


RUN git clone -b ${TARGET_BRANCH} --recursive https://github.com/ctlab/phantasus /root/phantasus

RUN R -e 'source("https://bioconductor.org/biocLite.R")'
RUN R -e 'devtools::install("/root/phantasus", build_vignettes=T)'

RUN cp -r /root/phantasus/inst/www/phantasus.js /var/www/html/phantasus

RUN a2enmod proxy_http

EXPOSE 80

RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN mkdir -p /var/phantasus/cache
RUN mkdir -p /var/phantasus/preloaded
RUN mkdir -p /var/phantasus/ocpu-root

CMD service apache2 start && \
   R -e 'library(phantasus); servePhantasus("0.0.0.0", 8001, cacheDir="/var/phantasus/cache", preloadedDir="/var/phantasus/preloaded")'

