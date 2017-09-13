FROM ubuntu

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
        locales && \
    apt-add-repository -y ppa:opencpu/opencpu-1.6 && \
    apt-get update && \
    apt-get install -y \
        opencpu-lib


RUN touch /etc/apache2/sites-available/opencpu2.conf
RUN printf "ProxyPass /ocpu/ http://localhost:8001/ocpu/\nProxyPassReverse /ocpu/ http://localhost:8001/ocpu\n" >> /etc/apache2/sites-available/opencpu2.conf
RUN a2ensite opencpu2

RUN sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -
RUN apt-get -y update && apt-get -y install \
    r-base \
    libprotobuf-dev \
    protobuf-compiler \
    r-cran-xml

RUN git clone -b develop --recursive https://github.com/ctlab/phantasus /root/phantasus

RUN R -e 'source("https://bioconductor.org/biocLite.R"); install.packages("XML", repo = "http://cran.gis-lab.info"); biocLite("Biobase"); biocLite("limma"); biocLite("org.Mm.eg.db")'
RUN R -e 'install.packages("devtools", repo = "http://cran.gis-lab.info"); library(devtools); install_github("hadley/scales"); install_github("assaron/GEOquery"); install("/root/phantasus")'

RUN cp -r /root/phantasus/inst/www/phantasus.js /var/www/html/phantasus

RUN a2enmod proxy_http

EXPOSE 80
EXPOSE 443

RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN mkdir -p /var/phantasus/cache
VOLUME ["/var/phantasus/cache"]

RUN mkdir -p /var/phantasus/preloaded
VOLUME ["/var/phantasus/preloaded"]

CMD service apache2 start && \
    R -e 'library(phantasus); servePhantasus("0.0.0.0", 8001, cacheDir="/var/phantasus/cache", preloadedDir="/var/phantasus/preloaded")'
