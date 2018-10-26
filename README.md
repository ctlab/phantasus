[![Travis-CI Build Status](https://travis-ci.org/ctlab/phantasus.svg?branch=master)](https://travis-ci.org/ctlab/phantasus)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ctlab/phantasus?branch=master&svg=true)](https://ci.appveyor.com/project/ctlab/phantasus)
[![codecov](https://codecov.io/gh/ctlab/phantasus/branch/master/graph/badge.svg)](https://codecov.io/gh/ctlab/phantasus)


This project contains several tools for gene expression analysis with help of R/Bioconductor`

## Installation

The package requires installing `devtools` package and Bioconductor library set up.

```{r}
devtools::install_github("ctlab/phantasus")
```

A warning could appear that the repository contain submodules. This warning 
can be safely ignored.


## Running

In R:

```{r}
library(phantasus)
servePhantasus('0.0.0.0', 8000, cacheDir='cache')
```

Open `http://localhost:8000` in your browser.

## System dependencies

There are several system packages that have to be installed on the system. The
names of these packages will be displayed during installation. On Ubuntu can
install them beforehand and all together using command:

```{bash}
sudo apt-get install libapparmor-dev libprotobuf-dev protobuf-compiler libcurl4-openssl-dev libssl-dev libxml2-dev
```

## Docker 

To simplify deployment phantasus Docker image can be used. It is build regularly and is available at https://hub.docker.com/r/dzenkova/phantasus 
You can run with the following commands:

```{bash}
docker pull dzenkova/phantasus
docker run -t -d -p 80:80 dzenkova/phantasus
```

Phantasus will be available at http://localhost

## Docker compose

Docker compose is available [here](docker-compose.yml) 
You can run with the following command:

```{bash}
docker-compose up -d
```

Phantasus will be available at http://localhost

## ARCHS4 

Phantasus tries to find ARCHS4 .h5 files at cacheDir/archs4
If you are running phantasus using our docker-compose then you can easily download ARCHS4 (~10GB) using with the following commands:

```{bash}
docker-compose run phantasus wget -P /var/phantasus/cache/archs4 https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5
docker-compose run phantasus wget -P /var/phantasus/cache/archs4 https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5
```

Beware that downloading ARCHS4 without volumes makes no sense.
