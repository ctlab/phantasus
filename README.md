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
and can be run with the following commands:

```{bash}
docker pull dzenkova/phantasus
docker run -t -d dzenkova/phantasus
```
