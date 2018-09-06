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


To install `phantasus` on R 3.4, use r-3.4 branch:

```{r}
devtools::install_github("ctlab/phantasus", ref="r-3.4")
```

## Running

In R:

```{r}
library(phantasus)
servePhantasus('0.0.0.0', 8000, cacheDir='cache')
```

Open `http://localhost:8000` in your browser.

## System dependencies

There are several system packages that have to be installed on the system. The names of these packages will be displayed during installation. On Ubuntu can install them beforehand and alltogether using command:

```{bash}
sudo apt-get install libapparmor-dev libprotobuf-dev protobuf-compiler libcurl4-openssl-dev libssl-dev libxml2-dev
```

## Installation example (CentOS 7)

Update packages and install `R` and dependencies:
```
yum update
yum install openssl-devel protobuf-compiler R R-Rcpp R-Rcpp-devel libcurl-devel libxml2-devel protobuf-devel git screen
```

Install `devtools` and `bioclite`:
```
R -e 'install.packages("devtools", repos="https://cloud.r-project.org")'
R -e 'if (!requireNamespace("BiocManager", quietly=TRUE)); BiocManager::install()'
     'install.packages("BiocManager"); BiocManager::install()'
```

Install `phantasus`:
```
R -e 'devtools::install_github("ctlab/phantasus", ref="r-3.4")'
```

Run `phantasus` at 80 port in screen:
```
screen -S phantasus-server
R -e 'library(phantasus); servePhantasus("0.0.0.0", 80, cacheDir="/var/phantasus/cache")'
```

Check `http://0.0.0.0`
