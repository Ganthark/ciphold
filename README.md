# ciphold
Adaptation of the SCAFFoLD cytometry tool.
The original project can be found on the [SCAFFoLD GitHub page](https://github.com/nolanlab/scaffold).
>Please note that this is the short version of the Readme. Please refer to the USER_MANUAL.html file for more details (download and open it with any Web Browser).

## Requirements
  * software: R(Version 3.4.3 to 3.5), Rstudio(optional)
  * R packages: Rcpp, cluster, flowCore, ggplot2, igraph, plyr, reshape, shiny, scales, grDevices, parallel, jsonlite, colourpicker, shinythemes, tcltk *(only on linux)*.

## Quick installation guide

  1. Run the following command in R/RStudio:
```
install.packages(c("Rcpp", "cluster", "ggplot2", "igraph", "reshape", "shiny", "jsonlite", "colourpicker", "shinythemes", "tcltk"))
```
  >You may be asked to reload your environment, if so, accept.
  
  2. Run the next commands:
```
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
```
  
  3. Download and extract the App in the desired folder from [GitHub](https://github.com/cipheLab/).
  
## Launching the shiny application

  1. Double click on the ciphold.Rproj file to open it. It will run a new R session with the working directory set to launch the shiny application.
  2. Run the following commands in R/RStudio:
```
library("shiny")
runApp("App/")
```
  3. The application should be running, if not, do all the steps from the "Setting up your environment" section in the good order. If it does not resolves the issue, please go to the "Known Issues" section.
  
### Notes: 
> You could also use **shiny::runApp("App/")**.
