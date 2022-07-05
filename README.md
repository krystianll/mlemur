# mlemur
## Installation instructions (Windows)
### Requirements
1. Install newest version of R (currently 4.2.1) in the default location from https://cran.r-project.org/bin/windows/base/

2. Install newest version of RTools (currently 4.2.0) in the default location from https://cran.r-project.org/bin/windows/Rtools/

3. Recommended: install newest version of RStudio from https://www.rstudio.com/products/rstudio/download/
### Package installation
1. Locate and run R.exe or RStudio.exe if installed

2. Install `devtools` package: copy the following command to R console and press Enter:
```
install.packages("devtools")
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

3. Install mlemur with all dependencies using the following command:
```
devtools::install_github("krystianll/mlemur", ref = "master", auth_token = "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", dependencies = TRUE)
```
You might be asked if you want to update the already installed packages. Type the number corresponding to the option `None` (usually `3`) and press Enter
## Running mlemur (Windows)
Mlemur can be run in both text and graphics mode. The graphics mode is initiated in a browser window.

To initiate the graphics mode, use the following command in R:
```
mlemur::mlemur()
```
Alternatively, you can download the following clickable script: https://github.com/krystianll/mlemur/blob/master/run_mlemur.bat (Please right-click on the link and click "Save as")
