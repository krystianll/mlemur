# mlemur
## Instructions (Windows)
### Installation using pre-compiled binary
1. Install newest version of R (currently 4.2.1) in the default location from https://cran.r-project.org/bin/windows/base/

2. Download mlemur zip file from https://github.com/krystianll/mlemur/blob/master/binaries/Windows/mlemur_0.9.zip

3. Locate and run R.exe

4. In the R main window, install required packages: copy the following command to R console and press Enter:
```
install.packages(c("devtools", "Rcpp", "reactable", "readxl", "writexl", "shiny", "shinyFeedback", "shinyWidgets", "shinyjs", "stats", "boot", ""BH, "hypergeo", "Ryacas"), dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

5. In the R main window, select Packages in the menu bar at the top of the window and then "Install package(s) from local zip files"

6. Locate the downloaded zip file and click Open

### Installation using source files (requires compilation)
1. Install newest version of R (currently 4.2.1) in the default location from https://cran.r-project.org/bin/windows/base/

2. Install newest version of RTools (currently 4.2.0) in the default location from https://cran.r-project.org/bin/windows/Rtools/

3. Install `devtools` package: copy the following command to R console and press Enter:
```
install.packages("devtools", dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

4. Install mlemur with all dependencies using the following command:
```
devtools::install_github("krystianll/mlemur", ref = "master", auth_token = "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", dependencies = TRUE)
```
### Running mlemur
Mlemur can be run in both text and graphical mode. The graphical mode is initiated in a browser window.

To initiate the graphical mode, use the following command in R:
```
mlemur::mlemur()
```
Alternatively, you can download the following clickable script: https://github.com/krystianll/mlemur/blob/master/run_mlemur.bat (Please right-click on the link and click "Save as". You can save the file in any location you choose). If Windows blocks the execution of the file, click "More info" and proceed
