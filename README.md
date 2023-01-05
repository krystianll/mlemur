# mlemur
## Instructions (Windows)
1. Install newest version of R (currently 4.2.2) in the default location from <https://cran.r-project.org/bin/windows/base/>

2. Install newest version of RTools (currently 4.2) in the default location from <https://cran.r-project.org/bin/windows/Rtools/>

3. Proceed to the next part. You can either follow instruction for installing the pre-compiled binary or compiling from the source file.
### Installation using pre-compiled binary
1. Download mlemur zip file from <https://github.com/krystianll/mlemur/blob/master/binaries/Windows/mlemur_0.9.6.zip>

2. Locate and run R.exe

3. In the R main window, install required packages: copy the following command to R console and press Enter:
```
install.packages(c("devtools", "Rcpp", "reactable", "readxl", "writexl", "shiny", "shinyFeedback", "shinyWidgets", "shinyjs", "boot", "BH", "hypergeo", "Ryacas"), dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

4. In the R main window, select "Packages" in the menu bar at the top of the window and then "Install package(s) from local zip files"

5. In the new window, locate the downloaded zip file and click Open
### Installation using source files (requires compilation)
1. Locate and run R.exe

2. Install `devtools` package: copy the following command to R console and press Enter:
```
install.packages("devtools", dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

3. Install mlemur with all dependencies using the following command:
```
devtools::install_github("krystianll/mlemur", ref = "master", auth_token = "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", dependencies = TRUE)
```
### Running mlemur
Mlemur can be run in both text and graphical mode. The graphical mode is initiated in a browser window.

To initiate the graphical mode, use the following command in R:
```
mlemur::mlemur()
```
Alternatively, you can download the following clickable script: <https://github.com/krystianll/mlemur/blob/master/run_mlemur.bat> (You can save the file in any location you choose). If Windows blocks the execution of the file, click "More info" and proceed
## Instructions (macOS)
1. Install newest version of R (currently 4.2.2) in the default location.
- For Intel Macs use the following link: <https://cran.r-project.org/bin/macosx/base/>
- For Apple Silicon (M1, M2) Macs use the following link: <https://cran.r-project.org/bin/macosx/big-sur-arm64/base/>
- If you're not sure which version to choose, click Apple logo in the top-left corner of the screen, choose About This Mac and inspect the Chip model

2. Don't open R yet. Open Terminal.app. To do this, press and hold Command ⌘ key and then press Space key. Type Terminal and press Enter.

3. In the Terminal window, copy the following command and press Enter:
```
xcode-select --install
```
Then follow the instructions in the new window.

4. Proceed to the next part. You can either follow instruction for installing the pre-compiled binary or compiling from the source file.
### Installation using pre-compiled binary
1. Locate and run R.app

2. Download mlemur tgz file corresponding to the verion of R you have installed.
- For Intel Macs use the following link: <https://github.com/krystianll/mlemur/blob/master/binaries/macOS_intel/mlemur_0.9.6.tgz>
- For Apple Silicon (M1, M2) Macs use the following link: <https://github.com/krystianll/mlemur/blob/master/binaries/macOS_arm/mlemur_0.9.6.tgz>
- If you're not sure which version to choose, inspect the text in R main window:
- Platform: aarch64-apple-darwinXX (64-bit) (with XX some number) means you should download the arm version
- Platform: x86_64-apple-darwinXX (64-bit) (with XX some number) means you should download the intel version

**Make sure that the file has .tgz extension! If it is saved as .tar, please change the extension by hand.**

3. In the R main window, install required packages: copy the following command to R console and press Enter:
```
install.packages(c("devtools", "Rcpp", "reactable", "readxl", "writexl", "shiny", "shinyFeedback", "shinyWidgets", "shinyjs", "boot", "BH", "hypergeo", "Ryacas"), dependencies = TRUE)
```
You might be asked to confirm by clicking "Yes" or to select server from which the packages will be downloaded - select the first one

4. In the R main window, select "Packages & Data" in the menu bar at the top of the window and then "Package Installer"

5. In the top-left menu, in the "Packages repository" section change CRAN (binaries) to Local Binary Package and then click Install in the bottom-right part of the window

6. Locate the mlemur_0.9.6.tgz file and press Open

Alternative 4-6. If the above doesn't work, in the R main window, select "Misc" in the menu bar at the top of the window and then "Change Working Directory…". Locate the folder containing mlemur_0.9.6.tgz file and press Open. Then in the R console use the following command and press Enter:
```
install.packages("mlemur_0.9.6.tgz")
```
### Installation using source files (requires compilation)
1. Open Terminal.app. To do this, press and hold Command ⌘ key and then press Space key. Type Terminal and press Enter.

2. In the Terminal window, copy the following command and press Enter to install homebrew:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
You might be prompted to type your password and press Enter. Press Enter again when asked to accept the creation of new folders. Installation might take several minutes.

3. In the Terminal window, copy the following command and press Enter to install boost libraries:
```
brew install boost
```
4. In the Terminal window, copy the following command and press Enter to install mpfr:
```
brew install mpfr
```
5. Locate and run R.app
6. Install `devtools` package: copy the following command to R console and press Enter:
```
install.packages("devtools", dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

7. Install mlemur with all dependencies using the following command:
```
devtools::install_github("krystianll/mlemur", ref = "master", auth_token = "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", dependencies = TRUE)
```
### Running mlemur
Mlemur can be run in both text and graphical mode. The graphical mode is initiated in a browser window.

To initiate the graphical mode, use the following command in R:
```
mlemur::mlemur()
```
Alternatively, you can download the following clickable script: <https://github.com/krystianll/mlemur/blob/master/run_mlemur.command> (You can save the file in any location you choose). Then open the Terminal window and copy the following command:
```
chmod u+x 
```
Make sure there is a space after "x". Drag and drop the run_mlemur.command file to Terminal window and then press Enter. Open the run_mlemur.command file. macOS will still try to block the execution of the file. Click on the Apple logo in the top-left corner of the screen, select System Preferences, then Security & Privacy. Click "Open Anyway" and then confirm by clicking "Open". The file will work from now on.
