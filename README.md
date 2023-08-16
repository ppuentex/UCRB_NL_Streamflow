## Understanding the Temporal Variability and Predictability of Streamflow Signatures in the Colorado River Basin

### Description
This repository contains the nonlinear analysis performed on streamflow time series in the Colorado River Basin and all output figures for manuscript. This workflow is intended for RStudio. The jupyter notebook included is for data processing purpose only. Use if intention is to give an updated Recorded Streamflow from [Bureau of Reclamation](https://www.usbr.gov/lc/region/g4000/NaturalFlow/current.html). 

### Getting Started 
1. If not installed, install RStudio by following these instructions [RStudio Installation](https://rstudio-education.github.io/hopr/starting.html)

2. Clone repo by pasting the following in the terminal
    ```git clone https://github.com/ppuentex/UCRB_NL_Streamflow.git```
    
3. Open RStudio, go to `File` > `Open Project...` navigate to cloned directory and choose `UCRB_NL_Streamflow.Rproj`
    - Now all project files are located within RStudio and should show up in the bottom right panel. 
    
4. Package installations - navigate to the console in RStudio and install the following packages before running any scripts 
    - Run `install.packages("[name of package]")` replace `[name of package]` with 
  ` zoo `,` fractal `, `quantmod`, `scatterplot3d` 
    - For example `install.packages("zoo")` , `install.packages("fractal")`, `install.packages("quadmod")`, etc 

  
### Folder Breakdown
Folders have a separate readme to further explain the content. 

#### Input_Data
These are the streamflow time series files from the [Bureau of Reclamation](https://www.usbr.gov/lc/region/g4000/NaturalFlow/current.html) and [Treeflow](https://www.treeflow.info/upper-colorado-basin). 

#### Output_Data 
These are the output analysis time series such as the streamflow signal, average Local Lyapunov Exponent (LLE), rolling average, and rolling variance. Including simple statistics for each gauge and time delay, dimension embedding parameters, and global average LLE for each gauge. 
