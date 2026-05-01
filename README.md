# Scripts to reproduce the figures in the CW3E-NCAR SMYLE paper 

The sub-directory ```DATA_SORT``` contains the scripts that can be used to pre-process the data that was used to make the figures.  The ```SIGNIF``` directory in here contains scripts used to calculate significance tests for the figures.  These scripts work on raw model data that would have to be downloaded by the user from the NCAR Geoscientific Data Exchange, but the output from these scripts can be obtained from https://gdex.ucar.edu/datasets/d651080/ (DOI: 10.5065/PPTK-R124) and can then be used to reproduce the figures using the scripts in the ```FIGURES``` sub-directory.

In the following, it is assumed that this package has been downloaded to $DIR

### Installing functions

To install the functions (located on smyle_utils) that are required for most of the jupyter notebooks and python scripts to process the data and generate the figures

```bash
cd $DIR
pip install -e . --user 
```
