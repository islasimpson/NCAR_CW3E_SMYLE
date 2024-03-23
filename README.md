# Scripts to reproduce the figures in Deflorio et al (XXXX)

The sub-directory DATA_SORT contains the scripts that can be used to pre-process the data that was used to make the figures.  These scrips work on raw model data that would have to be downloaded by the user, but the output from these scripts can be obtained from XXXXX.

In the following, it is assumed that this package has been downloaded to $DIR

### Installing functions

To install the functions (located on smyle_utils) that are required for most of the jupyter notebooks and python scripts to process the data and generate the figures

```bash
cd $DIR
pip install -e . --user 
```

### Data pre-processing

All the fields used for the analysis have been pre-processed and are provided in netcdf format at XXXXX.  These processed fields are derived from the raw model or observation-based output using the scripts that are located in $DIR/DATA_SORT


