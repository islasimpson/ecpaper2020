=================================================================================================
=== Package to perform the analysis described in Simpson et al (2010) Emergent constraints on ===
=== the large scale atmospheric circulation and regional hydroclimate: do they still work in  ===
=== CMIP6 and how much can they actually constrain the future?, submitted to J. Clim          ===
================================================================================================= 

In the following, it is assumed that this package has been downloaded to $DIR.

To install the sub-routines (located in ecpaper_utils) required for running the following scripts, do

cd $DIR
pip install -e . --user

============= DATA pre-processing ================

All the fields used for the analysis have been pre-processed and are provided in netcdf format.

The scripts used to process the data are found in $DIR/DATASORT 

Three sub-directories contain processing scripts for each of the emergent constraints considered

$DIR/DATASORT/JLAT
$DIR/DATASORT/VWIND
$DIR/DATASORT/CALP

CMIP5, CMIP6 and observational data are read and processed from a local archive in the Climate and Global Dynamics Laboratory at NCAR.  If a user wishes to start the data processing from scratch, they will have to change the file paths to point to the relevant location in their archives.

The multi-model large ensemble data were read and processed from the USClivar Multi-Model Large Ensembles archived on the cheyenne supercomputer at /glade/collections/cdg/data/CLIVAR-LE or availabale through the NCAR Climate Data Gateway (see cesm.ucar.edu/projects/community-projects/MMLEA/).

Processed data are provided within sub-directories for each emergent constraint and are used in the codes provided in $DIR/ANALYSIS





