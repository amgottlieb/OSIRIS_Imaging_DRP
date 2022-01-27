# OSIRIS_Imaging_DRP
This is a data reduction pipeline for OSIRIS imaging. It performs bias subtraction, flat division, bad pixel mask, cosmic ray cleaning, sky subtraction, astrometric correction, and stacking.

Pipeline requirements:
    python 3
    astropy
    astroscrappy
    astroquery
    ccdproc
    matplotlib
    numpy
    photutils
    scipy
    *local file- OSIRIS_imaging_gtcsetup.py


The pipeline can be run from any folder (where the gtcsetup.py is in the same folder) and is used with the following syntax:

python OSIRIS_phot_pipeline.py OBJNAME --workdir  "C:\Users\path\to\raw_files" --outputdir "C:\Users\test\path\to\output" --dobias --doflat --domask  --dowcs --dooverwrite False --filter g,i,z

OBJNAME is a required input and it is the name of the object that is given in the header under the keyword 'OBJECT'.
'workdir' is the full path to the raw files and 'outputdir' is the full path to where you want to store the final files (it does not have to be in the same folder as the raw files).

Or if the master bias, flat, and bad pixel mask files already exist:

python OSIRIS_phot_pipeline.py OBJNAME --workdir "C:\Users\path\to\raw_files" --outputdir "C:\Users\path\to\output" --bias MasterBias --flat MasterFlat --mask MasterBPM --dooverwrite False --filter g,i,z


Output files are:

    MasterBias.fits
    MasterFlatSloan_g.fits (and other filters)
    MasterBPMSloan_g.fits (and other filters)
    final images combined by filter split into ccd1 and ccd2 for both the target and standard star:
        objname_OSIRIS_filter_final_ccd1.fits and objname_OSIRIS_filter_final_ccd2.fits
        Standard_Star_OSIRIS_filter_final_ccd1.fits and Standard_Star_OSIRIS_filter_final_ccd2.fits
    folder containing diagnostic files:
        original images split into ccd1 and ccd2 for both the target and standard star
        CRmask images split into ccd1 and ccd2 for both the target and standard star
