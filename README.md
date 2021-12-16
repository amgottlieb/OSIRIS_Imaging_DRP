# OSIRIS_Imaging_DRP
This is a data reduction pipeline for OSIRIS imaging. It performs bias subtraction, flat division, bad pixel mask, cosmic ray cleaning, sky subtraction, astrometric correction, and stacking.

Pipeline requirements:

    astropy
    astroscrappy
    astroquery
    ccdproc
    matplotlib
    numpy
    photutils
    scipy
    *local file- OSIRIS_imaging_gtcsetup.py


The ipeline can be run from its folder and is used with the following syntax.

python OSIRIS_phot_pipeline.py --workdir  "C:\Users\path\to\GRB170428A-giz" --outputdir "C:\Users\test\path\to\GRB170428A-giz\output" GRB170428A-giz --dobias --doflat --domask  --dowcs --dooverwrite False --filter g,i,z

where GRB170428A-giz is the name of the object in the header.
Or if the master bias, flat, and bad pixel mask files already exist:

python OSIRIS_phot_pipeline.py --workdir "C:\Users\path\to\GRB170428A-giz" --outputdir "C:\Users\path\to\GRB170428A-giz\output" GRB170428A-giz --bias MasterBias --flat MasterFlat --mask MasterBPM --dooverwrite False --filter g,i,z


Output files are:

    MasterBias.fits
    MasterFlatSloan_g.fits (and other filters)
    MasterBPMSloan_g.fits (and other filters)
    original files split into ccd1 and ccd2
    final images combined by filter split into ccd1 and ccd2 for both the target and standard star:
        objname_OSIRIS_filter_final_ccd1.fits and objname_OSIRIS_filter_final_ccd2.fits
        Standard_Star_OSIRIS_filter_final_ccd1.fits and Standard_Star_OSIRIS_filter_final_ccd2.fits
