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


Pipeline can be run from its folder and is used with the following syntax.

python OSIRIS_phot_pipeline.py --workdir  "C:\Users\test\datafolder" --outputdir "C:\Users\test\datafolder\output" object_name --dobias --doflat --domask --dowcs --dooverwrite False --filter g,i,z

Or if files already exist:

python OSIRIS_phot_pipeline.py --workdir "C:\Users\test\data_folder" --outputdir "C:\Users\test\data_folder\output" object_name --bias MasterBias --flat MasterFlat --mask MasterBPM --dooverwrite False --filter g,i,z


Output files are:

    MasterBias.fits
    MasterFlatSloan_g.fits (and other filters)
    MasterBPMSloan_g.fits (and other filters)
    original files split into ccd1 and ccd2
    final images combined by filter split into ccd1 and ccd2 for both the target and standard star:
        objname_OSIRIS_filter_final_ccd1.fits and objname_OSIRIS_filter_final_ccd2.fits
        Standard_Star_OSIRIS_filter_final_ccd1.fits and Standard_Star_OSIRIS_filter_final_ccd2.fits
