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

If you don't have python installed already, install it for your os with Anaconda https://www.anaconda.com/products/individual

Store the pipeline and setup file in your favorite folder and add the folder to your PATH and PYTHON_PATH environment. 
On LINUX/mac in a terminal:

    export PATH=$PATH:/new/path
    export PYTHON_PATH=$PYTHON_PATH:/new/path

Then make the file an executable:  
    
    chmod u+x OSIRIS_imaging_pipeline.py

On WINDOWS, follow the instructions here: https://correlated.kayako.com/article/40-running-python-scripts-from-anywhere-under-windows 

Now the pipeline can be run from any folder (where the gtcsetup.py is in the same folder). 
You can run the pipeline from any folder using the full path to the data with the following syntax:

    python OSIRIS_phot_pipeline.py OBJNAME --workdir  "C:\Users\path\to\raw_files" --outputdir "C:\Users\test\path\to\output" --dobias --doflat --domask  --dowcs --dooverwrite False --filter g,i,z

OR you can move into the folder where your data is stored and then --workdir and --outputdir are not needed:

    python OSIRIS_phot_pipeline.py OBJNAME --dobias --doflat --domask  --dowcs --dooverwrite False --filter g,i,z

OBJNAME is a required input and it is the name of the object that is given in the header under the keyword 'OBJECT'.
'workdir' is the full path to the raw files and 'outputdir' is the full path to where you want to store the final files (it does not have to be in the same folder as the raw files).

Or if the master bias, flat, and bad pixel mask files already exist, give the names of their files with the following syntax (default file names are shown):

    python OSIRIS_phot_pipeline.py OBJNAME --bias MasterBias --flat MasterFlat --mask MasterBPM --dooverwrite False --filter g,i,z


The output files are:

    Calibration files:
        MasterBias.fits
        MasterFlatSloan_g.fits (and other filters)
        MasterBPMSloan_g.fits (and other filters)
        skymap for each filter split into ccd1 and ccd2 for both the target and standard star 
        stacked image after aligning with each other for each filter split into ccd1 and ccd2 for both the target and standard star 
    
    Final images combined by filter split into ccd1 and ccd2 for both the target and standard star:
        objname_OSIRIS_filter_final_ccd1.fits and objname_OSIRIS_filter_final_ccd2.fits
        Standard_Star_OSIRIS_filter_final_ccd1.fits and Standard_Star_OSIRIS_filter_final_ccd2.fits
    
    Folder containing diagnostic files:
        original images split into ccd1 and ccd2 for both the target and standard star
        CRmask images split into ccd1 and ccd2 for both the target and standard star
        images after aligning with eachother split into ccd1 and ccd2 for both the target and standard star
        
See *** for examples of a good master bias, flat, sky map, and final image.
