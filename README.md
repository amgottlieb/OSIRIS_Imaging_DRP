# OSIRIS_Imaging_DRP
This is a data reduction pipeline for OSIRIS imaging. It performs bias subtraction, flat division, bad pixel mask, cosmic ray cleaning, sky subtraction, astrometric correction, and stacking.

Pipeline requirements (install via 'pip install pkgname' unless otherwise noted):

    python 3 - If you don't have Python 3 installed already, install it for your os with Anaconda https://www.anaconda.com/products/individual
    
    astroalign
    astropy
    astroscrappy    conda install -c conda-forge astroscrappy
    astroquery
    ccdproc
    matplotlib
    numpy
    photutils
    scipy
    sep
    
NOTE: These two files come with the download on Github and must be in the same directory as the main script:   

    *local file- OSIRIS_imaging_gtcsetup.py
    *local file- OSIRIS_imaging_functions.py

Store the pipeline and setup/functions files in your favorite folder and add the folder to your PATH and PYTHON_PATH environment. 
To do this on LINUX/MAC in a terminal:

    export PATH=$PATH:/new/path
    export PYTHON_PATH=$PYTHON_PATH:/new/path

Then make the file an executable:  
    
    chmod u+x OSIRIS_imaging_pipeline.py

On WINDOWS, you can either use the AnacondaPrompt3 or you can install Ubuntu from the Microsoft Store. If you use Ubuntu, follow the instructions above. If you use Anaconda prompt, follow the instructions here: https://correlated.kayako.com/article/40-running-python-scripts-from-anywhere-under-windows 
If you are using the AnacondaPrompt3 and you get an error saying 'the following arguments are required: objectid', check for the first print statement showing the arguments that the program read in (ex: INPUT TEXT HERE). If it is an empty list, search for the Registry Editor in the Windows search bar and open it. Click the down arrow on each of the following folders: HKEY_CLASSES_ROOT, Applications, python.exe, shell, open. Then click on 'command'. In the window on the right, double click on the name (Default) , then add 
    
    %*"
    
to the end of the text (the whole string should like like:
    
    "C:\Users\username\anaconda3\python.exe" "%1" %*"
This allows Windows to pass all arguments to python.exe.

Now the pipeline can be run from any folder (where the the setup/function files are in the same folder). 
You can run the pipeline from any folder using the full path to the data with the following syntax:

    OSIRIS_phot_pipeline.py OBJNAME --workdir  "C:\Users\path\to\raw_files" --outputdir "C:\Users\test\path\to\output" --reduce_obj 2 --dobias --doflat --dobpm --docalib --docrmask --doskysub --dostack --dowcs --dointeractive --dooverwrite False --filter r,z

OR you can move into the folder where your data is stored and then --workdir and --outputdir are not needed:

    OSIRIS_phot_pipeline.py OBJNAME -reduce_obj 2 --dobias --doflat --dobpm --docalib --docrmask --doskysub --dostack --dowcs --dointeractive --dooverwrite False --filter r,z

OBJNAME is a required input and it is the name of the object that is given in the header under the keyword 'OBJECT'.
'workdir' is the full path to the raw files and 'outputdir' is the full path to where you want to store the final files (it does not have to be in the same folder as the raw files).

Or if the master bias, flat, and bad pixel mask files already exist, give the names of their files with the following syntax (default file names are shown):

    OSIRIS_phot_pipeline.py OBJNAME --bias MasterBias --flat MasterFlat --bpm MasterBPM --reduce_obj 2 --docalib --docrmask --doskysub --dostacking --dowcs --dointeractive --dooverwrite False --filter r,z

If they have the default file names, you can omit those parameters entirely:
    
    OSIRIS_phot_pipeline.py OBJNAME --reduce_obj 2 --docalib --docrmask --doskysub --dostacking --dowcs --dointeractive --dooverwrite False --filter r,z

If you have already done some of the steps with the pipeline and would like to skip them, omit the '--do'. For example, if you already did calibrations and crmask (in addition to bias, flat, bpm) the call
    
    OSIRIS_phot_pipeline.py OBJNAME --reduce_obj 2 --doskysub --dostacking --dowcs --dointeractive --dooverwrite False --filter r,z
    
The output files are:

    Main files:
        MasterBias.fits
        MasterFlatSloan_g.fits (and other filters)
        MasterBPMSloan_g.fits (and other filters)
        skymap for each filter split into ccd1 and ccd2 for both the target and standard star 
        stacked image after aligning with each other for each filter split into ccd1 and ccd2 for both the target and standard star 
    
    Final images combined by filter split into ccd1 and ccd2 for both the target and standard star:
        objname_OSIRIS_filter_final_ccd1.fits and objname_OSIRIS_filter_final_ccd2.fits
        Standard_Star_OSIRIS_filter_final_ccd1.fits and Standard_Star_OSIRIS_filter_final_ccd2.fits
    
    Diagnostic folder which contains the following:
        1-calib
            images with bias, flat applied, split into ccd1 and ccd2 for both the target and standard star
        2-bpm
            images with the bad pixel mask applied, split into ccd1 and ccd2 for both the target and standard star
        3-crmask
            CRmask: Just the cosmic ray mask images split into ccd1 and ccd2 for both the target and standard star
            CRmask_applied: Science images with the cosmic ray mask applied split into ccd1 and ccd2 for both the target and standard star
        4-skymap
            bkgimg - the background image as determined by source extractor (sep.Background())
            mask1 - the mask after the first round of sigma clipping to remove stars
            mask2 - the mask after the second round of sigma clipping to remove stars
            sciN_final - 
        5-astrometry
            images after aligning with eachother split into ccd1 and ccd2 for both the target and standard star
        
See the output folder in INSERT HERE for examples of a good master bias, flat, sky map, and final image as well as the raw data (which is also in the copy folder).

Some parameters you may want to change are:

    -INSERT HERE

