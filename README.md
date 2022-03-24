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
    
to the end of the text so that the whole string should look like:
    
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
    aligned - stacked image after aligning with each other for each filter split into ccd1 and ccd2 for both the target and standard star 
    Final images combined by filter split into ccd1 and ccd2 for both the target and standard star:
        final_objname_OSIRIS_filter_ccd1.fits and final_objname_OSIRIS_filter_ccd2.fits
        Standard_Star_OSIRIS_filter_ccd1.fits and final_Standard_Star_OSIRIS_filter_ccd2.fits

Diagnostic folder which contains the following folders:
    1-calib
        images with bias, flat applied, split into ccd1 and ccd2 for both the target and standard star
    2-bpm
        images with the bias, flat and bad pixel mask applied, split into ccd1 and ccd2 for both the target and standard star
    3-crmask
        CRmask: The cosmic ray mask images split into ccd1 and ccd2 for both the target and standard star
        CRmask_applied: Science images with the cosmic ray mask applied split into ccd1 and ccd2 for both the target and standard star
    4-skymap
        sciN_final - Science images with the bkg (and sky) applied split into ccd1 and ccd2 for both the target and standard star
        final_mask - final mask used before creating 
    5-astrometry
        pad_aligned - images after aligning with eachother split into ccd1 and ccd2 for both the target and standard star
        footprint- footprints that show which direction the images were shifted; NOTE: if you're using ds9 to view these files, change the color to aips0
        
See the output folder in https://drive.google.com/drive/folders/12mIXEg5bA2IgjTOwK6NLanjZS-1Ls407?usp=sharing for examples of a good master bias, flat, sky map, and final image as well as the raw data (which is also in the copy folder).

If the default values don't work for you, some parameters you may want to change (with explanations) can be found at the top of the OSIRIS_imaging_functions.py file as well as below:

- cutoff = 2 arcsec; cutoff is used in calc_astrometry_accuracy around line 960. Separations between OSIRIS objects and GAIA stars that are less than this are considered good
- pad = 200 pixels; pad is used in do_stacking around line 1050. It is the number of rows/columns to add onto each side of the image before aligning and combining (in pixels)
-ellipse_lim = 1. is used in get_osiris_obj_cat around line 700; objects with ellipticities less than this will be included in the OSIRIS detected object catalog/for astrometry
-patch_radius = 5 pixels; it is used in plot_images and plot_check_images to set the radius of the circles around the detected stars
- parameters in detect_cosmics

Note: the remaining varibles are used in do_interactive_astrometry around line 800

-cross_match_r = patch_radius; as long as you click within the circle of this radius, this should be fine
- pix_limit = 0 pixels; it determines how close to the edge of the image stars can be; stars at pixel values less than this will not be included in the detected object catalog/for astrometry
- n_stars_init = 40; n_stars determines how many of the brightest stars (n_stars) will be displayed when doing interactive photometry. If there are less than 40 detected sources, n_stars will default to the number of detected sources.
- select_n_stars = 6 is the minimum number of stars you must choose when doing interactive astrometry. Absolute minimum should be 3, 6 is sufficient, but you can choose as many as you want
- sip_deg = 2 determines the distortion correction; this should not be more than 2

