# OSIRIS_Imaging_DRP
This is a data reduction pipeline for OSIRIS imaging. It performs bias subtraction, flat division, bad pixel mask, cosmic ray cleaning, sky subtraction, astrometric correction, and stacking.

Pipeline requirements:

    ccdproc

    astropy
    
    photutils
    
    astroscrappy
    
    astroquery
    
    scipy
    
    numpy
    
    matplotlib
    
    *local file- OSIRIS_imaging_gtcsetup.py


Pipeline can be run from its folder and is used with the following syntax.


 python OSIRIS_phot_pipeline.py
 
 --workdir  "C:\Users\test\datafolder"
 
 --outputdir "C:\Users\test\datafolder\output"
 
 object_name
 
 --dobias
 
 --doflat
 
 --domask
 
 --dowcs
 
 --dooverwrite False
 
 --filter g,i,z
