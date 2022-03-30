#!/usr/bin/env python

r"""
Created on Wed Nov  3 11:39:52 2021.

@author: amy

Pipeline requirements:
    ccdproc
    astropy
    photutils
    astroscrappy
    astroquery
    scipy
    numpy
    matplotlib
    spalipy
    *local file- OSIRIS_imaging_setup.py

The pipeline can be run from any folder and is used with the following syntax:

python OSIRIS_phot_pipeline.py
OBJNAME
--workdir "C:\Users\path\to\raw_files"
--outputdir "C:\Users\test\path\to\output"
--dooverwrite True
--dobias
--doflat
--domask
--docrmask
--doskysub
--dowcs
--dointeractive
--filter g,i,z

OBJNAME is a required input and it is the name of the object that is given in
the header under the keyword 'OBJECT'. 'workdir' is the full path to the raw
files and 'outputdir' is the full path to where you want to store the final
files (it does not have to be in the same folder as the raw files).

Or if the master bias, flat, and bad pixel mask files already exist:

python OSIRIS_phot_pipeline.py
OBJNAME
--workdir "C:\Users\path\to\raw_files"
--outputdir "C:\Users\path\to\output"
--bias MasterBias
--flat MasterFlat
--mask MasterBPM
--dooverwrite False
--filter g,i,z

Things to update in the future:
     handle .gz files
     create configuration file for parameters that can be changed
         or just make note of them up here
     do something else if no gaia stars are found

NOTE: If --dooverwrite is True but you don't have --dobias --doflat and --domask,
all files will be deleted so there will be no master files and the program
will crash.
"""
# # Local dependencies
import OSIRIS_imaging_setup as gtcsetup
import OSIRIS_imaging_functions as gtcdo
####
from time import time
from glob import glob
import sys
import copy
import ccdproc
import numpy as np
from astropy.nddata import CCDData
import astropy.units as u
from astropy.io import fits
from datetime import datetime
import os
import warnings
###
# from matplotlib import use
# use('Agg')


__author__ = "A. Gottlieb"
__created__ = "2021-11-03"


def main(argv):
    """Imaging reduction pipeline."""
    tstart = time()
    times = []
    actions = []

    warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')

    # Instrument specific constants
    saturation = 65535
    gain = 0.95
    rdnoise = 4.5  # electrons
    nccd = 2
    # ccdgap  = 37   # pixels

    # ############################## SETUP #################################
    # Get input arguments
    args, use_slash = gtcsetup.read_args()
    if args.doall:
        args = gtcsetup.change_args(args)

    now = datetime.now()  # current date and time
    date_time = now.strftime("%m-%d-%Y_%H-%M-%S")

    # Create log file name based on current date/time so you have a new file
    # each time you run
    new_logfile = args.workdir+args.logfile[:-4]+date_time+'.txt'

    log_fname = open(new_logfile, 'a')
    gtcsetup.print_both(log_fname, 'Writing output to file', new_logfile)
    gtcsetup.print_both(log_fname, 'ARGUMENTS:', argv[1:])

    # ########################### ADMIN ####################################
    msg = 'Step 1: Administration'
    gtcsetup.print_both(log_fname, msg)

    # Get working directory
    topdir = args.workdir
    gtcsetup.print_both(log_fname, topdir)

    # Move into the working directory
    os.chdir(topdir)
    gtcsetup.print_both(log_fname, '    Changing to directory', topdir)
    raw_path = topdir

    if not os.path.exists(raw_path):
        sys.exit('*** FATAL ERROR *** Data directory not found!')
    gtcsetup.print_both(log_fname, '    Data are in directory', raw_path)

    gtcsetup.print_both(log_fname, '    Working on object ',
                        args.objectid, 'in filter(s)', args.filt)

    # get all fits files in directory
    gtcsetup.print_both(log_fname, '    Getting list of all files')
    raw_list = [i.replace(os.sep, '/') for i in glob(raw_path+'0*.fits*')]
    if len(raw_list) == 0:
        raw_list = [i.replace(os.sep, '/') for i in glob(raw_path+'0*.fits.gz*')]

    gtcsetup.print_both(log_fname, np.array(raw_list))
    gtcsetup.print_both(log_fname, '    Overwriting:', args.dooverwrite)

    # directory for all diagnostic images such as sky frames or cosmic ray mask
    diag_path = args.outputdir+'diagnostic/'
    gtcsetup.print_both(log_fname, 'Diagnostic path:', diag_path)

    # sort files into different folders based on header information
    # fb, ff etc are lists of files containing the full path
    fb, ff, fsci, fstd, all_paths = gtcsetup.sort_files(args.objectid,
                                                        raw_path, raw_list,
                                                        use_slash, args.outputdir,
                                                        args.dooverwrite,
                                                        diag_path, log_fname)
    calib_path, bpm_path, crmask_path, skymap_path, astrom_path = all_paths

    # print the observing log;
    # ****NOTE: THIS CAN BE USED TO CHECK HEADER INFORMATION FOR ALL FILES****
    # I created a script called update header.py to change obstype from focus
    # to object
    gtcsetup.print_obslog(args.objectid, args.filt, fb, ff, fsci, fstd, log_fname)

    # if doobslog is true (e.g. only print the observing log), stop here
    if args.doobslog:
        sys.exit()

    gtcsetup.print_both(log_fname, '--------------------------------------')

    if len(fsci) == 0:
        sys.exit('No science images found; check your OBJNAME.')

    # #######################################################
    # ##################### MASTER BIAS #####################
    # #######################################################

    msg = 'Step 2: Master Bias'
    gtcsetup.print_both(log_fname, msg)

    gtcsetup.print_both(log_fname, '    Dobias = ', args.dobias)
    if args.dobias:
        gtcsetup.print_both(log_fname, '    Creating master bias')
        mbias = [None, None]
        mbias, bias_fname = gtcsetup.createMaster(fb, 'bias', nccd, mbias,
                                                  gain, rdnoise, args.outputdir,
                                                  args.biasfile, '', log_fname)
    else:
        # Just noting file must be named MasterBias
        gtcsetup.print_both(log_fname, '    Reading in bias',
                            args.outputdir+args.biasfile+'.fits')
        mbias = [CCDData.read(args.outputdir+args.biasfile+'.fits',
                              hdu=x+1, unit=u.electron) for x in range(nccd)]
        if not os.path.exists(args.outputdir+args.biasfile+'.fits'):
            sys.exit('*** FATAL ERROR *** Bias file not found!')

    times, actions = gtcsetup.save_times(
        tstart, time(), 'MasterBias', 'creating master bias',
        times, actions, log_fname)

    #########################################################################

    # Get all filters taken in observations
    all_obs_filt = gtcsetup.get_filters(args.filt.split(','))
    gtcsetup.print_both(log_fname, '-------------------------------------------')
    gtcsetup.print_both(log_fname, 'Now looping through filters:', all_obs_filt)

    # Loop through filters
    for filt in all_obs_filt:
        gtcsetup.print_both(log_fname, 'Running through filter', filt)

        # #######################################################
        # ##################### MASTER FLAT #####################
        # #######################################################
        msg = 'Step 3: Master Flat; filter '+filt
        gtcsetup.print_both(log_fname, msg)

        gtcsetup.print_both(log_fname, '    Doflat', args.doflat)
        if args.doflat:

            mflat, flatname = gtcsetup.createMaster(ff, 'flat', nccd, mbias,
                                                    gain, rdnoise,
                                                    args.outputdir,
                                                    args.flatfile, filt, log_fname)

        else:
            # Again noting the file must be called MasterFlatSloan_g.fits etc
            flatname = args.outputdir+args.flatfile+filt+'.fits'
            gtcsetup.print_both(log_fname, '    Reading in Flat', flatname)
            mflat = [CCDData.read(flatname,
                                  hdu=x+1, unit=u.electron) for x in range(nccd)]
            if not os.path.exists(flatname):
                sys.exit('*** FATAL ERROR *** Flat file not found!')

        times, actions = gtcsetup.save_times(
            tstart, time(), 'MasterFlat:'+filt, 'creating master flat',
            times, actions, log_fname)
        gtcsetup.print_both(log_fname, '---------------------')

        # #######################################################
        # ################### BAD PIXEL MASK ####################
        # #######################################################

        msg = 'Step 4: Mask Bad Pixels; filter '+filt
        gtcsetup.print_both(log_fname, msg)

        gtcsetup.print_both(log_fname, '    Domask', args.domask)
        if args.domask:
            with fits.open(flatname) as fits_open:
                hdr = fits_open[0].header
            mask = [ccdproc.ccdmask(mflat[k])  # , findbadcolumns=False)
                    for k in range(nccd)]
            mask_ccd = [CCDData(data=x.astype('uint8'),
                                unit=u.dimensionless_unscaled) for x in mask]
            mask_hdu = fits.HDUList([fits.PrimaryHDU(header=hdr)])
            for x in mask_ccd:
                mask_hdu.append(fits.ImageHDU(x.data))
            # Bad pixels have value of 1, good pixels have a value of 0
            gtcsetup.print_both(log_fname, '    Writing bad pixel mask to ',
                                args.outputdir+'MasterBPM'+filt+'.fits')
            mask_hdu.writeto(args.outputdir+'MasterBPM'+filt+'.fits', overwrite=True)
        else:
            # Again noting the file must be called MasterBPMSloan_g.fits etc
            bpmask_name = args.outputdir+args.maskfile+filt+'.fits'
            gtcsetup.print_both(
                log_fname, '    Reading in bad pixel mask', bpmask_name)
            mask_ccd = [CCDData.read(bpmask_name,
                                     hdu=x+1, unit=u.dimensionless_unscaled)
                        for x in range(nccd)]

            if not os.path.exists(bpmask_name):
                sys.exit('*** FATAL ERROR *** Master bad pixel mask file not found!')

            # NOTE: sigma clipping is getting rid of vertical streaks in the science
            # images that are not in the master flat which is what the bad pixel mask
            # is based on

        times, actions = gtcsetup.save_times(
            tstart, time(), 'MasterBPM: '+filt, 'creating master bad pixel mask',
            times, actions, log_fname)

        gtcsetup.print_both(log_fname, '---------------------')

        ############################################################
        if int(args.reduce_obj) == 0:
            gtcsetup.print_both(log_fname, 'Only reducing the science target')
            objs = [fsci]
            roots = [args.objectid.replace(' ', '')+'_OSIRIS_']
        elif int(args.reduce_obj) == 1:
            gtcsetup.print_both(log_fname, 'Only reducing the standard star')
            objs = [fstd]
            roots = ['Standard_Star_OSIRIS_']
        else:
            gtcsetup.print_both(
                log_fname, 'Reducing the science target AND the standard star')
            roots = [args.objectid.replace(
                ' ', '')+'_OSIRIS_', 'Standard_Star_OSIRIS_']
            objs = [fsci, fstd]

        # Loop through target object and standard star
        for z, obj in enumerate(objs):

            obj = np.array(obj)
            root = roots[z]

            gtcsetup.print_both(log_fname, '---------------------')

            msg = 'Step 5: Science Frames; filter '+filt+'; object:'+root
            gtcsetup.print_both(log_fname, msg)
            gtcsetup.print_both(
                log_fname, 'Applying master bias, flat, bpm, and trimming')
            gtcsetup.print_both(log_fname, 'Found', len(obj), 'files')

            # #######################################################
            # ########### CALIBRATIONS (BIAS, FLAT, BPM) ############
            # #######################################################

            # Either do calibrations or read in images from diagnostic folder
            # where calibrations have already been applied
            if args.docalib:
                all_sci_calib, all_headers, all_filts = gtcdo.do_calib(
                    obj, filt, log_fname, nccd, mbias, mflat, mask_ccd, gain,
                    rdnoise, calib_path, bpm_path, root)
            else:
                all_sci_calib = gtcsetup.read_in_files(
                    bpm_path, root+filt, log_fname)
                all_headers, all_filts = gtcsetup.get_header_info(obj, filt)

            times, actions = gtcsetup.save_times(
                tstart, time(), 'Calib '+root,
                'apply master bias, flat, bpm, and trim object '+root,
                times, actions, log_fname)

            # #######################################################
            # ################## COSMIC RAY REMOVAL #################
            # #######################################################
            msg = 'Step 6: Cosmic Ray Removal'
            gtcsetup.print_both(log_fname, msg)
            # cosmic ray cleaning must be done before sky subtraction
            # NOTE: it's still finding saturated stars as crs
            if args.docrmask is True:

                all_sci_proc = gtcdo.do_crmask(
                    log_fname, all_sci_calib, mask_ccd, gain, rdnoise, saturation,
                    obj[all_filts], all_headers, root, filt, nccd, crmask_path)

            else:
                # Check for diagnostic folder with files; if files exist, read
                # them in, otherwise, skip
                check = 'CRmask_applied_'+root+filt
                use_path = crmask_path
                gtcsetup.print_both(log_fname, 'Checking for files named ',
                                    check, ' in directory', use_path)
                check_fnames = [f for f in os.listdir(use_path) if check in f]
                gtcsetup.print_both(log_fname, 'Found', len(check_fnames),
                                    'files in directory')

                if os.path.isdir(use_path) and len(check_fnames) > 0:

                    all_sci_proc = gtcsetup.read_in_files(
                        use_path, check, log_fname)
                else:
                    all_sci_proc = copy.deepcopy(all_sci_calib)
                    gtcsetup.print_both(
                        log_fname, '     Skipping cosmic ray rejection')

            times, actions = gtcsetup.save_times(
                tstart, time(), 'Cosmic ray mask',
                'creating and applying a cosmic ray mask',
                times, actions, log_fname)
            gtcsetup.print_both(log_fname, '---------------------')

            # #######################################################
            # ################### SKY SUBTRACTION ###################
            # #######################################################

            msg = 'Step 7: Background and sky subtraction'
            gtcsetup.print_both(log_fname, msg)

            # If there is only one
            if len(all_sci_proc) == 1:
                gtcsetup.print_both(
                    log_fname, 'Only 1 image; not doing sky subtraction')
                doskysub = False
            else:
                doskysub = args.doskysub

            gtcsetup.print_both(log_fname, '    doskysub', doskysub)

            # If the user wants to subtract the sky, do it, otherwise just subtract
            # the background
            if doskysub is True:

                gtcsetup.print_both(log_fname, 'Doing sky subtraction')

                sci_final, sci_skymap = gtcdo.do_bkg_sky_subtraction(
                    all_sci_proc, all_headers, log_fname, skymap_path,
                    obj[all_filts], root, filt)

                # Write out the sky image
                gtcsetup.write_ccd(all_headers[0][0], all_headers[0][1:],
                                   sci_skymap, args.outputdir,
                                   'skymap.fits', root, filt, log_fname)
                to_write = True

            else:
                # Check for diagnostic folder with files; if files exist, read
                # them in, otherwise, skip
                check = 'final_'+root+filt
                use_path = skymap_path
                gtcsetup.print_both(log_fname, 'Checking for files named ',
                                    check, ' in directory', use_path)
                check_fnames = [f for f in os.listdir(use_path) if check in f]
                gtcsetup.print_both(log_fname, 'Found', len(check_fnames),
                                    'files in directory')

                if os.path.isdir(use_path) and len(check_fnames) > 0:
                    gtcsetup.print_both(log_fname, 'Looking for files with ',
                                        check, 'in the file name')
                    sci_final = gtcsetup.read_in_files(
                        use_path, check, log_fname)
                    to_write = False
                else:
                    gtcsetup.print_both(log_fname, 'Doing bkg subtraction')
                    sci_final = gtcdo.do_bkg_subtraction(
                        all_sci_proc, all_headers, log_fname, skymap_path)
                    to_write = True

            if to_write is True:
                # Write out each calibrated science image
                for i, f in enumerate(obj[all_filts]):
                    gtcsetup.write_ccd(all_headers[i][0], all_headers[i][1:],
                                       sci_final[i], skymap_path,
                                       f[:-5]+'_final.fits', root, filt, log_fname)

            times, actions = gtcsetup.save_times(
                tstart, time(), 'Bkg/Sky subtraction',
                'doing sky/bkg subtraction',
                times, actions, log_fname)
            gtcsetup.print_both(log_fname, '---------------------')

            # #################################s########################
            # ################## ASTROMETRY PT 1 #######################
            # ##########################################################

            # Align images with eachother and then median combine
            msg = 'Step 8: Astrometry'
            gtcsetup.print_both(log_fname, msg)

            if args.dostack:

                gtcsetup.print_both(
                    log_fname, '     Part 1- Aligning images with eachother')

                final_aligned_image = gtcdo.do_stacking(
                    sci_final, all_headers, args, root, filt, astrom_path, log_fname)

            else:

                # Check for diagnostic folder with files; if files exist, read
                # them in, otherwise, skip
                if os.path.isdir(args.outputdir) and len(os.listdir(
                        args.outputdir)) > 0:

                    final_aligned_image = gtcsetup.read_in_files(
                        args.outputdir, 'aligned_'+root+filt, log_fname)[0]

                else:
                    # ###### CHECK THIS ############
                    final_aligned_image = sci_final

            times, actions = gtcsetup.save_times(
                tstart, time(), 'Align pt 1', 'aligning images with eachother',
                times, actions, log_fname)

            # #################################s########################
            # ################## ASTROMETRY PT 2 #######################
            # ##########################################################
            if args.dowcs:

                gtcsetup.print_both(
                    log_fname, '     Part 2- Aligning median combined image',
                    ' with GAIA')

                # Do the interactive astrometry or automatic astrometry
                if args.dointeractive:

                    # Final image is written to fits file inside this function
                    gtcsetup.print_both(log_fname, 'Working on CCD 1')

                    final_name = args.outputdir+'final_'+root+filt+'_ccd1.fits'
                    fname = args.outputdir+'aligned_'+root+filt+'_ccd1.fits'
                    gtcdo.do_interactive_astrometry(
                        final_name, fname, filt, '1', astrom_path, log_fname)
                    ####
                    gtcsetup.print_both(log_fname, 'Working on CCD 2')

                    final_name = args.outputdir+'final_'+root+filt+'_ccd2.fits'
                    fname = args.outputdir+'aligned_'+root+filt+'_ccd2.fits'
                    gtcdo.do_interactive_astrometry(
                        final_name, fname, filt, '2', astrom_path, log_fname)

                    times, actions = gtcsetup.save_times(
                        tstart, time(), 'Align pt 2',
                        'aligning images with gaia interactively',
                        times, actions, log_fname)

                else:
                    final_name = args.outputdir+'final_auto_'+root+filt+'_ccd'
                    updated_headers = gtcdo.do_auto_astrometry(
                        final_aligned_image, filt, args.seeing, astrom_path,
                        final_name, log_fname)

                    times, actions = gtcsetup.save_times(
                        tstart, time(), 'Align pt 2',
                        'aligning images with gaia automatically',
                        times, actions, log_fname)

                    # Write out the final image to the output folder
                    gtcsetup.write_ccd(all_headers[0][0], updated_headers,
                                       final_aligned_image, args.outputdir,
                                       'final_auto.fits', root, filt, log_fname)

            times, actions = gtcsetup.save_times(
                tstart, time(), 'End filter loop',
                'executing one object in one filter', times, actions, log_fname)

            gtcsetup.print_both(
                log_fname, '--------------------------------------------------')

        times, actions = gtcsetup.save_times(tstart, time(), 'End object loop',
                                             'executing all objects in one filter',
                                             times, actions, log_fname)

    times, actions = gtcsetup.save_times(
        tstart, time(), 'End all loops',
        'executing all objects in all filters', times, actions, log_fname)

    print(times)
    print(actions)

    gtcsetup.print_pipeline_times(times, actions, log_fname)

    gtcsetup.print_both(log_fname, '*** Done ***')

    log_fname.close()


if __name__ == "__main__":

    main(sys.argv)
