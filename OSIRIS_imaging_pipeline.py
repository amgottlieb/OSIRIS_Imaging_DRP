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
    *local file- OSIRIS_imaging_gtcsetup.py


Pipeline can be run from its folder and is used with the following syntax.

  python OSIRIS_phot_pipeline.py
  --workdir "C:\Users\test\data_folder"
  --outputdir "C:\Users\test\data_folder\output"
  object_name
  --dobias
  --doflat
  --domask
  --dooverwrite False
  --filter g,i,z

If you already have master bias, flats, bpmask they must be named
MasterBias.fits or MasterFlatSloan_g.fits or MasterBPMSloan_g.fits OR be
specified via --bias etc below.
To run after creating master bias+flat+mask:

  python OSIRIS_phot_pipeline.py
  --workdir "C:\Users\test\data_folder"
  --outputdir "C:\Users\test\data_folder\output"
  object_name
  --bias MasterBias
  --flat MasterFlat
  --mask MasterBPM
  --dooverwrite False
  --filter g,i,z

Things to update in the future:
     add check on astrometry
     handle .gz files
     create configuration file for parameters that can be changed
         or just make note of them up here
     save diagnostic figures
     do something else if no gaia stars are found

NOTE: If --dooverwrite is True but you don't have --dobias --doflat and --domask,
all files will be deleted and the program will crash
"""

__author__ = "A. Gottlieb"
__created__ = "2021-11-03"


from time import time
from glob import glob
import os
import sys
import copy
import ccdproc
import numpy as np
from astropy.nddata import CCDData
import astropy.units as u
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.io import fits
from astropy.modeling import models
from astropy.coordinates import SkyCoord
from astropy.wcs import wcs
from astroquery.ipac.irsa import Irsa
from photutils import DAOStarFinder
from photutils.background import MeanBackground
from astroscrappy import detect_cosmics
import matplotlib.pyplot as plt
from datetime import datetime

# # Local dependencies
import OSIRIS_imaging_setup as gtcsetup


def main(argv):
    r"""Pipeline can be run from its folder and is used with the following syntax.

    python OSIRIS_phot_pipeline.py
    --workdir "C:\Users\test\data_folder"
    --outputdir "C:\Users\test\data_folder\output"
    object_name
    --dobias
    --doflat
    --domask
    --dooverwrite False
    --filter g,i,z

    If you already have master bias, flats, bpmask they must be named
    MasterBias.fits or MasterFlatSloan_g.fits or MasterBPMSloan_g.fits OR be
    specified via --bias etc below.
    To run after creating master bias+flat+mask:

    python OSIRIS_phot_pipeline.py
    --workdir "C:\Users\test\data_folder"
    --outputdir "C:\Users\test\data_folder\output"
    object_name
    --bias MasterBias
    --flat MasterFlat
    --mask MasterBPM
    --dooverwrite False
    --filter g,i,z
    """
    tstart = time()

    # Instrument specific constants
    saturation = 65535
    gain = 0.95
    rdnoise = 4.5  # electrons
    nccd = 2
    # ccdgap  = 37   # pixels

    to_plot = False

    # Get input arguments
    args, use_slash = gtcsetup.read_args()

    now = datetime.now()  # current date and time
    date_time = now.strftime("%m-%d-%Y_%H-%M-%S")
    new_logfile = args.outputdir+args.logfile[:-4]+date_time+'.txt'

    log_fname = open(new_logfile, 'a')
    gtcsetup.print_both(log_fname, 'Writing output to file', new_logfile)

    msg = 'Step 1: Administration'
    # gtcsetup.print_both(log_fname, bcolors.HEADER + bcolors.BOLD +
    # '\n{}\n'.format(msg) + bcolors.ENDC)
    gtcsetup.print_both(log_fname, msg)

    # Get working directory
    topdir = args.workdir
    gtcsetup.print_both(log_fname, topdir)

    # move into the working directory
    os.chdir(topdir)
    gtcsetup.print_both(log_fname, '    Changing to directory', topdir)
    raw_path = topdir

    if not os.path.exists(raw_path):
        sys.exit('*** FATAL ERROR *** Data directory not found!')
    gtcsetup.print_both(log_fname, '    Data are in directory', raw_path)

    gtcsetup.print_both(log_fname, '    Working on object ',
                        args.objectid, 'in filters', args.filt)

    # get all fits files in directory
    gtcsetup.print_both(log_fname, '    Getting list of all files')
    raw_list = [i.replace(os.sep, '/') for i in glob(raw_path+'0*.fits*')]
    gtcsetup.print_both(log_fname, np.array(raw_list))
    gtcsetup.print_both(log_fname, '    Overwriting:', args.dooverwrite)

    # directory for all diagnostic images such as sky frames or cosmic ray mask
    diag_path = args.outputdir+'diagnostic/'
    gtcsetup.print_both(log_fname, 'Diagnostic path:', diag_path)
    # sort files into different folders based on header information
    fb, ff, fsci, fstd = gtcsetup.sort_files(args.objectid,
                                             raw_path, raw_list,
                                             use_slash, args.outputdir,
                                             args.dooverwrite,
                                             diag_path, log_fname)

    # print the observing log;
    # ****NOTE: THIS CAN BE USED TO CHECK HEADER INFORMATION FOR ALL FILES****
    # I created a script called update header.py to change obstype from focus
    # to object
    gtcsetup.print_obslog(args.objectid, args.filt, fb, ff, fsci, fstd)

    # TO DO: include update header script?

    # if dolog is true (e.g. only print the observing log), stop here
    if args.dolog:
        sys.exit()

    gtcsetup.print_both(log_fname, '--------------------------------------')
    msg = 'Step 2: Master Bias'
    # gtcsetup.print_both(log_fname, bcolors.HEADER + bcolors.BOLD +
    # '\n{}\n'.format(msg) + bcolors.ENDC)
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

    # Get all filters taken in observations
    all_filt = gtcsetup.get_filters(args.filt.split(','))
    gtcsetup.print_both(log_fname, '-------------------------------------------')
    gtcsetup.print_both(log_fname, 'Now looping through filters:', all_filt)

    for filt in all_filt:
        gtcsetup.print_both(log_fname, 'Running through filter', filt)

        msg = 'Step 3: Master Flat; filter '+filt
        # gtcsetup.print_both(log_fname, bcolors.HEADER + bcolors.BOLD +
        # '\n{}\n'.format(msg) + bcolors.ENDC)
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

        gtcsetup.print_both(log_fname, '---------------------')
        msg = 'Step 4: Mask Bad Pixels; filter '+filt
        # gtcsetup.print_both(log_fname, bcolors.HEADER + bcolors.BOLD +
        # '\n{}\n'.format(msg) + bcolors.ENDC)
        gtcsetup.print_both(log_fname, msg)

        gtcsetup.print_both(log_fname, '    Domask', args.domask)
        if args.domask:
            with fits.open(flatname) as fits_open:
                hdr = fits_open[0].header
            mask = [ccdproc.ccdmask(mflat[k], findbadcolumns=False)
                    for k in range(nccd)]
            mask_ccd = [CCDData(data=x.astype('uint8'),
                                unit=u.dimensionless_unscaled) for x in mask]
            mask_hdu = fits.HDUList([fits.PrimaryHDU(header=hdr)])
            for x in mask_ccd:
                mask_hdu.append(fits.ImageHDU(x.data))
            # Bad pixels have value of 1, good pixels have a value of 0
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
                sys.exit('*** FATAL ERROR *** Bad pixel mask file not found!')

        objs = [fsci, fstd]

        for z, obj in enumerate(objs):
            if z == 0:
                root = args.objectid.replace(' ', '')+'_OSIRIS_'
            else:
                root = 'Standard_Star_OSIRIS_'
            gtcsetup.print_both(log_fname, '---------------------')
            msg = 'Step 5: Science Frames; filter '+filt+'; object:'+root
            gtcsetup.print_both(log_fname, msg)
            # gtcsetup.print_both(log_fname, bcolors.HEADER + bcolors.BOLD +
            # '\n{}\n'.format(msg) + bcolors.ENDC)
            for j, f in enumerate(obj):
                with fits.open(f) as fits_open:
                    # Note: because images have been trimmed, the wcs information is
                    # no longer correct; subtract off the trimmed amount to get close
                    hdr = fits_open[0].header

                    hdr1 = fits_open[1].header
                    trim1 = hdr1['TRIMSEC'].split(':')[0][1:]
                    hdr1['CRPIX1'] = float(hdr1['CRPIX1'])-float(trim1)

                    hdr2 = fits_open[2].header
                    trim2 = hdr2['TRIMSEC'].split(':')[0][1:]
                    hdr2['CRPIX1'] = float(hdr2['CRPIX1'])-float(trim2)

                    hdrs = [hdr1, hdr2]
                    hdr_filt = hdr['FILTER2']

                    # get initial wcs info from header
                    if j == 0:
                        # NOTE: somehow updating hdr1 and 2 above also updates
                        # this header too...
                        wcs_ref = [wcs.WCS(fits_open[k+1]) for k in range(nccd)]

                    # Skip over this image if its filter does not match the
                    # current one
                    if hdr_filt != filt:
                        continue

                    gtcsetup.print_both(log_fname, '    Working on file', f)
                    # get raw frame
                    sci_raw = [CCDData.read(f, hdu=k+1, unit='adu')
                               for k in range(nccd)]

                    # apply bias, flat, bpm corrections
                    # again noting that Nora had the oscan, oscan_model and trim;
                    # I didn't change this
                    sci_proc_init = [ccdproc.ccd_process(
                        x, oscan='[3:22,:]',
                        # x.header['BIASSEC'],
                        oscan_model=models.Chebyshev1D(3),
                        trim=x.header['TRIMSEC'],
                        master_bias=mbias[k],
                        master_flat=mflat[k],
                        bad_pixel_mask=np.array(mask_ccd[k]),
                        gain=gain*u.electron/u.adu,
                        readnoise=rdnoise*u.electron,
                        gain_corrected=True)
                        for k, x in enumerate(sci_raw)]

                    sci_proc_no_bpm = [ccdproc.ccd_process(
                        x, oscan='[3:22,:]',
                        # x.header['BIASSEC'],
                        oscan_model=models.Chebyshev1D(3),
                        trim=x.header['TRIMSEC'],
                        master_bias=mbias[k],
                        master_flat=mflat[k],
                        bad_pixel_mask=None,
                        gain=gain*u.electron/u.adu,
                        readnoise=rdnoise*u.electron,
                        gain_corrected=True)
                        for k, x in enumerate(sci_raw)]

                    plt.figure()
                    plt.imshow(sci_proc_no_bpm[1].data-sci_proc_init[1].data,
                               interpolation=None, origin='lower')
                    gtcsetup.print_both(
                        log_fname, ((sci_proc_no_bpm[1].data-sci_proc_init[1].data)))
                    sys.exit()
                    # cosmic ray cleaning must be done before sky subtraction
                    # NOTE: it's still finding saturated stars as crs

                    gtcsetup.print_both(log_fname, '    Removing cosmic rays')
                    cleaned_sci = []
                    cr_mask = []
                    for x in range(len(sci_proc_init)):

                        test_mask, _clean = detect_cosmics(np.array(
                            sci_proc_init[x].data),
                            # NOTE: want None, otherwise bad pixels aren't
                            # masked out?
                            inmask=None,
                            sigclip=4.5,  # lower values flag more pixels as
                            # cosmic rays
                            sigfrac=0.3,
                            objlim=10.0,  # increase if centers of stars flagged
                            # as cosmic rays
                            gain=gain, readnoise=rdnoise, satlevel=saturation,
                            niter=4, sepmed=True,
                            cleantype='meanmask', fsmode='median',
                            psfmodel='gauss', psffwhm=2.5, psfsize=7,
                            psfk=None, psfbeta=4.765,
                            verbose=True)

                        cleaned_sci.append(_clean)
                        cr_mask.append(np.multiply(test_mask, 1))

                    if to_plot is True:
                        gtcsetup.plot_cr(sci_proc_init, cr_mask, cleaned_sci[1],
                                         mask_ccd, diag_path, f[:-5]+'_CRmask.png',
                                         root, filt, log_fname)

                    # check how cosmic ray mask is affecting the image
                    # plt.figure()
                    # plt.imshow(cleaned_sci[1]-sci_proc_init[1].data,
                    #            interpolation=None, origin='lower')

                    sci_proc = copy.deepcopy(sci_proc_init)
                    sci_proc[0].data = cleaned_sci[0]
                    sci_proc[1].data = cleaned_sci[1]

                    # Write out cosmic ray masks to diagnostic folders
                    gtcsetup.write_ccd(hdr, hdrs, cr_mask, diag_path,
                                       f[:-5]+'_CRmask.fits', root, filt, log_fname)

                    # Get the mean sky of the current image; Nora set these vals
                    mean_sky = MeanBackground(SigmaClip(sigma=3., maxiters=10))

                    # Create the skymap for the current image using the mean sky?
                    sci_skymap = [
                        mean_sky.calc_background(x)*np.ones(np.shape(x))*u.electron
                        for x in sci_proc]

                    # Write out cosmic ray masks to diagnostic folders
                    gtcsetup.write_ccd(hdr, hdrs, sci_skymap, diag_path,
                                       f[:-5]+'_skymap.fits', root, filt, log_fname)

                    # Subtract the skymap from the current image in both ccds
                    sky_sub = [x.subtract(sci_skymap[k],
                                          propagate_uncertainties=True,
                                          handle_meta='first_found')
                               for k, x in enumerate(sci_proc)]

                    # Divide sky subtracted image by the exposure time
                    sci_final = [x.divide(hdr['EXPTIME']*u.second,
                                          propagate_uncertainties=True,
                                          handle_meta='first_found')
                                 for x in sky_sub]

                    # if args.dowcs:
                    for n in range(len(sci_final)):
                        # if n == 0:
                        #     continue
                        # gtcsetup.print_both(log_fname, '    -----------------')
                        gtcsetup.print_both(
                            log_fname, '    Working on wcs info for CCD', n)
                        ima = sci_final[n]
                        imafile = f

                        # get approximate center of ccd
                        sky = wcs_ref[n].pixel_to_world(
                            ima.shape[1]/2, ima.shape[0]/2)

                        # Get gaia comparison catalog centered on ccd
                        # using cone with radius 6 arcmin (Nora chose this)
                        # NOTE: radius of 6 gives circle, larger radius gives
                        # weird shape that doesn't cover the ccd...
                        gaia = Irsa.query_region(SkyCoord(sky.ra, sky.dec,
                                                          frame='fk5'),
                                                 catalog="gaia_dr2_source",
                                                 spatial="Cone", radius=6*u.arcmin)

                        if len(gaia) == 0:
                            gtcsetup.print_both(log_fname, 'No GAIA stars found \
                                                within search radius for filter.')
                            # TO DO: do something else? exclude this image from
                            # further processing
                            break

                        # again Nora chose these values, I didn't change them
                        mean, median, std = sigma_clipped_stats(ima, sigma=3.0)
                        daofind = DAOStarFinder(fwhm=8., threshold=8.*std)
                        sources = daofind(np.asarray(ima))  # sources in image
                        if len(sources) == 0:
                            sys.exit('No sources detected within image.')

                        # First pass at fixing astrometry (cross reference gaia
                        # catalog with sources in image to get matched list,
                        w, good = gtcsetup.do_astrometry(gaia, sources,
                                                         wcs_ref[n], ima, log_fname)

                        # Second pass
                        w2, good2 = gtcsetup.do_astrometry(
                            gaia, sources, w, ima, log_fname)

                        if to_plot is True:
                            title = 'image ' + str(j)+' ccd ' + \
                                str(n) + ' filter '+filt

                            gtcsetup.plot_sources(wcs_ref[n], wcs_ref[n], ima,
                                                  'before correcting'+title,
                                                  sources, gaia,
                                                  np.arange(len(sources)))
                            #     gtcsetup.plot_sources(
                            #         w, w, ima, 'first pass'+title,
                            #         sources, gaia, good)

                            gtcsetup.plot_sources(w2, w2, ima,
                                                  'second pass'+title,
                                                  sources, gaia, good2)

                        fits_open[n+1].header.update(w2.to_header(relax=True))
                        hdrs[n].update(w2.to_header(relax=True))

                    # Write out 2 ccds to 2 separate files with correct headers
                    gtcsetup.write_ccd(hdr, hdrs, sci_final, args.outputdir,
                                       imafile, root, filt, log_fname)

                    gtcsetup.print_both(log_fname, '--------------------------')

            # Get all files needed to combine
            gtcsetup.combine_files(args.outputdir, root, filt, diag_path,
                                   use_slash, log_fname)

            gtcsetup.print_both(
                log_fname, '--------------------------------------------------')
    gtcsetup.print_both(log_fname, 'Total execution time',
                        ((time()-tstart))/60., 'min')
    # gtcsetup.print_both(log_fname, bcolors.OKBLUE +
    #       '\nTotal execution time {0:.1f} min'.format(((time()-tstart) / 60)) +
    #       bcolors.ENDC)
    gtcsetup.print_both(log_fname, '*** Done ***')
    log_fname.close()


if __name__ == "__main__":
    main(sys.argv[1:])
