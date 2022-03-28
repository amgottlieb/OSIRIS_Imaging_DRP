"""Contains functions used in OSIRIS_imaging_pipeline.py."""

from astropy.io import fits
import ccdproc
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
from astropy.modeling import models
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord
from astropy.wcs import wcs, WCS
import sep
from photutils.segmentation import make_source_mask
from astroquery.ipac.irsa import Irsa
from astropy import coordinates
from astroquery.skyview import SkyView
import time
from matplotlib.backend_bases import MouseButton
import matplotlib.font_manager as font_manager
import astroalign as aa
from astroscrappy import detect_cosmics
import copy

import OSIRIS_imaging_setup as gtcsetup

# Variables the user may want to change

# parameters in detect_cosmics

# cutoff is used in calc_astrometry_accuracy around line 1190
# separations between OSIRIS objects and GAIA stars that are less than this
# are considered good
cutoff = 2  # arcsec

# pad is used in do_stacking around line 1325
# it is the number of rows/columns to add onto each side of the image before
# aligning and combining (in pixels)
pad = 200

# ellipse_lim is used in get_osiris_obj_cat around line 874
# objects with ellipticities less than this will be included in the
# OSIRIS detected object catalog/for astrometry
ellipse_lim = 1.

# patch_radius is used in plot_images and plot_check_images
# has units of pixels
patch_radius = 5

# the remaining varibles are used in do_interactive_astrometry around line 980

# default is same as patch radius;
# as long as you click within the circle, this should be fine
cross_match_r = patch_radius

# pix_limit determines how close to the edge of the image stars can be;
# stars at pixel values less than this will not be included in the catalog/for
# astrometry
pix_limit = 0

#  n_stars determines how many of the brightest stars (n_stars) will
# be displayed when doing interactive photometry.
n_stars_init = 40

# select_n_stars is the minimum number of stars you must choose when doing
# interactive astrometry. Absolute minimum should be 3, 6 is sufficient,
# but you can choose as many as you want
select_n_stars = 6

# sip_deg determines the distortion correction
sip_deg = 2

# tolerance  is used in auto_astrometry to set how far large of an area the
# program can look for matching Gaia stars. If your astrometry is good,
# this can be smaller.
tolerance = 15


# do not change this:
pix_scale = 0.254  # arcsec per pixel


def do_calib(obj, filt, log_fname, nccd, mbias, mflat, mask_ccd, gain,
             rdnoise, calib_path, bpm_path, root):
    """Apply bias, flat and bad pixel mask to science frames.

    Parameters
    ----------
    obj (list) : list of file names
    filt (str) : current filter of the images (ex: Sloan_z)
    log_fname (str) : full path to log file
    nccd (int) : 2 ccds for OSIRIS
    mbias (array) : master bias [ccd1, ccd2]
    mflat (array) : master flat [ccd1, ccd2]
    mask_ccd (array) : master bad pixel mask [ccd1, ccd2]
    gain (float) : instrument constant
    rdnoise (float) : instrument constant
    calib_path (str) : directory to write to bias/flat corrected images to
    bpm_path (str) : directory to write to bias/flat/bpm corrected images to
    root (str) : contains object name

    Returns
    -------
    all_sci_proc (list): list of images after applying bias, flat, bpm and trimming
    all_headers (list) : list of the headers from each of the 3 extensions
    all_filts (list)   : list of True/False; True if the filter matches the
                         current reduction filter
    """
    all_sci_proc = []
    all_headers = []
    all_filts = []

    # loop over each file in fsci or fstd
    for j, f in enumerate(obj):

        with fits.open(f) as fits_open:
            # Note: because images have been trimmed, the wcs information
            # (crpix) is no longer correct; subtract off the number of pixels
            # trimmed to get close to the original wcs
            hdr = fits_open[0].header

            hdr1 = fits_open[1].header
            trim1 = hdr1['TRIMSEC'].split(':')[0][1:]
            hdr1['CRPIX1'] = float(hdr1['CRPIX1'])-float(trim1)
            trimsec1 = gtcsetup.update_trimsec(hdr1['TRIMSEC'])

            hdr2 = fits_open[2].header
            trim2 = hdr2['TRIMSEC'].split(':')[0][1:]
            hdr2['CRPIX1'] = float(hdr2['CRPIX1'])-float(trim2)
            trimsec2 = hdr2['TRIMSEC']

            trimsec = [trimsec1, trimsec2]
            hdrs = [hdr1, hdr2]
            hdr_filt = hdr['FILTER2']

            all_headers.append([hdr, hdr1, hdr2])

            # # get initial wcs info from header
            # if j == 0:
            #     # NOTE: somehow updating hdr1 and 2 above also updates
            #     # this header too...
            #     wcs_ref = [wcs.WCS(fits_open[k+1]) for k in range(nccd)]

        # Skip over this image if its filter does not match the
        # current one
        if hdr_filt != filt:
            all_filts.append(False)
            continue
        else:
            all_filts.append(True)

        gtcsetup.print_both(log_fname, '    Working on file', f)

        # get raw frame
        sci_raw = [CCDData.read(f, hdu=k+1, unit='adu')
                   for k in range(nccd)]

        # #######################################################
        # ################# APPLY CALIBRATIONS ##################
        # #######################################################

        # apply bias, flat, bpm corrections
        # again noting that Nora had the oscan, oscan_model and trim;
        # I didn't change this
        sci_proc_init = [ccdproc.ccd_process(
            x, oscan='[3:22,:]',
            oscan_model=models.Chebyshev1D(3),
            trim=trimsec[k],
            master_bias=mbias[k],
            master_flat=mflat[k],
            bad_pixel_mask=np.array(mask_ccd[k]),
            # ^I don't remember why I formatted the bpm this way
            gain=gain*u.electron/u.adu,
            readnoise=rdnoise*u.electron,
            gain_corrected=True)
            for k, x in enumerate(sci_raw)]

        # apply bias and flat but not bad pixel mask just for
        # diagnostic purposes
        sci_proc_no_bpm = [ccdproc.ccd_process(
            x, oscan='[3:22,:]',
            oscan_model=models.Chebyshev1D(3),
            trim=trimsec[k],
            master_bias=mbias[k],
            master_flat=mflat[k],
            bad_pixel_mask=None,
            gain=gain*u.electron/u.adu,
            readnoise=rdnoise*u.electron,
            gain_corrected=True)
            for k, x in enumerate(sci_raw)]

        # Write out these files to individual images
        # NOTE: this function creates a 'master' combined image
        # ('combining' 1 image with itself)
        # otherwise the bad pixel mask will not be applied
        gtcsetup.write_one_image(hdr, hdrs, sci_proc_no_bpm, calib_path,
                                 f[:-5]+'_calib.fits', root, filt,
                                 log_fname, nccd)
        gtcsetup.write_one_image(hdr, hdrs, sci_proc_init, bpm_path,
                                 f[:-5]+'_bpm.fits', root, filt,
                                 log_fname, nccd)

        all_sci_proc.append(sci_proc_init)

    return all_sci_proc, all_headers, np.array(all_filts)


def do_crmask(log_fname, all_sci_calib, mask_ccd, gain, rdnoise, saturation,
              obj, all_headers, root, filt, nccd, crmask_path):
    """Detect cosmic rays and remove them from the science frames.

    Parameters
    ----------
    log_fname (str) : full path to log file
    all_sci_calib (list): list of images after applying bias, flat, bpm and trimming
    mask_ccd (array) : master bad pixel mask [ccd1, ccd2]
    gain (float) : instrument constant
    rdnoise (float) : instrument constant
    saturation (float) : instrument constant
    obj (list) : list of file names
    all_headers (list) : list of the headers from each of the 3 extensions
    root (str) : contains object name
    filt (str) : current filter of the images (ex: Sloan_z)
    nccd (int) : 2 ccds for OSIRIS
    crmask_path (str) : directory to write to cosmic ray corrected images to

    Returns
    -------
    all_sci_proc (list): list of images after removing cosmic rays
    """
    gtcsetup.print_both(log_fname, '    Removing cosmic rays')
    all_sci_proc = []

    for j, f in enumerate(obj):
        sci_proc_init = all_sci_calib[j]
        cleaned_sci = []
        cr_mask = []
        for ccd in range(len(sci_proc_init)):
            # cr_cleaned = ccdproc.cosmicray_lacosmic(x, gain=1.0,
            # sig_clip=5)

            # NOTE: some large saturated stars may be flagged as crs
            # but I'm going to assume we don't care about them
            test_mask, _clean = detect_cosmics(np.array(
                sci_proc_init[ccd].data),
                inmask=np.array(mask_ccd[ccd].data),  # None,
                sigclip=4.5,  # lower values flag more pixels as
                # cosmic rays
                sigfrac=0.3,
                objlim=5.0,  # increase if CENTERS of stars flagged
                # as cosmic rays
                # NOTE: Tried 5000, it didn't find as many saturated
                # stars but it also barely found any cosmic rays
                gain=gain, readnoise=rdnoise, satlevel=saturation,
                niter=4, sepmed=True,
                cleantype='meanmask', fsmode='median',
                psfmodel='gauss', psffwhm=2.5, psfsize=7,
                psfk=None, psfbeta=4.765,  # these are all defaults
                verbose=True)
            gtcsetup.print_both(log_fname, '------')
            cleaned_sci.append(_clean)
            cr_mask.append(np.multiply(test_mask, 1))

        # assign the cosmic ray cleaned image to a CCDData object
        sci_proc = copy.deepcopy(sci_proc_init)
        sci_proc[0].data = cleaned_sci[0]
        sci_proc[1].data = cleaned_sci[1]

        all_sci_proc.append(sci_proc)

        # create the cosmic ray mask as a CCDData object
        cr_mask_ccd = [CCDData(data=x.astype('uint8'),
                               unit=u.dimensionless_unscaled)
                       for x in cr_mask]

        # Write out the cosmic ray mask to its diagnostic folder
        gtcsetup.write_one_image(all_headers[j][0], all_headers[j][1:],
                                 cr_mask_ccd, crmask_path,
                                 f[:-5]+'_CRmask.fits', root, filt,
                                 log_fname, nccd)

        # Write out the image with the cosmic ray mask applied to its
        # diagnostic folder
        gtcsetup.write_one_image(all_headers[j][0], all_headers[j][1:],
                                 sci_proc, crmask_path,
                                 f[:-5]+'_CRmask_applied.fits', root,
                                 filt, log_fname, nccd)

    return all_sci_proc


#############################################################################


def do_bkg_sky_subtraction(all_sci_proc, all_headers, log_fname,
                           skymap_path, fnames, root, filt):
    """Subtract background and sky image.

    Parameters
    ----------
    all_sci_proc (list) : list of all science images; they have have been
        bias/flat/bpm/crmask corrected
    all_headers (list) : list of headers for whole fits file, ccd 1 and ccd 2
    log_fname (str) : full path to log file
    skymap_path (str) : path to the diagnostic folder containing intermediate files
    fnames (str) : name of the object/target
    root (str) : contains object name
    filt (str) : filter of the images

    Returns
    -------
    sci_final (arr) : list of all images with both ccds after sky and bkg subtraction
    sci_skymap (arr) : sky image with both ccds (median combined science image
                                                 with stars masked out)
    """
    sci_skymap = []
    all_bkg_sub_images = []
    all_bkg_images = []

    # Focus on one ccd at a time
    for ccd in range(len(all_sci_proc[0])):
        gtcsetup.print_both(log_fname, '     CCD', ccd+1)

        # do one ccd for all images
        all_masked_bkg_sub_images = []
        bkg_sub_images = []
        bkg_images = []
        for image in range(len(all_sci_proc)):

            gtcsetup.print_both(log_fname, '         Image', image+1)

            # use source extractor to calculate bakground and subtract it off
            try:
                bkg = sep.Background(all_sci_proc[image][ccd].data)
            except ValueError:
                bkg = sep.Background(
                    all_sci_proc[image][ccd].data.byteswap().newbyteorder())
            bkg_image1 = bkg.back()
            bkg_image = CCDData(bkg_image1, unit='electron')
            bkg_images.append(bkg_image)

            # write out the background image to the diagnostic folder
            new_hdul = fits.HDUList()
            new_hdul.append(fits.ImageHDU(header=all_headers[image][0]))
            new_hdul.append(fits.ImageHDU(data=bkg_image.data,
                                          header=all_headers[image][ccd+1]))
            new_hdul.writeto(skymap_path+root+fnames[image].split('/')[-1][:-5] +
                             '_'+filt+'_bkg_img_ccd'+str(ccd+1)+'.fits',
                             overwrite=True)

            # subtract the background from the current image and append it
            # to the list
            bkg_sub_img = all_sci_proc[image][ccd].subtract(
                bkg_image, propagate_uncertainties=True,
                handle_meta='first_found')
            bkg_sub_images.append(bkg_sub_img)

            # np.save(skymap_path+'test_data'+str(image)+'_ccd'+str(ccd+1) +
            #         '.npy', all_sci_proc[image][ccd].data)
            # np.save(skymap_path+'test_mask'+str(image)+'_ccd'+str(ccd+1) +
            #         '.npy', all_sci_proc[image][ccd].mask)

            # Start creating sky image
            # 2 iterations of sigma clipping to remove stars; can uncomment
            # below to add a third iteration
            mask1 = make_source_mask(bkg_sub_img.data, nsigma=3, npixels=3,
                                     mask=all_sci_proc[image][ccd].mask)

            # new_hdul = fits.HDUList()
            # new_hdul.append(fits.ImageHDU(header=all_headers[image][0]))
            # new_hdul.append(fits.ImageHDU(data=mask1.astype(int),
            #                               header=all_headers[image][ccd+1]))
            # new_hdul.writeto(skymap_path+root+fnames[image].split('/')[-1][:-5] +
            #                  '_'+filt+'_mask1_ccd'+str(ccd+1)+'.fits',
            #                  overwrite=True)

            intermediate = CCDData(
                bkg_sub_img.data, unit='electron', mask=mask1)

            mask2 = make_source_mask(intermediate, nsigma=3, npixels=3)

            # new_hdul = fits.HDUList()
            # new_hdul.append(fits.ImageHDU(header=all_headers[image][0]))
            # new_hdul.append(fits.ImageHDU(data=mask2.astype(int),
            #                               header=all_headers[image][ccd+1]))
            # new_hdul.writeto(skymap_path+root+fnames[image].split('/')[-1][:-5] +
            #                  '_'+filt+'_mask2_ccd'+str(ccd+1)+'.fits',
            #                  overwrite=True)

            # intermediate2 = CCDData(intermediate.data, unit='electron',
            #                         mask=mask2)

            # mask3 = make_source_mask(intermediate2, nsigma=3, npixels=3)

            # new_hdul = fits.HDUList()
            # new_hdul.append(fits.ImageHDU(header=hdr))
            # new_hdul.append(fits.ImageHDU(data=mask3.astype(int), header=hdrs[0]))
            # new_hdul.writeto(skymap_path+root+fnames[image].split('/')[-1][:-5] +
            #                  '_'+filt+'_mask3_ccd'+str(ccd+1)+'.fits',
            #                  overwrite=True)

            # combine all masks; don't actually need to sum here, mask 3 = mask
            mask = mask1+mask2  # +mask3

            # Write out star mask image to the diagnostic file
            new_hdul = fits.HDUList()
            new_hdul.append(fits.ImageHDU(header=all_headers[image][0]))
            new_hdul.append(fits.ImageHDU(data=mask.astype(int),
                                          header=all_headers[image][ccd+1]))
            new_hdul.writeto(skymap_path+root+fnames[image].split('/')[-1][:-5] +
                             '_'+filt+'_finalmask_ccd'+str(ccd+1)+'.fits',
                             overwrite=True)

            # mask out the stars in the background subtracted image to get sky image
            bkg_sub_img_masked = CCDData(bkg_sub_img.data, unit='electron',
                                         mask=mask)

            all_masked_bkg_sub_images.append(bkg_sub_img_masked)

        # median combine the images where stars have been masked out
        # and background has been subtracted (sky image) for 1 ccd

        masked_median_img = ccdproc.combine(all_masked_bkg_sub_images,
                                            method='median')

        sci_skymap.append(masked_median_img)
        all_bkg_sub_images.append(bkg_sub_images)
        all_bkg_images.append(bkg_images)

    new_bkg_sub_images = [[all_bkg_sub_images[0][i], all_bkg_sub_images[1][i]]
                          for i in range(len(all_bkg_sub_images[0]))]

    # subtract the sky from all images
    sky_sub = [[x.subtract(sci_skymap[k], propagate_uncertainties=True,
                           handle_meta='first_found')
                for k, x in enumerate(j)] for j in new_bkg_sub_images]

    # Divide sky subtracted image by the exposure time
    sci_final = [[x.divide(all_headers[i][0]['EXPTIME']*u.second,
                           propagate_uncertainties=True,
                           handle_meta='first_found')
                  for i, x in enumerate(imgs)] for imgs in sky_sub]

    return sci_final, sci_skymap


############################################################################


def do_bkg_subtraction(all_sci_proc, all_headers, log_fname, skymap_path):
    """Only do background subtraction.

    Parameters
    ----------
    all_sci_proc (list) : list of all science images; they have have been
        bias/flat/bpm/crmask corrected
    hdr (?) : header for whole fits file
    log_fname (str) : full path to log file
    skymap_path (str) : path to the diagnostic folder containing intermediate files
    fnames (str) : name of the object/target

    Returns
    -------
    sci_final (arr) : list of all images with both ccds after bkg subtraction
    """
    all_bkg_sub_images = []
    # Focus on one ccd at a time
    for ccd in range(len(all_sci_proc[0])):
        gtcsetup.print_both(log_fname, '     CCD', ccd+1)

        # Do one ccd for all images
        bkg_sub_images = []
        for image in range(len(all_sci_proc)):

            gtcsetup.print_both(log_fname, '         Image', image+1)

            # Use source extractor to calculate bakground and subtract it off
            try:
                bkg = sep.Background(all_sci_proc[image][ccd].data)
            except ValueError:
                bkg = sep.Background(
                    all_sci_proc[image][ccd].data.byteswap().newbyteorder())

            bkg_image1 = bkg.back()
            bkg_image = CCDData(bkg_image1, unit='electron')

            bkg_sub_img = all_sci_proc[image][ccd].subtract(
                bkg_image, propagate_uncertainties=True,
                handle_meta='first_found')
            bkg_sub_images.append(bkg_sub_img)

        all_bkg_sub_images.append(bkg_sub_images)

    new_bkg_sub_images = [[all_bkg_sub_images[0][i], all_bkg_sub_images[1][i]]
                          for i in range(len(all_bkg_sub_images[0]))]

    # Divide sky subtracted image by the exposure time
    sci_final = [[x.divide(all_headers[i][0]['EXPTIME']*u.second,
                           propagate_uncertainties=True,
                           handle_meta='first_found')
                  for i, x in enumerate(imgs)] for imgs in new_bkg_sub_images]

    return sci_final


# ##########################################################
# ################ Astrometry functions ####################
# ##########################################################


def tellme(s):
    """Print text and change the title of a plot to that text.

    Parameters
    ----------
    s (str) : text to print/set as title

    Returns
    -------
    None
    """
    print(s)
    plt.suptitle(s, fontsize=16)
    plt.draw()


def plot_images(img, ref_img, w, ref_wcs, sky, tbl_crds, obj_ra,
                obj_dec, n_stars, pts, select, patch_radius=patch_radius):
    """Plot the OSIRIS image and a reference SDSS image w/detected stars in each.

    Parameters
    ----------
    img (array) : stacked science image; 1 ccd
    ref_img (array) : SDSS image for reference
    w (?) : wcs information for the science image
    ref_wcs (?) : wcs information for the reference image
    sky (SkyCoord) : coordinates of the center of the ccd
    tbl_crds (list) : list of Gaia star coordinates
    obj_ra (list) : list of RAs of detected OSIRIS stars
    obj_dec (list) : list of Decs of detected OSIRIS stars
    n_stars (int) : number of stars to display on the plot
    pts (list) : pairs of points corresponding to stars the user selected
    select (str) : either gaia or osiris to designate which objects the user
                   selected

    Returns
    -------
    fig (figure) : figure containing the OSIRIS image, SDSS image, and detected
                    stars overplotted
    all_patches (list) : list of all patch objects (circle with x,y,r) for either
                        Gaia or OSIRIS depending on the variable select
    """
    font_prop = font_manager.FontProperties(size=14)

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(121, projection=w)
    ax.set_xlabel("RA", fontproperties=font_prop)
    ax.set_ylabel("Dec", fontproperties=font_prop)

    ax.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))

    if select == 'gaia':
        for j in range(len(tbl_crds.ra)):
            x, y = w.all_world2pix(tbl_crds.ra.deg[j], tbl_crds.dec.deg[j], 0)
            ax.add_patch(patches.Circle((x, y), radius=patch_radius, ec='b',
                                        fc='none', lw=2, alpha=0.8,
                                        label='Reference: RA = %f, Dec = %f' %
                                        (tbl_crds.ra.deg[j], tbl_crds.dec.deg[j]),
                                        picker=True))
    else:
        ax.plot(tbl_crds.ra, tbl_crds.dec, 'o', transform=ax.get_transform('fk5'),
                mec='b', mfc='none', label='GAIA')

    if select == 'osiris':
        for j in range(len(obj_ra[:n_stars])):
            x, y = w.all_world2pix(obj_ra[j], obj_dec[j], 0)
            ax.add_patch(patches.Circle((x, y), radius=patch_radius, ec='k',
                                        fc='none', lw=2, alpha=0.8,
                                        label='Position: x = %f, y = %f' %
                                        (x, y), picker=True))
    else:
        ax.plot(obj_ra[:n_stars], obj_dec[:n_stars], 'o',
                transform=ax.get_transform('fk5'), mec='k', mfc='none',
                label='OSIRIS')

    for i, p in enumerate(pts):
        ax.text(p[0], p[1], str(i+1), c='magenta', fontsize='large')

    ax.set_title('Your OSIRIS image', fontproperties=font_prop)
    all_patches = ax.patches

    #########################

    ax2 = fig.add_subplot(122, projection=ref_wcs)
    ax2.set_xlabel("RA", fontproperties=font_prop)
    ax2.set_ylabel("Dec", fontproperties=font_prop)

    ref = ref_img[0].data
    ax2.imshow(ref, origin='lower', cmap='gray', vmin=np.nanmean(
        ref)-np.nanstd(ref), vmax=np.nanmean(ref)+np.nanstd(ref))

    ax2.plot(tbl_crds.ra, tbl_crds.dec, 'o', transform=ax2.get_transform('fk5'),
             mec='b', mfc='none', label='GAIA')
    ax2.plot(obj_ra[:n_stars], obj_dec[:n_stars], 'o',
             transform=ax2.get_transform('fk5'), mec='k', mfc='none',
             label='OSIRIS')

    ax2.set_title('SDSS image for reference',
                  fontproperties=font_prop)

    plt.legend()

    return fig, all_patches


def select_points(img, select_n_stars, color):
    """Set up interactive point selection.

    Parameters
    ----------
    img (array) : stacked science image; 1 ccd
    select_n_stars (int) : minimum number of stars the user must select
    color (str) : color of cirlce the user must select

    Returns
    -------
    pts (list) : pairs of points corresponding to stars the user selected
    """
    while True:
        pts = []
        while len(pts) < select_n_stars:
            tellme('Right click inside a '+color+' circle to select at least ' +
                   str(select_n_stars) +
                   ' stars on the LEFT plot. \n'
                   ' Use the magnifying glass to zoom in.\n'
                   'To remove a selected star, '
                   'click the back button on your mouse. \n'
                   'Press enter when done.')

            pts = np.asarray(plt.ginput(-1, timeout=-1,
                                        mouse_add=MouseButton(3),
                                        mouse_pop=MouseButton(8)))
            print(pts)

            if len(pts) < select_n_stars:
                tellme('Too few points, starting over')
                time.sleep(1)  # Wait a second

        check_fig = plt.figure()
        plt.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
            img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))
        plt.plot(np.array(pts).T[0], np.array(pts).T[1], 'o', c='orange')

        tellme('Happy? Key click for yes, mouse click for no')

        if plt.waitforbuttonpress():
            break
            check_fig.close()

    return pts


def select_yes_no(img, select_n_stars):
    """Set up interactive point selection.

    Parameters
    ----------
    img (array) : stacked science image; 1 ccd
    select_n_stars (int) : minimum number of stars the user must select

    Returns
    -------
    pts (list) : one pairs of points corresponding to the user's selection
    """
    while True:
        pts = []
        while len(pts) < select_n_stars:
            tellme('Right click to select yes or no\n'
                   'To remove your selection, '
                   'click the back button on your mouse. \n'
                   'Press enter when done.')

            pts = np.asarray(plt.ginput(-1, timeout=-1,
                                        mouse_add=MouseButton(3),
                                        mouse_pop=MouseButton(8)))
            print(pts)

            if len(pts) < select_n_stars:
                tellme('Too few points, starting over')
                time.sleep(1)  # Wait a second

        tellme('Happy? Key click for yes, mouse click for no')

        if plt.waitforbuttonpress():
            break

    return pts


def plot_check_img(img, final_wcs, tbl_crds):
    """Plot the OSIRIS image w/Gaia stars to check updated astrometry.

    Parameters
    ----------
    img (array) : stacked science image with corrected astrometry; 1 ccd
    final_wcs (?) : updated wcs information for the science image
    tbl_crds (list) : list of Gaia star coordinates

    Returns
    -------
    fig3 (figure) : figure containing the OSIRIS image in x/y and ra/dec, and
                    detected stars overplotted
    all_patches (list) : list of all 2 patch objects (circle with x,y,r) for
                        yes/no
    """
    font_prop = font_manager.FontProperties(size=14)

    fig3 = plt.figure(figsize=(12, 10))

    ax = fig3.add_subplot(121)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    ax.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))

    x1 = 100
    y = img.shape[0]-100
    x2 = img.shape[1]-100

    ax.text(x1, y, 'YES', c='g')
    ax.text(x2, y, 'NO', c='g')

    ax.add_patch(patches.Circle((x1, y), radius=patch_radius*2, ec='k', fc='none',
                                lw=2, alpha=0.8,
                                label='Position: x = %f, y = %f' %
                                (x1, y), picker=True))
    ax.add_patch(patches.Circle((x2, y), radius=patch_radius*2, ec='k', fc='none',
                                lw=2, alpha=0.8,
                                label='Position: x = %f, y = %f' %
                                (x2, y), picker=True))
    all_patches = ax.patches

    #########################

    ax2 = fig3.add_subplot(122, projection=final_wcs)
    ax2.set_xlabel("RA", fontproperties=font_prop)
    ax2.set_ylabel("Dec", fontproperties=font_prop)

    ax2.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))
    ax2.plot(tbl_crds.ra, tbl_crds.dec, '*', transform=ax2.get_transform('fk5'),
             mec='b', mfc='none', label='GAIA')

    # ax2.set_title('SDSS zband image for reference',
    #               fontproperties=font_prop)

    # plt.legend()

    return fig3, all_patches


def get_gaia_img_cat(img, pix_limit, w, sky, log_fname):
    """Get a catalog of stars that are within the limits of the OSIRSI image.

    Parameters
    ----------
    img (array) : stacked science image; 1 ccd
    pix_limit (int) :
    w (?) : wcs information for the science image
    sky (SkyCoord) : coordinates of the center of the ccd
    log_fname (?) : opened log file

    Returns
    -------
    tbl_crds (list) : list of Gaia star coordinates
    """
    # Getting corners of the ccd to select only stars that are actually in the
    # image
    pairs = [(pix_limit, pix_limit), (pix_limit, img.shape[0]-pix_limit),
             (img.shape[1]-pix_limit, pix_limit), (img.shape[1]-pix_limit,
                                                   img.shape[0]-pix_limit)]
    polygon = []
    for pair in pairs:
        sky_pair = w.all_pix2world(pair[0], pair[1], 0)
        polygon.append((sky_pair[0], sky_pair[1]))

    # Get gaia comparison catalog of stars using corners of OSIRIS image as region
    gaia = Irsa.query_region(SkyCoord(sky.ra, sky.dec, frame='fk5'),
                             catalog="gaia_dr2_source", spatial="Polygon",
                             polygon=polygon)
    gtcsetup.print_both(log_fname, 'Found ', len(gaia), ' stars in the GAIA catalog')

    #######################################
    gaia_flux_init = np.array(gaia['phot_rp_mean_flux'])

    # Sort gaia stars based on flux from ** to ** (bright/faint)
    ind_init = np.argsort(-gaia_flux_init)
    # print(ind_init)
    # print(gaia_flux_init[ind_init])

    # Remove all nan flux entries
    gaia_ind = []
    for i in ind_init:
        if not np.isnan(gaia_flux_init[i]):  # == False:
            gaia_ind.append(i)
    gaia_ind = np.array(gaia_ind)

    # gaia_fluxes = gaia_flux_init[gaia_ind]
    # print(gaia_fluxes)

    tbl_crds = coordinates.SkyCoord(gaia['ra'][gaia_ind], gaia['dec'][gaia_ind],
                                    unit=(u.deg, u.deg), frame='fk5')

    return tbl_crds


def get_osiris_obj_cat(img, fname, w, pix_limit, log_fname):
    """Detect sources in the OSIRIS image using source extractor.

    Detect sources
    Eliminate objects with source extractor flags (saturated etc)
    Elliminate objects that are too close to the edge
    Eliminate objects that are too elliptical
    Convert x/y to RA/Dec

    Parameters
    ----------
    img (array) : stacked science image; 1 ccd
    pix_limit (int) : value in pixels that determines how far away from the edge
                      of the image a star must be
    w (?) : wcs information for the science image
    sky (SkyCoord) : coordinates of the center of the ccd
    log_fname (?) : opened log file

    Returns
    -------
    obj_ra (list) : list of RAs of detected OSIRIS stars
    obj_dec (list) : list of Decs of detected OSIRIS stars
    final_x (list) : list of x coordinates in pixels of detect OSIRIS stars
    final_y (list) : list of y coordinates in pixels of detect OSIRIS stars
    """
    # detect sources in OSIRIS image using source extractor
    try:
        bkg = sep.Background(img)
    except ValueError:
        bkg = sep.Background(img.byteswap().newbyteorder())
    data_sub = img  # -bkg #bkg already subtracted
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)
    gtcsetup.print_both(log_fname, 'Found ', len(objects), 'objects in ', fname)

    objects = objects[objects['flag'] == 0]
    gtcsetup.print_both(log_fname, len(objects),
                        'objects remain after removing those with flags')

    # check how circular/ellipsoidal;
    axes_diff = np.array([round(abs(objects['a'][i]-objects['b'][i]), 3)
                          for i in range(len(objects))])
    # ellipse_lim = 1.
    sorted_axes = np.array([x for _, x in sorted(
        zip(objects['flux'], axes_diff))])[::-1]
    # print(len(good))

    # Sort x and y based on the flux; we want to take the n brightest stars
    # sorted_flux = np.array(sorted(objects['flux'])[::-1])
    sorted_x = np.array([x for _, x in sorted(
        zip(objects['flux'], objects['x']))])[::-1]
    sorted_y = np.array([x for _, x in sorted(
        zip(objects['flux'], objects['y']))])[::-1]

    # Eliminate stars that are close to the edge AND too elliptical
    # pix_limit is how much of the edge of the ccd to exclude
    final_x = []
    final_y = []
    for i in range(len(sorted_x)):
        x = sorted_x[i]
        y = sorted_y[i]
        e = sorted_axes[i]
        if x > pix_limit and (x < img.shape[1]-pix_limit) and (
                y > pix_limit) and (y < img.shape[0]-pix_limit):
            if e < ellipse_lim:
                final_x.append(x)
                final_y.append(y)
    final_x = np.array(final_x)
    final_y = np.array(final_y)

    # NOTE: Pixel to world conversion differs from DS9 because DS9 does not take
    # into account SIP distortions while all_world2pix() does
    # (however, wcs_world2pix() does not).

    # Convert detected objects in OSIRIS image from x,y to ra,dec
    obj_ra = []
    obj_dec = []
    for i in range(len(final_x)):
        pos = w.all_pix2world(final_x[i], final_y[i], 0)
        obj_ra.append(pos[0])
        obj_dec.append(pos[1])
    obj_ra = np.array(obj_ra)
    obj_dec = np.array(obj_dec)

    return obj_ra, obj_dec, final_x, final_y


def in_circle(center_x, center_y, radius, x, y):
    """Determine if an (x,y) point is within a circle.

    Parameters
    ----------
    center_x (float) : x coordinate of the circle in pixels
    center_y (float) : y coordinate of the circle in pixels
    radius (float) : radius of the circle in pixels
    x (float) : x coordinate of the point in pixels
    y (float) : y coordinate of the circle in pixels

    Returns
    -------
    True or False if the squared distance is less than the radius of the circle
    squared
    """
    square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
    return square_dist <= radius ** 2


def cross_match(all_patches, selected_pts, radius):
    """Cross match the user selected points w/the coordinates of detected sources.

    Parameters
    ----------
    all_patches (list)
    selected_pts (list)
    radius (float) : radius of the circle in pixels

    Returns
    -------
    final_pts (array) :

    """
    final_pts = []
    for pt in selected_pts:
        # print('Pt', pt)
        for patch in all_patches:
            p_cen = patch.center
            incirc = in_circle(p_cen[0], p_cen[1], radius, pt[0], pt[1])
            if incirc:
                # print('found')
                final_pts.append([patch.center[0], patch.center[1]])
                break
    return np.array(final_pts)


def final_check(img, final_wcs, tbl_crds, log_fname):
    """Plot the OSIRIS image w/Gaia stars to check updated astrometry.

    Parameters
    ----------
    img (array) : stacked science image with corrected astrometry; 1 ccd
    final_wcs (?) : updated wcs information for the science image
    tbl_crds (list) : list of Gaia star coordinates
    log_fname (?) : opened log file

    Returns
    -------
    yes (astrometry is good), OR no (astrometry is not good) OR
    final_check (if user did not select yes or no)

    """
    # Plot the OSIRIS image with the new wcs and overplot Gaia points to make
    # sure astrometry is good; user selects yes or no
    fig3, all_patches = plot_check_img(img, final_wcs, tbl_crds)
    text = 'Right click inside a GREEN circle to select yes or no.'
    tellme(text)
    plt.waitforbuttonpress()

    selected_pts = select_points(img, 1, 'GREEN')
    # print(selected_pts)

    # CROSS MATCH SELECTED POINTS WITH ALL PATCHES
    final_pts = cross_match(all_patches, selected_pts, patch_radius*2)

    plt.close(fig3)

    if len(final_pts) == 0:
        gtcsetup.print_both(log_fname, 'Invalid selection; try again')
        return final_check(img, final_wcs, tbl_crds)

    elif final_pts[-1][0] == 100:
        gtcsetup.print_both(log_fname, 'You selected yes')
        return 'yes'

    elif final_pts[-1][0] == img.shape[1]-100:
        gtcsetup.print_both(log_fname, 'You selected no')
        return 'no'


def do_interactive_astrometry(final_name, fname, full_filt, ccd, astrom_path,
                              log_fname):
    """Allow the user to select stars in an image to perform astrometry on.

    This script will display the current OSIRIS image next to an SDSS image
    with markers on the brightest OSIRIS target stars and markers on GAIA reference
    stars.
    The user will then select at least 6 stars in the OSIRIS image that have
    good GAIA matches. Then the same plot will show and the user will select
    the same stars in the GAIA catalog in the same order. Numbers are shown
    next to the selected stars so you do not have to remember the order.
    Once you select target and corresponding reference stars, the astrometry
    correction is performed and applied. The corrected image is then displayed
    so that the user can inspect it and decide if it is okay or if they want to
    redo it by selecting yes or no in the same way the stars were selected.

    Note: I am using GAIA ra/dec because it has better accuracy, but I am
    using SDSS image because I could not find a GAIA image

    Note: the stars you choose will usually have good astrometry correction,
    but other stars may not. Therefore, it seems best to choose as many stars
    as possible

    Note: Things the user may want to change at the very top
        - pix_limit
        - n_stars
        - select_n_stars
        - sip_deg

    Parameters
    ----------
    final_name (str) : the name to write the corrected image to
    fname (str) : the name of the image to perform astrometry on (will be 1 ccd)
    full_filt (str) : the current filter (ex: Sloan_z)
    ccd (str) : 1 or 2
    astrom_path (str) : directory to write to astrometry related images to
    log_fname (?) : log file to append print statements to

    Returns
    -------
    None
    """
    gtcsetup.print_both(log_fname, 'Doing interactive astrometry')

    filt = full_filt.split('_')[-1]  # 'z'

    img = fits.open(fname)[1].data.byteswap().newbyteorder()
    hdr = fits.open(fname)[1].header
    w = WCS(hdr)

    # Get the approximate center of the current image in ra/dec
    sky = w.pixel_to_world(img.shape[1]/2, img.shape[0]/2)
    gtcsetup.print_both(log_fname, 'Center of ccd:', sky.ra.deg, sky.dec.deg)

    # pix_limit = 0
    # Detect sources in the OSIRIS image and convert x/y to ra/dec using the
    # original wcs information
    obj_ra, obj_dec, final_x, final_y = get_osiris_obj_cat(img, fname, w,
                                                           pix_limit, log_fname)

    # Get ra/dec of gaia stars that cover the OSIRIS image
    tbl_crds = get_gaia_img_cat(img, pix_limit, w, sky, log_fname)

    # Set how many OSIRIS stars should be displayed for astrometry
    # n_stars = 40
    n_stars = min(n_stars_init, len(obj_ra))

    # Get sdss image to compare to OSIRIS image
    ref_img = SkyView.get_images(position=str(sky.ra.deg)+','+str(sky.dec.deg),
                                 # survey='SDSSz',
                                 survey='SDSS'+filt,
                                 coordinates='J2000',
                                 pixels=str(img.shape[1])+','+str(img.shape[0]))
    if len(ref_img) == 0:
        ref_img = SkyView.get_images(position=str(sky.ra.deg)+','+str(sky.dec.deg),
                                     survey='DSS',
                                     coordinates='J2000',
                                     pixels=str(img.shape[1])+','+str(img.shape[0]))
    print(ref_img)
    ref_img = ref_img[0]
    # Get wcs information from SDSS image
    ref_wcs = wcs.WCS(ref_img[0].header)

    # Plot the OSIRIS image and the SDSS image with detected objects and
    # Gaia sources overlayed; select OSIRIS stars
    fig, all_patches = plot_images(img, ref_img, w, ref_wcs, sky, tbl_crds,
                                   obj_ra, obj_dec, n_stars, [], 'osiris')

    # Wait for user to select at least n points
    # select_n_stars = 6
    text = 'You will select '+str(select_n_stars)+' BLACK (OSIRIS) stars first.\n' \
        'Then you will select the corresponding Gaia stars (in blue) in the' \
        'same order. \n To remove a selected star, click the back button' \
        ' on your mouse. Click anywhere to begin.'
    tellme(text)
    plt.waitforbuttonpress()
    # Points selected by user
    pts = select_points(img, select_n_stars, 'BLACK')

    # Cross match where the user selected with the actual positions of the stars
    source_xy = cross_match(all_patches, pts, cross_match_r)

    plt.close(fig)

    # Plot the OSIRIS image and the SDSS image with detected objects and
    # Gaia sources overlayed; select GAIA stars
    fig2, all_patches2 = plot_images(img, ref_img, w, ref_wcs, sky, tbl_crds,
                                     obj_ra, obj_dec, n_stars, pts, 'gaia')

    # Wait for user to select at least n points
    tellme('You will now select 6 BLUE (gaia) stars in the same order'
           ' (pink numbers). Click anywhere to begin.')
    plt.waitforbuttonpress()
    pts2 = select_points(img, select_n_stars, 'BLUE')

    # Cross match where the user selected with the actual positions of the stars
    final_pts2 = cross_match(all_patches2, pts2, cross_match_r)
    plt.close(fig2)

    # Convert x/y to ra/dec using OSIRIS image wcs info
    ref_radec = np.array([w.all_pix2world(final_pts2[i][0], final_pts2[i][1], 0)
                          for i in range(len(final_pts2))])

    # Create a SkyCoord object for the GAIA stars to be inputted on next line
    coords = SkyCoord(np.array(ref_radec), frame='fk5', unit='deg')

    # Calculate new wcs
    # sip_deg = 2
    new_wcs = fit_wcs_from_points(xy=source_xy.T, world_coords=coords,
                                  projection='TAN', sip_degree=sip_deg)

    # proj_point=SkyCoord( sky.ra.deg, sky.dec.deg, unit='deg'),
    gtcsetup.print_both(log_fname, 'Old wcs:')
    gtcsetup.print_both(log_fname, w)
    gtcsetup.print_both(log_fname, '-------')
    gtcsetup.print_both(log_fname, 'New wcs:')
    gtcsetup.print_both(log_fname, new_wcs)

    # Update the header with the new wcs
    hdr.update(new_wcs.to_header(relax=True))
    final_wcs = WCS(hdr)

    answer = final_check(img, final_wcs, tbl_crds, log_fname)

    calc_astrometry_accuracy(img, final_wcs, final_x, final_y, tbl_crds, astrom_path,
                             ccd, final_name, log_fname)

    if answer == 'no':
        gtcsetup.print_both(log_fname, 'REDOING ASTROMETRY ')
        do_interactive_astrometry(final_name, fname, log_fname)
    else:
        gtcsetup.print_both(log_fname, 'Done with astrometry')

        new_hdul = fits.HDUList()
        new_hdul.append(fits.ImageHDU(header=hdr))
        new_hdul.append(fits.ImageHDU(data=img, header=hdr))
        new_hdul.writeto(final_name, overwrite=True)
        gtcsetup.print_both(log_fname, 'Writing final astrometry corrected image to',
                            final_name)


def calc_astrometry_accuracy(img, final_wcs, final_x, final_y, tbl_crds, astrom_path,
                             ccd, final_name, log_fname):
    """Compare OSIRIS stars with corrected astrometry to Gaia stars.

    On the left show the OSIRIS image with all Gaia stars overplotted.
    On the right, show a histogram of the distances between matched OSIRIS/Gaia
        stars that are

    Parameters
    ----------
    img (array) : stacked astrometry corrected science image; 1 ccd
    final_wcs (?) : updated wcs information for the science image
    final_x (list) : list of x coordinates in pixels of detect OSIRIS stars
    final_y (list) : list of y coordinates in pixels of detect OSIRIS stars
    tbl_crds (list) : list of Gaia star coordinates
    astrom_path (str) : directory to write to astrometry related images to
    ccd (str) : 1 or 2
    final_name (str) : the name to write the corrected image to
    log_fname (str) : full path to log file

    Returns
    -------
    None
    """
    font_prop = font_manager.FontProperties(size=14)

    obj_ra = []
    obj_dec = []
    for i in range(len(final_x)):
        pos = final_wcs.all_pix2world(final_x[i], final_y[i], 0)
        obj_ra.append(pos[0])
        obj_dec.append(pos[1])
    obj_ra = np.array(obj_ra)
    obj_dec = np.array(obj_dec)

    obj_crds = coordinates.SkyCoord(
        obj_ra, obj_dec, unit=(u.deg, u.deg), frame='fk5')

    idx, d2d, d3d = obj_crds.match_to_catalog_sky(tbl_crds)

    # cutoff = 2  # arcsec
    good = np.where(d2d.arcsec < cutoff)[0]
    # good = np.arange(len(d2d))

    fig = plt.figure(figsize=(12, 14))
    ax = fig.add_subplot(121, projection=final_wcs)
    ax.set_xlabel("RA", fontproperties=font_prop)
    ax.set_ylabel("Dec", fontproperties=font_prop)
    ax.set_title("CCD"+ccd, fontproperties=font_prop)

    ax.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))
    ax.plot(obj_ra[good], obj_dec[good], '.', transform=ax.get_transform('fk5'),
            mec='k', mfc='none', label='OSIRIS')
    ax.plot(obj_ra[good], obj_dec[good], 'o', transform=ax.get_transform('fk5'),
            mec='r', mfc='none', label='OSIRIS matched')
    ax.plot(tbl_crds.ra, tbl_crds.dec, '*', transform=ax.get_transform('fk5'),
            mec='b', mfc='none', label='GAIA')
    plt.legend()

    ax2 = fig.add_subplot(122)
    ax2.hist(d2d[good].arcsec)
    ax2.set_xlabel('Separation (arcsec)', fontproperties=font_prop)
    ax2.set_title("CCD"+ccd+': Median = '+str(round(
        np.median(d2d[good].arcsec), 6))+'arcsec \n Stdv: '+str(round(
            np.std(d2d[good].arcsec), 6)),
        fontproperties=font_prop)

    plt.savefig(astrom_path+final_name.split('/')
                [-1][:-5]+'_astrometry_accuracy_ccd'+ccd+'.png')
    plt.suptitle(
        'Click q four times or close all plots to continue with pipeline reduction.')
    plt.show()


def crop_padding(img):
    """Crop out all padding in an image (rows/columns with all zeros).

    Parameters
    ----------
    img (arr) : one ccd image

    Returns
    -------
    cropped_img (arr) : same image after removing columns/rows with all zeros
    """
    ymax, xmax = img.shape

    # test_img = copy.deepcopy(img)
    # test_img[np.isnan(test_img)] = -10000
    # plt.figure()
    # plt.imshow(img, origin='lower', interpolation='none', vmin=np.nanmean(
    #     img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))
    # plt.show()

    # print(img.shape)

    # Calculate the sum of all rows
    xsums = np.nansum(img, axis=0)
    # Get the indices where the sum equals zero
    xinds = np.where(xsums == 0)[0]
    # print('xinds', len(xinds), xinds)
    # Subtract the indices from eachother (ex: [1,2,3,5] -> [1,1,2])
    xdiffs = np.diff(xinds)
    # print('xdiffs', xdiffs)
    # Find where there is a large jump in index (ex: [1,2,3,1000,1001,1002])
    final_xind = np.where(xdiffs > 100)[0][0]
    # print(final_xind)
    # print(xinds[final_xind])
    # This is the beginning/end of the actual data ex: x1=4, x2=999
    x1 = xinds[final_xind]
    x2 = xinds[final_xind+1]

    # Repeat for columns
    ysums = np.nansum(img, axis=1)
    yinds = np.where(ysums == 0)[0]
    # print(yinds)
    ydiffs = np.diff(yinds)
    # print(ydiffs)
    final_yind = np.where(ydiffs > 100)[0][0]
    # print(final_yind)
    y1 = yinds[final_yind]
    y2 = yinds[final_yind+1]

    # x1 = 198
    # x2 = 1046
    # y1 = 122
    # y2 = 2050
    # print('x/y', x1, x2, y1, y2)
    # Select the part of the image that contains data and not all zeros
    cropped_img = img[y1:y2, x1:x2]
    # print(cropped_img.shape)

    # plt.figure()
    # plt.imshow(cropped_img, origin='lower', interpolation='none',
    #             vmin=np.nanmean(img)-np.nanstd(img),
    #             vmax=np.nanmean(img)+np.nanstd(img))
    # plt.show()

    return cropped_img, x1, y1


def do_stacking(sci_final, all_headers, args, root, filt, astrom_path, log_fname):
    """Align images with eachother and stack/median combine them into one image.

    Parameters
    ----------
    sci_final (list): list of images after applying all calibrations/corrections
    all_headers (list) : list of the headers from each of the 3 extensions
    args (?) : user inputted arguments from the command line
    root (str) : contains object name
    filt (str) : current filter of the images (ex: Sloan_z)
    astrom_path (str) : directory to write to astrometry related images to
    log_fname (str) : full path to log file

    Returns
    -------
    final_aligned_image (list): [ccd1,ccd2] image after aligning/combining
                                individual images
    """
    if len(sci_final) == 1:
        gtcsetup.print_both(log_fname, 'ONLY 1 IMAGE; NOT COMBINING')

        # Write out the aligned/median combined images
        gtcsetup.write_ccd(all_headers[0][0], all_headers[0][1:],
                           sci_final[0], args.outputdir,
                           'aligned.fits', root, filt, log_fname)

        return sci_final

    # Loop through each ccd so that at the end of the loop we can combine
    # all the images in one ccd
    padded_final_aligned_image = []
    final_aligned_image = []
    all_aligned_images = []
    footprints = []
    padding = ((pad, pad), (pad, pad))
    constant_vals = ((np.nan, np.nan), (np.nan, np.nan))
    for ccd in range(len(sci_final[0])):

        gtcsetup.print_both(log_fname, '     CCD', ccd+1)

        footprint_1ccd = []
        aligned_images = []
        # Loop through each image for a particular ccd
        for i, image in enumerate(sci_final):

            gtcsetup.print_both(log_fname, '          Image', i+1)

            # If this is the first image, use it as the reference
            # Also add some padding of zeros around the edges so that
            # you don't lose data when you combine the images
            if i == 0:
                ref_imgs_init = np.array(image[ccd].data)
                # ################# MAY NEED TO UPDATE ##############
                # pad = 200  # pixels; checked GRB210704A, dithered by ~80 pix
                ######################################
                ref_imgs = np.pad(ref_imgs_init, padding,
                                  'constant', constant_values=constant_vals)
                aligned_images.append(ref_imgs)

                # MOVE ON TO THE NEXT IMAGE
                continue

            # Pad the image
            use_image_init = image[ccd].data

            use_image = np.pad(use_image_init, padding,
                               'constant', constant_values=constant_vals)
            # constant_values=((np.nan, np.nan), (np.nan, np.nan)))
            # Align the current image with the first reference image
            try:
                img_aligned, footprint = aa.register(
                    use_image, ref_imgs,
                    propagate_mask=True, detection_sigma=3.0)
            except ValueError:
                img_aligned, footprint = aa.register(
                    use_image.byteswap().newbyteorder(), ref_imgs,
                    propagate_mask=True, detection_sigma=3.0)

            # plt.figure()
            # plt.imshow(footprint*1, origin='lower')
            # plt.show()

            # add the footprint image to an array
            footprint_1ccd.append(np.array(footprint*10))

            # add the aligned image to an array
            aligned_images.append(np.array(img_aligned))

        # Median combine the aligned images for one ccd
        aligned_images = np.array(aligned_images)
        all_aligned_images.append(aligned_images)
        median_aligned_img = np.nanmedian(aligned_images, axis=0)
        padded_final_aligned_image.append(median_aligned_img)

        # Too many nans for ccdproc.combine? so I used numpy instead
        # median_aligned_img = ccdproc.combine(aligned_images,
        #                       method='median')

        # Crop out the padding in the median combined image
        cropped_img1, x1, y1 = crop_padding(median_aligned_img)

        all_headers[0][ccd + 1]['CRPIX1'] = float(
            all_headers[0][ccd+1]['CRPIX1'])+float(x1)
        all_headers[0][ccd + 1]['CRPIX2'] = float(
            all_headers[0][ccd+1]['CRPIX2'])+float(y1)

        cropped_img = median_aligned_img
        # gtcsetup.write_ccd(all_headers[i][0], all_headers[i][1:],
        #                    [CCDData(cropped_img, unit='electron'),
        #                     CCDData(cropped_img, unit='electron')], astrom_path,
        #                    'cropped_aligned'+str(i+1)+'_ccd'+str(ccd+1)+'.fits',
        #                    root, filt, log_fname)
        # print(cropped_img.shape)
        # Append the final aligned image as to be written later

        # CCDData(cropped_img, unit='electron'))
        final_aligned_image.append(cropped_img)

        footprints.append(np.array(footprint_1ccd))

    # Write out the aligned/median combined images
    gtcsetup.write_ccd(all_headers[0][0], all_headers[0][1:],
                       [final_aligned_image[0], final_aligned_image[1]],
                       args.outputdir,
                       'aligned.fits', root, filt, log_fname)

    # gtcsetup.write_ccd(all_headers[0][0], all_headers[0][1:],
    #                    padded_final_aligned_image, args.outputdir,
    #                    'padded_aligned.fits', root, filt, log_fname)
    # print(len(footprints), len(footprints[0]), len(footprints[0][0]))

    for i in range(len(footprint_1ccd)):
        gtcsetup.write_ccd(all_headers[i][0], all_headers[i][1:],
                           [footprints[0][i], footprints[1][i]], astrom_path,
                           'footprint'+str(i+1)+'.fits', root, filt, log_fname)

    for i in range(len(footprint_1ccd)):
        gtcsetup.write_ccd(all_headers[i][0], all_headers[i][1:],
                           [all_aligned_images[0][i], all_aligned_images[1][i]],
                           astrom_path,
                           'pad_aligned'+str(i+1)+'.fits', root, filt, log_fname)

    return final_aligned_image


def get_matches(all_gtc_inds, all_gaia_inds):
    """Test."""
    flat_list = [item for sublist in all_gaia_inds for item in sublist]
    dup_inds = [i for i, x in enumerate(flat_list) if flat_list.count(x) > 1]
    dupes = list(set(np.array(flat_list)[dup_inds]))

    flags = np.zeros(len(all_gaia_inds))
    # 0 = good, no duplicates or confusion
    # 1 = >1 gaia source found within search radius
    # 2 = confused- gaia source found more than once

    for i in range(len(all_gaia_inds)):
        flag = ''
        if len(all_gaia_inds[i]) > 1:
            flag = '1'
        if all_gaia_inds[i][0] in dupes:
            flag = '2'
        if flag == '':
            flag = '0'
        flags[i] = float(flag)

    good_inds = np.where(flags == 0.)[0]
    bad_inds = np.where(flags != 0.)[0]

    gtc = np.array(all_gtc_inds)[good_inds]
    gaia = np.array(all_gaia_inds)[good_inds]
    gaia_flat = np.array([item for sublist in gaia for item in sublist])

    gtc_bad = np.array(all_gtc_inds)[bad_inds]
    gaia_bad = np.array(all_gaia_inds)[bad_inds]
    # print(gaia_bad)
    gaia_flat_bad = np.array([item for sublist in gaia_bad for item in sublist])

    return gtc, gaia_flat, gtc_bad, gaia_flat_bad

# final_gtc_inds,final_gaia_inds,bad_gtc,bad_gaia=get_matches(all_gtc_inds,all_gaia_inds)


def calc_dist(xcen, ycen, radius, xpts, ypts):
    """Test."""
    inds = []
    distances = []
    for i in range(len(xpts)):
        dist = (xpts[i]-xcen)**2+(ypts[i]-ycen)**2
        if dist <= radius**2:
            inds.append(i)
            distances.append(dist)
    return np.array(inds), np.array(distances)


def plot_auto_astrometry(img, w, img_xy, img_radec, ref_img, ref_wcs, ref_xy,
                         ref_radec, all_gtc_inds, all_gaia_inds, bad_gtc_inds,
                         bad_gaia_inds):
    """Test."""
    colors = ['r', 'g', 'b', 'y', 'cyan', 'k', 'm']*10
    m = np.nanmean(img)
    s = np.nanstd(img)

    fig = plt.figure(figsize=(10, 10))
    # left bottom width height
    ax = fig.add_axes([0.1, 0.1, 0.45, 0.8])  # , projection=w)
    ax.imshow(img, cmap='gray', interpolation='none',
              origin='lower', vmin=m-s, vmax=m+s)
    # ax.set_title("Source Image")

    for i in range(len(all_gaia_inds)):
        ax.plot(img_xy[0][all_gtc_inds[i]], img_xy[1][all_gtc_inds[i]],
                'o', c=colors[i])  # , transform=ax.get_transform('fk5'))
        ax.plot(ref_xy[0][all_gaia_inds[i]], ref_xy[1][all_gaia_inds[i]],
                'o', c=colors[i], alpha=0.5)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    ax2 = fig.add_axes([0.5, 0.1, 0.45, 0.8], projection=w)
    ax2.imshow(img, cmap='gray', interpolation='none',
               origin='lower', vmin=m-s, vmax=m+s)

    for i in range(len(all_gaia_inds)):
        ax2.plot(img_radec.ra.deg[all_gtc_inds[i]],
                 img_radec.dec.deg[all_gtc_inds[i]], 'o', c=colors[i],
                 transform=ax2.get_transform('fk5'))

    for i in range(len(bad_gtc_inds)):
        ax2.plot(img_radec.ra.deg[bad_gtc_inds[i]],
                 img_radec.dec.deg[bad_gtc_inds[i]], 'x',
                 c='orange', transform=ax2.get_transform('fk5'))
    ax2.set_xlabel('RA')
    ax2.set_ylabel('Dec')
    plt.suptitle('OSIRIS image; press q to close and continue data reduction.')


def do_auto_astrometry(final_aligned_image, full_filt, seeing, astrom_path,
                       final_name, log_fname):
    """Automatically correct the astrometry of an image.

    Parameters
    ----------
    gaia (table?) : table containing ra/dec of reference stars
    source (table?) : table containing x/y of sources in image
    wcs_ref (?) : wcs projection information from original image
    ima (array) : image to fix wcs information for


    Returns
    -------
    w (?) : wcs projection information? updated
    good (array) : list of indices corresponding to sources in the image that
        were within a certain distance in arcsec of a gaia catalog star
        (aka matched)**
    """
    gtcsetup.print_both(log_fname, 'Doing AUTOMATIC astrometry')

    filt = full_filt.split('_')[-1]  # 'z'

    all_headers = []
    for ccd in range(len(final_aligned_image)):

        gtcsetup.print_both(log_fname, 'Working on CCD', ccd+1)

        img = final_aligned_image[ccd].data
        hdr = final_aligned_image[ccd].header

        w = final_aligned_image[ccd].wcs

        # Get the approximate center of the current image in ra/dec
        sky = w.pixel_to_world(img.shape[1]/2, img.shape[0]/2)
        gtcsetup.print_both(log_fname, 'Center of ccd:', sky.ra.deg, sky.dec.deg)

        # Get sdss image to compare to OSIRIS image
        ref_img = SkyView.get_images(position=str(sky.ra.deg)+','+str(sky.dec.deg),
                                     # survey='SDSSz',
                                     survey='SDSS'+filt,
                                     coordinates='J2000',
                                     pixels=str(img.shape[1])+','+str(img.shape[0]))
        if len(ref_img) == 0:
            ref_img = SkyView.get_images(
                position=str(sky.ra.deg)+','+str(sky.dec.deg), survey='DSS',
                coordinates='J2000', pixels=str(img.shape[1])+','+str(img.shape[0]))
        # print(ref_img)
        ref_img = ref_img[0]
        # print(ref_img[0].data)
        # Get wcs information from SDSS image
        ref_wcs = wcs.WCS(ref_img[0].header)
        # print(ref_wcs)

        # pix_limit = 0
        # Detect sources in the OSIRIS image and convert x/y to ra/dec using the
        # original wcs information
        obj_ra, obj_dec, final_x, final_y = get_osiris_obj_cat(
            img.byteswap().newbyteorder(), '', w, pix_limit, log_fname)

        # Get ra/dec of gaia stars that cover the OSIRIS image
        ref_radec = get_gaia_img_cat(img.data, pix_limit, w, sky, log_fname)

        print('number of gaia sources: ', len(ref_radec))
        print('number of sources in image: ', len(obj_ra))

        try:
            ref_x, ref_y = w.all_world2pix(ref_radec.ra, ref_radec.dec, 0,
                                           maxiter=100,
                                           tolerance=5e-4,
                                           detect_divergence=True,
                                           adaptive=True, quiet=True)
        except wcs.NoConvergence as e:
            gtcsetup.print_both(log_fname, "Indices of diverging points: {0}"
                                .format(e.divergent))
            gtcsetup.print_both(log_fname, "Indices of poorly converging points: {0}"
                                .format(e.slow_conv))
            gtcsetup.print_both(log_fname,
                                [i for i in range(len(ref_radec.ra))
                                 if i not in e.divergent and i not in e.slow_conv])
            gtcsetup.print_both(log_fname, len(ref_radec.ra))
            # gtcsetup.print_both(log_fname,
            #                     "Best solution:\n{0}".format(e.best_solution))
            # gtcsetup.print_both(log_fname,
            #                     "Achieved accuracy:\n{0}".format(e.accuracy))
            raise e

        search_radius = seeing/pix_scale*tolerance
        print(search_radius)
        n_stars = min(n_stars_init, len(ref_x), len(final_x))
        print(n_stars)
        fig, all_patches = plot_images(img, ref_img, w, ref_wcs, sky, ref_radec,
                                       obj_ra, obj_dec, n_stars, [], 'osiris',
                                       patch_radius=search_radius)

        # m = np.nanmean(img.data)
        # s = np.nanstd(img.data)
        # fig = plt.figure()
        # ax = fig.add_axes([0.1, 0.1, 0.45, 0.8])  # , projection=)
        # ax.imshow(img.data, origin='lower', interpolation='none',
        #           vmin=m-s, vmax=m+s, cmap='gray')
        # # ax.plot(final_x, final_y, 'o', label='image sources')
        # for i in range(len(final_x)):
        #     ax.add_patch(patches.Circle((final_x[i], final_y[i]),
        #                                 radius=cutoff/pix_scale, ec='b',
        #                                 fc='none', lw=2, alpha=0.8))
        # ax.plot(ref_x, ref_y, '.', label='gaia sources', c='r')
        # ax.set_title('before correcting or trimming')
        # ax.legend()
        # ax2 = fig.add_axes([0.5, 0.1, 0.45, 0.8], projection=w)
        # ax2.imshow(img.data, origin='lower', interpolation='none',
        #            vmin=m-s, vmax=m+s, cmap='gray')
        # ax2.plot(obj_ra, obj_dec, 'o', c='b', label='image sources',
        #          transform=ax2.get_transform('fk5'))
        # ax2.plot(ref_radec.ra, ref_radec.dec, '.', c='r',
        #          label='gaia sources', transform=ax2.get_transform('fk5'))

        ref_radec = ref_radec
        ref_xy = np.array([ref_x, ref_y])
        img_xy = np.array([final_x, final_y])
        # np.array([obj_ra[good],obj_dec[good]])
        img_radec = SkyCoord(list(zip(obj_ra, obj_dec)), frame='fk5', unit='deg')

        all_gaia_inds = []
        all_gtc_inds = []
        for i in range(n_stars):
            xcen = final_x[i]  # use_img_xy[0][i]
            ycen = final_y[i]  # use_img_xy[1][i]
            inds, distances = calc_dist(xcen, ycen, search_radius, ref_x, ref_y)
            if len(inds) != 0:
                print(i, xcen, ycen, inds, distances)
                all_gaia_inds.append(inds)
                all_gtc_inds.append(i)

        final_gtc_inds, final_gaia_inds, bad_gtc, bad_gaia = get_matches(
            all_gtc_inds, all_gaia_inds)

        plot_auto_astrometry(img, w, img_xy, img_radec, ref_img[0].data,
                             ref_wcs, ref_xy, ref_radec,
                             all_gtc_inds, all_gaia_inds, bad_gtc, bad_gaia)

        stars_xy = img_xy.T[final_gtc_inds]
        coords = SkyCoord(np.array([ref_radec[final_gaia_inds].ra.deg,
                                    ref_radec[final_gaia_inds].dec.deg]).T,
                          frame='fk5', unit='deg')

        new_wcs = fit_wcs_from_points(xy=stars_xy.T, world_coords=coords.T,
                                      projection='TAN', sip_degree=2)

        hdr.update(new_wcs.to_header(relax=True))
        final_wcs = WCS(hdr)

        calc_astrometry_accuracy(img, final_wcs, final_x, final_y, ref_radec,
                                 astrom_path, str(ccd+1), final_name +
                                 str(ccd+1)+'.fits',
                                 log_fname)
        all_headers.append(hdr)

    return all_headers
