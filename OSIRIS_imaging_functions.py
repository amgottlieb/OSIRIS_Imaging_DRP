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


def do_calib(obj, filt, log_fname, nccd, mbias, mflat, mask_ccd, gain,
             rdnoise, calib_path, bpm_path, root):
    """Test."""
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
    """Test."""
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

    return all_sci_proc  # , cr_mask_ccd


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

        masked_median_img = ccdproc.combine(
            all_masked_bkg_sub_images, method='median')

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
                obj_dec, n_stars, pts, select):
    """Test."""
    font_prop = font_manager.FontProperties(size=14)

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(121, projection=w)
    # ax = plt.subplot(projection=w)
    ax.set_xlabel("RA", fontproperties=font_prop)  # INCREASE SIZE
    ax.set_ylabel("Dec", fontproperties=font_prop)  # INCREASE SIZE

    # ax2 = fig.add_axes([0.5, 0.1, 0.45, 0.8], projection=w)
    ax.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))

    if select == 'gaia':
        for j in range(len(tbl_crds.ra)):
            x, y = w.all_world2pix(tbl_crds.ra.deg[j], tbl_crds.dec.deg[j], 0)
            ax.add_patch(patches.Circle((x, y), radius=5, ec='b', fc='none',
                                        lw=2, alpha=0.8,
                                        label='Reference: RA = %f, Dec = %f' %
                                        (tbl_crds.ra.deg[j], tbl_crds.dec.deg[j]),
                                        picker=True))
    else:
        ax.plot(tbl_crds.ra, tbl_crds.dec, 'o', transform=ax.get_transform('fk5'),
                mec='b', mfc='none', label='GAIA')

    if select == 'osiris':
        for j in range(len(obj_ra[:n_stars])):
            x, y = w.all_world2pix(obj_ra[j], obj_dec[j], 0)
            ax.add_patch(patches.Circle((x, y), radius=5, ec='k', fc='none',
                                        lw=2, alpha=0.8,
                                        label='Position: x = %f, y = %f' %
                                        (x, y), picker=True))
    else:
        ax.plot(obj_ra[:n_stars], obj_dec[:n_stars], 'o',
                transform=ax.get_transform('fk5'), mec='k', mfc='none',
                label='OSIRIS')

    for i, p in enumerate(pts):
        ax.text(p[0], p[1], str(i+1), c='magenta', fontsize='large')

    ax.set_title('Your OSIRIS image', fontproperties=font_prop)  # INCREASE SIZE
    all_patches = ax.patches

    #########################

    ax2 = fig.add_subplot(122, projection=ref_wcs)
    ax2.set_xlabel("RA", fontproperties=font_prop)  # INCREASE SIZE
    ax2.set_ylabel("Dec", fontproperties=font_prop)  # INCREASE SIZE

    ref = ref_img[0].data
    ax2.imshow(ref, origin='lower', cmap='gray', vmin=np.nanmean(
        ref)-np.nanstd(ref), vmax=np.nanmean(ref)+np.nanstd(ref))

    ax2.plot(tbl_crds.ra, tbl_crds.dec, 'o', transform=ax2.get_transform('fk5'),
             mec='b', mfc='none', label='GAIA')
    ax2.plot(obj_ra[:n_stars], obj_dec[:n_stars], 'o',
             transform=ax2.get_transform('fk5'), mec='k', mfc='none',
             label='OSIRIS')

    ax2.set_title('SDSS zband image for reference',
                  fontproperties=font_prop)  # INCREASE SIZE

    plt.legend()

    return fig, ref_wcs, all_patches


def select_points(img, select_n_stars, color):
    """Test."""
    while True:
        pts = []
        while len(pts) < select_n_stars:
            tellme('Right click to select 6 '+color +
                   ' stars on the LEFT plot;\n'
                   'To remove a selected star, '
                   'click the back button on your mouse. \n'
                   'press enter when done.')

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

        # plt.plot(np.array(pts2).T[0],np.array(pts2).T[1],'o',label='GAIA')
        # print(pts)

        # tellme('Happy? Press y for yes, n for no.')
        tellme('Happy? Key click for yes, mouse click for no')

        # val = check_if_done()

        if plt.waitforbuttonpress():
            break
            check_fig.close()

    return pts


def select_yes_no(img, select_n_stars):
    """Test."""
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

        # check_fig = plt.figure()
        # plt.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        #     img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))
        # plt.plot(np.array(pts).T[0], np.array(pts).T[1], 'o', c='orange')

        # # plt.plot(np.array(pts2).T[0],np.array(pts2).T[1],'o',label='GAIA')
        # # print(pts)

        # # tellme('Happy? Press y for yes, n for no.')
        # tellme('Happy? Key click for yes, mouse click for no')

        # val = check_if_done()

        if plt.waitforbuttonpress():
            break
            # check_fig.close()

    return pts


def plot_check_img(img, final_wcs, tbl_crds):
    """Test."""
    font_prop = font_manager.FontProperties(size=14)

    fig3 = plt.figure(figsize=(12, 10))
    # letter = fig3.canvas.mpl_connect('key_press_event', press)

    ax = fig3.add_subplot(121)  # , projection=final_wcs)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    ax.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))

    x1 = 100
    y = img.shape[0]-100
    x2 = img.shape[1]-100

    ax.text(x1, y, 'YES', c='g')
    ax.text(x2, y, 'NO', c='g')
    # xs = [x1, x2]

    ax.add_patch(patches.Circle((x1, y), radius=5, ec='k', fc='none',
                                lw=2, alpha=0.8,
                                label='Position: x = %f, y = %f' %
                                (x1, y), picker=True))
    ax.add_patch(patches.Circle((x2, y), radius=5, ec='k', fc='none',
                                lw=2, alpha=0.8,
                                label='Position: x = %f, y = %f' %
                                (x2, y), picker=True))
    all_patches = ax.patches

    #########################

    ax2 = fig3.add_subplot(122, projection=final_wcs)
    ax2.set_xlabel("RA", fontproperties=font_prop)  # INCREASE SIZE
    ax2.set_ylabel("Dec", fontproperties=font_prop)  # INCREASE SIZE

    ax2.imshow(img, origin='lower', cmap='gray', vmin=np.nanmean(
        img)-np.nanstd(img), vmax=np.nanmean(img)+np.nanstd(img))
    ax2.plot(tbl_crds.ra, tbl_crds.dec, '*', transform=ax2.get_transform('fk5'),
             mec='b', mfc='none', label='GAIA')

    # ax2.set_title('SDSS zband image for reference',
    #               fontproperties=font_prop)  # INCREASE SIZE

    # plt.legend()

    return fig3, all_patches


def get_gaia_img_cat(img, pix_limit, w, sky, log_fname):
    """Test."""
    ######################################
    # Get gaia comparison catalog centered on ccd
    # using cone with radius 6 arcmin (Nora chose this)
    # NOTE: radius of 6 gives circle, larger radius gives
    # weird shape that doesn't cover the ccd...
    # gaia = Irsa.query_region(SkyCoord(sky.ra, sky.dec, frame='fk5'),
    #                          catalog="gaia_dr2_source", spatial="Cone",
    #                          radius=6*u.arcmin)

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
    # print(gaia_ind)
    # print(gaia_flux_init)

    # gaia_fluxes = gaia_flux_init[gaia_ind]
    # print(gaia_fluxes)
    tbl_crds = coordinates.SkyCoord(gaia['ra'][gaia_ind], gaia['dec'][gaia_ind],
                                    unit=(u.deg, u.deg), frame='fk5')

    return tbl_crds


def get_osiris_obj_cat(img, fname, w, pix_limit, log_fname):
    """Test."""
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
    # ****not currently used to remove objects
    axes_diff = np.array([round(abs(objects['a'][i]-objects['b'][i]), 3)
                          for i in range(len(objects))])
    # print(axes_diff)

    # Sort x and y based on the flux; we want to take the n brightest stars
    # sorted_flux = np.array(sorted(objects['flux'])[::-1])
    sorted_x = np.array([x for _, x in sorted(
        zip(objects['flux'], objects['x']))])[::-1]
    sorted_y = np.array([x for _, x in sorted(
        zip(objects['flux'], objects['y']))])[::-1]

    # Eliminate stars that are close to the edge
    # pix_limit is how much of the edge of the ccd to exclude
    final_x = []
    final_y = []
    for i in range(len(sorted_x)):
        x = sorted_x[i]
        y = sorted_y[i]
        if x > pix_limit and (x < img.shape[1]-pix_limit) and (
                y > pix_limit) and (y < img.shape[0]-pix_limit):
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
    # print(obj_ra)
    # print(obj_dec)

    return obj_ra, obj_dec, final_x, final_y


def in_circle(center_x, center_y, radius, x, y):
    """Test."""
    square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
    return square_dist <= radius ** 2


def cross_match(all_patches, selected_pts, radius):
    """Test."""
    final_pts = []
    for pt in selected_pts:
        print('Pt', pt)
        for patch in all_patches:
            p_cen = patch.center
            incirc = in_circle(p_cen[0], p_cen[1], radius, pt[0], pt[1])
            if incirc:
                print('found')
                final_pts.append([patch.center[0], patch.center[1]])
                break
    return np.array(final_pts)


def final_check(img, final_wcs, tbl_crds, log_fname):
    """Test."""
    # Plot the OSIRIS image with the new wcs and overplot Gaia points to make
    # sure astrometry is good; user selects yes or no
    fig3, all_patches = plot_check_img(img, final_wcs, tbl_crds)
    text = 'Right click to select yes or no'
    tellme(text)
    plt.waitforbuttonpress()

    select_n_stars = 1
    selected_pts = select_yes_no(img, select_n_stars)
    # print(selected_pts)

    # CROSS MATCH SELECTED POINTS WITH ALL PATCHES
    final_pts = cross_match(all_patches, selected_pts, 5)

    # print(final_pts)
    # print(img.shape)
    plt.close(fig3)
    # print('------------')

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

    Note: Things the user may want to change
        - pix_limit
        - obj_frac
        - n_stars
        - select_n_stars
        - sip_deg


    NOTE TO SELF*** ADD IN FILTER**** ADD ELLIPTICITY CHECK ***

    Parameters
    ----------
    final_name (str) : the name to write the corrected image to
    fname (str) : the name of the image to perform astrometry on (will be 1 ccd)
    full_filt (str) : the current filter (ex: Sloan_z)
    ccd 
    astrom_path
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

    pix_limit = 0
    # Detect sources in the OSIRIS image and convert x/y to ra/dec using the
    # original wcs information
    obj_ra, obj_dec, final_x, final_y = get_osiris_obj_cat(img, fname, w,
                                                           pix_limit, log_fname)

    # Get ra/dec of gaia stars that cover the OSIRIS image
    tbl_crds = get_gaia_img_cat(img, pix_limit, w, sky, log_fname)

    # Set how many OSIRIS stars should be displayed for astrometry
    obj_frac = 4
    n_stars = round(len(obj_ra)/obj_frac)
    # print(n_stars)
    n_stars = 40

    # Get sdss image for testing purposes
    ref_img = SkyView.get_images(position=str(sky.ra.deg)+','+str(sky.dec.deg),
                                 # survey='SDSSz',
                                 survey='SDSS'+filt,
                                 coordinates='J2000',
                                 pixels=str(img.shape[1])+','+str(img.shape[0]))[0]
    # Get wcs information from SDSS image
    ref_wcs = wcs.WCS(ref_img[0].header)

    # Plot the OSIRIS image and the SDSS image with detected objects and
    # Gaia sources overlayed; select OSIRIS stars
    fig, ref_wcs, all_patches = plot_images(img, ref_img, w, ref_wcs, sky, tbl_crds,
                                            obj_ra, obj_dec, n_stars, [], 'osiris')

    # Wait for user to select at least n points
    select_n_stars = 6
    text = 'You will select '+str(select_n_stars)+' BLACK (OSIRIS) stars first.\n' \
        'Then you will select the corresponding Gaia stars (in blue) in the' \
        'same order. \n To remove a selected star, click the back button' \
        ' on your mouse.'
    tellme(text)
    plt.waitforbuttonpress()
    # Points selected by user
    pts = select_points(img, select_n_stars, 'BLACK')

    # Cross match where the user selected with the actual positions of the stars
    source_xy = cross_match(all_patches, pts, 5)

    plt.close(fig)

    # Plot the OSIRIS image and the SDSS image with detected objects and
    # Gaia sources overlayed; select GAIA stars
    fig2, ref_wcs2, all_patches2 = plot_images(img, ref_img, w, ref_wcs, sky,
                                               tbl_crds, obj_ra, obj_dec,
                                               n_stars, pts, 'gaia')

    # Wait for user to select at least n points
    tellme('You will select 6 BLUE (gaia) stars in the same order (green numbers). ')
    plt.waitforbuttonpress()
    pts2 = select_points(img, select_n_stars, 'BLUE')

    # Cross match where the user selected with the actual positions of the stars
    final_pts2 = cross_match(all_patches2, pts2, 5)
    plt.close(fig2)

    # Convert x/y to ra/dec using OSIRIS image wcs info
    ref_radec = np.array([w.all_pix2world(final_pts2[i][0], final_pts2[i][1], 0)
                          for i in range(len(final_pts2))])

    # Create a SkyCoord object for the GAIA stars to be inputted on next line
    coords = SkyCoord(np.array(ref_radec), frame='fk5', unit='deg')

    # Calculate new wcs
    sip_deg = 2
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
                             ccd, log_fname)

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
                             ccd, log_fname):
    """Test."""
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

    cutoff = 2  # arcsec
    good = np.where(d2d.arcsec < cutoff)[0]

    fig = plt.figure(figsize=(12, 14))
    # letter = fig3.canvas.mpl_connect('key_press_event', press)

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
        np.median(d2d[good].arcsec), 6))+'arcsec',
        fontproperties=font_prop)

    plt.savefig(astrom_path+'astrometry_accuracy_ccd'+ccd+'.png')
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

    # Calculate the sum of all rows
    xsums = np.sum(img, axis=0)
    # Get the indices where the sum equals zero
    xinds = np.where(xsums == 0)[0]
    # Subtract the indices from eachother (ex: [1,2,3,5] -> [1,1,2])
    xdiffs = np.diff(xinds)
    # Find where there is a large jump in index (ex: [1,2,3,1000,1001,1002])
    final_xind = np.where(xdiffs > 1)[0]
    # This is the beginning/end of the actual data
    x1 = xinds[final_xind][0]
    x2 = xinds[final_xind+1][0]

    # Repeat for columns
    ysums = np.sum(img, axis=1)
    yinds = np.where(ysums == 0)[0]
    ydiffs = np.diff(yinds)
    final_yind = np.where(ydiffs > 1)[0]
    y1 = yinds[final_yind][0]
    y2 = yinds[final_yind+1][0]

    # Select the part of the image that contains data and not all zeros
    cropped_img = img[y1:y2, x1:x2]

    return cropped_img


def do_stacking(sci_final, all_headers, args, root, filt, astrom_path, log_fname):
    """Test."""
    if len(sci_final) == 1:
        gtcsetup.print_both(log_fname, 'ONLY 1 IMAGE; NOT COMBINING')

        # Write out the aligned/median combined images
        gtcsetup.write_ccd(all_headers[0][0], all_headers[0][1:],
                           sci_final[0], args.outputdir,
                           'aligned.fits', root, filt, log_fname)

        return sci_final

    # Loop through each ccd so that at the end of the loop we can combine
    # all the images in one ccd
    final_aligned_image = []
    footprints = []
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
                pad = 200  # pixels checked GRB210704A, dithered by ~80 pix
                ######################################
                padding = ((pad, pad), (pad, pad))
                ref_imgs = np.pad(ref_imgs_init, padding,
                                  'constant', constant_values=0)
                aligned_images.append(ref_imgs)

                # MOVE ON TO THE NEXT IMAGE
                continue

            # Pad the image
            use_image_init = image[ccd].data
            use_image = np.pad(use_image_init, padding,
                               'constant', constant_values=(
                                   (np.nan, np.nan), (np.nan, np.nan)))

            # Align the current image with the first reference image
            try:
                img_aligned, footprint = aa.register(
                    use_image, ref_imgs,
                    propagate_mask=True, detection_sigma=3.0)
            except ValueError:
                img_aligned, footprint = aa.register(
                    use_image.byteswap().newbyteorder(), ref_imgs,
                    propagate_mask=True, detection_sigma=3.0)

            # add the footprint image to an array
            footprint_1ccd.append(footprint)

            # add the aligned image to an array
            aligned_images.append(np.array(img_aligned))

        # Median combine the aligned images for one ccd
        aligned_images = np.array(aligned_images)
        median_aligned_img = np.nanmedian(aligned_images, axis=0)

        # Too many nans for ccdproc.combine? so I used numpy instead
        # median_aligned_img = ccdproc.combine(aligned_images,
        #                       method='median')

        # Crop out the padding in the median combined image
        cropped_img = crop_padding(median_aligned_img)

        # Append the final aligned image as a CCDData object to be
        # written later
        final_aligned_image.append(CCDData(cropped_img, unit='electron'))
        footprints.append(footprint_1ccd)

    # Write out the aligned/median combined images
    gtcsetup.write_ccd(all_headers[0][0], all_headers[0][1:],
                       final_aligned_image, args.outputdir,
                       'aligned.fits', root, filt, log_fname)

    for i in range(len(footprint_1ccd)):
        gtcsetup.write_ccd(all_headers[i][0], all_headers[i][1:],
                           footprints.T[i], astrom_path,
                           'footprint.fits', root, filt, log_fname)

    return final_aligned_image


def do_auto_astrometry(gaia, sources, wcs_ref, ima, log_fname):
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
    # Get skycoord formatted ra/dec of refrence gaia stars
    sky_ref_radec_init = SkyCoord(list(zip(np.array(gaia['ra']),
                                           np.array(gaia['dec']))),
                                  frame='fk5', unit='deg')

    # Convert gaia ra/dec into pixels using the original wcs info
    sky_ref_x = []
    sky_ref_y = []
    good_ref = []
    print('number of gaia sources: ', len(gaia['ra']))
    print('number of sources in image: ', len(sources))

    try:
        x, y = wcs_ref.all_world2pix(gaia['ra'], gaia['dec'], 1, maxiter=100,
                                     tolerance=5e-4, detect_divergence=True,
                                     adaptive=True, quiet=True)
    except wcs.NoConvergence as e:
        gtcsetup.print_both(log_fname, "Indices of diverging points: {0}"
                            .format(e.divergent))
        gtcsetup.print_both(log_fname, "Indices of poorly converging points: {0}"
                            .format(e.slow_conv))
        gtcsetup.print_both(log_fname,
                            [i for i in range(len(gaia['ra']))
                             if i not in e.divergent and i not in e.slow_conv])
        gtcsetup.print_both(log_fname, len(gaia['ra']))
        # gtcsetup.print_both(log_fname,
        #                     "Best solution:\n{0}".format(e.best_solution))
        # gtcsetup.print_both(log_fname,
        #                     "Achieved accuracy:\n{0}".format(e.accuracy))
        raise e

    plt.figure()
    plt.plot(sources['xcentroid'], sources['ycentroid'], 'o', label='image sources')
    plt.plot(x, y, '.', label='gaia sources')
    plt.title('before correcting or trimming')
    plt.legend()

    for j in range(len(gaia['ra'])):

        # exclude sources that are too close to the edge;
        # I randomly chose 20 pixels as the limit
        if x[j] > 20 and y[j] > 20 and x[j] < ima.shape[1]-20. and (
                y[j] < ima.shape[0]-20.):
            sky_ref_x.append(x)
            sky_ref_y.append(y)
            good_ref.append(j)
    sky_ref_x = np.array(sky_ref_x)
    sky_ref_y = np.array(sky_ref_y)

    # Get stars that are not close to the edge
    sky_ref_radec = sky_ref_radec_init[good_ref]

    # Get / format x/y coordinates of sources found in image
    sky_img_xy = np.array([list(sources['xcentroid']),
                           list(sources['ycentroid'])])

    # Convert all x/y coordinates of sources in image to ra/dec
    sky_img_radec_init = wcs_ref.wcs_pix2world(
        sky_img_xy[0], sky_img_xy[1], 1)

    # get skycoord formatted ra/dec of sources (same format as reference stars)
    sky_img_radec = SkyCoord(
        list(zip(sky_img_radec_init[0], sky_img_radec_init[1])),
        frame='fk5', unit='deg')

    # match reference catalog with source list
    idx, d2d, d3d = sky_img_radec.match_to_catalog_sky(sky_ref_radec)
    # Only take sources that are within a certain distance
    # in this case, the median distance;
    cutoff = np.median(d2d.arcsec)*2.  # TO DO: UPDATE THIS #############
    good = np.where(d2d.arcsec < cutoff)[0]
    use_ref_radec = sky_ref_radec[idx][good]
    use_img_xy = np.array([sky_img_xy[0][good], sky_img_xy[1][good]])

    gaia_sky_coords = []
    for i in range(len(use_ref_radec)):
        gaia_sky_coords.append(
            [use_ref_radec[i].ra.degree, use_ref_radec[i].dec.degree])
    gaia_sky_coords = np.array(gaia_sky_coords)
    print(gaia_sky_coords.shape)

    plt.figure()
    plt.hist(d2d.arcsec, bins=20)
    plt.axvline(cutoff, c='k')

    # Get the new wcs information
    # format of lists required by fit_wcs_from_points:
    # use_img_xy = np.array([[xlist],[ylist]])
    # use_ref_radec=SkyCoord([(ra,dec),(ra,dec)...],frame='fk5',unit='deg')
    # NOTE on sip_degree=2: not sure why I picked this but it seems to work?
    w = fit_wcs_from_points(xy=use_img_xy, world_coords=use_ref_radec,
                            projection='TAN', sip_degree=2)
    print(w)

    sources_sky_coords = w.wcs_pix2world(
        list(zip(use_img_xy[0], use_img_xy[1])), 0)
    # print(sources_sky_coords)
    print(sources_sky_coords.shape)
    # check how good new coordinate system is
    plt.figure()
    plt.plot(sources_sky_coords[:, 0],
             sources_sky_coords[:, 1], 'o', label='sources in image')
    plt.plot(gaia_sky_coords[:, 0], gaia_sky_coords[:, 1], '.', label='gaia sources')
    plt.legend()
    print('-----------------------------------------')
    return w, good
