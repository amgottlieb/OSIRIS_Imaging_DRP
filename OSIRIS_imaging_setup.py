"""Test."""

from __future__ import print_function
from glob import glob
import argparse
import shutil
import os
# import sys
from astropy.io import fits
# from astropy.io.fits import getheader
import ccdproc
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
from astropy.modeling import models
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize
import matplotlib.patches as patches
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord
from astropy.wcs import wcs


def read_args():
    """Create wrapper to parse and customize command-line arguments."""
    parser = parse_args()
    args, use_slash = setup_args(parser)
    return args, use_slash


def parse_args():
    """Define and read arguments expected from command line."""
    parser = argparse.ArgumentParser(
        description='A Python data reduction pipeline'
        'for GTC/OSIRIS imaging')
    add = parser.add_argument

    # mandatory arguments
    add('objectid',
        help='Object name as given in the FITS file header. If the'
        ' object name contains spaces, enclose the name with'
        ' quotation marks.')

    # optional arguments
    add('--workdir', dest='workdir', default='./',
        help='Set the directory ')

    add('--outputdir', dest='outputdir', default='./output/',
        help='Set the directory to put results in')

    add('--filter', dest='filt', default='r',
        help='Filter of the images to be reduced. Can be multiple ex: g,i,z \
        Default: r')

    add('--dobias', dest='dobias', action='store_true',
        help='Make the master bias instead of reading it from input file.')

    add('--bias', dest='biasfile', default='MasterBias',
        help='Name of the input master bias file.')

    add('--doflat', dest='doflat', action='store_true',
        help='Make the master flat instead of reading it from input file.')

    add('--flat', dest='flatfile', default='MasterFlat',
        help='Name of the input master flat file.')

    add('--domask', dest='domask', action='store_true',
        help='Make mask to remove bad pixels instead of reading it from input file.')

    add('--mask', dest='maskfile', default='BPmask',
        help='Name of the input mask file.')

    add('--dowcs', dest='dowcs', action='store_true',
        help='Improve the astrometric solution.')

    add('--dolog', dest='dolog', action='store_true',
        help='Create an observing log without going through the full reduction.')

    add('--rawdir', dest='rawdir', default='raw/',
        help='Directory of the raw data.')

    add('--dooverwrite', dest='dooverwrite', default=False,
        help='Overwrites any existing files.')

    return parser


def setup_args(parser):
    """Any manipulation of the arguments that may be required."""
    args = parser.parse_args()
    # check for \ or / depending on windows or linux/mac
    # print('system',sys.platform)
    print('    CHECKING for slashes in ', args.workdir)
    if '\\\\' in args.workdir:
        use_slash = '\\\\'
    elif '/' in args.workdir:
        use_slash = '/'
    elif "\\" in args.workdir:
        use_slash = '\\'

    if not args.workdir[-1] == use_slash:
        args.workdir += use_slash
    if not args.outputdir[-1] == use_slash:
        args.outputdir += use_slash

    return args, use_slash


def sort_files(obj, raw_path, raw_list, use_slash, outputdir, dooverwrite):
    r"""Sort/move files into different folders based on file names.

    Can be sorted based on header information (commented code below, may need
    to be modified). I moved away from this because the header information kept
    getting messed up (obstype = focus for science, etc).

    Parameters
    ----------
    obj (str) : name of main science target
    raw_path (str) : full path to raw data (C:\User etc)
    raw_list (array) : array of file names for raw data; includes full path to
                        file
    use_slash (str) : \ or \\ depending on os, windows or linux
    outputdir (str) : full path to folder where output files are stored
                        (master files, processed data)
    dooverwrite (bool) : if True, overwrite current existing directories;
                        WILL LOSE OUTPUT FOLDER IF TRUE

    Returns
    -------
    bias_list, flat_list, sci_list, std_list (array) : lists of all files in
    each directory; each file contains full path
    """
    # NOTE: raw_list already has file path in it
    print('    Sorting files...')
    # this creates the directories if they don't already exist
    bias_path, flat_path, sci_path, std_path = setup_directories(
        raw_path, use_slash, outputdir, dooverwrite)

    for f in raw_list:
        if 'Bias' in f:
            move_path = bias_path
        elif 'Flat' in f:
            move_path = flat_path
        elif 'Image' in f:
            with fits.open(f) as file_open:
                hdr = file_open[0].header
                target = hdr['OBJECT'].replace(' ', '')
                print('        Objectid', obj, 'check against target', target)
                if target == obj:
                    move_path = sci_path
                else:
                    move_path = std_path
        else:
            move_path = ''

        if move_path != '':
            shutil.copy2(f, move_path)

    bias_list = [i.replace(os.sep, '/') for i in glob(bias_path+'0*.fits*')]
    flat_list = [i.replace(os.sep, '/') for i in glob(flat_path+'0*.fits*')]
    std_list = [i.replace(os.sep, '/') for i in glob(std_path+'0*.fits*')]
    sci_list = [i.replace(os.sep, '/') for i in glob(sci_path+'0*.fits*')]

    return bias_list, flat_list, sci_list, std_list


def setup_directories(raw_path, use_slash, outputdir, dooverwrite):
    r"""Create necessary folders and copy all files into the copy folder.

    Parameters
    ----------
    raw_path (str) : full path to raw data (C:\User etc)
    use_slash (str) : \ or \\ depending on os, windows or linux
    outputdir (str) : full path to folder where output files are stored
                        (master files, processed data)
    dooverwrite (bool) : if True, overwrite current existing directories;
                        WILL LOSE OUTPUT FOLDER IF TRUE

    Returns
    -------
    bias_path (str) : full path to bias folder (C:\User etc)
    flat_path (str) : full path to flat folder (C:\User etc)
    sci_path (str) : full path to science folder (C:\User etc)
    std_path (str) : full path to standard folder (C:\User etc)
    """
    bias_path = raw_path+'bias'+use_slash
    flat_path = raw_path+'flat'+use_slash
    sci_path = raw_path+'object'+use_slash
    std_path = raw_path+'standard'+use_slash
    # Copy all data to this directory in case you need to start over
    copy_path = raw_path+'copy'+use_slash

    paths = [bias_path, flat_path, sci_path, std_path, outputdir]

    print('        Overwriting folders?', dooverwrite)

    # if there isn't a copy directory already, make it and copy all files into it
    if not os.path.exists(copy_path):
        print('        Creating copy directory')
        os.makedirs(copy_path)
        print('        Copying all fits files into copy directory')
        copy_list = [i.replace(os.sep, '/') for i in glob(raw_path+'0*.fits*')]
        for f in copy_list:
            shutil.copy2(f, copy_path)

    # go through the remaining folders and create them if they don't already exist
    # or you want to overwrite them
    for i, path1 in enumerate(paths):
        if not os.path.exists(path1):
            os.makedirs(path1)
            print('        Creating directory', path1)
        elif dooverwrite == 'True' or dooverwrite is True:
            shutil.rmtree(path1)
            os.makedirs(path1)
            print('        Overwriting directory', path1)
        else:
            print('        ', path1, ': Directory already exists; not overwriting.')

    return bias_path, flat_path, sci_path, std_path


def print_obslog_folder(folder, out, verbose):
    """Print observing log for a specific folder.

    Parameters
    ----------
    folder (str) : full path to folder that you want to print the observing log for
    out (str) : observing log output file

    Returns
    -------
    None; Print statements

    """
    for i, f in enumerate(folder):
        with fits.open(f) as file_open:
            hdr = file_open[0].header
            obsid = os.path.basename(f).split('-')
            filt = hdr['FILTER2']
            grism = hdr['FILTER4']
            target = hdr['OBJECT'].replace(' ', '')
            expo = hdr['EXPTIME']
            date = hdr['DATE-OBS']
            airmass = hdr['AIRMASS']
            msg = '{0:2d} {1:14s} {2:<10s} {3:<10s} {4:<23s}  {5:<10f}  {6:5f}  {7:<14s}'.format(
                i+1, obsid[0], filt, grism, date, airmass, expo, target)
        print(msg, file=out)
        if verbose:
            print(msg)


def print_obslog(obj, flt, bias, flat, sci, std, verbose=True):
    """Print header + observing logs for all folders. Does not return anything.

    Parameters
    ----------
    obj (str) : name of main science target
    flt (str) : list of filters, used for output file name; formatted as 'g,r'
    bias (array) : list of all files in the bias directory
    flat (array) : list of all files in the flat directory
    sci (array) : list of all files in the science directory
    std (array) : list of all files in the standard directory
    Note: all arrays above contain full path for each file

    Returns
    -------
    None
    """
    # Create output obs log file name
    output = '{0}_{1}.obslog'.format(obj.replace(' ', '_'), flt.strip(','))
    print('{0}\n-----'.format(output))

    # Open observing log output file in current working directory
    out = open(output, 'w')
    head = '{0:<16s}  {1:<10s} {2:<10s} {3:<23s}  {4:<10s}  {5:<14s}  {6:<14s}'.format(
        'ObservationID', 'Filter', 'Grism', 'Start Time', 'Airmass', 'Exposure', 'Target')
    print(head, file=out)
    if verbose:
        print(head)

    print_obslog_folder(bias, out, verbose)
    print_obslog_folder(flat, out, verbose)
    print_obslog_folder(sci, out, verbose)
    print_obslog_folder(std, out, verbose)

    return


def createMaster(full_flist, frametype, nccd, mbias, gain, rdnoise, outputdir,
                 name, filt):
    """Given list of files, create Master bias or flat.

    Usage:
    mflat, flatname = gtcsetup.createMaster(ff, 'flat', nccd, mbias,
                                        gain, rdnoise,
                                        args.outputdir,
                                        args.flatfile, filt)

    Parameters
    ----------
    full_list (array) : list of all files in flat or bias folder
    frametype (str) : 'bias' or 'flat'
    nccd (int) : 2 for osiris
    mbias (array) : [None, None] for frametype = bias; mbias = master bias if
            creating master flat
    gain (float) : instrument constant
    rdnoise (float) : instrument constants
    outputdir (str) : directory to write to
    name (str) : name of the output file to write to (INCLUDES FULL PATH)
    filt (str) : for example, Sloan_g; doesn't mastter for bias;
                filters are looped over in main code

    Returns
    -------
    master (array) :[ccd1, ccd2] processed data
    fname (str) : file name (including full path) that master was written to
    """
    print(full_flist)

    # Get all filters present (really only used for flat)
    all_filters = get_filters_from_header(full_flist)

    # if creating a master bias, don't add a filter to the filename
    if frametype == 'bias':
        add = ''
        flist = full_flist
    # if creating a master flat, add a filter to the filename and select only
    # the relevant filters
    else:
        add = filt
        use = np.where(np.array(all_filters) == filt)[0]
        flist = np.array(full_flist)[use]

    # add filter to filename if necessary
    fname = name+add+'.fits'

    # go through each file, process it, write out as temp file
    for j, f in enumerate(flist):
        # os.system("tar -xzf {}".format(f))
        # f=f.replace('.gz','')
        with fits.open(f) as fits_open:
            hdr = fits_open[0].header

        # Read in raw data for both ccds
        raw = [CCDData.read(f, hdu=x+1, unit='adu') for x in range(nccd)]
        print('Processing', frametype, f)

        # Process raw data for each ccd
        # notes: Nora had oscan, oscan_model, and trim set to these values;
        #        I haven't tried changing them
        proc = [ccdproc.ccd_process(x, oscan='[3:22,:]',  # x.header['BIASSEC'],
                                    oscan_model=models.Chebyshev1D(3),
                                    trim=x.header['TRIMSEC'],
                                    master_bias=mbias[k],
                                    gain=gain*u.electron/u.adu,
                                    readnoise=rdnoise*u.electron)
                for k, x in enumerate(raw)]

        # Write out processed data to temporary files
        hdu = fits.HDUList([fits.PrimaryHDU(header=hdr)])
        for x in proc:
            hdu.append(fits.ImageHDU(x.data, header=x.header))
        hdu.writeto('tmpfile'+str(j)+'.fits', overwrite=True)
    tmp = glob('tmpfile*fits')

    print('Combining processed', frametype, ' files')
    # Note: Nora also had these values, I didn't change them
    master = [ccdproc.combine(tmp, hdu=x+1, unit=u.electron, method='median',
                              sigma_clip=True, sigma_clip_low_thresh=5,
                              sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median) for x in range(nccd)]

    # Normalize the flat if creating master flat
    if frametype == 'flat':
        print('Normalizing flat')
        for i, x in enumerate(master):
            master[i].data /= master[i].data.max()

    # Write out the master flat or bias file
    master_hdu = fits.HDUList([fits.PrimaryHDU(header=hdr)])
    for x in master:
        master_hdu.append(fits.ImageHDU(x.data, header=x.header))
    print('Writing master ', frametype, ' to ', outputdir+fname)
    master_hdu.writeto(outputdir+fname, overwrite=True)

    # Remove all temporary files
    for f in tmp:
        os.remove(f)

    return master, outputdir+fname


def get_filters_from_header(flist):
    """Loop through all files in flist and get filter from header.

    Parameters
    ----------
    flist (array) : list of files including full path for each

    Returns
    -------
    filt_list (array) : list of filters for each file in flist
    """
    filt_list = []
    for i, f in enumerate(flist):
        with fits.open(f) as file_open:
            hdr = file_open[0].header
            filt = hdr['FILTER2']
            filt_list.append(filt)
    return filt_list


def get_filters(filt_list):
    """Get correct filter names.

    Parameters
    ----------
    filt_list (array) : list of filters (most likely 'g' or 'r' etc)

    Returns
    -------
    new_filt (array) : list of updated filters (now 'Sloan_g' etc)
    """
    new_filt = []
    for filt in filt_list:
        if 'Sloan' not in filt:
            new_filt.append('Sloan_'+filt)
        else:
            new_filt.append(filt)
    return new_filt


def do_astrometry(gaia, sources, wcs_ref, ima):
    """Correct the astrometry of an image.

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
    # get skycoord formatted ra/dec of refrence gaia stars
    sky_ref_radec_init = SkyCoord(list(zip(np.array(gaia['ra']),
                                           np.array(gaia['dec']))),
                                  frame='fk5', unit='deg')

    # Convert gaia ra/dec into pixels using the original wcs info
    sky_ref_x = []
    sky_ref_y = []
    good_ref = []

    try:
        x, y = wcs_ref.all_world2pix(gaia['ra'], gaia['dec'], 1, maxiter=100,
                                     tolerance=5e-4, detect_divergence=True,
                                     adaptive=True, quiet=True)
    except wcs.NoConvergence as e:
        print("Indices of diverging points: {0}"
              .format(e.divergent))
        print("Indices of poorly converging points: {0}"
              .format(e.slow_conv))
        print([i for i in range(len(gaia['ra']))
               if i not in e.divergent and i not in e.slow_conv])
        print(len(gaia['ra']))
        # print("Best solution:\n{0}".format(e.best_solution))
        # print("Achieved accuracy:\n{0}".format(e.accuracy))
        raise e

    for j in range(len(gaia['ra'])):

        # exclude sources that are too close to the edge;
        # I randomly chose 20 as the limit
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
    cutoff = np.median(d2d.arcsec)  # TO DO: UPDATE THIS #############
    good = np.where(d2d.arcsec < cutoff)[0]
    use_ref_radec = sky_ref_radec[idx][good]
    use_img_xy = np.array([sky_img_xy[0][good], sky_img_xy[1][good]])

    # plt.figure()
    # plt.hist(d2d.arcsec, bins=20)

    # Get the new wcs information
    # format of lists required by fit_wcs_from_points:
    # use_img_xy = np.array([[xlist],[ylist]])
    # use_ref_radec=SkyCoord([(ra,dec),(ra,dec)...],frame='fk5',unit='deg')
    # NOTE on sip_degree=2: not sure why I picked this but it seems to work?
    w = fit_wcs_from_points(xy=use_img_xy, world_coords=use_ref_radec,
                            projection='TAN', sip_degree=2)

    return w, good


def write_ccd(hdr, hdrs, sci_final, outputdir, imafile, root, filt):
    """Write out 2 ccds to separate files."""
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(header=hdr))
    new_hdul.append(fits.ImageHDU(data=sci_final[0].data, header=hdrs[0]))
    new_hdul.writeto(outputdir+imafile.split('/')[-1][:-5] +
                     '_'+root+filt+'_ccd1.fits', overwrite=True)

    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(header=hdr))
    new_hdul.append(fits.ImageHDU(data=sci_final[1].data, header=hdrs[1]))
    new_hdul.writeto(outputdir+imafile.split('/')[-1][:-5] +
                     '_'+root+filt+'_ccd2.fits', overwrite=True)


def combine_files(outputdir, root, filt):
    """Median combine a list of images.

    Parameters
    ----------
    outputdir (str) : directory files will be written to
    root (str) : contains object name
    filt (str) : filter of the images

    Returns
    -------
    None
    """
    combine_list = [outputdir+o for o in os.listdir(outputdir)
                    if o.endswith(root+filt+'_ccd1.fits')]
    print('    combining files for ccd1:', combine_list)
    ima = ccdproc.combine(combine_list, method='median', sigma_clip=True,
                          sigma_clip_func=np.ma.median, unit='adu')
    ima.write(outputdir+root+filt+'_final_ccd1.fits', overwrite=True)

    combine_list = [outputdir+o for o in os.listdir(outputdir)
                    if o.endswith(root+filt+'_ccd2.fits')]
    print('    combining files for ccd2', combine_list)
    ima = ccdproc.combine(combine_list, method='median', sigma_clip=True,
                          sigma_clip_func=np.ma.median, unit='adu')
    ima.write(outputdir+root+filt+'_final_ccd2.fits', overwrite=True)


#############################################################################


def plot_sources(wcs_sources, wcs_gaia, ima, title, sources, gaia, good):
    """Plot circles around detected sources in image and gaia sources.

    Parameters
    ----------
    wcs_sources (?) : wcs projection information/object for original image
    wcs_gaia (?) : wcs projection information/object for gaia sources?
        Note: wcs gaia is actually the same
    ima (array) : image data to plot
    title (str) : title for the plot
    sources (array) : table containing x and y positions of detected stars in ima
    gaia (array) : table containing ra and dec of stars in gaia catalog
    good (array) : indices for 'sources' that were matched with the gaia catalog

    Returns
    -------
    None
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    ax = plt.subplot(projection=wcs_sources)
    ax.imshow(ima, cmap='gray', norm=ImageNormalize(ima, vmin=0,
                                                    vmax=5, clip=True))

    all_c = []
    all_label = []
    for j in range(len(gaia)):
        # Convert gaia ra/dec to pixel space using wcs_gaia projection
        x, y = (wcs_gaia.all_world2pix(gaia['ra'][j], gaia['dec'][j], 1,
                                       maxiter=100, tolerance=5e-4, quiet=True,
                                       detect_divergence=True, adaptive=True))
        c = patches.Circle((x, y), radius=5, ec='r', fc='none',
                           lw=2, alpha=0.8,
                           label='Reference: RA = %f, Dec = %f' %
                           (gaia['ra'][j], gaia['dec'][j]),
                           picker=False)
        ax.add_patch(c)
        if j == 0:
            all_c.append(c)
            all_label.append('Gaia sources')

    test = []
    for k in range(len(sources)):
        if k in good:
            color = 'g'  # matches with gaia
            leg = 'matched sources'
        else:
            color = 'b'  # not matched with gaia
            leg = 'unmatched sources'
        c = patches.Circle((sources['xcentroid'][k],
                            sources['ycentroid'][k]),
                           radius=5, ec=color,
                           fc='none', lw=2, alpha=0.8,
                           label='Position: x = %f, y = %f' %
                           (sources['xcentroid'][k],
                            sources['ycentroid'][k]), picker=False)
        ax.add_patch(c)
        if color not in test:
            all_c.append(c)
            all_label.append(leg)
            test.append(color)

        ax.legend(all_c, all_label)
        ax.set_title(title)


def plot_cr(sci, cr_mask, sci_clean, bpm_mask):
    """Plot cosmic ray cleaned image vs original and compare BPM to CR mask.

    Parameters
    ----------
    sci (2xNxN array) : contains 2 ccds; original image before cosmic ray cleaning
    cr_mask (2xNxN array) : contains 2 ccds; shows pixels that have been masked
        as crs
    sci_clean (2xNxN array) : containes 2 ccds; image after cosmic ray cleaning
    bpm_mask (2xNxN array) : contains 2 ccds: shows pixels that have been masked
        as bad

    Returns
    -------
    None
    """
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.ravel()
    axes[0].imshow(sci[1].data, cmap='gray', interpolation='none', origin='lower')
    axes[0].set_title("Original Image")

    axes[1].imshow(sci_clean[1], cmap='gray', interpolation='none', origin='lower')
    axes[1].set_title("Cosmic ray cleaned image")

    axes[2].imshow(bpm_mask[1].data, cmap='gray',
                   interpolation='none', origin='lower')
    axes[2].set_title("Bad pixel mask")

    axes[3].imshow(cr_mask[1], cmap='gray', interpolation='none', origin='lower')
    axes[3].set_title("Cosmic ray mask")
