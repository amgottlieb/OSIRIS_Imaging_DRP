"""Contains functions used in OSIRIS_imaging_pipeline.py."""

from __future__ import print_function
from glob import glob
import argparse
import shutil
import os
from astropy.io import fits
import ccdproc
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
from astropy.modeling import models
from photutils.segmentation import make_source_mask


def print_both(file, *args):
    """Print a statement (args) to the command line and to a file."""
    toprint = ' '.join([str(arg) for arg in args])
    print(toprint)
    file.write(toprint+' \n')


def read_args():
    """Create wrapper to parse and customize command-line arguments."""
    # print('READ ARGS')
    parser = parse_args()
    args, use_slash = setup_args(parser)
    # print('ARGS', args)
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
        help='Set the working directory. Default: ./ ')

    add('--outputdir', dest='outputdir', default='./output/',
        help='Set the directory to put results in. Default: ./output/')

    add('--filter', dest='filt', default='r',
        help='Filter of the images to be reduced. Can be multiple ex: g,i,z \
        Default: r')

    add('--reduce_obj', dest='reduce_obj', default=2,
        help='Choose which object to reduce: 0 for science, 1 for standard, \
            2 for both. Default: 2')

    add('--doobslog', dest='doobslog', action='store_true',
        help='Create an observing log without going through the full reduction.')

    add('--dobias', dest='dobias', action='store_true',
        help='Make the master bias instead of reading it from input file.')

    add('--bias', dest='biasfile', default='MasterBias',
        help='Name of the input master bias file. Default: MasterBias')

    add('--doflat', dest='doflat', action='store_true',
        help='Make the master flat instead of reading it from input file.')

    add('--flat', dest='flatfile', default='MasterFlat',
        help='Name of the input master flat file. Default: MasterFlat')

    add('--maskflat', dest='maskflat', default=False, help='If there are stars \
        in your flat that are not removed with the normal method, set this value \
            to true to mask out the stars.')

    add('--dobpm', dest='domask', action='store_true',
        help='Create the master bad pixel mask instead of reading it from the \
        input file.')

    add('--bpm', dest='maskfile', default='MasterBPM',
        help='Name of the input mask file.  Default: MasterBPM')

    add('--docalib', dest='docalib', action='store_true',
        help='Apply calibrations (bias, flat, bpm) to science images.')

    add('--docrmask', dest='docrmask', action='store_true',
        help='Run cosmic ray rejection.')

    add('--doskysub', dest='doskysub', action='store_true',
        help='Perform sky subtraction.')

    add('--dostack', dest='dostack', action='store_true',
        help='Align all images with eachother and median combine them into one \
            image.')

    add('--dowcs', dest='dowcs', action='store_true',
        help='Improve the astrometric solution.')

    add('--dointeractive', dest='dointeractive', action='store_true',
        help='Manually choose stars to use for performing astrometry.')

    add('--seeing', dest='seeing', default=1.0,
        help='Seeing of the images. Check the header or the QA \
            (quality assurance) file.')

    add('--logfile', dest='logfile', default='log.txt',
        help='Name of the file which contains all print statements from the \
            pipeline.')

    add('--dooverwrite', dest='dooverwrite', default=False,
        help='Overwrites any existing files. Default: False')

    return parser


def setup_args(parser):
    """Any manipulation of the arguments that may be required."""
    args = parser.parse_args()
    # check for \ or / depending on windows or linux/mac
    # gtcsetup.print_both(log_fname, 'system',sys.platform)
    # print('    CHECKING for slashes in ', args.workdir)
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


def sort_files(obj, raw_path, raw_list, use_slash, outputdir, dooverwrite,
               diagnostic_dir, log_fname):
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
    log_fname (str) : full path to log file

    Returns
    -------
    bias_list, flat_list, sci_list, std_list (array) : lists of all files in
    each directory; each file contains full path
    """
    # NOTE: raw_list already has file path in it
    print_both(log_fname, '    Sorting files...')

    # this function creates the directories if they don't already exist
    paths = setup_directories(raw_path, use_slash, outputdir, diagnostic_dir,
                              dooverwrite, log_fname)
    bias_path, flat_path, sci_path, std_path = paths[0], paths[1], paths[2], paths[3]

    # Go through each file and determine what folder it should go into
    # based on its file name which contains 'Bias', 'Flat' etc
    for f in raw_list:
        if 'Bias' in f:
            move_path = bias_path
        elif 'Flat' in f:
            move_path = flat_path
        elif 'Image' in f:
            # Check which object is in the image (target vs standard)
            with fits.open(f) as file_open:
                hdr = file_open[0].header
                target = hdr['OBJECT'].replace(' ', '')
                print_both(log_fname, '        Objectid', obj,
                           'check against target', target)
                if target == obj:
                    move_path = sci_path
                else:
                    move_path = std_path
        else:
            move_path = ''
        # Move the file into its ** path
        if move_path != '':
            shutil.copy2(f, move_path)

    # Get lists of all files in each folder
    bias_list1 = [i.replace(os.sep, '/') for i in glob(bias_path+'0*.fits*')]
    flat_list1 = [i.replace(os.sep, '/') for i in glob(flat_path+'0*.fits*')]
    std_list1 = [i.replace(os.sep, '/') for i in glob(std_path+'0*.fits*')]
    sci_list1 = [i.replace(os.sep, '/') for i in glob(sci_path+'0*.fits*')]

    # Rename each file to be 'bias1.fits, bias2.fits' etc
    bias_list = rename(bias_list1, 'bias')
    flat_list = rename(flat_list1, 'flat')
    sci_list = rename(sci_list1, 'sci')
    std_list = rename(std_list1, 'std')

    return bias_list, flat_list, sci_list, std_list, paths[6:]


def rename(files, name):
    """Rename a list of files to have a given name and index number.

    Parameters
    ----------
    files (list) : list of strings containing the full path of files
    name (str)   : name to rename the file to (ex, 'bias')

    Returns
    -------
    new_files (list) : list of strings containing the full path of the renmaed files
    """
    new_files = []
    for i, f in enumerate(files):
        # Example f: /user/documents/0003029174-20210706-OSIRIS-OsirisBias.fits

        # Get the file path as a list ['user','documents']
        split_path = f.split('/')[:-1]

        # Join the file path with / at the end:  /user/documents/
        path = ('/').join(split_path)+'/'

        new_fname = path+name+str(i+1)+'.fits'  # ex bias1.fits

        # Rename the file
        shutil.move(f, new_fname)

        # Save the new name to a list
        new_files.append(new_fname)

    return new_files


def setup_directories(raw_path, use_slash, outputdir, diagnostic_dir,
                      dooverwrite, log_fname):
    r"""Create necessary folders and copy all files into the copy folder.

    Parameters
    ----------
    raw_path (str)     : full path to raw data (C:\User etc)
    use_slash (str)    : \ or \\ depending on os, windows or linux
    outputdir (str)    : full path to folder where output files are stored
                         (master files, processed data)
    dooverwrite (bool) : if True, overwrite current existing directories;
                         WILL LOSE OUTPUT FOLDER IF TRUE
    log_fname (str)    : full path to log file

    Returns
    -------
    bias_path (str) : full path to bias folder (C:\User etc)
    flat_path (str) : full path to flat folder (C:\User etc)
    sci_path (str)  : full path to science folder (C:\User etc)
    std_path (str)  : full path to standard folder (C:\User etc)
    """
    bias_path = raw_path+'bias'+use_slash
    flat_path = raw_path+'flat'+use_slash
    sci_path = raw_path+'object'+use_slash
    std_path = raw_path+'standard'+use_slash
    # Copy all data to this 'copy' directory in case you need to start over
    copy_path = raw_path+'copy'+use_slash
    # more diagnostic folders
    calib_path = diagnostic_dir+'1-calib'+use_slash
    bpm_path = diagnostic_dir+'2-bpm'+use_slash
    crmask_path = diagnostic_dir+'3-crmask'+use_slash
    skymap_path = diagnostic_dir+'4-skymap'+use_slash
    astrom_path = diagnostic_dir+'5-astrometry'+use_slash

    paths = [bias_path, flat_path, sci_path, std_path, outputdir,
             diagnostic_dir, calib_path, bpm_path, crmask_path,
             skymap_path, astrom_path]

    print_both(log_fname, '        Overwriting folders?', dooverwrite)

    # If there isn't a copy directory already, make it and copy all files into it
    if not os.path.exists(copy_path):
        print_both(log_fname, '        Creating copy directory')
        os.makedirs(copy_path)
        print_both(log_fname, '        Copying all fits files into copy directory')
        copy_list = [i.replace(os.sep, '/') for i in glob(raw_path+'0*.fits*')]
        for f in copy_list:
            shutil.copy2(f, copy_path)

    # Go through the remaining folders and create them if they don't already exist
    # or you want to overwrite them
    for i, path1 in enumerate(paths):
        if not os.path.exists(path1):
            os.makedirs(path1)
            print_both(log_fname, '        Creating directory', path1)
        elif dooverwrite == 'True' or dooverwrite is True:
            shutil.rmtree(path1)
            os.makedirs(path1)
            print_both(log_fname, '        Overwriting directory', path1)
        else:
            print_both(log_fname, '        ', path1,
                       ': Directory already exists; not overwriting.')

    return paths


def print_obslog_folder(folder, out, verbose, log_fname):
    """Print observing log for a specific folder.

    Parameters
    ----------
    folder (str) : full path to folder that you want to print the observing log for
    out (str) : observing log output file
    log_fname (str) : full path to log file

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
            print_both(log_fname, msg)


def print_obslog(obj, flt, bias, flat, sci, std, log_fname, verbose=True):
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
    log_fname (str) : full path to log file

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
        print_both(log_fname, head)

    print_obslog_folder(bias, out, verbose, log_fname)
    print_obslog_folder(flat, out, verbose, log_fname)
    print_obslog_folder(sci, out, verbose, log_fname)
    print_obslog_folder(std, out, verbose, log_fname)

    return


def update_trimsec(trimsec):
    """
    There are more bad columns in ccd1 than the header says.

    This function changes one of the trimsec x values by 3 pixels so that
    these columns will be removed.

    Parameters
    ----------
    trimsec (str) : string from the header of a fits file,
                    looks like: '[26:1050,1:2000]'

    Returns
    -------
    join2 (str) : the updated string in the same format as trimsec
                    ex: '[26:1047,1:2000]'
    """
    split1 = trimsec.split(',')  # ['[26:1050','1:2000']
    split2 = [x.split(':') for x in split1]  # [['[26', '1050'],['1','2000]']]

    if int(split2[0][1]) > 1046:
        split2[0][1] = str(int(split2[0][1])-3)  # '1047'

    join1 = [':'.join(x) for x in split2]  # ['[26:1047','1:2000]']
    join2 = ','.join(join1)  # '[26:1047, 1:2000]'

    return join2


def createMaster(full_flist, frametype, nccd, mbias, gain, rdnoise, outputdir,
                 name, filt, log_fname):
    """Given list of files, median combine them to create a master file.

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
    log_fname (str) : full path to log file

    Returns
    -------
    master (array) :[ccd1, ccd2] processed data
    fname (str) : file name (including full path) that master was written to
    """
    print_both(log_fname, full_flist)

    # Get all filters present (really only used for flat)
    all_filters = get_filters_from_header(full_flist)

    # If creating a master bias, don't add a filter to the filename
    if frametype == 'bias':
        add = ''
        flist = full_flist
    # If creating a master flat, add a filter to the filename and select only
    # the relevant filters
    else:
        add = filt
        use = np.where(np.array(all_filters) == filt)[0]
        flist = np.array(full_flist)[use]

    # Add filter to filename if necessary
    fname = name+add+'.fits'

    all_exps = []
    # Go through each file, process it, write out as temp file
    for j, f in enumerate(flist):
        # os.system("tar -xzf {}".format(f))
        # f=f.replace('.gz','')
        with fits.open(f) as fits_open:
            hdr = fits_open[0].header
            exp = hdr['EXPTIME']
            all_exps.append(exp)

        # Read in raw data for both ccds
        raw = [CCDData.read(f, hdu=x+1, unit='adu') for x in range(nccd)]
        print_both(log_fname, '     Processing', frametype, f)

        # Modifying the trimsec for ccd1 b/c there are 3 columns on the right
        # that are bad and they mess up the background estimation
        trimsec = [x.header['TRIMSEC'] for x in raw]
        trimsec[0] = update_trimsec(trimsec[0])

        # Process raw data for each ccd
        # Notes: Nora had oscan, oscan_model, and trim set to these values;
        #        I haven't tried changing them
        proc = [ccdproc.ccd_process(x, oscan='[3:22,:]',  # x.header['BIASSEC'],
                                    oscan_model=models.Chebyshev1D(3),
                                    trim=trimsec[k],  # x.header['TRIMSEC'],
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

    print_both(log_fname, '   Combining processed', frametype, ' files')

    # Normalize the flat if creating master flat
    if frametype == 'flat':
        # if maskflat is True:
        #     master = []
        #     for ccd in range(nccd):
        #         all_masked_images = []

        #         for i, f in enumerate(tmp):
        #             tmp_fits = [CCDData.read(f, hdu=x+1, unit='adu', memmap=False)
        #                         for x in range(nccd)]
        #             tmp_fits[ccd].data /= all_exps[i]

        #             # tmp_fits = [
        #             #     x.data.divide(all_exps[i]*u.second, handle_meta='first found') for x in tmp_fits_init]
        #             # sci_final = [[x.divide(all_headers[i][0]['EXPTIME']*u.second,
        #             #                        propagate_uncertainties=True,
        #             #                        handle_meta='first_found')
        #             # for i, x in enumerate(imgs)] for imgs in sky_sub]

        #             mask = make_source_mask(tmp_fits[ccd].data, nsigma=3, npixels=3)

        #             new_hdul = fits.HDUList(fits.PrimaryHDU(header=hdr))
        #             new_hdul.append(fits.ImageHDU(data=mask.astype(int), header=hdr))
        #             new_hdul.writeto(outputdir+fname[:-5]+'_mask' +
        #                              '_ccd'+str(ccd)+'.fits', overwrite=True)

        #             # mask out the stars in the flats before median combining
        #             # sometimes they're very close to the background value that
        #             # sigma clipping does not remove them
        #             masked_img = CCDData(tmp_fits[ccd].data,  # mask=mask,
        #                                  header=tmp_fits[ccd].header, unit='adu')
        #             all_masked_images.append(masked_img)

        #         masked_median_img = ccdproc.combine(all_masked_images,
        #                                             method='median')
        #         master.append(masked_median_img)

        #     print(master[0].header)
        #     print(master[1].header)

        master = [ccdproc.combine(tmp, hdu=x+1, unit=u.electron, method='median',
                                  sigma_clip=True, sigma_clip_low_thresh=5,
                                  sigma_clip_high_thresh=5,
                                  sigma_clip_func=np.ma.median) for x in range(nccd)]

        print_both(log_fname, 'Normalizing flat')
        for i, x in enumerate(master):
            # print(master[i].data.max)
            # print(np.nanmax(master[i].data))
            # print(master[i].data)
            master[i].data /= np.nanmax(master[i].data)  # master[i].data.max()

    else:
        # Note: Nora also had these values, I didn't change them
        # everything above/below 5 sigma will be clipped
        # if you increase to 10 sigma, then less will be clipped
        master = [ccdproc.combine(tmp, hdu=x+1, unit=u.electron, method='median',
                                  sigma_clip=True, sigma_clip_low_thresh=5,
                                  sigma_clip_high_thresh=5,
                                  sigma_clip_func=np.ma.median) for x in range(nccd)]

    # Write out the master flat or bias file
    master_hdu = fits.HDUList([fits.PrimaryHDU(header=hdr)])
    for x in master:
        master_hdu.append(fits.ImageHDU(x.data, header=x.header))
    print_both(log_fname, 'Writing master ', frametype, ' to ', outputdir+fname)
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


def write_one_image(hdr, hdrs, sci_final, outputdir, imafile, root, filt,
                    log_fname, nccd):
    """Write out images from 2 ccds to separate files.

    Need to 'combine' one file into one master file
    so that the bad pixel mask is actually applied.

    Parameters
    ----------
    hdr (?) : header for whole fits file
    hdrs (?) : header for ccd 1 and 2
    sci_final (CCDData arr) : the 2 ccd images
    outputdir (str) : directory files will be written to
    imafile (str) : base name that the images will be written to
    root (str) : contains object name
    filt (str) : filter of the images
    log_fname (str) : full path to log file
    nccd (int) : number of ccds; normally 2

    Returns
    -------
    None
    """
    # Not sure why I chose these values
    sig_clip_low = 5
    sig_clip_high = 5

    # NOTE: sigma clipping is getting rid of vertical streaks in the science
    # images that are not in the master flat which is what the bad pixel mask
    # is based on
    master = [ccdproc.combine([sci_final[x]], hdu=x+1, unit=u.electron,
                              method='median',
                              sigma_clip=True, sigma_clip_low_thresh=sig_clip_low,
                              sigma_clip_high_thresh=sig_clip_high,
                              sigma_clip_func=np.ma.median) for x in range(nccd)]

    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(header=hdr))
    new_hdul.append(fits.ImageHDU(data=master[0].data, header=hdrs[0]))
    new_hdul.writeto(outputdir+imafile.split('/')[-1][:-5] +
                     '_'+root+filt+'_ccd1.fits', overwrite=True)

    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(header=hdr))
    new_hdul.append(fits.ImageHDU(data=master[1].data, header=hdrs[1]))
    new_hdul.writeto(outputdir+imafile.split('/')[-1][:-5] +
                     '_'+root+filt+'_ccd2.fits', overwrite=True)

    print_both(log_fname, '      Writing files to',
               outputdir+imafile.split('/')[-1][:-5] +
               '_'+root+filt+'_ccd 1 and 2.fits')


def write_ccd(hdr, hdrs, sci_final, outputdir, imafile, root, filt, log_fname):
    """Write out 2 ccds to separate files.

    Parameters
    ----------
    hdr (?) : header for whole fits file
    hdrs (?) : header for ccd 1 and 2
    sci_final (CCDData arr) : the 2 ccd images
    outputdir (str) : directory files will be written to
    imafile (str) : base name that the images will be written to
    root (str) : contains object name
    filt (str) : filter of the images
    log_fname (str) : full path to log file

    Returns
    -------
    None

    """
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
    print_both(log_fname, '      Writing files to',
               outputdir+imafile.split('/')[-1][:-5] +
               '_'+root+filt+'_ccd 1 and 2 .fits')


def get_header_info(obj, filt):
    """Get all headers and filters, and adjust the CRPIX keywords.

    Also check which files match the current reduction filter.

    Parameters
    ----------
    obj (list) : list of file names
    filt (str) : current filter of the images (ex: Sloan_z)

    Returns
    -------
    all_headers (list) : list of the headers from each of the 3 extensions
    all_filts (list)   : list of True/False; True if the filter matches the
                         current reduction filter

    """
    all_filts = []
    all_headers = []
    for j, f in enumerate(obj):

        with fits.open(f) as fits_open:
            # Note: because images have been trimmed, the wcs information
            # (crpix) is no longer correct; subtract off the number of pixels
            # trimmed to get close to the original wcs
            hdr = fits_open[0].header

            hdr1 = fits_open[1].header
            trim1 = hdr1['TRIMSEC'].split(':')[0][1:]
            hdr1['CRPIX1'] = float(hdr1['CRPIX1'])-float(trim1)

            hdr2 = fits_open[2].header
            trim2 = hdr2['TRIMSEC'].split(':')[0][1:]
            hdr2['CRPIX1'] = float(hdr2['CRPIX1'])-float(trim2)

            hdr_filt = hdr['FILTER2']

            all_headers.append([hdr, hdr1, hdr2])

        if hdr_filt != filt:
            all_filts.append(False)
        else:
            all_filts.append(True)

    return all_headers, all_filts


def read_in_files(path, root, log_fname):
    """Read in the separate ccd1/2 files as CCDData objects.

    Parameters
    ----------
    path (str) : full path to files
    root (str) : contains object name
    log_fname (str) : opened log file

    Returns
    -------
    all_images (list) : list of pairs of images [[ccd1,ccd2], [ccd1,ccd2]...]
    """
    all_images = []
    if os.path.isdir(path):
        ccd1 = [f for f in os.listdir(path) if root in f and 'ccd1' in f]
        ccd2 = [f for f in os.listdir(path) if root in f and 'ccd2' in f]

        for i in range(len(ccd1)):
            print_both(log_fname, 'Reading in ', ccd1[i])
            print_both(log_fname, 'Reading in ', ccd2[i])
            img1 = CCDData.read(path+ccd1[i], hdu=1, unit=u.electron)
            img2 = CCDData.read(path+ccd2[i], hdu=1, unit=u.electron)
            all_images.append([img1, img2])

        return all_images
    else:
        print_both(log_fname, 'FOLDER DOES NOT EXIST:', path)
        return []


def save_times(tstart, tnow, action, text, times, actions, log_fname):
    """Print out the current runtime; save it and the current reduction step.

    Parameters
    ----------
    tstart (time/float) : time at the start of the reduction
    tnow (time/float) : current time after completeing the current reduction step
    action (str) : current reduction step
    text (str) : statement about the current reduction step to print
    times (list) : list of completion times of all previous reduction steps
    actions (list) : list of completed reduction step names
    log_fname (?) : opened log file

    Returns
    -------
    times (list) : updated list of completion times of all previous reduction steps
    actions (list) : updated list of completed reduction step names
    """
    print_both(log_fname, 'Time after '+text,
               (tnow-tstart)/60., 'min')
    times.append((tnow-tstart)/60.)
    actions.append(action)

    return times, actions


def print_pipeline_times(times, actions, log_fname):
    """Print out each step and how long it took to run.

    Parameters
    ----------
    times (list) : list of completion times of all reduction steps
    actions (list) : list of completed reduction step names
    log_fname (?) : opened log file

    Returns
    -------
    None
    """
    print(times, actions)

    print_both(log_fname, 'Action   Time to complete Action    Total pipeline time')
    for i in range(len(actions)):
        if i == 0:
            print_both(log_fname, actions[i], '-', round(
                times[i], 4))
        else:
            print_both(log_fname, actions[i], round(times[i]-times[i-1], 4),
                       round(times[i], 4))
