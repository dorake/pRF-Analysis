from __future__ import division, print_function


def convert_edf_2_hdf5(edf_file, low_pass_pupil_f=6.0, high_pass_pupil_f=0.01):
    """converts the edf_file to hdf5 using hedfpy

    Requires hedfpy

    Parameters
    ----------
    edf_file : str
        absolute path to edf file.
    low_pass_pupil_f : float
        low pass cutoff frequency for band-pass filtering of the pupil signal
    high_pass_pupil_f : float
        high pass cutoff frequency for band-pass filtering of the pupil signal
    Returns
    -------
    hdf5_file : str
        absolute path to hdf5 file.
    """

    import os
    import os.path as op
    from hedfpy import HDFEyeOperator
    import tempfile

    tempdir = tempfile.mkdtemp()
    temp_edf = op.join(tempdir, op.split(edf_file)[-1])

    os.system('cp ' + edf_file + ' ' + temp_edf)

    hdf5_file = op.join(tempdir, op.split(
        op.splitext(edf_file)[0] + '.h5')[-1])
    alias = op.splitext(op.split(edf_file)[-1])[0]

    ho = HDFEyeOperator(hdf5_file)
    ho.add_edf_file(temp_edf)
    ho.edf_message_data_to_hdf(alias=alias)
    ho.edf_gaze_data_to_hdf(
        alias=alias, pupil_hp=high_pass_pupil_f, pupil_lp=low_pass_pupil_f)

    return hdf5_file


def combine_eye_hdfs_to_nii_hdf(nii_hdf5_file, eye_hdf_filelist, new_alias='eye'):
    import os.path as op
    import tables as tb

    nii_hf = tb.open_file(nii_hdf5_file, 'a')
    eye_hfs = [tb.open_file(ef, 'r') for ef in eye_hdf_filelist]
    eye_group_names = [op.splitext(op.split(ef)[-1])[0]
                       for ef in eye_hdf_filelist]

    try:
        nii_group = nii_hf.get_node("/", name=new_alias, classname='Group')
    except tb.NoSuchNodeError:
        print('Adding group ' + new_alias + ' to ' + nii_hdf5_file)
        nii_group = nii_hf.create_group("/", new_alias, new_alias)

    for ef, en in zip(eye_hfs, eye_group_names):
        ef.copy_node(where='/' + en, newparent=nii_group,
                     newname=en, overwrite=True, recursive=True)
        ef.close()

    nii_hf.close()

    return nii_hdf5_file


def mask_nii_2_hdf5(in_files, mask_files, hdf5_file, folder_alias):
    """masks data in in_files with masks in mask_files,
    to be stored in an hdf5 file

    Takes a list of 3D or 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    in_files : list
        list of absolute path to functional nifti-files.
        all nifti files are assumed to have the same ndim
    mask_files : list
        list of absolute path to mask nifti-files.
        mask_files are assumed to be 3D
    hdf5_file : str
        absolute path to hdf5 file.
        folder_alias : str
                name of the to-be-created folder in the hdf5 file.

    Returns
    -------
    hdf5_file : str
        absolute path to hdf5 file.
    """

    import nibabel as nib
    import os.path as op
    import numpy as np
    import tables

    success = True

    mask_data = [np.array(nib.load(mf).get_data(), dtype=bool)
                 for mf in mask_files]
    nifti_data = [nib.load(nf).get_data() for nf in in_files]

    mask_names = [op.split(mf)[-1].split('_vol.nii.gz')[0]
                  for mf in mask_files]
    nifti_names = [op.split(nf)[-1].split('.nii.gz')[0] for nf in in_files]

    h5file = tables.open_file(hdf5_file, mode="a", title=hdf5_file)
    # get or make group for alias folder
    try:
        folder_alias_run_group = h5file.get_node(
            "/", name=folder_alias, classname='Group')
    except tables.NoSuchNodeError:
        print('Adding group ' + folder_alias + ' to this file')
        folder_alias_run_group = h5file.create_group(
            "/", folder_alias, folder_alias)

    for (roi, roi_name) in zip(mask_data, mask_names):
        # get or make group for alias/roi
        try:
            run_group = h5file.get_node(
                where="/" + folder_alias, name=roi_name, classname='Group')
        except tables.NoSuchNodeError:
            print('Adding group ' + folder_alias +
                  '_' + roi_name + ' to this file')
            run_group = h5file.create_group(
                "/" + folder_alias, roi_name, folder_alias + '_' + roi_name)

        h5file.create_array(run_group, roi_name, roi, roi_name +
                            ' mask file for reconstituting nii data from masked data')

        for (nii_d, nii_name) in zip(nifti_data, nifti_names):
            print('roi: %s, nifti: %s' % (roi_name, nii_name))
            n_dims = len(nii_d.shape)
            if n_dims == 3:
                these_roi_data = nii_d[roi]
            elif n_dims == 4:   # timeseries data, last dimension is time.
                these_roi_data = nii_d[roi, :]
            else:
                print(
                    "n_dims in data {nifti} do not fit with mask".format(nii_name))
                success = False

            h5file.create_array(run_group, nii_name, these_roi_data,
                                roi_name + ' data from ' + nii_name)

    h5file.close()

    return hdf5_file


def roi_data_from_hdf(data_types_wildcards, roi_name_wildcard, hdf5_file, folder_alias):
    """takes data_type data from masks stored in hdf5_file

    Takes a list of 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    data_types_wildcards : list
        list of data types to be loaded.
        correspond to nifti_names in mask_2_hdf5
    : str
        wildcard for masks. 
        corresponds to mask_name in mask_2_hdf5.
    hdf5_file : str
        absolute path to hdf5 file.
    folder_alias : str
        name of the folder in the hdf5 file from which data
        should be loaded.

    Returns
    -------
    output_data : list
        list of numpy arrays corresponding to data_types and roi_name_wildcards
    """
    import tables
    import itertools
    import fnmatch
    import numpy as np

    h5file = tables.open_file(hdf5_file, mode="r")

    try:
        folder_alias_run_group = h5file.get_node(
            where='/', name=folder_alias, classname='Group')
    except NoSuchNodeError:
        # import actual data
        print('No group ' + folder_alias + ' in this file')
        # return None

    all_roi_names = h5file.list_nodes(
        where='/' + folder_alias, classname='Group')
    roi_names = [
        rn._v_name for rn in all_roi_names if roi_name_wildcard in rn._v_name]
    if len(roi_names) == 0:
        print('No rois corresponding to ' +
              roi_name_wildcard + ' in group ' + folder_alias)
        # return None

    data_arrays = []
    for roi_name in roi_names:
        try:
            roi_node = h5file.get_node(
                where='/' + folder_alias, name=roi_name, classname='Group')
        except tables.NoSuchNodeError:
            print('No data corresponding to ' +
                  roi_name + ' in group ' + folder_alias)
            pass
        all_data_array_names = h5file.list_nodes(
            where='/' + folder_alias + '/' + roi_name)
        data_array_names = [adan._v_name for adan in all_data_array_names]
        selected_data_array_names = list(itertools.chain(
            *[fnmatch.filter(data_array_names, dtwc) for dtwc in data_types_wildcards]))

        # if sort_data_types:
        selected_data_array_names = sorted(selected_data_array_names)
        if len(data_array_names) == 0:
            print('No data corresponding to ' + str(selected_data_array_names) +
                  ' in group /' + folder_alias + '/' + roi_name)
            pass
        else:
            print('Taking data corresponding to ' + str(selected_data_array_names) +
                  ' from group /' + folder_alias + '/' + roi_name)
            data_arrays.append([])
            for dan in selected_data_array_names:
                data_arrays[-1].append(
                    eval('roi_node.__getattr__("' + dan + '").read()'))

            # stack across timepoints or other values per voxel
            data_arrays[-1] = np.hstack(data_arrays[-1])
    # stack across regions to create a single array of voxels by values (i.e. timepoints)
    all_roi_data_np = np.vstack(data_arrays)

    h5file.close()

    return all_roi_data_np

def get_scaninfo(in_file):
    """ Extracts info from nifti file.

    Extracts affine, shape (x, y, z, t), dynamics (t), voxel-size and TR from
    a given nifti-file.

    Parameters
    ----------
    in_file : nifti-file (*.nii.gz or *.nii)
        Nifti file to extract info from.

    Returns
    -------
    TR : float
        Temporal resolution of functional scan
    shape : tuple
        Dimensions of scan (x, y, z, time)
    dyns : int
        Number of dynamics (time)
    voxsize : tuple
        Size of voxels (x, y, z = slice thickness)
    affine : np.ndarray
        Affine matrix.
    """
    import nibabel as nib

    nifti = nib.load(in_file)
    affine = nifti.affine
    shape = nifti.shape
    dyns = nifti.shape[-1]
    voxsize = nifti.header['pixdim'][1:4]
    TR = float(nifti.header['pixdim'][4])

    return TR, shape, dyns, voxsize, affine

def leave_one_out_lists(input_list):
    """leave_one_out_lists takes creates a list of lists, with each element
    of the input_list left out of the returned lists once, in order.


    Parameters
    ----------
    input_list : list
        list of items, for instance absolute paths to nii files

    Returns
    -------
    output_data : list
        list of lists
    """

    out_lists = []
    for x in input_list:
        out_lists.append([y for y in input_list if y != x])

    return out_lists


def suff_file_name(in_file, suff='_av_loo', extension='.nii.gz'):
    out_file = in_file[:-len(extension)] + suff + extension

    return out_file


def combine_cv_prf_fit_results_all_runs(basedir, fit_dir, avg_folder):
    """combine_fit_results_one_fold combines a per-slice

    Parameters
    ----------
    basedir : str
        absolute path to directory in which to search for the input
        files to the fits that happened on cartesius.
    fit_dir : str
        absolute path to directory in which to find the fit result files.
    avg_folder: 'loo' or 'all'

    Returns
    -------
    output_files : list
        absolute paths to output nifti files.
    """

    import os
    import glob

    avg_files = sorted(glob.glob(os.path.join(basedir, avg_folder, '*.nii.gz')))
    #import pdb ; pdb.set_trace()
    avg_files_no_ext = [os.path.join(
        fit_dir, os.path.split(afl)[-1][:-7]) for afl in avg_files]

    output_files = []
    for basefilename in avg_files_no_ext:
        output_files.append(
            combine_cv_prf_fit_results_one_fold(
                basefilename=basefilename)
        )

    return output_files

def combine_cv_prf_fit_results_one_fold(basefilename):
    """combine_fit_results_one_fold combines a per-slice

    Parameters
    ----------
    basefilename : str
        absolute path to stem of per-slice files.
 
    Returns
    -------
    output_file : str
        absolute path to output nifti file.
    """
    import os
    import nibabel as nb
    import glob
    import numpy as np

    output_file = os.path.join(basefilename + '_est_all.nii.gz')

    if os.path.isfile(output_file):
        print(output_file + ' already exists. deleting....')
        os.remove(output_file)

    in_files = sorted(glob.glob(basefilename + '_est_*.nii.gz'))

    nif = nb.load(in_files[0])
    all_fits = np.zeros(nif.get_data().shape)+np.nan

    for i, aff in enumerate(in_files):

        print(aff)
        index_z = int(aff.split('_')[-1][:-7])

        data_to_put = nb.load(aff).get_data()
        
        all_fits[:,:,index_z,:] = data_to_put[:,:,index_z,:]

    print('creating ' + output_file)
    img = nb.Nifti1Image(all_fits, affine=nif.affine, header=nif.header)
    img.to_filename(output_file)

    return output_file

def select_condition_files(in_files, condition):
    return [f for f in in_files if condition in f]

def pickfirst(files):
    if isinstance(files, list):
        if len(files) > 0:
            return files[0]
        else:
            return files
    else:
        return files

def FS_T1_file(freesurfer_subject_ID, freesurfer_subject_dir):
    # housekeeping function for finding T1 file in FS directory

    import os.path as op
    return op.join(freesurfer_subject_dir, freesurfer_subject_ID, 'mri', 'T1.mgz')

def average_over_runs(in_files, func = 'mean', output_filename = None):
    """Converts data in a nifti-file to percent signal change.

    Takes a list of 4D fMRI nifti-files and averages them. 
    That is, the resulting file has the same dimensions as each
    of the input files. average_over_runs assumes that all nifti
    files to be averaged have identical dimensions.

    Parameters
    ----------
    in_files : list
        Absolute paths to nifti-files.
    func : string ['mean', 'median'] (default: 'mean')
        the function used to calculate the 'average'
    output_filename : str
        path to output filename

    Returns
    -------
    out_file : str
        Absolute path to average nifti-file.
    """

    import nibabel as nib
    import numpy as np
    import os
    import bottleneck as bn

    template_data = nib.load(in_files[0])
    dims = template_data.shape
    affine = template_data.affine
    header = template_data.header
    all_data = np.zeros([len(in_files)]+list(dims))
    
    for i in range(len(in_files)):
        d = nib.load(in_files[i])
        all_data[i] = d.get_data()

    if func == 'mean':
        av_data = all_data.mean(axis = 0)
    elif func == 'median':
        # weird reshape operation which hopeully fixes an issue in which
        # np.median hogs memory and lasts amazingly long
        all_data = all_data.reshape((len(in_files),-1))
        av_data = bn.nanmedian(all_data, axis = 0)
        av_data = av_data.reshape(dims)

    img = nib.Nifti1Image(av_data, affine=affine, header=header)

    if output_filename == None:
        new_name = os.path.basename(in_files[0]).split('.')[:-2][0] + '_av.nii.gz'
        out_file = os.path.abspath(new_name)
    else:
        out_file = os.path.abspath(output_filename)
    nib.save(img, out_file)

    return out_file

def rescale_timeseries(psc_file, rsq_file, output_file):
    import nibabel as nb
    import numpy as np
    from IPython import embed as shell

    npsc = nb.load(psc_file)
    psc_data = npsc.get_data()
    rsq_data = np.abs(nb.load(rsq_file).get_data())

    op_image = nb.Nifti1Image(
        psc_data * rsq_data[..., np.newaxis], header=npsc.header, affine=npsc.affine)
    op_image.to_filename(output_file)

    return output_file


def natural_sort(l):
    import re

    def convert(text): return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): return [convert(c)
                                   for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def add_subject_2_cortex(fs_ID,fs_dir):
    import os
    import shutil
    from cortex.freesurfer import import_subj

    sub_cortex_path = os.path.join(fs_dir,'cortex','db',fs_ID)
    try:
        shutil.rmtree(sub_cortex_path)
    except OSError:
        pass
    import_subj(subject=fs_ID, sname=fs_ID, freesurfer_subject_dir=fs_dir)
    return fs_ID

def add_flat_2_cortex(fs_ID,fs_dir):
    from cortex.freesurfer import import_flat
    import_flat(subject = fs_ID, patch='full', freesurfer_subject_dir=fs_dir)
    return fs_ID

def add_session_2_cortex(epi_2_T1, epi, T1, subject, xfmname, transpose=False):
    from cortex.xfm import Transform

    if transpose:
        epi_file = transpose_epi(epi)
    else:
        epi_file = epi

    xfm = Transform.from_fsl(xfm=epi_2_T1, func_nii=epi_file, anat_nii=T1)
    # Save as pycortex 'coord' transform
    xfm.save(subject, xfmname, 'coord')

    return xfmname

def create_cortical_mask(subject, xfmname, type, epi, temp_dir_cortex):
    import cortex
    import nibabel as nb
    import os#, tempfile

    #tempdir = tempfile.gettempdir()
    try:
        os.makedirs(temp_dir_cortex)
    except OSError:
        pass

    out_file = os.path.join(temp_dir_cortex, 'cortex_%s.nii.gz'%type)

    mask = cortex.get_cortical_mask(subject, xfmname, type=type)
    epi_nii = nb.load(epi)
    dims = epi_nii.shape

    mask_img = nb.Nifti1Image(dataobj=mask.transpose((2,1,0)), affine=epi_nii.affine, header=epi_nii.header)
    
    mask_img.to_filename(out_file)

    return out_file

def create_roi_masks(subject, xfmname, epi, op_folder, types=['cortical', 'thick', 'thin', 'cortical-liberal', 'cortical-conservative'], threshold=None, **kwargs):
    import cortex
    import nibabel as nb
    import os

    # template file
    epi_nii = nb.load(epi)
    dims = epi_nii.shape

    all_masks = {}

    for t in types:
        try:
            os.makedirs(os.path.join(op_folder, t))
        except OSError:
            pass

        print('generating "%s" masks, outputting to %s'%(t, os.path.join(op_folder, t)))

        masks = cortex.get_roi_masks(subject, xfmname, gm_sampler=t, threshold=threshold, **kwargs)
        # subject, xfmname, roi_list=None, gm_sampler='cortical', split_lr=False, allow_overlap=False, fail_for_missing_rois=True, exclude_empty_rois=False, threshold=None, return_dict=True
        for mk in masks.keys():
            mask_img = nb.Nifti1Image(dataobj=masks[mk].transpose((2,1,0)), affine=epi_nii.affine, header=epi_nii.header)
            mask_img.to_filename(os.path.join(op_folder, t, mk+'.nii.gz'))

        all_masks.update({t: masks})

    return all_masks


def transpose_epi(epi):
    import nibabel as nb

    epi_nii = nb.load(epi)
    epi_T_img = nb.Nifti1Image(dataobj=epi_nii.get_data().transpose((2,1,0)), affine=epi_nii.affine, header=epi_nii.header)
    epi_T_fn = epi[:-7] + '_T.nii.gz'
    epi_T_img.to_filename(epi_T_fn)

    return epi_T_fn

def set_pycortex_config_file(project_folder):
    # Import necessary modules
    import os
    import cortex
    from pathlib import Path

    # Define the new database and colormaps folder
    pycortex_db_folder = project_folder + '/db/'
    pycortex_cm_folder = project_folder + '/colormaps/'

    # Get pycortex config file location
    pycortex_config_file  = cortex.options.usercfg

    # Create name of new config file that will be written
    new_pycortex_config_file = pycortex_config_file[:-4] + '_new.cfg'

    # Create the new config file
    Path(new_pycortex_config_file).touch()

    # Open the config file in read mode and the newly created one in write mode.
    # Loop over every line in the original file and copy it into the new one.
    # For the lines containing either 'filestore' or 'colormap', it will
    # change the saved folder path to the newly created one above (e.g. pycortex_db_folder)
    with open(pycortex_config_file, 'r') as fileIn:
        with open(new_pycortex_config_file, 'w') as fileOut:

            for line in fileIn:

                if 'filestore' in line:
                    newline = 'filestore=' + pycortex_db_folder
                    fileOut.write(newline)
                    newline = '\n'

                elif 'Colormaps' in line:
                    newline = 'Colormaps=' + pycortex_cm_folder
                    fileOut.write(newline)
                    newline = '\n'

                else:
                    newline = line

                fileOut.write(newline)

    # Renames the original config file als '_old' and the newly created one to the original name
    os.rename(pycortex_config_file, pycortex_config_file[:-4] + '_old.cfg')
    os.rename(new_pycortex_config_file, pycortex_config_file)
    return None
