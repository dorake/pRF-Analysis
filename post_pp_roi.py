"""
-----------------------------------------------------------------------------------------
post_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create roi-masks and save all data in hdf5 format
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name (e.g. 'sub-001')
-----------------------------------------------------------------------------------------
Output(s):
hdf5 files
-----------------------------------------------------------------------------------------
To run:
cd /home/szinte/projects/pRF_analysis/
python post_pp_roi.py sub-002
-----------------------------------------------------------------------------------------
"""

# General imports
# ---------------
import os
import sys
import json
import glob
import nibabel as nb
import numpy as np
import pdb
opj = os.path.join

# Bokeh imports
# ---------------
from bokeh.io import output_notebook, show, save, output_file, export_png, export_svgs
from bokeh.layouts import row, column, gridplot

# Function import
# ---------------
from pRF_gazeMod.utils.utils import create_roi_masks, roi_data_from_hdf, mask_nii_2_hdf5
from pRF_gazeMod.utils.plot_class import PlotOperator

# Define analysis parameters
# --------------------------
with open('analysis_settings.json') as f:                                                                   # get main analysis settings
    json_s                          =   f.read()
    analysis_info                   =   json.loads(json_s)

sub_id                          =   str(sys.argv[1])                                                        # subject name (e.g. sub-001)
xfmname                         =   'fmriprep'                                                              # pycortex transform name
rois                            =   analysis_info['rois']                                                   # define regions of interest

# Check if it received the 'save' argument to save the figures as svg files
try:saving_figs                 =   True if str(sys.argv[2]).lower() == 'save' else False
except IndexError: saving_figs  =   False


# Define folders
# --------------
base_dir                        =   opj(analysis_info['aeneas_project_directory'],'pp', sub_id)
reg_file                        =   opj(base_dir, 'reg', 'epi.nii.gz')                                       # epi registration files
masks_folder                    =   opj(base_dir, 'masks')                                                   # mask folder
h5_folder                       =   opj(base_dir, 'h5')                                                      # define h5 folder
pRFbasic_folder                 =   opj(base_dir, 'figs','prf_basic')                                        # pRF analysis main folder
pRFbasic_folder_L               =   opj(base_dir, 'figs','prf_basic','LH')                                   # pRF analysis left folder
pRFbasic_folder_R               =   opj(base_dir, 'figs','prf_basic','RH')                                   # pRF analysis right folder
pRFbasic_svg                    =   opj(base_dir, 'figs','svg_files')                                        # pRF analysis svg folder

# Create new folders
# ------------------
for folder in [h5_folder, pRFbasic_folder_L, pRFbasic_folder_R, pRFbasic_svg]:
    try: os.makedirs(folder)                                                                                 # create folder
    except OSError: pass


# Python 3 Warning because of Pycortex
sys.exit('This script works with Pycortex, which requires Python 2.x to work without issues. The current environment seems to have a higher version. Aborting.') if sys.version_info[0] > 2 else None


# # Create ROI masks
# # ----------------
# # create left/right hemisphere roi masks
# create_roi_masks(                               subject         =   sub_id,                                           # subject names
#                                                 xfmname         =   xfmname,                                          # xfmname for pycortex db
#                                                 epi             =   reg_file,                                         # registration file
#                                                 op_folder       =   masks_folder,                                     # output folder
#                                                 types           =   analysis_info['pycortex_roi_mask_types'],         # specify how to sample the cortical gray matter
#                                                 threshold       =   None,                                             # value used to convert probablistic ROI values to a boolean mask
#                                                 split_lr        =   True)                                             # specify whether to separate ROIs in to left and right hemispheres

# Create HDF5 files
# -----------------

# define file name of data to include in hdf5
loo_avg_files                   =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id, 'av_gaze*','loo','*.nii.gz')))                    # get averaged loo files

loo_prf_files                   =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id, 'deriv_gazeAll','all','loo_*.nii.gz')))           # get loo pRF analysis including pos/neg pRF  

loo_prf_files_left              =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id, 'deriv_gazeLeft','all','loo_*.nii.gz')))          # get loo pRF analysis including pos/neg pRF  

loo_prf_files_right             =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
                                                        'pp', sub_id, 'deriv_gazeRight','all','loo_*.nii.gz')))         # get loo pRF analysis including pos/neg pRF  

files                           =   [ (loo_avg_files,       'loo_avg'),      (loo_prf_files,        'loo_prf'), 
                                      (loo_prf_files_left,  'loo_prf_left'), (loo_prf_files_right,  'loo_prf_right')]

for roi_mask in analysis_info['pycortex_roi_mask_types']:
    for roi in analysis_info['rois']:
        mask_files                      =   sorted(glob.glob(opj(masks_folder,roi_mask,\
                                                        '{roi}*.nii.gz'.format(roi = roi))))                            # get roi mask file name
        h5_file                         =   opj(h5_folder,'{roi}.h5'.format(\
                                                        roi = roi))                                                     # define h5 file

        try: os.system('rm '+ h5_file)                                                                                  # delete old file
        except: pass
        
        for file, fName in files:
            _                           =   mask_nii_2_hdf5(                                                            # add masked loo data
                                            in_files        =   file,                                                   # absolute path to functional nifti-files.
                                            mask_files      =   mask_files,                                             # absolute path to mask nifti-files
                                            hdf5_file       =   h5_file,                                                # absolute path to hdf5 file
                                            folder_alias    =   fName)                                                  # name of the to-be-created folder


# Draw main analysis figure
# -------------------------
# Initialize data dictionary that will save all data arrays
data = {}
for roi_text in rois:
    # Get data
    # --------
    for hemisphere in ['_L', '_R']:

        pRF_mean                          =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_mean_all'],                 # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder

        pRF_stim_ratio                    =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_stim_ratio_all'],           # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder
        pRF_ecc                           =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_ecc_all'],                  # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder
        pRF_retinal_x_gain                =   roi_data_from_hdf(                                                                 
                                                        data_types_wildcards =    ['loo_prf_retinal_x_gain_all'],           # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder
        pRF_screen_x_gain                 =   roi_data_from_hdf(                                                                 
                                                        data_types_wildcards =    ['loo_prf_screen_x_gain_all'],            # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder
        pRF_retinal_gain_index            =   roi_data_from_hdf(                                                                 
                                                        data_types_wildcards =    ['loo_prf_retinal_gain_index_all'],       # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder
        
        pRF_amp_change                    =   roi_data_from_hdf(                                                                 
                                                        data_types_wildcards =    ['loo_prf_amplitude_change_all'],         # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder
        pRF_y_change                      =   roi_data_from_hdf(                                                                 
                                                        data_types_wildcards =    ['loo_prf_y_change_all'],                 # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemisphere,                # roi name
                                                        hdf5_file            =    opj(h5_folder,roi_text+'.h5'),        # file location
                                                        folder_alias         =    'loo_prf')                            # alias folder

        var_name                             =   'pRF_fit{0}'.format(hemisphere)

        data.update({var_name                :   np.zeros((pRF_mean.shape[0],14)) })                                         # create empty ndarray                                        
        data[var_name][:,:-7]                =   pRF_mean[:,:7]                                                              # first 6 columns (exclude 2 empty column)
        data[var_name][:, -7]                =   pRF_stim_ratio                                                              # add last columns
        data[var_name][:, -6]                =   pRF_ecc                                                                     # add last columns
        data[var_name][:, -5]                =   pRF_retinal_x_gain                                                          # add last columns
        data[var_name][:, -4]                =   pRF_screen_x_gain                                                           # add last columns
        data[var_name][:, -3]                =   pRF_retinal_gain_index                                                      # add last columns
        data[var_name][:, -2]                =   pRF_amp_change                                                              # add last columns
        data[var_name][:, -1]                =   pRF_y_change                                                                # add last columns

        data_fit_mask                        =   data[var_name][:,6] >= analysis_info['rsq_threshold_roi']                            # define rsq mask for left hemi
        data[var_name]                       =   data[var_name][data_fit_mask]                                                   # mask left hemi data


    # get both hemisphere        
    data.update({'pRF_fit_LR'                :   np.zeros((data['pRF_fit_L'].shape[0] + 
                                                           data['pRF_fit_R'].shape[0], data['pRF_fit_R'].shape[1]))})    # create empty ndarray
    data['pRF_fit_LR'][0:data['pRF_fit_L'].shape[0], :] =  data['pRF_fit_L']                                                                  # fill it with left hemi.
    data['pRF_fit_LR'][  data['pRF_fit_L'].shape[0]:,:] =  data['pRF_fit_R']                                                                  # and fill it with right hemi.



    
    # define main parameters
    # ----------------------
    param_all                       =   dict(
                                            saving_figs             =   saving_figs,
                                            svg_folder              =   pRFbasic_svg,                               # folder to save .svg files
                                            roi_t                   =   roi_text,                                   # current region of interest
                                            p_width                 =   400,                                        # individual plot width
                                            p_height                =   400,                                        # individual plot height
                                            min_border_large        =   10,                                         # large border between figure and axis
                                            min_border_small        =   5,                                          # small border between figure and axis
                                            bg_color                =   tuple([229,229,229]),                       # background color         
                                            stim_color              =   tuple([250,250,250]),                       # stimuli color         
                                            hist_fill_color         =   tuple([255,255,255]),                       # histogram color
                                            hist_line_color         =   tuple([0,0,0]),                             # histogram color
                                            stim_radius             =   analysis_info['stim_radius'],               # stimulus radius
                                            rsq_threshold           =   analysis_info['rsq_threshold'],             # rsquare threshold
                                            cmap                    =   'Spectral',                                 # colormap to use
                                            cmap_steps              =   8,                                          # colormap steps
                                            col_offset              =   0,                                          # colormap offset
                                            vmin                    =   0,                                          # val min to draw with colormap
                                            vmax                    =   4,                                          # val max to draw with colormap
                                            leg_xy_max_ratio        =   1.8                                         # xy axis max based on vmax
                                            )
    # Initialize plotting class
    # -------------------------
    plotter                         =   PlotOperator(**param_all)



    # pRF map
    # -------
    f_pRFmap                        =   {'lh' : [], 'rh' : [], 'lrh' : []} 
    
    params_pRFmap                   =   param_all                                                                   # get main params
    _                               =   params_pRFmap.update(                                                       # update parameters
                                            dict(
                                            x_range                 =   (-8, 8),                                    # x axis limits
                                            y_range                 =   (-8, 8),                                    # y axis limits
                                            x_label                 =   'Horizontal coordinate (dva)',              # x axis label
                                            y_label                 =   'Vertical coordinate (dva)',                # y axis label
                                            x_source_label          =   'x',
                                            y_source_label          =   'y',
                                            x_tick_steps            =   2,                                          # x axis major tick steps
                                            y_tick_steps            =   2,                                          # y axis major tick steps
                                            hist_range              =   (0,0.5),                                    # histogram x/y range
                                            hist_steps              =   0.5,                                        # x axis major tick steps
                                            h_hist_bins             =   16,                                         # hor. histogram bins
                                            v_hist_bins             =   16))                                         # ver. histogram bins
    
    for hemi, list_name, ext in [('-left', 'lh', '_L'), ('-right', 'rh', '_R'), ('', 'lrh', '_LR')]:
        title                           =   '{roi}{hemi} hemisphere: pRF position and size (n={voxel})'.\
                                             format(roi = roi_text, hemi = hemi, voxel = data['pRF_fit%s'%ext].shape[0])            # define titlte
        
        _                               =   params_pRFmap.update(                                                       # update parameters
                                                dict(              
                                                dataMat                 =   data['pRF_fit%s' % ext],                                  # main data matrix
                                                main_fig_title          =   title,
                                                hemisphere              =   list_name,
                                                typeData                =   None))                                     # title
        f_pRFmap[list_name], _          =   plotter.draw_figure(parameters = params_pRFmap, plot = 'map')               # draw pRF map figure

    
    # pRF ecc
    # -------
    f_pRFecc_basic                      =   {'lh' : [], 'rh' : [], 'lrh' : []}
    first_pRFecc_fig_basic              =   []

    for numData, typeData in enumerate(['Size','CV R2','R2','Amplitude','Stim-ratio']):

        params_pRFecc                   =   param_all                                                               # get main params     
        _                               =   params_pRFecc.update(                                                   # update parameters
                                                dict(
                                                x_range                 =   (0, 10),                                # x axis limits
                                                x_label                 =   'Eccentricity (dva)',                   # x axis label
                                                x_tick_steps            =   2,                                      # x axis major tick steps
                                                x_source_label          =   'ecc',
                                                vmin                    =   0,
                                                vmax                    =   4,
                                                hist_range              =   (0,0.5),                                # histogram x/y range
                                                hist_steps              =   0.5,                                    # x axis major tick steps
                                                h_hist_bins             =   20),                                     # hor. histogram bins
                                                condition               =   'ecc'
                                                )

        if typeData == 'Size':
            _                           =   params_pRFecc.update(                                                   # update parameters
                                                dict(
                                                y_range                 =   (0, 4),                                 # y axis limits (adapted for size)
                                                y_label                 =   'Size (dva)',                           # y axis label
                                                y_source_label          =   'sigma',                                # y axis source val in data source
                                                y_tick_steps            =   1,                                      # y axis major tick steps
                                                v_hist_bins             =   16)                                      # ver. histogram bins
                                                )
        elif typeData == 'CV R2':
            _                           =   params_pRFecc.update(                                                   # update parameters
                                                dict(
                                                y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
                                                y_label                 =   'Cross-validated R2',                   # y axis label
                                                y_source_label          =   'cv_rsq',                               # y axis source val in data source
                                                y_tick_steps            =   0.2,                                    # y axis major tick steps
                                                v_hist_bins             =   10)                                     # ver. histogram bins
                                                )
        elif typeData == 'R2':
            _                           =   params_pRFecc.update(                                                   # update parameters
                                                dict(
                                                y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
                                                y_label                 =   'R2',                                   # y axis label
                                                y_source_label          =   'rsq',                                  # y axis source val in data source
                                                y_tick_steps            =   0.2,                                    # y axis major tick steps
                                                v_hist_bins             =   10)                                     # ver. histogram bins
                                                )
        elif typeData == 'Amplitude':
            _                           =   params_pRFecc.update(                                                   # update parameters
                                                dict(
                                                y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
                                                y_label                 =   'Amplitude (z-score)',                  # y axis label
                                                y_source_label          =   'beta',                                 # y axis source val in data source
                                                y_tick_steps            =   0.25,                                   # y axis major tick steps
                                                v_hist_bins             =   10)                                     # ver. histogram bins   
                                                )

        elif typeData == 'Stim-ratio':
            _                           =   params_pRFecc.update(                                                   # update parameters
                                                dict(
                                                y_range                 =   (-10, 110),                             # y axis limits (adapted for size)
                                                y_label                 =   'Stimulus ratio (%)',                   # y axis label
                                                y_source_label          =   'stim_ratio',                           # y axis source val in data source
                                                y_tick_steps            =   10,                                     # y axis major tick steps
                                                v_hist_bins             =   12)                                     # ver. histogram bins
                                                )
        

        for hemi, list_name, ext in [('-left', 'lh', '_L'), ('-right', 'rh', '_R'), ('', 'lrh', '_LR')]:
            title                           =   '{roi}{hemi} hemisphere: Eccentricity vs. {typeData}'.\
                                                 format(roi = roi_text, hemi = hemi, typeData = typeData)           # define title
            _                               =   params_pRFecc.update(                                                   # update parameters
                                                    dict(
                                                    dataMat                 =   data['pRF_fit%s'%ext],                              # main data matrix
                                                    ecc_gainRatio_counter   =   0,
                                                    main_fig_title          =   title,
                                                    hemisphere              =   list_name,
                                                    typeData                =   typeData))                                  # title
            
            old_fig = [] if list_name == 'lh' else first_pRFecc_fig_basic[numData]

            out1, out2                      =   plotter.draw_figure(parameters = params_pRFecc,
                                                                    # old_main_fig = old_fig,                                                                    
                                                                    plot = 'ecc')
            _                               =   f_pRFecc_basic[list_name].append(out1)
            _                               =   first_pRFecc_fig_basic.append(out2)


    # pRF cov
    # -------
    f_pRFcov                        =   {'lh' : [], 'rh' : [], 'lrh' : []}
    params_pRFcov                   =   param_all                                                                   # get main params
    _                               =   params_pRFcov.update(                                                       # update parameters
                                            dict(
                                            x_range                 =   (-8, 8),                                    # x axis limits
                                            y_range                 =   (-8, 8),                                    # y axis limits
                                            x_label                 =   'Horizontal coordinate (dva)',              # x axis label
                                            y_label                 =   'Vertical coordinate (dva)',                # y axis label
                                            x_tick_steps            =   2,                                          # x axis major tick steps
                                            y_tick_steps            =   2,                                          # y axis major tick steps
                                            dataMat                 =   data['pRF_fit_LR'],
                                            smooth_factor           =   15,                                         # pixel per degree in image
                                            cmap                    =   'viridis',                                  # colormap to use
                                            cmap_steps              =   10,                                         # colormap steps
                                            col_offset              =   0,                                          # colormap offset
                                            vmin                    =   0,                                          # val min to draw with colormap
                                            vmax                    =   1,                                          # val max to draw with colormap
                                            colorbar_tick_steps     =   0.2,                                        # colorbar axis major tick steps
                                            condition               =   'cov',
                                            colorbar_label          =   'pRF coverage (norm.)')                     # colorbar label
                                            )
    
    for hemi, list_name, ext in [('-left', 'lh', '_L'), ('-right', 'rh', '_R'), ('', 'lrh', '_LR')]:
        title                           =   '{roi}{hemi} hemisphere: pRF coverage (n={voxel})'.\
                                             format(roi = roi_text, hemi = hemi, voxel = data['pRF_fit_L'].shape[0])            # define titlte
        _                               =   params_pRFcov.update(dict(              
                                                dataMat                 =   data['pRF_fit%s'%ext],                                  # main data matrix
                                                main_fig_title          =   title,
                                                hemisphere              =   list_name,
                                                typeData                =   None))
        f_pRFcov[list_name],_           =   plotter.draw_figure(parameters = params_pRFcov, plot = 'cov')



    # Combining Figures
    # -----------------
    for list_name in ['lh', 'rh', 'lrh']:
        all_fL                          =     gridplot([
                                                [f_pRFmap[list_name],                f_pRFecc_basic[list_name][0],     f_pRFcov[list_name]],                # specify figure 1st row
                                                [f_pRFecc_basic[list_name][1],       f_pRFecc_basic[list_name][2],     None],                               # specify figure 2nd row
                                                [f_pRFecc_basic[list_name][3],       f_pRFecc_basic[list_name][4],     None]])
        folder_name                     =   opj(pRFbasic_folder, list_name.upper()) if list_name != 'lrh' else pRFbasic_folder
        output_file_html                =   opj(folder_name,'%s_pRF.html'%roi_text)                                                          # define html file
        _                               =   output_file(output_file_html, title='%s pRF analysis'%roi_text)                                   # get html saved
        _                               =   save(all_fL)
    

  
