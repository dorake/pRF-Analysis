"""
-----------------------------------------------------------------------------------------
main_analysis.py
-----------------------------------------------------------------------------------------
Goal of the script:
Per voxel main analysis launcher
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name (e.g. 'sub-001')
-----------------------------------------------------------------------------------------
Output(s):
hdf5 files
-----------------------------------------------------------------------------------------
To run:
cd /home/szinte/projects/pRF_analysis/
python main_analysis.py sub-002
-----------------------------------------------------------------------------------------
To do:
- eror in pRFcor with second plots not being Y axis
- too many plot in pRFcor and in general
- put colorbar in pRFshift?
- fix colorbar to have white central
-----------------------------------------------------------------------------------------
"""

# General imports
# ---------------
import os, sys, json, glob
import nibabel as nb
import numpy as np
opj = os.path.join

# Bokeh imports
# ---------------
from bokeh.io import output_notebook, show, save, output_file, export_png, export_svgs
from bokeh.layouts import row, column, gridplot

# Function import
# ---------------
from pRF_gazeMod.utils.utils import roi_data_from_hdf, mask_nii_2_hdf5
from pRF_gazeMod.utils.prf import weighted_avg_std
from pRF_gazeMod.utils.plot_class import PlotOperator

# Define analysis parameters
# --------------------------
with open('analysis_settings.json') as f:                                                                           # get main analysis settings
    json_s                          =   f.read()
    analysis_info                   =   json.loads(json_s)

sub_id                          =   str(sys.argv[1])                                                                # subject name (e.g. sub-001)
rois                            =   analysis_info['rois']                                                           # define regions of interest

# Check if it received the 'save' argument to save the figures as svg files
try:saving_figs                 =   True if str(sys.argv[2]).lower() == 'save' else False
except IndexError: saving_figs  =   False


# Python 3 Warning because of Pycortex
sys.exit('This script works with Pycortex, which requires Python 2.x to work without issues. The current environment seems to have a higher version. Aborting.') if sys.version_info[0] > 2 else None

# Define folders
# --------------
base_folder =   opj(analysis_info['aeneas_project_directory'], 'pp', sub_id)
folders     =   {}
for folder in ['figs', 'masks', 'h5']:
    folders.update({folder : opj(base_folder, folder)})


for folder in ['ecc_x_gain','ecc_amp', 'fit_cor', 'prf_shift', 'prf_spatial_shift', 'prf_gain', 'prf_roi_basic', 'prf_roi_main', 'svg_files']:
    folders.update({folder : opj(folders['figs'], folder)})
    
    # Save svg seperately because we dont need LH and RH
    if folder == 'svg_files':
        try: os.makedirs(folders[folder])
        except OSError: pass

    # Create new folders
    # ------------------
    for fileExt in ['LH', 'RH']:
        fileName = opj(folders[folder], fileExt)
        try: os.makedirs(fileName)
        except OSError: pass


# Main analysis
# ------------- 
# Create all necessary arrays
data = {}                                                                          
for hemisphere in ['L', 'R', 'LR']:
    for gaze, columnLength in zip(['gazeAll', 'gazeLeft', 'gazeRight', 'gazeGain'], [14, 14, 14, 1]):
        for statistic in ['avg', 'std', 'num']:
            colLen = columnLength if statistic != 'num' else 1 
            variable_name = 'pRF_{gaze}_{statistic}_{hemisphere}'.format(gaze = gaze, statistic = statistic, hemisphere = hemisphere)
            data.update({variable_name : np.zeros([len(rois), colLen])})


# Get the data of both hemispheres (left, right), both conditions (gazeAll, gazeLeft, gazeRight), save them and calculate averages
for numroi, roi_text in enumerate(rois):
    for hemi in ['_L', '_R']:
        pRF_mean              =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_mean_all'],         # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                        hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                        folder_alias         =    'loo_prf' + '_left')              # alias folder

        pRF_stim_ratio        =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_stim_ratio_all'],   # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                        hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                        folder_alias         =    'loo_prf' + '_left')              # alias folder
        pRF_ecc               =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_ecc_all'],          # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                        hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                        folder_alias         =    'loo_prf' + '_left')              # alias folders
        var_name = 'pRF_gazeLeft_fit{hemi}'.format(hemi=hemi)
        data_gazeLeft                          =   np.zeros((pRF_mean.shape[0],9))
        data_gazeLeft[:,:-2]         =   pRF_mean[:,:7]                                                        # first 8 columns
        data_gazeLeft[:, -2]         =   pRF_stim_ratio                                                        # add last columns
        data_gazeLeft[:, -1]         =   pRF_ecc                                                               # add last columns
        pRF_Left_mask         =   (data_gazeLeft[:, 6] >= analysis_info['rsq_threshold']) & (data_gazeLeft[:, 7] >= analysis_info['in_stim_threshold'])        # define rsq mask for left hemi



        pRF_mean              =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_mean_all'],         # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                        hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                        folder_alias         =    'loo_prf' + '_right')              # alias folder

        pRF_stim_ratio        =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_stim_ratio_all'],   # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                        hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                        folder_alias         =    'loo_prf' + '_right')              # alias folder
        pRF_ecc               =   roi_data_from_hdf(
                                                        data_types_wildcards =    ['loo_prf_ecc_all'],          # data contained in the file
                                                        roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                        hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                        folder_alias         =    'loo_prf' + '_right')              # alias folders
        var_name = 'pRF_gazeRight_fit{hemi}'.format(hemi=hemi)
        data_gazeRight                          =   np.zeros((pRF_mean.shape[0],9))
        data_gazeRight[:,:-2]         =   pRF_mean[:,:7]                                                        # first 8 columns
        data_gazeRight[:, -2]         =   pRF_stim_ratio                                                        # add last columns
        data_gazeRight[:, -1]         =   pRF_ecc                                                               # add last columns

        pRF_Right_mask       =   (data_gazeRight[:, 6] >= analysis_info['rsq_threshold']) & (data_gazeRight[:, 7] >= analysis_info['in_stim_threshold'] )         # define rsq mask for left hemi



        for idx, gaze, ext in ([(0, 'gazeAll', ''), (1, 'gazeLeft', '_left'), (2, 'gazeRight', '_right')]):
            pRF_mean              =   roi_data_from_hdf(
                                                            data_types_wildcards =    ['loo_prf_mean_all'],         # data contained in the file
                                                            roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                            folder_alias         =    'loo_prf' + ext)              # alias folder

            pRF_stim_ratio        =   roi_data_from_hdf(
                                                            data_types_wildcards =    ['loo_prf_stim_ratio_all'],   # data contained in the file
                                                            roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                            folder_alias         =    'loo_prf' + ext)              # alias folder

            pRF_ecc               =   roi_data_from_hdf(
                                                            data_types_wildcards =    ['loo_prf_ecc_all'],          # data contained in the file
                                                            roi_name_wildcard    =    roi_text + hemi,              # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),# file location
                                                            folder_alias         =    'loo_prf' + ext)              # alias folders


            var_name = 'pRF_{gaze}_fit{hemi}'.format(gaze=gaze, hemi=hemi)
            data.update({var_name : np.zeros((pRF_mean.shape[0],14))})
            data[var_name][:,:-7]         =   pRF_mean[:,:7]                                                        # first 8 columns
            data[var_name][:, -7]         =   pRF_stim_ratio                                                        # add last columns
            data[var_name][:, -6]         =   pRF_ecc                                                               # add last columns
            
            # import ipdb ; ipdb.set_trace()
            if 'pRF_gazeAll' in var_name:
                pRF_retinal_x_gain     =   roi_data_from_hdf(                                                                 
                                                            data_types_wildcards =    ['loo_prf_retinal_x_gain_all'],            # data contained in the file
                                                            roi_name_wildcard    =    roi_text+ hemi,                        # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),     # file location
                                                            folder_alias         =    'loo_prf')                             # alias folder
                pRF_screen_x_gain      =   roi_data_from_hdf(                                                                 
                                                            data_types_wildcards =    ['loo_prf_screen_x_gain_all'],             # data contained in the file
                                                            roi_name_wildcard    =    roi_text+ hemi,                        # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),     # file location
                                                            folder_alias         =    'loo_prf')                             # alias folder
                pRF_retinal_gain_index =   roi_data_from_hdf(                                                                 
                                                            data_types_wildcards =    ['loo_prf_retinal_gain_index_all'],        # data contained in the file
                                                            roi_name_wildcard    =    roi_text+ hemi,                        # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),     # file location
                                                            folder_alias         =    'loo_prf')                             # alias folder
                                
                pRF_amp_change         =   roi_data_from_hdf(                                                                 
                                                            data_types_wildcards =    ['loo_prf_amplitude_change_all'],          # data contained in the file
                                                            roi_name_wildcard    =    roi_text+ hemi,                        # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),     # file location
                                                            folder_alias         =    'loo_prf')                             # alias folder

                pRF_y_change           =   roi_data_from_hdf(                                                                 
                                                            data_types_wildcards =    ['loo_prf_y_change_all'],                  # data contained in the file
                                                            roi_name_wildcard    =    roi_text+ hemi,                        # roi name
                                                            hdf5_file            =    opj(folders['h5'],roi_text+'.h5'),     # file location
                                                            folder_alias         =    'loo_prf')                             # alias folder

                # import ipdb ; ipdb.set_trace()
                data[var_name][:, -5]         =   pRF_retinal_x_gain                                                         # add last columns
                data[var_name][:, -4]         =   pRF_screen_x_gain                                                          # add last columns
                data[var_name][:, -3]         =   pRF_retinal_gain_index                                                     # add last columns
                data[var_name][:, -2]         =   pRF_amp_change
                data[var_name][:, -1]         =   pRF_y_change                                                               # add last columns
            

            data[var_name]                =   data[var_name][pRF_Left_mask & pRF_Right_mask]                                 # mask left hemi data

            data['pRF_{0}_avg{1}'.format(gaze,hemi)][numroi,:],\
            data['pRF_{0}_std{1}'.format(gaze,hemi)][numroi,:],\
            data['pRF_{0}_num{1}'.format(gaze,hemi)][numroi,:]                  =    weighted_avg_std(
                                                            matrix_input        =    data[var_name],                 # input matrix to average
                                                            matrix_weight       =    data['pRF_gazeAll_fit%s' % hemi][:,7])               # input matrix to average
    # # get both hemisphere
    # # -------------------
    for gaze in ['gazeAll', 'gazeLeft', 'gazeRight']:
        # gazeAll
        data.update({'pRF_%s_fit_LR' % gaze : np.zeros((data['pRF_%s_fit_L' % gaze].shape[0] +
                                                      data['pRF_%s_fit_R' % gaze].shape[0], 
                                                      data['pRF_%s_fit_R' % gaze].shape[1]))})                       # create empty ndarray
        data['pRF_%s_fit_LR' % gaze][0: data['pRF_%s_fit_L' % gaze].shape[0], :] = data['pRF_%s_fit_L' % gaze]       # fill it with left hemi.
        data['pRF_%s_fit_LR' % gaze][   data['pRF_%s_fit_L' % gaze].shape[0]:,:] = data['pRF_%s_fit_R' % gaze]       # and fill it with right hemi.
        
        data['pRF_%s_avg_LR' % gaze][numroi,:], data['pRF_%s_std_LR' % gaze][numroi,:],\
        data['pRF_%s_num_LR' % gaze][numroi,:]                                  =    weighted_avg_std(
                                                           matrix_input         =    data['pRF_%s_fit_LR' % gaze],      # input matrix to average
                                                           matrix_weight        =    data['pRF_gazeAll_fit_LR'][:,6])   # input matrix to average


    # ---------------------------------
    # Across condition correlation plot
    # ---------------------------------
    # define main parameters
    param_all                           =   dict(
                                            saving_figs             =   saving_figs,
                                            svg_folder              =   folders['svg_files'],                       # folder to save .svg files
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
                                            leg_xy_max_ratio        =   1.8,                                         # xy axis max based on vmax
                                            hist_range              =   (0,0.5),                                    # histogram x/y range
                                            hist_steps              =   0.5,
                                            across_subjects         =   False
                                            )
    # Initialize plotting class
    # -------------------------
    plotter                             =   PlotOperator(**param_all)


    # # pRF ecc with x gain as color
    # # -----------------------------
    # f_pRFecc_x_gain                     = {'lh': [], 'rh': [], 'lrh': [], 'first_fig': []}
    # first_pRFecc_fig_x_gain             = [] 

    # for numData, typeData in enumerate(['Size','CV R2','R2','Amplitude','Stim-ratio']):

    #     params_pRFecc                   =   param_all                                                               # get main params     
    #     _                               =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             x_range                 =   (0, 10),                                # x axis limits
    #                                             x_label                 =   'Eccentricity (dva)',                   # x axis label
    #                                             x_tick_steps            =   2,                                      # x axis major tick steps
    #                                             x_source_label          =   'ecc',
    #                                             vmin                    =   -0.5,
    #                                             vmax                    =   0.5,
    #                                             cmap_steps              =   255,
    #                                             hist_range              =   (0,0.5),                                # histogram x/y range
    #                                             hist_steps              =   0.5,                                    # x axis major tick steps
    #                                             h_hist_bins             =   20),                                     # hor. histogram bins
    #                                             condition               =   'ecc_gainRatio',
    #                                             gain                    =   'x',
    #                                             cmap                    =   'RdBu_r',                                 # colormap to use
    #                                             )

    #     if typeData == 'Size':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 4),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'Size (dva)',                           # y axis label
    #                                             y_source_label          =   'sigma',                                # y axis source val in data source
    #                                             y_tick_steps            =   1,                                      # y axis major tick steps
    #                                             v_hist_bins             =   16)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'CV R2':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'Cross-validated R2',                   # y axis label
    #                                             y_source_label          =   'cv_rsq',                               # y axis source val in data source
    #                                             y_tick_steps            =   0.2,                                    # y axis major tick steps
    #                                             v_hist_bins             =   20)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'R2':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'R2',                                   # y axis label
    #                                             y_source_label          =   'rsq',                                  # y axis source val in data source
    #                                             y_tick_steps            =   0.2,                                    # y axis major tick steps
    #                                             v_hist_bins             =   20)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'Amplitude':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (-0.75, 0.75),                          # y axis limits (adapted for size)
    #                                             y_label                 =   'Amplitude (z-score)',                  # y axis label
    #                                             y_source_label          =   'beta',                                 # y axis source val in data source
    #                                             y_tick_steps            =   0.25,                                   # y axis major tick steps
    #                                             v_hist_bins             =   24)                                     # ver. histogram bins   
    #                                             )

    #     elif typeData == 'Stim-ratio':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (-10, 110),                             # y axis limits (adapted for size)
    #                                             y_label                 =   'Stimulus ratio (%)',                   # y axis label
    #                                             y_source_label          =   'stim_ratio',                           # y axis source val in data source
    #                                             y_tick_steps            =   10,                                     # y axis major tick steps
    #                                             v_hist_bins             =   24))                                    # ver. histogram bins
        
    #     for hemi, list_name, ext in [('-left', 'lh', '_L'), ('-right', 'rh', '_R'), ('', 'lrh', '_LR')]:
        
    #         title                           =   '{roi}{hemi} hemisphere: Eccentricity vs. {typeData}'.\
    #                                              format(roi = roi_text, hemi = hemi, typeData = typeData)           # define title
    #         _                               =   params_pRFecc.update(                                                   # update parameters
    #                                                 dict(
    #                                                 dataMat                 =   data['pRF_gazeAll_fit%s'%ext],                              # main data matrix
    #                                                 ecc_gainRatio_counter   =   0,
    #                                                 main_fig_title          =   title,
    #                                                 hemisphere              =   list_name,
    #                                                 typeData                =   typeData))                                     # title
       
    #         out1, out2                      =   plotter.draw_figure(parameters = params_pRFecc, 
    #                                                                 plot = 'ecc')
    #         _                               =   f_pRFecc_x_gain[list_name].append(out1)
    #         _                               =   first_pRFecc_fig_x_gain.append(out2)



    # # pRF ecc with amplitude gain as color
    # # -----------------------------
    # f_pRFecc_amp                            = {'lh': [], 'rh': [], 'lrh': [], 'first_fig': []}
    # first_pRFecc_fig_amp                    = [] 
    # numData_amp                             = 0

    # for numData, typeData in enumerate(['Size','CV R2','R2','Amplitude','Stim-ratio']):

    #     params_pRFecc                   =   param_all                                                               # get main params     
    #     _                               =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             x_range                 =   (0, 10),                                # x axis limits
    #                                             x_label                 =   'Eccentricity (dva)',                   # x axis label
    #                                             x_tick_steps            =   2,                                      # x axis major tick steps
    #                                             x_source_label          =   'ecc',
    #                                             vmin                    =   -0.5,
    #                                             vmax                    =   0.5,
    #                                             cmap_steps              =   255,
    #                                             hist_range              =   (0,0.5),                                # histogram x/y range
    #                                             hist_steps              =   0.5,                                    # x axis major tick steps
    #                                             h_hist_bins             =   20),                                     # hor. histogram bins
    #                                             condition               =   'ecc_gainRatio',
    #                                             gain                    =   'amplitude',
    #                                             cmap                    =   'RdBu_r'                                # colormap to use
    #                                             )

    #     if typeData == 'Size':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 4),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'Size (dva)',                           # y axis label
    #                                             y_source_label          =   'sigma',                                # y axis source val in data source
    #                                             y_tick_steps            =   1,                                      # y axis major tick steps
    #                                             v_hist_bins             =   16)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'CV R2':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'Cross-validated R2',                   # y axis label
    #                                             y_source_label          =   'cv_rsq',                               # y axis source val in data source
    #                                             y_tick_steps            =   0.2,                                    # y axis major tick steps
    #                                             v_hist_bins             =   20)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'R2':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits (adapted for size)
    #                                             y_label                 =   'R2',                                   # y axis label
    #                                             y_source_label          =   'rsq',                                  # y axis source val in data source
    #                                             y_tick_steps            =   0.2,                                    # y axis major tick steps
    #                                             v_hist_bins             =   20)                                     # ver. histogram bins
    #                                             )            
    #     elif typeData == 'Amplitude':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (-0.75, 0.75),                          # y axis limits (adapted for size)
    #                                             y_label                 =   'Amplitude (z-score)',                  # y axis label
    #                                             y_source_label          =   'beta',                                 # y axis source val in data source
    #                                             y_tick_steps            =   0.25,                                   # y axis major tick steps
    #                                             v_hist_bins             =   24)                                     # ver. histogram bins   
    #                                             )

    #     elif typeData == 'Stim-ratio':
    #         _                           =   params_pRFecc.update(                                                   # update parameters
    #                                             dict(
    #                                             y_range                 =   (-10, 110),                             # y axis limits (adapted for size)
    #                                             y_label                 =   'Stimulus ratio (%)',                   # y axis label
    #                                             y_source_label          =   'stim_ratio',                           # y axis source val in data source
    #                                             y_tick_steps            =   10,                                     # y axis major tick steps
    #                                             v_hist_bins             =   12)                                     # ver. histogram bins
    #                                             )
        
    #     for hemi, list_name, ext in [('-left', 'lh', '_L'), ('-right', 'rh', '_R'), ('', 'lrh', '_LR')]:
        
    #         title                           =   '{roi}{hemi} hemisphere: Eccentricity vs. {typeData}'.\
    #                                              format(roi = roi_text, hemi = hemi, typeData = typeData)           # define title
    #         _                               =   params_pRFecc.update(                                                   # update parameters
    #                                                 dict(
    #                                                 dataMat                 =   data['pRF_gazeAll_fit%s'%ext],                              # main data matrix
    #                                                 ecc_gainRatio_counter   =   0,
    #                                                 main_fig_title          =   title,                                  # title
    #                                                 hemisphere              =   list_name,
    #                                                 typeData                =   typeData))                                     # title
       
    #         out1, out2                      =   plotter.draw_figure(parameters = params_pRFecc, 
    #                                                                 plot = 'ecc')
    #         _                               =   f_pRFecc_amp[list_name].append(out1)
    #         _                               =   first_pRFecc_fig_amp.append(out2)


    # # Combine figures
    # # ---------------
    # for hemi, letter in [('lh', 'L'), ('rh', 'R'), ('lrh', '')]:

    #     all_f                           =   gridplot([[f_pRFecc_x_gain[hemi][0], f_pRFecc_x_gain[hemi][1], f_pRFecc_x_gain[hemi][2]],
    #                                                   [f_pRFecc_x_gain[hemi][3], f_pRFecc_x_gain[hemi][4],  None           ]])                                    
    #     figs_folder                     =   opj(folders['ecc_x_gain'], hemi.upper()) if hemi != 'lrh' else folders['ecc_x_gain']
    #     output_file_html                =   opj(figs_folder,'{0}{1}_prfEcc_x_gain.html'.format(roi_text, letter))                             # define html file
    #     _                               =   output_file(output_file_html, title='{0} {1} left/right ecc.'.format(roi_text, letter))      # get html saved
    #     _                               =   save(all_f)

    
    # # Combine figures
    # # ---------------
    # for hemi, letter in [('lh', 'L'), ('rh', 'R'), ('lrh', '')]:

    #     all_f                           =   gridplot([[f_pRFecc_amp[hemi][0], f_pRFecc_amp[hemi][1], f_pRFecc_amp[hemi][2]],
    #                                                   [f_pRFecc_amp[hemi][3], f_pRFecc_amp[hemi][4],  None               ]])                                    
    #     figs_folder                     =   opj(folders['ecc_amp'], hemi.upper()) if hemi != 'lrh' else folders['ecc_amp']
    #     output_file_html                =   opj(figs_folder,'{0}{1}_prfEccGain_amplitude.html'.format(roi_text, letter))                             # define html file
    #     _                               =   output_file(output_file_html, title='{0} {1} left/right ecc.'.format(roi_text, letter))      # get html saved
    #     _                               =   save(all_f)



    # # Fit Cor
    # # -------
    # f_pRF_cor = {'lh' : [], 'rh' : [], 'lrh' : []}
    
    # for typeData in ['Hor. coord.','Ver. coord.','Eccentricity','Size','CV R2', 'R2', 'Amplitude','Stim-ratio']:

    #     params_pRFcor                   =   param_all                                                               # get main params

    #     if (typeData == 'Hor. coord.') or (typeData == 'Ver. coord.'):
    #         ax = 'x' if typeData == 'Hor. coord.' else 'y'
    #         _                           =   params_pRFcor.update(dict(
    #                                         x_range                 =   (-8, 8),                                    # x axis limits
    #                                         y_range                 =   (-8, 8),                                    # y axis limits
    #                                         x_label                 =   'Gaze Left - %s (dva)' % typeData,            # x axis label
    #                                         y_label                 =   'Gaze Right - %s (dva)' % typeData,           # y axis label
    #                                         x_source_label          =   '%s_left'%ax,                                   # x axis source val in data source
    #                                         y_source_label          =   '%s_right'%ax,                                  # y axis source val in data source
    #                                         xy_source_label         =   '%s_all'%ax,                                    # x/y axis source val in data source
    #                                         x_tick_steps            =   2,                                          # x axis major tick steps
    #                                         h_hist_bins             =   20,                                         # hor. histogram bins
    #                                         y_tick_steps            =   2,                                          # y axis major tick steps
    #                                         v_hist_bins             =   20                                         # ver. histogram bins
    #                                         ))
    #     elif typeData == 'Eccentricity':
    #         _                           =   params_pRFcor.update(                                                   # update parameters
    #                                         dict(
    #                                         x_range                 =   (0, 10),                                    # x axis limits
    #                                         y_range                 =   (0, 10),                                    # y axis limits
    #                                         x_label                 =   'Gaze Left - Eccentricity (dva)',           # x axis label
    #                                         y_label                 =   'Gaze Right - Eccentricity (dva)',          # y axis label
    #                                         x_source_label          =   'y_left',                                   # x axis source val in data source
    #                                         y_source_label          =   'y_right',                                  # y axis source val in data source
    #                                         xy_source_label         =   'y_all',                                    # x/y axis source val in data source
    #                                         x_tick_steps            =   2,                                          # x axis major tick steps
    #                                         h_hist_bins             =   20,                                         # hor. histogram bins
    #                                         y_tick_steps            =   2,                                          # y axis major tick steps
    #                                         v_hist_bins             =   20,                                         # ver. histogram bins
    #                                         ))

    #     elif typeData == 'Size':
    #         _                           =   params_pRFcor.update(                                                   # update parameters
    #                                         dict(
    #                                         x_range                 =   (0, 4),                                     # x axis limits
    #                                         y_range                 =   (0, 4),                                     # y axis limits
    #                                         x_label                 =   'Gaze Left - Size (dva)',                   # x axis label
    #                                         y_label                 =   'Gaze Right - Size (dva)',                  # y axis label
    #                                         x_source_label          =   'sigma_left',                               # x axis source val in data source
    #                                         y_source_label          =   'sigma_right',                              # y axis source val in data source
    #                                         xy_source_label         =   'sigma_all',                                # x/y axis source val in data source
    #                                         x_tick_steps            =   1,                                          # x axis major tick steps
    #                                         h_hist_bins             =   16,                                         # hor. histogram bins
    #                                         y_tick_steps            =   1,                                          # y axis major tick steps
    #                                         v_hist_bins             =   16,                                         # ver. histogram bins
    #                                         ))
    #     elif typeData == 'CV R2':
    #         _                           =   params_pRFcor.update(                                                   # update parameters
    #                                         dict(
    #                                         x_range                 =   (0, 1),                                     # x axis limits
    #                                         y_range                 =   (0, 1),                                     # y axis limits
    #                                         x_label                 =   'Gaze Left - CV R2',                        # x axis label
    #                                         y_label                 =   'Gaze Right - CV R2',                       # y axis label
    #                                         x_source_label          =   'cv_rsq_left',                              # x axis source val in data source
    #                                         y_source_label          =   'cv_rsq_right',                             # y axis source val in data source
    #                                         xy_source_label         =   'cv_rsq_all',                               # x/y axis source val in data source
    #                                         x_tick_steps            =   0.2,                                        # x axis major tick steps
    #                                         h_hist_bins             =   20,                                         # hor. histogram bins
    #                                         y_tick_steps            =   0.2,                                        # y axis major tick steps
    #                                         v_hist_bins             =   20,                                         # ver. histogram bins
    #                                         ))
    #     elif typeData == 'R2':
    #         _                           =   params_pRFcor.update(                                                   # update parameters
    #                                         dict(
    #                                         x_range                 =   (0, 1),                                     # x axis limits
    #                                         y_range                 =   (0, 1),                                     # y axis limits
    #                                         x_label                 =   'Gaze Left - R2',                           # x axis label
    #                                         y_label                 =   'Gaze Right - R2',                          # y axis label
    #                                         x_source_label          =   'rsq_left',                                 # x axis source val in data source
    #                                         y_source_label          =   'rsq_right',                                # y axis source val in data source
    #                                         xy_source_label         =   'rsq_all',                                  # x/y axis source val in data source
    #                                         x_tick_steps            =   0.2,                                        # x axis major tick steps
    #                                         h_hist_bins             =   20,                                         # hor. histogram bins
    #                                         y_tick_steps            =   0.2,                                        # y axis major tick steps
    #                                         v_hist_bins             =   20,                                         # ver. histogram bins
    #                                         ))

    #     elif typeData == 'Amplitude':
    #         _                           =   params_pRFcor.update(                                                   # update parameters
    #                                         dict(
    #                                         x_range                 =   (-0.75, 0.75),                              # x axis limits
    #                                         y_range                 =   (-0.75, 0.75),                              # y axis limits
    #                                         x_label                 =   'Gaze Left - Amplitude (z)',                # x axis label
    #                                         y_label                 =   'Gaze Right - Amplitude (z)',               # y axis label
    #                                         x_source_label          =   'beta_left',                                # x axis source val in data source
    #                                         y_source_label          =   'beta_right',                               # y axis source val in data source
    #                                         xy_source_label         =   'beta_all',                                 # x/y axis source val in data source
    #                                         x_tick_steps            =   0.25,                                       # x axis major tick steps
    #                                         h_hist_bins             =   24,                                         # hor. histogram bins
    #                                         y_tick_steps            =   0.25,                                       # y axis major tick steps
    #                                         v_hist_bins             =   24,                                         # ver. histogram bins
    #                                         ))

    #     elif typeData == 'Stim-ratio':
    #         _                           =   params_pRFcor.update(                                                   # update parameters
    #                                         dict(
    #                                         x_range                 =   (-10, 110),                                 # x axis limits
    #                                         y_range                 =   (-10, 110),                                 # y axis limits
    #                                         x_label                 =   'Gaze Left - Stim ratio (%)',               # x axis label
    #                                         y_label                 =   'Gaze Right - Stim ratio (%)',              # y axis label
    #                                         x_source_label          =   'stim_ratio_left',                          # x axis source val in data source
    #                                         y_source_label          =   'stim_ratio_right',                         # y axis source val in data source
    #                                         xy_source_label         =   'stim_ratio_all',                           # x/y axis source val in data source
    #                                         x_tick_steps            =   10,                                         # x axis major tick steps
    #                                         h_hist_bins             =   12,                                         # hor. histogram bins
    #                                         y_tick_steps            =   10,                                         # y axis major tick steps
    #                                         v_hist_bins             =   12,                                         # ver. histogram bins
    #                                         ))
        
        
        
    #     for hemi, ext, list_name in [('_L', '-left', 'lh'), ('_R', '-right', 'rh'), ('_LR', '', 'lrh')]:
        
    #         title                           =   '{0}{1} hemi: {2} (n={3})'.\
    #                                             format(roi_text, ext, typeData, data['pRF_gazeAll_fit{0}'.format(hemi)].shape[0])
    #         _                               =   params_pRFcor.update(dict(
    #                                                 main_fig_title          =    title,                            
    #                                                 dataMat                 =   [data['pRF_gazeLeft_fit{0}'.format(hemi)], 
    #                                                                              data['pRF_gazeRight_fit{0}'.format(hemi)],
    #                                                                              data['pRF_gazeAll_fit{0}'.format(hemi)]],
    #                                                 hemisphere              =   list_name,
    #                                                 condition               =   'cor',
    #                                                 typeData                =   typeData))                                     # title
    #                                                                                      # title

    #         out1, _                         =   plotter.draw_figure(parameters = params_pRFcor, plot = 'cor')
    #         _                               =   f_pRF_cor[list_name].append(out1)
 
    # for hemi, letter in [('lh', 'L'), ('rh', 'R'), ('lrh', '')]:
    #     # Combine figures
    #     # ---------------
    #     all_f                           =   gridplot([
    #                                                  [f_pRF_cor[hemi][0],    f_pRF_cor[hemi][1]],                                    # specify figure 1st row
    #                                                  [f_pRF_cor[hemi][2],    f_pRF_cor[hemi][3]],                                    # specify figure 2nd row
    #                                                  [f_pRF_cor[hemi][4],    f_pRF_cor[hemi][5]],                                    # specify figure 3rd row
    #                                                  [f_pRF_cor[hemi][6],    f_pRF_cor[hemi][7]]])                                    
    #     figs_folder                     =   opj(folders['fit_cor'], hemi.upper()) if hemi != 'lrh' else folders['fit_cor']
    #     output_file_html                =   opj(figs_folder,'{0}{1}_fitCor.html'.format(roi_text, letter))                             # define html file
    #     _                               =   output_file(output_file_html, title='{0} {1} left/right cor.'.format(roi_text, letter))      # get html saved
    #     _                               =   save(all_f)




    # # --------------
    # # pRF gain plots
    # # --------------

    # params_pRFgain                      =   param_all                                                               # get main params     
    # _                                   =   params_pRFgain.update(                                                  # update parameters
    #                                             dict(
    #                                             gaze_shift              =   analysis_info['gazeShift'],             # Right gaze position - Left gaze position
    #                                             x_range                 =   (-50, 50),                              # x axis limits
    #                                             x_label                 =   'Gaze gain ratio (%)',                  # x axis label
    #                                             x_tick_steps            =   10,                                     # x axis major tick steps
    #                                             hist_range              =   (0,0.5),                                # histogram x/y range
    #                                             hist_steps              =   0.5,                                    # x axis major tick steps
    #                                             h_hist_bins             =   50,                                     # hor. histogram bins
    #                                             x_source_label          =   'x_gain',                               # source label for x axis
    #                                             reg_weight_source_label =   'cv_rsq_all',                           # regression weight source label
    #                                             cmap                    =   'RdBu',                                 # colormap to use
    #                                             cmap_steps              =   200,                                    # colormap steps
    #                                             col_offset              =   0,                                      # colormap offset
    #                                             vmin                    =   -50,                                    # val min to draw with colormap
    #                                             vmax                    =   50,
    #                                             condition               =   'gain')                                     # val max to draw with colormap
    #                                             )


    # f_pRFgain = {'lh': [], 'rh': [], 'lrh': []}
    # first_pRFgain_fig = {'L': [], 'R': [], 'LR': []}

    # for typeData in ['Ecc.','Size','CV R2','Amplitude','Stim-ratio']:

    #     if typeData == 'Ecc.':
    #         _                           =   params_pRFgain.update(                                                  # update parameters
    #                                         dict(
    #                                             y_range                 =   (0, 10),                                # y axis limits
    #                                             y_label                 =   'Eccentricity (dva)',                   # y axis label
    #                                             y_tick_steps            =   2,                                      # y axis major tick steps
    #                                             y_source_label          =   'ecc_all',                              # y source label
    #                                             v_hist_bins             =   20)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'Size':
    #         _                           =   params_pRFgain.update(                                                  # update parameters
    #                                         dict(
    #                                             y_range                 =   (0, 4),                                 # y axis limits
    #                                             y_label                 =   'Size (dva)',                           # y axis label
    #                                             y_tick_steps            =   1,                                      # y axis major tick steps
    #                                             y_source_label          =   'sigma_all',                            # y source label
    #                                             v_hist_bins             =   16)                                     # ver. histogram bins
    #                                             )

    #     elif typeData == 'CV R2':
    #         _                           =   params_pRFgain.update(                                                  # update parameters
    #                                         dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits
    #                                             y_label                 =   'Cross-validated R2',                   # y axis label
    #                                             y_tick_steps            =   0.2,                                    # y axis major tick steps
    #                                             y_source_label          =   'cv_rsq_all',                           # y source label
    #                                             v_hist_bins             =   20)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'Amplitude':
    #         _                           =   params_pRFgain.update(                                                  # update parameters
    #                                         dict(
    #                                             y_range                 =   (0, 1),                                 # y axis limits
    #                                             y_label                 =   'Amplitude (z)',                        # y axis label
    #                                             y_tick_steps            =   0.25,                                   # y axis major tick steps
    #                                             y_source_label          =   'beta_all',                             # y source label
    #                                             v_hist_bins             =   16)                                     # ver. histogram bins
    #                                             )
    #     elif typeData == 'Stim-ratio':
    #         _                           =   params_pRFgain.update(                                                  # update parameters
    #                                         dict( 
    #                                             y_range                 =   (-10, 110),                             # y axis limits
    #                                             y_label                 =   'Stim ratio (%)',                       # y axis label
    #                                             y_tick_steps            =   10,                                     # y axis major tick steps
    #                                             y_source_label          =   'stim_ratio_all',                       # y source label
    #                                             v_hist_bins             =   12)                                     # ver. histogram bins
    #                                             )
        
    #     for hemi, ext, list_name in [('L', '-left', 'lh'), ('R', '-right', 'rh'), ('LR', '', 'lrh')]:
    #         title                           =   '{0}{1} hemi: Gain vs. {2} (n={3})'.\
    #                                             format(roi_text, ext, typeData, data['pRF_gazeAll_fit_{0}'.format(hemi)].shape[0])

    #         _                               =   params_pRFgain.update(dict(
    #                                             main_fig_title          =    title,
    #                                             dataMat                 =   [data['pRF_gazeLeft_fit_{0}'.format(hemi)],
    #                                                                          data['pRF_gazeRight_fit_{0}'.format(hemi)],
    #                                                                          data['pRF_gazeAll_fit_{0}'.format(hemi)]],
    #                                             hemisphere              =   list_name,
    #                                             typeData                =   typeData))                                     # title


    #         out1, out2, out3, out4          =   plotter.draw_figure(parameters   = params_pRFgain,
    #                                                                 old_main_fig = first_pRFgain_fig[hemi],
    #                                                                 plot         = 'gain')
    #         _                               =   f_pRFgain[list_name].append(out1)
    #         _                               =   first_pRFgain_fig[hemi].append(out2)
    #         data['pRF_gazeGain_avg_{0}'.format(hemi)][numroi,:]    =   out3
    #         data['pRF_gazeGain_std_{0}'.format(hemi)][numroi,:]    =   out4


    # # Combine figures
    # # ---------------
    # for hemi, letter in [('lh', 'L'), ('rh', 'R'), ('lrh', '')]:
    #     all_f                          =   gridplot([
    #                                                  [f_pRFgain[hemi][0],  f_pRFgain[hemi][1],  f_pRFgain[hemi][2]],
    #                                                  [f_pRFgain[hemi][3],  f_pRFgain[hemi][4],  None],
    #                                                  ])
    #     figs_folder                     =   opj(folders['prf_gain'], hemi.upper()) if hemi != 'lrh' else folders['prf_gain']
    #     output_file_html                =   opj(figs_folder,'{0}{1}_pRFGain.html'.format(roi_text, letter))                             # define html file
    #     _                               =   output_file(output_file_html, title='{0} {1} left/right cor.'.format(roi_text, letter))      # get html saved
    #     _                               =   save(all_f)



    # # --------------
    # # pRF shift maps
    # # --------------
    # params_pRFshift             =   param_all                                                                       # get main params
    # _                           =   params_pRFshift.update(dict(
    #                                         x_range                 =   (-8, 8),                                    # x axis limits
    #                                         y_range                 =   (-5, 5),                                    # y axis limits
    #                                         x_label                 =   'Horizontal coord. (dva)',                  # x axis label
    #                                         y_label                 =   'Vertical coord. (dva)',                    # y axis label
    #                                         x_tick_steps            =   2,                                          # x axis major tick steps
    #                                         y_tick_steps            =   2,                                          # y axis major tick steps
    #                                         stim_radius             =   2.65,                                       # 2.65 instead of 3 (seems that the circle radius gets somehow distorted
    #                                                                                                                 # by the rectangle shape of the figure. 2.65 is chosen from testing)
    #                                         p_width                 =   int(400 * 1.6),                             # individual plot width
    #                                         p_height                =   400,
    #                                         condition               =   'shift'
    #                                         ))


    # f_pRFshift = {'lh': [], 'rh': [], 'lrh':[]}
    # first_pRFshift_fig = {'lh': [], 'rh': [], 'lrh': []}

    # for typeData in ['xy pRF shift','x pRF shift','y pRF shift']:
    #     if typeData == 'xy pRF shift':
    #         _                           =   params_pRFshift.update(dict(
    #                                         xs                      =   'xs',                                       # define datasource list of list of x coordinates to plot for lines
    #                                         ys                      =   'ys',                                       # define datasource list of list of y coordinates to plot for lines
    #                                         x_left                  =   'x_left',                                   # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_left                  =   'y_left',                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         x_right                 =   'x_right',                                  # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_right                 =   'y_right'                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         ))

            
    #     elif typeData == 'x pRF shift':
    #         _                           =   params_pRFshift.update(dict(
    #                                         xs                      =   'xs',                                       # define list of list of x coordinates to plot for lines
    #                                         ys                      =   'ys_mean',                                  # define list of list of y coordinates to plot for lines
    #                                         x_left                  =   'x_left',                                   # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_left                  =   'y_mean',                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         x_right                 =   'x_right',                                  # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_right                 =   'y_mean'                                    # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         ))
    #     elif typeData == 'y pRF shift':
    #         _                           =   params_pRFshift.update(dict(
    #                                         xs                      =   'xs_mean',                                  # define list of list of x coordinates to plot for lines
    #                                         ys                      =   'ys',                                       # define list of list of y coordinates to plot for lines
    #                                         x_left                  =   'x_mean',                                   # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_left                  =   'y_left',                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         x_right                 =   'x_mean',                                   # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_right                 =   'y_right'                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         ))

        
    #     for hemi, ext, list_name in [('_L', '-left', 'lh'), ('_R', '-right', 'rh'), ('_LR', '', 'lrh')]:
    #         title                           =   '{0}{1} hemi: {2} (n={3})'.\
    #                                             format(roi_text, ext, typeData, data['pRF_gazeAll_fit{0}'.format(hemi)].shape[0])
            
    #         _                               =   params_pRFshift.update(dict(
    #                                                 main_fig_title          =    title,                            
    #                                                 dataMat                 =   [data['pRF_gazeLeft_fit{0}'.format(hemi)], 
    #                                                                              data['pRF_gazeRight_fit{0}'.format(hemi)],
    #                                                                              data['pRF_gazeAll_fit{0}'.format(hemi)]],
    #                                                 hemisphere              =   list_name,
    #                                                 typeData                =   typeData))                                     # title


    #         out1, out2                      =   plotter.draw_figure(parameters = params_pRFshift, 
    #                                                                 old_main_fig = first_pRFshift_fig[list_name],
    #                                                                 plot = 'shift')
    #         _                               =   f_pRFshift[list_name].append(out1)
    #         _                               =   first_pRFshift_fig[list_name].append(out2)


    # # Combine figures
    # # ---------------
    # for hemi, letter in [('lh', 'L'), ('rh', 'R'), ('lrh', '')]:

    #     all_f                           =   gridplot([
    #                                                  [f_pRFshift[hemi][0], f_pRFshift[hemi][1]], [f_pRFshift[hemi][2]]
    #                                                  ])                                    
    #     figs_folder                     =   opj(folders['prf_shift'], hemi.upper()) if hemi != 'lrh' else folders['prf_shift']
    #     output_file_html                =   opj(figs_folder,'{0}{1}_prfShift.html'.format(roi_text, letter))                             # define html file
    #     _                               =   output_file(output_file_html, title='{0} {1} left/right cor.'.format(roi_text, letter))      # get html saved
    #     _                               =   save(all_f)




    # # --------------------------
    # # pRF spatiotopic shift maps
    # # --------------------------
    # # X-Axis
    # dataMat_gazeLeft                         =   data['pRF_gazeLeft_fit_LR']
    # dataMat_gazeLeft[:,0]                    =   data['pRF_gazeLeft_fit_LR'][:,0] -4
    # dataMat_gazeRight                        =   data['pRF_gazeRight_fit_LR']
    # dataMat_gazeRight[:,0]                   =   data['pRF_gazeRight_fit_LR'][:,0] +4

    # # Y-Axis
    # dataMat_gazeLeft                         =   data['pRF_gazeLeft_fit_LR']
    # dataMat_gazeLeft[:,1]                    =   data['pRF_gazeLeft_fit_LR'][:,1] +2
    # dataMat_gazeRight                        =   data['pRF_gazeRight_fit_LR']
    # dataMat_gazeRight[:,1]                   =   data['pRF_gazeRight_fit_LR'][:,1] +2

    # params_pRFshift_spatio      =   param_all                                                                       # get main params
    # _                           =   params_pRFshift_spatio.update(dict(
    #                                         x_range                 =   (-8, 8),                                    # x axis limits
    #                                         y_range                 =   (-5, 5),                                    # y axis limits
    #                                         x_label                 =   'Horizontal coord. (dva)',                  # x axis label
    #                                         y_label                 =   'Vertical coord. (dva)',                    # y axis label
    #                                         x_tick_steps            =   2,                                          # x axis major tick steps
    #                                         y_tick_steps            =   1,                                          # y axis major tick steps
    #                                         p_width                 =   int(400 * 1.6),                                        # individual plot width
    #                                         p_height                =   400,
    #                                         stim_radius             =   2.65,                                       # 2.65 instead of 3 (seems that the circle radius gets somehow distorted
    #                                                                                                                 # by the rectangle shape of the figure. 2.65 is chosen from testing)
    #                                         dataMat                 =   [dataMat_gazeLeft,
    #                                                                      dataMat_gazeRight,
    #                                                                      data['pRF_gazeAll_fit_LR']],
    #                                         condition               =   'shift_spatio'
    #                                         ))


    # f_pRFshift                  = {'lh': [], 'rh': [], 'lrh': []}
    # first_pRFshift_fig          = {'lh': [], 'rh': [], 'lrh': []}

    # for typeData in ['xy pRF shift','x pRF shift']:
    #     if typeData == 'xy pRF shift':
    #         _                           =   params_pRFshift_spatio.update(dict(
    #                                         xs                      =   'xs',                                       # define datasource list of list of x coordinates to plot for lines
    #                                         ys                      =   'ys',                                       # define datasource list of list of y coordinates to plot for lines
    #                                         x_left                  =   'x_left',                                   # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_left                  =   'y_left',                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         x_right                 =   'x_right',                                  # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_right                 =   'y_right'                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         ))

    #     elif typeData == 'x pRF shift':
    #         _                           =   params_pRFshift_spatio.update(dict(
    #                                         xs                      =   'xs',                                       # define list of list of x coordinates to plot for lines
    #                                         ys                      =   'ys_mean',                                  # define list of list of y coordinates to plot for lines
    #                                         x_left                  =   'x_left',                                   # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_left                  =   'y_mean',                                   # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         x_right                 =   'x_right',                                  # define datasource x gazeLeft coordinates to plot for left pRF dots
    #                                         y_right                 =   'y_mean'                                    # define datasource y gazeLeft coordinates to plot for left pRF dots
    #                                         ))

    #     for hemi, ext, list_name in [('_L', '-left', 'lh'), ('_R', '-right', 'rh'), ('_LR', '', 'lrh')]:
    #         title                           =   '{0}{1} hemi: {2} (n={3})'.\
    #                                             format(roi_text, ext, typeData, data['pRF_gazeAll_fit{0}'.format(hemi)].shape[0])
            
    #         _                               =   params_pRFshift_spatio.update(dict(
    #                                                 main_fig_title          =    title,
    #                                                 hemisphere              =   list_name,
    #                                                 typeData                =   typeData))                                     # title


    #         out1, out2                      =   plotter.draw_figure(parameters = params_pRFshift_spatio, 
    #                                                                 old_main_fig = first_pRFshift_fig[list_name],
    #                                                                 plot = 'shift_spatiotopic')
    #         _                               =   f_pRFshift[list_name].append(out1)
    #         _                               =   first_pRFshift_fig[list_name].append(out2)


    # # Combine figures
    # # ---------------
    # for hemi, letter in [('lh', 'L'), ('rh', 'R'), ('lrh', '')]:

    #     all_f                           =   gridplot([
    #                                                  [f_pRFshift[hemi][0], f_pRFshift[hemi][1]]
    #                                                  ])                                    
    #     figs_folder                     =   opj(folders['prf_spatial_shift'], hemi.upper()) if hemi != 'lrh' else folders['prf_spatial_shift']
    #     output_file_html                =   opj(figs_folder,'{0}{1}_prfShift_spatial.html'.format(roi_text, letter))                             # define html file
    #     _                               =   output_file(output_file_html, title='{0} {1} left/right shift.'.format(roi_text, letter))      # get html saved
    #     _                               =   save(all_f)




# -------------------
# Basic ROI analysis
# -------------------
params_pRFroi                       =   param_all                                                                   # get main params     
_                                   =   params_pRFroi.update(                                                       # update parameters
                                            dict(
                                            gaze_shift              =   analysis_info['gazeShift'],                 # Right gaze position - Left gaze position
                                            x_range                 =   (0, len(rois)),                             # x axis limits
                                            x_label                 =   'Region of interest',                       # x axis label
                                            rois                    =   rois,                                       # rois x label
                                            bar_width               =   0.25,                                       # bar width
                                            bar_width_gain          =   0.5,                                     
                                            condition               =   'basic_roi')
                                            )

f_pRFroi = {'lh': [], 'rh': [], 'lrh': []}

for typeData in ['Ecc.','Size','CV R2','Amplitude','Voxels', 'Stim-ratio']:                     
    if typeData == 'Ecc.':
        _                           =   params_pRFroi.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (0, 5),                                    # y axis limits
                                            y_label                 =   'Eccentricity (dva)',                       # y axis label
                                            y_tick_steps            =   1,                                          # y axis major tick steps
                                            y_source_label          =   'ecc',                                      # y source label
                                            tootips_end             =   'dva',                                      # tooltip end
                                            tootips_txt             =   'Ecc.')                                    # tooltip text
                                            )
    elif typeData == 'Size':
        _                           =   params_pRFroi.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (0, 4),                                     # y axis limits
                                            y_label                 =   'Size (dva)',                               # y axis label
                                            y_tick_steps            =   1,                                          # y axis major tick steps
                                            y_source_label          =   'sigma',                                    # y source label
                                            tootips_end             =   'dva',                                      # tooltip end
                                            tootips_txt             =   'Size')                                    # tooltip text
                                            )

    elif typeData == 'CV R2':
        _                           =   params_pRFroi.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (0, 1),                                     # y axis limits
                                            y_label                 =   'Cross-validated R2',                       # y axis label
                                            y_tick_steps            =   0.2,                                        # y axis major tick steps
                                            y_source_label          =   'cv_rsq',                                   # y source label
                                            tootips_end             =   '',                                         # tooltip end
                                            tootips_txt             =   'CV R2')                                   # tooltip text
                                            )
    elif typeData == 'Amplitude':
        _                           =   params_pRFroi.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (0, 1),                                     # y axis limits
                                            y_label                 =   'Amplitude (z)',                            # y axis label
                                            y_tick_steps            =   0.25,                                       # y axis major tick steps
                                            y_source_label          =   'beta',                                     # y source label
                                            tootips_end             =   '',                                         # tooltip end
                                            tootips_txt             =   'Amp')                                     # tooltip text
                                            )
    elif typeData == 'Stim-ratio':
        _                           =   params_pRFroi.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (0, 100),                                   # y axis limits
                                            y_label                 =   'Stim ratio (%)',                           # y axis label
                                            y_tick_steps            =   10,                                         # y axis major tick steps
                                            y_source_label          =   'stim_ratio',                               # y source label
                                            tootips_end             =   '%',                                        # tooltip end
                                            tootips_txt             =   'Stim-ratio')                              # tooltip text
                                            )
    elif typeData == 'Voxels':
        _                           =   params_pRFroi.update(dict( 
                                            y_range                 =   (0, 500),                                   # y axis limits
                                            y_label                 =   'Number of voxels',                         # y axis label
                                            y_tick_steps            =   100,                                        # y axis major tick steps
                                            y_source_label          =   'voxel',                                    # y source label
                                            tootips_end             =   '',                                         # tooltip end
                                            tootips_txt             =   '# vox'))
    

    for hemi, title_ext, list_name in [('L', 'Left', 'lh'), ('R', 'Right', 'rh'), ('LR', '', 'lrh')]:
        title                           =   '{0} hemisphere: {1}'.format(title_ext, typeData)
        _                               =   params_pRFroi.update(dict(
                                                main_fig_title          =    title,
                                                dataMat                 =   [data['pRF_gazeAll_avg_{0}'.format(hemi)],   data['pRF_gazeAll_std_{0}'.format(hemi)],
                                                                             data['pRF_gazeAll_num_{0}'.format(hemi)],   data['pRF_gazeLeft_avg_{0}'.format(hemi)],
                                                                             data['pRF_gazeLeft_std_{0}'.format(hemi)],  data['pRF_gazeLeft_num_{0}'.format(hemi)],
                                                                             data['pRF_gazeRight_avg_{0}'.format(hemi)], data['pRF_gazeRight_std_{0}'.format(hemi)],
                                                                             data['pRF_gazeRight_num_{0}'.format(hemi)], data['pRF_gazeGain_avg_{0}'.format(hemi)],
                                                                             data['pRF_gazeGain_std_{0}'.format(hemi)],  rois,
                                                                             data['pRF_gazeLeft_fit_LR'],                data['pRF_gazeRight_fit_LR'],
                                                                             data['pRF_gazeAll_fit_LR']],
                                                hemisphere              =   list_name,
                                                typeData                =   typeData))                                     # title

        out1, out2                      =   plotter.draw_figure(parameters = params_pRFroi, plot = 'roi')
        _                               =   f_pRFroi[list_name].append(out1)


#import ipdb ; ipdb.set_trace()

# Combine figures
# ---------------
for hemi, letter in [('lh', 'L_'), ('rh', 'R_'), ('lrh', '')]:

    all_fL                          =   gridplot([
                                            [f_pRFroi[hemi][0],   f_pRFroi[hemi][1],  f_pRFroi[hemi][2],  f_pRFroi[hemi][3]],       # specify figure 1st row
                                            [f_pRFroi[hemi][4],   f_pRFroi[hemi][5],  None,  None],                 # specify figure 2nd row
                                            ])
    figs_folder                     =   opj(folders['prf_roi_basic'], hemi.upper()) if hemi != 'lrh' else folders['prf_roi_basic']
    output_L_file_html              =   opj(figs_folder,'{0}pRFroi_basic.html'.format(letter))                                            # define html file
    _                               =   output_file(output_L_file_html, title='ALL ROI{0} pRF roi'.format(' ' + letter))                      # get html saved
    _                               =   save(all_fL)




# -------------------
# Main ROI analysis
# -------------------
params_pRFroi_main                  =   param_all                                                                   # get main params     
_                                   =   params_pRFroi_main.update(                                                       # update parameters
                                            dict(
                                            gaze_shift              =   analysis_info['gazeShift'],                 # Right gaze position - Left gaze position
                                            x_range                 =   (0, len(rois)),                             # x axis limits
                                            x_label                 =   'Region of interest',                       # x axis label
                                            rois                    =   rois,                                       # rois x label
                                            bar_width               =   0.25,                                       # bar width
                                            bar_width_gain          =   0.5,
                                            condition               =   'main_roi')                                     
                                            )

f_pRFroi = {'lh': [], 'rh': [], 'lrh': []}

for typeData in ['Retinal X Gain', 'Screen X Gain', 'Retinal Gain Index', 'Amplitude Change', 'Y Change']:                     
    if typeData == 'Retinal X Gain':
        _                           =   params_pRFroi_main.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (-20, 110),                                  # y axis limits
                                            y_label                 =   'Retinal Gain (%)',                      # y axis label
                                            y_tick_steps            =   20,                                         # y axis major tick steps
                                            y_source_label          =   'retinal_x_gain',                                # y source label
                                            tootips_end             =   '%',                                      # tooltip end
                                            tootips_txt             =   'retinal gain ratio')                         # tooltip text
                                            )
    

    elif typeData == 'Screen X Gain':
        _                           =   params_pRFroi_main.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (-20, 110),                                  # y axis limits
                                            y_label                 =   'Screen Gain (%)',                      # y axis label
                                            y_tick_steps            =   20,                                         # y axis major tick steps
                                            y_source_label          =   'screen_x_gain',                                # y source label
                                            tootips_end             =   '%',                                      # tooltip end
                                            tootips_txt             =   'screen gain ratio')                         # tooltip text
                                            )
    elif typeData == 'Retinal Gain Index':
        _                           =   params_pRFroi_main.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (-0.5, 1.5),                                  # y axis limits
                                            y_label                 =   'Retinal Gain Index',                      # y axis label
                                            y_tick_steps            =   0.5,                                         # y axis major tick steps
                                            y_source_label          =   'retinal_gain_index',                                # y source label
                                            tootips_end             =   '',                                      # tooltip end
                                            tootips_txt             =   'Retinal Gain Index')                         # tooltip text
                                            )
    elif typeData == 'Amplitude Change':
        _                           =   params_pRFroi_main.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (-20, 110),                                  # y axis limits
                                            y_label                 =   'Amplitude Change',                      # y axis label
                                            y_tick_steps            =   20,                                         # y axis major tick steps
                                            y_source_label          =   'amplitude_change',                                # y source label
                                            tootips_end             =   '%',                                      # tooltip end
                                            tootips_txt             =   'Amplitude Change')                         # tooltip text
                                            )
                                            
    elif typeData == 'Y Change':
        _                           =   params_pRFroi_main.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (-20, 110),                                  # y axis limits
                                            y_label                 =   'Y Change',                      # y axis label
                                            y_tick_steps            =   20,                                         # y axis major tick steps
                                            y_source_label          =   'y_change',                                # y source label
                                            tootips_end             =   '%',                                      # tooltip end
                                            tootips_txt             =   'Y Change')                         # tooltip text
                                            )


    for hemi, title_ext, list_name in [('LR', '', 'lrh')]:
        title                           =   '{0} hemisphere: {1}'.format(title_ext, typeData)
        _                               =   params_pRFroi_main.update(dict(
                                                main_fig_title          =    title,
                                                dataMat                 =   [data['pRF_gazeAll_avg_{0}'.format(hemi)],   data['pRF_gazeAll_std_{0}'.format(hemi)],
                                                                             data['pRF_gazeAll_num_{0}'.format(hemi)],   data['pRF_gazeLeft_avg_{0}'.format(hemi)],
                                                                             data['pRF_gazeLeft_std_{0}'.format(hemi)],  data['pRF_gazeLeft_num_{0}'.format(hemi)],
                                                                             data['pRF_gazeRight_avg_{0}'.format(hemi)], data['pRF_gazeRight_std_{0}'.format(hemi)],
                                                                             data['pRF_gazeRight_num_{0}'.format(hemi)], data['pRF_gazeGain_avg_{0}'.format(hemi)],
                                                                             data['pRF_gazeGain_std_{0}'.format(hemi)],  rois,
                                                                             data['pRF_gazeLeft_fit_LR'],                data['pRF_gazeRight_fit_LR'],
                                                                             data['pRF_gazeAll_fit_LR']],
                                                hemisphere              =   list_name,
                                                typeData                =   typeData))                                     # title

        out1, out2                      =   plotter.draw_figure(parameters = params_pRFroi_main, plot = 'roi')
        _                               =   f_pRFroi[list_name].append(out1)


# Combine figures
# ---------------
for hemi, letter in [('lrh', '')]:

    all_fL                          =   gridplot([
                                            [f_pRFroi[hemi][0],   f_pRFroi[hemi][1]],       # specify figure 1st row
                                            [f_pRFroi[hemi][2],   f_pRFroi[hemi][3]],
                                            [f_pRFroi[hemi][4],   None            ]])
    figs_folder                     =   opj(folders['prf_roi_main'], hemi.upper()) if hemi != 'lrh' else folders['prf_roi_main']
    output_L_file_html              =   opj(figs_folder,'{0}pRFroi_main.html'.format(letter))                                            # define html file
    _                               =   output_file(output_L_file_html, title='ALL ROI{0} pRF roi'.format(' ' + letter))                      # get html saved
    _                               =   save(all_fL)

