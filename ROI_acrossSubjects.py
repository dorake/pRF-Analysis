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
from operator import itemgetter
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

rois                            =   analysis_info['rois']                                                           # define regions of interest

# Check if it received the 'save' argument to save the figures as svg files
try:saving_figs                 =   True if str(sys.argv[1]).lower() == 'save' else False
except IndexError: saving_figs  =   False


# Python 3 Warning because of Pycortex
sys.exit('This script works with Pycortex, which requires Python 2.x to work without issues. The current environment seems to have a higher version. Aborting.') if sys.version_info[0] > 2 else None

subjects_data = {}
for subj in ['sub-004', 'sub-003', 'sub-002']:
    sub_id = subj
    
    # Define folders
    # --------------
    base_folder =   opj(analysis_info['aeneas_project_directory'], 'pp', sub_id)
    folders     =   {}
    for folder in ['figs', 'masks', 'h5']:
        folders.update({folder : opj(base_folder, folder)})

    try: os.makedirs(opj(base_folder, 'figs', 'across_subjects'))
    except OSError: pass
    folders.update({'across_subjects' : opj(base_folder, 'figs', 'across_subjects')})

    for folder in ['ecc_x_gain','ecc_amp', 'fit_cor', 'prf_shift', 'prf_spatial_shift', 'prf_gain', 'prf_roi_basic', 'prf_roi_main', 'svg_files']:
        folders.update({folder : opj(folders['across_subjects'], folder)})
        
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

    subjects_data.update({sub_id : data})    


def across_subjects_values(subjects_data):
    # Average lists across subjects
    subjects            = list(subjects_data.keys())
    arrays              = list(subjects_data[subjects[0]].keys())

    average_data        = {}
    std_data            = {}

    for idx_a, array in enumerate(arrays):
        
        # Get subject with largest array
        array_list      = []
        for sub in subjects:
            array_list.append((subjects_data[sub][array].shape[0], sub))
        
        # Sort array sizes to get subject name with largest array
        array_list.sort(key = itemgetter(1), reverse = True)
        largest_array_subject = array_list[0][1]

        # Create nan array with shape of largest array
        ar = np.zeros([len(subjects)] + list(subjects_data[largest_array_subject][array].shape))

        # copy subject's array into just created array to mean across subjects
        for idx_s, sub in enumerate(subjects):
            subject_shape = subjects_data[sub][array].shape
            ar[idx_s, : subject_shape[0], : subject_shape[1]] = subjects_data[sub][array]
        
        print('Across Subjects Averaging of Array: {0}'.format(array))
        mean_array = np.nanmean(ar, axis = 0)
        std_array  = np.nanstd(ar, axis = 0)

        average_data.update({array : mean_array})
        std_data.update({array : std_array})

    return average_data, std_data


average_data, std_data = across_subjects_values(subjects_data = subjects_data)


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
                                        across_subjects         =   True
                                        )
# Initialize plotting class
# -------------------------
plotter                             =   PlotOperator(**param_all)




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
        title                           =   'Both Hemispheres: {0}'.format(typeData)
        _                               =   params_pRFroi.update(dict(
                                                main_fig_title          =    title,
                                                dataMat                 =   [average_data['pRF_gazeAll_avg_{0}'.format(hemi)],   average_data['pRF_gazeAll_std_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeAll_num_{0}'.format(hemi)],   average_data['pRF_gazeLeft_avg_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeLeft_std_{0}'.format(hemi)],  average_data['pRF_gazeLeft_num_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeRight_avg_{0}'.format(hemi)], average_data['pRF_gazeRight_std_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeRight_num_{0}'.format(hemi)], average_data['pRF_gazeGain_avg_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeGain_std_{0}'.format(hemi)],  rois,
                                                                             average_data['pRF_gazeLeft_fit_LR'],                average_data['pRF_gazeRight_fit_LR'],
                                                                             average_data['pRF_gazeAll_fit_LR']],
                                                hemisphere              =   list_name,
                                                typeData                =   typeData))                                     # title

        out1, out2                      =   plotter.draw_figure(parameters = params_pRFroi, plot = 'roi')
        _                               =   f_pRFroi[list_name].append(out1)

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

# for typeData in ['Retinal X Gain', 'Screen X Gain', 'Retinal Gain Index', 'Amplitude Change', 'Y Change']:
for typeData in ['Retinotopy Index']:                                          
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
    elif typeData == 'Retinotopy Index':
        _                           =   params_pRFroi_main.update(                                                       # update parameters
                                        dict(
                                            y_range                 =   (-0.5, 1.5),                                  # y axis limits
                                            y_label                 =   'Retinotopy Index',                      # y axis label
                                            y_tick_steps            =   0.5,                                         # y axis major tick steps
                                            y_source_label          =   'retinal_gain_index',                                # y source label
                                            tootips_end             =   '',                                      # tooltip end
                                            tootips_txt             =   'Retinotopy Index')                         # tooltip text
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
        title                           =   'Both Hemispheres: {0}'.format(typeData)
        _                               =   params_pRFroi_main.update(dict(
                                                main_fig_title          =    title,
                                                dataMat                 =   [average_data['pRF_gazeAll_avg_{0}'.format(hemi)],   average_data['pRF_gazeAll_std_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeAll_num_{0}'.format(hemi)],   average_data['pRF_gazeLeft_avg_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeLeft_std_{0}'.format(hemi)],  average_data['pRF_gazeLeft_num_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeRight_avg_{0}'.format(hemi)], average_data['pRF_gazeRight_std_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeRight_num_{0}'.format(hemi)], average_data['pRF_gazeGain_avg_{0}'.format(hemi)],
                                                                             average_data['pRF_gazeGain_std_{0}'.format(hemi)],  rois,
                                                                             average_data['pRF_gazeLeft_fit_LR'],                average_data['pRF_gazeRight_fit_LR'],
                                                                             average_data['pRF_gazeAll_fit_LR']],
                                                hemisphere              =   list_name,
                                                typeData                =   typeData))                                     # title

        out1, out2                      =   plotter.draw_figure(parameters = params_pRFroi_main, plot = 'roi')
        _                               =   f_pRFroi[list_name].append(out1)


# Combine figures
# ---------------
for hemi, letter in [('lrh', '')]:

    # all_fL                          =   gridplot([
    #                                         [f_pRFroi[hemi][0],   f_pRFroi[hemi][1]],       # specify figure 1st row
    #                                         [f_pRFroi[hemi][2],   f_pRFroi[hemi][3]],
    #                                         [f_pRFroi[hemi][4],   None            ]])
    all_fL                          =   gridplot([[f_pRFroi[hemi][0]]])
    figs_folder                     =   opj(folders['prf_roi_main'], hemi.upper()) if hemi != 'lrh' else folders['prf_roi_main']
    output_L_file_html              =   opj(figs_folder,'{0}pRFroi_main.html'.format(letter))                                            # define html file
    _                               =   output_file(output_L_file_html, title='ALL ROI{0} pRF roi'.format(' ' + letter))                      # get html saved
    _                               =   save(all_fL)

