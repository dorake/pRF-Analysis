def create_post_fmriprep_workflow(name='post_fmriprep_workflow'):
    """
    -----------------------------------------------------------------------------------------
    create_post_fmriprep_workflow(name='post_fmriprep_workflow'):
    -----------------------------------------------------------------------------------------
    Goal of the function:
    Perform different motion correction registration
    -----------------------------------------------------------------------------------------
    Input(s):
    name: name of workflow
    -----------------------------------------------------------------------------------------
    Output(s):
    post_fmriprep_workflow: pre-processing nipype workflow
    -----------------------------------------------------------------------------------------
    Workflow input(s):
    inputspec.fmriprep_directory: fmriprep directory
    inputspec.output_directory: post-fmriprep directory
    inputspec.FS_subject_dir: freesurfer subject directory

    inputspec.sub_id: freesurfer subject id
    inputspec.av_func: averaging function
    inputspec.pycortex_mask_types: list of masks
    inputspec.temp_dir_cortex: tempory directory for pycortex masks
    -----------------------------------------------------------------------------------------
    Workflow output(s):  
    => data through datasink
    -----------------------------------------------------------------------------------------
    """

    # Importation
    # ===========
    # general import
    import os.path as op

    # Nipype import
    import nipype.pipeline as pe
    from nipype.interfaces import fsl, freesurfer
    from nipype.interfaces.io import SelectFiles, DataSink
    from nipype.interfaces.utility import Function, IdentityInterface, Rename

    # Functions import
    from ..utils.GLM import fit_BD_glm_CV_one_fold, fit_pp_GLM
    from ..utils.utils import leave_one_out_lists, suff_file_name, pickfirst, FS_T1_file, select_condition_files, average_over_runs, add_session_2_cortex, create_cortical_mask, add_subject_2_cortex, add_flat_2_cortex

    # WORKFLOW DEFINITION
    # ===================
    post_fmriprep_workflow      =   pe.Workflow(
                                        name            = name)                                             # main workflow: pre-processing after fmriprep

    # DEFINE NODES
    # ============

    # General nodes
    # -------------
    # workflow input
    input_node                  =   pe.Node( IdentityInterface(
                                        fields          = ['fmriprep_directory',                            # fmriprep directory
                                                           'output_directory',                              # post-fmriprep directory
                                                           'FS_subject_dir',                                # freesurfer subject directory
                                                           'sub_id',                                        # subject name
                                                           'av_func',                                       # function used to average across run
                                                           'pycortex_mask_types',                           # different types of mask types of pycortex 
                                                           'temp_dir_cortex']),                             # temporary directory of pycortex masks
                                        name            = 'inputspec')                                      # node name

    # data input node: find and define files on the disk
    datasource_templates        =   dict(   
                                        func            = '{sub_id}/func/*space-T1w_preproc.nii.gz',        # bold data preprocessed by fmriprep
                                        confounds       = '{sub_id}/func/*confounds.tsv')                   # fmriprep confouds parameters
    
    datasource                  =   pe.Node( SelectFiles(                                                   # Find files on the disk and feed them into the workflow
                                        datasource_templates,                                               # see above
                                        sort_filelist   = True,                                             # Sort_filelist option to control whether sorting should be applied.
                                        aise_on_empty   = False),                                           # Generate exception if list is empty for a given field
                                        name            = 'datasource')                                     # node name


    # Datasink (Nipype's standard output module)
    datasink                    =   pe.Node( DataSink(),                                                    # DataSink function
                                        name            = 'sinker')                                         # node name

    datasink.inputs.parameterization                    = False                                             # Parameter controling whether thea data stored with DataSink is in a folder structure that contains iterable information or not


    # Filtering etc. nodes
    # --------------------

    # TF GLM: Performs a Cosine GLM on nifti-file
    pp_tf_glm                   =   pe.MapNode( Function(
                                        input_names     = ['niifile',                                       # Absolute path to nifti-file.
                                                           'fmriprep_design_matrix_file'],                  # fmriprep confouds parameters
                                        output_names    = ['res_file',                                      # Absolute path to nifti-file filtered nifti file
                                                           'res_Z_file'],                                   # Absolute path to nifti-file filtered residual files
                                        function        = fit_pp_GLM),                                      # Performs a Cosine GLM on nifti-file
                                        name            = 'tf_GLM',                                         # map node name
                                        iterfield       = ['niifile',                                       # input that should be iterated over
                                                           'fmriprep_design_matrix_file'])                  # input that should be iterated over

    # condition selection
    scf_PRF_gazeLeft            =   pe.Node( Function(
                                        input_names     = ['in_files',                                      # list of files
                                                          'condition'],                                     # text of a specific condition written in the file name
                                        output_names    = ['out_files'],                                    # list of files including the condition name
                                        function        = select_condition_files),                          # select the file of a specific condition
                                        name            = 'scf_PRF_gazeLeft')                               # name of the node

    scf_PRF_gazeLeft.inputs.condition                   = '_task-gazeLeft_'                                 # condition name in the files

    scf_PRF_gazeRight           =   pe.Node( Function(
                                        input_names     = ['in_files',                                      # list of files
                                                           'condition'],                                    # text of a specific condition written in the file name
                                        output_names    = ['out_files'],                                    # list of files including the condition name
                                        function        = select_condition_files),                          # selection the file of a specific condition
                                        name            = 'scf_PRF_gazeRight')                              # name of the node

    scf_PRF_gazeRight.inputs.condition                  = '_task-gazeRight_'                                # condition name in the files

    # node to select the first of a list 
    pickfirst_N                 =   pe.Node( Function(
                                        input_names     = ['files'],                                        # list of files
                                        output_names    = ['out_file'],                                     # first picked file name
                                        function        = pickfirst),                                       # select first out of a list
                                        name            = 'pickfirst_N')                                    # name of the node

    # find T1 file in FS directory
    FS_T1_file_node             = pe.Node( Function(
                                        input_names     = ['freesurfer_subject_ID',                         # freesurfer subject name
                                                           'freesurfer_subject_dir'],                       # freesurfer subject directory
                                        output_names    = 'T1_mgz_path',                                    # T1 file path
                                        function        = FS_T1_file),                                      # find T1 file in FS directory
                                        name            = 'FS_T1_file_node')                                # name of the node

    # convert files to niigz
    mriConvert_N                = pe.Node(freesurfer.MRIConvert(                                            # mri_convert is a general purpose utility for converting between different file formats
                                        out_type        = 'niigz'),                                         # output format
                                        name            = 'mriConvert_N')                                   # name of the node

    # Average conditions nodes
    # ------------------------

    # leave one out list
    loo_gazeLeft                = pe.Node(Function(
                                        input_names     = ['input_list'],                                   # list of items, for instance absolute paths to nii files
                                        output_names    = ['out_lists'],                                    # list of lists
                                        function        = leave_one_out_lists),                             # Creates a list of lists, with each element of the input_list leave out of the returned lists once, in order.
                                        name            = 'leave_one_out_lists_gazeLeft')                   # name of the node

    # leave one out averaging
    av_loo_gazeLeft             = pe.MapNode(Function(
                                        input_names     = ['in_files',                                      # absolute paths to nifti-files
                                                           'func',                                          # the function used to calculate the 'average' ('mean','median')
                                                           'output_filename'],                              # path where to output filename
                                        output_names    = ['out_file'],                                     # absolute path to averaged nifti-file.
                                        function        = average_over_runs),                               # takes a list of 4D fMRI nifti-files and averages them
                                        name            = 'average_over_runs_loo_gazeLeft',                 # name of the mapnode
                                        iterfield       = ['in_files',                                      # input that should be iterated over
                                                           'output_filename'])                              # input that should be iterated over

    # averaging over all runs
    av_all_gazeLeft             = pe.Node(Function(
                                        input_names     = ['in_files',                                      # absolute paths to nifti-files
                                                           'func'],                                         # the function used to calculate the 'average' ('mean','median')
                                        output_names    = ['out_file'],                                     # absolute path to averaged nifti-file.
                                        function        = average_over_runs),                               # takes a list of 4D fMRI nifti-files and averages them
                                        name            = 'average_over_runs_all_gazeLeft')                 # name of the mapnode

    # rename loo file
    rename_loo_gazeLeft         = pe.MapNode(Function(
                                        input_names     = ['in_file'],                                      # path to a nifti file
                                        output_names    = ['out_file'],                                     # path to the suffixed file
                                        function        = suff_file_name),                                  # put suffix _avg_loo.nii.gz at the end of a file path
                                        name            = 'rename_loo_gazeLeft',                            # name of the mapnode
                                        iterfield       = ['in_file'])                                      # input that should be iterated over
    
    # leave one out list
    loo_gazeRight               = pe.Node(Function(
                                        input_names     = ['input_list'],                                   # list of items, for instance absolute paths to nii files
                                        output_names    = ['out_lists'],                                    # list of lists
                                        function        = leave_one_out_lists),                             # Creates a list of lists, with each element of the input_list leave out of the returned lists once, in order.
                                        name            = 'leave_one_out_lists_gazeRight')                  # name of the node

    # leave one out averaging
    av_loo_gazeRight            = pe.MapNode(Function(
                                        input_names     = ['in_files',                                      # absolute paths to nifti-files
                                                           'func',                                          # the function used to calculate the 'average' ('mean','median')
                                                           'output_filename'],                              # path where to output filename
                                        output_names    = ['out_file'],                                     # absolute path to averaged nifti-file.
                                        function        = average_over_runs),                               # takes a list of 4D fMRI nifti-files and averages them
                                        name            ='average_over_runs_loo_gazeRight',                 # name of the mapnode
                                        iterfield       = ['in_files',                                      # input that should be iterated over
                                                          'output_filename'])                               # input that should be iterated over

    # rename loo file
    rename_loo_gazeRight        = pe.MapNode(Function(
                                        input_names     = ['in_file'],                                      # path to a nifti file
                                        output_names    = ['out_file'],                                     # path to the suffixed file
                                        function        = suff_file_name),                                  # put suffix _avg_loo.nii.gz at the end of a file path
                                        name            = 'rename_loo_gazeRight',                           # name of the mapnode
                                        iterfield       = ['in_file'])                                      # input that should be iterated over


    # averaging over all runs
    av_all_gazeRight            = pe.Node(Function(
                                        input_names     = ['in_files',                                      # absolute paths to nifti-files
                                                           'func'],                                         # the function used to calculate the 'average' ('mean','median')
                                        output_names    = ['out_file'],                                     # absolute path to averaged nifti-file.
                                        function        = average_over_runs),                               # takes a list of 4D fMRI nifti-files and averages them
                                        name            = 'average_over_runs_all_gazeRight')                # name of the mapnode

    # Registration nodes
    # ------------------

    # get mean image
    mean_bold                   = pe.Node(fsl.maths.MeanImage(                                              # generate a mean image across a given dimension.
                                        dimension       = 'T'),                                             # dimension to mean across
                                        name            = 'mean_space')                                     # name of the node

    # register to freesurfer volume
    bbregister_N                = pe.Node(freesurfer.BBRegister(                                            # register a volume to a Freesurfer anatomical
                                        init            = 'header',                                         # initialize registration spm, fsl, header
                                        contrast_type   = 't2',                                             # contrast type of image t1 or t2
                                        out_fsl_file    = True ),                                           # write the transformation matrix in FSL FLIRT format
                                        name            = 'bbregister_N')                                   # name of the node
    # remame epi files
    rename_epi                  = pe.Node(Rename(                                                           # change the name of a file based on a mapped format string
                                        format_string   = 'epi.nii.gz',                                     # formatting string for output template
                                        keep_ext        = False),                                           # keep in_file extension, replace non-extension component of name
                                        name            = 'rename_epi')                                     # name of the node

    # rename bbregistered files
    bbrdat                      = pe.Node(Rename(                                                           # change the name of a file based on a mapped format string
                                        format_string   = 'bbr.dat',                                        # formatting string for output template
                                        keep_ext        = False),                                           # keep in_file extension, replace non-extension component of name
                                        name            = 'bbrdat')                                         # name of the node
    # rename bbregistered matrices
    bbrmat                      = pe.Node(Rename(                                                           # change the name of a file based on a mapped format string
                                        format_string   = 'bbr.mat',                                        # formatting string for output template
                                        keep_ext        = False),                                           # keep in_file extension, replace non-extension component of name
                                        name            = 'bbrmat')                                         # name of the node

    # add subject to cortex
    add_subject_2_cortex_N       = pe.Node(Function(
                                        input_names     = ['fs_ID',                                         # freesurfer id
                                                           'fs_dir'],                                       # freesurfer directory
                                        output_names    = ['fs_ID'],                                        # freesurfer id (just to get something)
                                        function        = add_subject_2_cortex),                            # transform fsl to cortex format
                                        name            = 'add_subject_2_cortex_N')                         # name of the node

    # add flat surface to cortex
    add_flat_2_cortex_N         =  pe.Node(Function(
                                        input_names     = ['fs_ID',                                         # freesurfer id
                                                           'fs_dir'],                                       # freesurfer directory
                                        output_names    = ['fs_ID'],                                        # freesurfer id (just to get something)
                                        function        = add_flat_2_cortex),                               # transform fsl flat to cortex flat format
                                        name            = 'add_flat_2_cortex_N')                            # name of the node

    # tranform to cortex format
    add_session_2_cortex_N      = pe.Node(Function(
                                        input_names     = ['epi_2_T1',                                      # epi to T1 affine transform
                                                           'epi',                                           # epi file
                                                           'T1',                                            # T1 weighted file
                                                           'subject',                                       # subject name
                                                           'xfmname'],                                      # transform name
                                        output_names    = ['xfmname'],                                      # transform name
                                        function        = add_session_2_cortex),                            # transform fsl to cortex format
                                        name            = 'add_session_2_cortex_N')                         # name of the node

    add_session_2_cortex_N.inputs.xfmname               = 'fmriprep'                                        # xfmname input

    # create cortical cortex mask
    create_cortical_mask_N      = pe.MapNode(Function(
                                        input_names     = ['subject',                                       # subject name
                                                           'xfmname',                                       # transform name
                                                           'type',                                          # mask type ['cortical','thin','thick','nearest']
                                                           'epi',                                           # epi file
                                                           'temp_dir_cortex'],                              # temporary directory for pycortex masks
                                        output_names    = ['out_file'],                                     # boolean mask array for cortical voxels in functional space
                                        function        = create_cortical_mask),                            # get cortical mask for a particular transform
                                        iterfield       = ['type'],                                         # input that should be iterated over
                                        name            = 'create_cortical_mask_N')                         # map node name


    # WORKFLOW CONNECTIONS
    # ====================

    #---------------------------------+-----------------------+---------------------------+-----------------------+-----------------------------------+-------------------------------------------------------------------------------+
    #                                 |      source node      |     source output         |   destination node    |   destination input               | Connection goal                                                               |
    #---------------------------------+-----------------------+---------------------------+-----------------------+-----------------------------------+-------------------------------------------------------------------------------+

    # data source
    post_fmriprep_workflow.connect(     input_node,             'fmriprep_directory',       datasource,             'base_directory'                ) # Give datasource the fmriprep data directory
    post_fmriprep_workflow.connect(     input_node,             'sub_id',                   datasource,             'sub_id'                        ) # Give datasource the subject name
    
    # data sink
    post_fmriprep_workflow.connect(     input_node,             'output_directory',         datasink,               'base_directory'                ) # Give datasink the output directory path

    # TF-GLM filtering
    post_fmriprep_workflow.connect(     datasource,             'func',                     pp_tf_glm,              'niifile'                       ) # Give tf_glm function the functional data
    post_fmriprep_workflow.connect(     datasource,             'confounds',                pp_tf_glm,              'fmriprep_design_matrix_file'   ) # Give tf_glm function the confounds of fmriprep
    post_fmriprep_workflow.connect(     pp_tf_glm,              'res_file',                 datasink,               'tf.res'                        ) # Give the tf_glm filtered file to the datasink
    post_fmriprep_workflow.connect(     pp_tf_glm,              'res_Z_file',               datasink,               'tf.Z'                          ) # Give the residuals values image to datasink 
    post_fmriprep_workflow.connect(     pp_tf_glm,              'res_Z_file',               scf_PRF_gazeLeft,       'in_files'                      ) # Give the residuals values image to condition selector gazeLeft (?) (why the residuals and not the filted file)
    post_fmriprep_workflow.connect(     pp_tf_glm,              'res_Z_file',               scf_PRF_gazeRight,      'in_files'                      ) # Give the residuals values image to condition selector gazeRight (?) (why the residuals and not the filted file)

    # all run average - gazeLeft
    post_fmriprep_workflow.connect(     scf_PRF_gazeLeft,       'out_files',                av_all_gazeLeft,        'in_files'                      ) # Give files to average
    post_fmriprep_workflow.connect(     input_node,             'av_func',                  av_all_gazeLeft,        'func'                          ) # Specified averaging function for loo averaging
    post_fmriprep_workflow.connect(     av_all_gazeLeft,        'out_file',                 datasink,               'av_gazeLeft.all'               ) # Give averaged files to datasink

    # Leave one out averaging - gazeLeft
    post_fmriprep_workflow.connect(     scf_PRF_gazeLeft,       'out_files',                rename_loo_gazeLeft,    'in_file'                       ) # Give selected files to rename
    post_fmriprep_workflow.connect(     scf_PRF_gazeLeft,       'out_files',                loo_gazeLeft,           'input_list'                    ) # Give selected files to leave one out selection
    post_fmriprep_workflow.connect(     input_node,             'av_func',                  av_loo_gazeLeft,        'func'                          ) # Specified averaging function for loo averaging
    post_fmriprep_workflow.connect(     loo_gazeLeft,           'out_lists',                av_loo_gazeLeft,        'in_files'                      ) # Give selected files to averaged
    post_fmriprep_workflow.connect(     rename_loo_gazeLeft,    'out_file',                 av_loo_gazeLeft,        'output_filename'               ) # Specified output filename
    post_fmriprep_workflow.connect(     av_loo_gazeLeft,        'out_file',                 datasink,               'av_gazeLeft.loo'               ) # Give averaged files to datasink

    # all run average - gazeRight
    post_fmriprep_workflow.connect(     scf_PRF_gazeRight,      'out_files',                av_all_gazeRight,       'in_files'                      ) # Give files to average
    post_fmriprep_workflow.connect(     input_node,             'av_func',                  av_all_gazeRight,       'func'                          ) # Specified averaging function for loo averaging
    post_fmriprep_workflow.connect(     av_all_gazeRight,       'out_file',                 datasink,               'av_gazeRight.all'               ) # Give averaged files to datasink

    # Leave one out averaging - gazeRight
    post_fmriprep_workflow.connect(     scf_PRF_gazeRight,      'out_files',                rename_loo_gazeRight,   'in_file'                       ) # Give selected files to rename
    post_fmriprep_workflow.connect(     scf_PRF_gazeRight,      'out_files',                loo_gazeRight,          'input_list'                    ) # Give selected files to leave one out selection
    post_fmriprep_workflow.connect(     input_node,             'av_func',                  av_loo_gazeRight,       'func'                          ) # Specified averaging function for loo averaging
    post_fmriprep_workflow.connect(     loo_gazeRight,          'out_lists',                av_loo_gazeRight,       'in_files'                      ) # Give selected files to averaged
    post_fmriprep_workflow.connect(     rename_loo_gazeRight,   'out_file',                 av_loo_gazeRight,       'output_filename'               ) # Specified output filename
    post_fmriprep_workflow.connect(     av_loo_gazeRight,       'out_file',                 datasink,               'av_gazeRight.loo'              ) # Give averaged files to datasink

    # Registration: get files from freesurfer
    post_fmriprep_workflow.connect(     datasource,             'func',                     pickfirst_N,            'files'                         ) # Give pre-processed func file to selec the first of the list (why the first?) 
    post_fmriprep_workflow.connect(     pickfirst_N,            'out_file',                 mean_bold,              'in_file'                       ) # Mean the selected func file over time
    post_fmriprep_workflow.connect(     input_node,             'sub_id',                   bbregister_N,           'subject_id'                    ) # Give subject name to bbregister
    post_fmriprep_workflow.connect(     input_node,             'FS_subject_dir',           bbregister_N,           'subjects_dir'                  ) # Give the FS subject directory to bbregister
    post_fmriprep_workflow.connect(     mean_bold,              'out_file',                 bbregister_N,           'source_file'                   ) # Give the meaned func file to bbregister
    post_fmriprep_workflow.connect(     bbregister_N,           'out_fsl_file',             bbrmat,                 'in_file'                       ) # Give the bbregistered FLIRT style registration file to rename as bbr.mat 
    post_fmriprep_workflow.connect(     bbregister_N,           'out_reg_file',             bbrdat,                 'in_file'                       ) # Give the bbregisterd registration file to rename as bbr.dat
    post_fmriprep_workflow.connect(     mean_bold,              'out_file',                 rename_epi,             'in_file'                       ) # give meaned bold to rename as epi
    post_fmriprep_workflow.connect(     bbrmat,                 'out_file',                 datasink,               'reg.@bbrmat'                   ) # give renamed BBRegister FLIRT style registration to datasink
    post_fmriprep_workflow.connect(     bbrdat,                 'out_file',                 datasink,               'reg.@bbrdat'                   ) # give renamed BBRegister registration to datasink
    post_fmriprep_workflow.connect(     rename_epi,             'out_file',                 datasink,               'reg.@epi'                      ) # give renamed mean bold epi to datasink

    # Registration: create pycortex database
    post_fmriprep_workflow.connect(     input_node,             'sub_id',                   add_subject_2_cortex_N, 'fs_ID'                         ) # give freesurfer subject id to create subject in database
    post_fmriprep_workflow.connect(     input_node,             'FS_subject_dir',           add_subject_2_cortex_N, 'fs_dir'                        ) # give freesurfer directory to create subject in database
    post_fmriprep_workflow.connect(     add_subject_2_cortex_N, 'fs_ID',                    add_flat_2_cortex_N,    'fs_ID'                         ) # give freesurfer subject id to create flat surfaces in subject database
    post_fmriprep_workflow.connect(     input_node,             'FS_subject_dir',           add_flat_2_cortex_N,    'fs_dir'                        ) # give freesurfer directory to create flat surfaces in subject database

    # Registration: create pycortex files
    post_fmriprep_workflow.connect(     add_flat_2_cortex_N,    'fs_ID',                    add_session_2_cortex_N, 'subject'                       ) # give subject name to create pycortex file
    post_fmriprep_workflow.connect(     bbrmat,                 'out_file',                 add_session_2_cortex_N, 'epi_2_T1'                      ) # give bbregistered FLIRT style registration to create pycortex file
    post_fmriprep_workflow.connect(     rename_epi,             'out_file',                 add_session_2_cortex_N, 'epi'                           ) # give mean bold file to create pycortex file
    post_fmriprep_workflow.connect(     input_node,             'sub_id',                   FS_T1_file_node,        'freesurfer_subject_ID'         ) # give freesurfer subject name to find the free surfer t1 file
    post_fmriprep_workflow.connect(     input_node,             'FS_subject_dir',           FS_T1_file_node,        'freesurfer_subject_dir'        ) # give freesurfer subject directory to find the free surfer t1 file
    post_fmriprep_workflow.connect(     FS_T1_file_node,        'T1_mgz_path',              mriConvert_N,           'in_file'                       ) # give the t1w freesurfer file to convert in niigz format
    post_fmriprep_workflow.connect(     mriConvert_N,           'out_file',                 add_session_2_cortex_N, 'T1'                            ) # give the t1w converted freesurfer file to create pycortex file
    post_fmriprep_workflow.connect(     add_session_2_cortex_N, 'xfmname',                  create_cortical_mask_N, 'xfmname'                       ) # give transform name to create the pycortex masks
    post_fmriprep_workflow.connect(     rename_epi,             'out_file',                 create_cortical_mask_N, 'epi'                           ) # give mean bold file to create the pycortex masks
    post_fmriprep_workflow.connect(     input_node,             'pycortex_mask_types',      create_cortical_mask_N, 'type'                          ) # give the type of masks to create the pycortex masks
    post_fmriprep_workflow.connect(     input_node,             'sub_id',                   create_cortical_mask_N, 'subject'                       ) # give subject name to create the pycortex masks
    post_fmriprep_workflow.connect(     input_node,             'temp_dir_cortex',          create_cortical_mask_N, 'temp_dir_cortex'               ) # give temporary directory to create the pycortex masks
    post_fmriprep_workflow.connect(     create_cortical_mask_N, 'out_file',                 datasink,               'masks.@cortex'                 ) # give the pycortex masks to datasink

    #---------------------------------+-----------------------+---------------------------+-----------------------+-----------------------------------+-------------------------------------------------------------------------------+

    # OUTPUT(S)
    # =========
    return post_fmriprep_workflow