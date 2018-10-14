from __future__ import division, print_function
from popeye.spinach import generate_og_receptive_field
from IPython import embed as shell


def to_polar(all_prf_data, scaling_array=None):
    import numpy as np
    # polar angle conversion, assuming that x and y are 0 and 1 in array
    complex_peak_polar = all_prf_data[..., 0] + 1j * all_prf_data[..., 1]
    if scaling_array != None:
        normed_complex_peak_polar = complex_peak_polar / \
            np.abs(complex_peak_polar)
        normed_scaled_complex_peak_polar = normed_complex_peak_polar * scaling_array
        return normed_scaled_complex_peak_polar
    else:
        return complex_peak_polar


def to_ecc(all_prf_data, min_ecc=0, max_ecc=30):
    import numpy as np
    # eccentricity conversion, assuming that x and y are 0 and 1 in array
    ecc = np.sqrt(all_prf_data[..., 0]**2 + all_prf_data[..., 1]**2)

    ecc[ecc < min_ecc] = min_ecc
    ecc[ecc > max_ecc] = max_ecc

    return ecc

def convert_fit_results(prf_filenames,
                        output_dir,
                        stim_radius     =   3.0,
                        typeData        =  'loo'):
    """
    Convert pRF fitting value in different parameters for following analysis
   
    Parameters
    ----------
    prf_filenames: absolute paths to prf result files. 
    output_dir: absolute path to directory into which to put the resulting files.
    stim_radius: stimulus radius in deg
    typeData: sting indicating the type of pRF analysis made ('loo': based on leave-one-out analysis; 'all': based on average of all )

    Returns
    -------
    None
    """

    # Imports
    # -------
    # General imports
    import os
    import nibabel as nb
    import glob
    import numpy as np

    # Popeye imports
    from popeye.spinach import generate_og_receptive_fields

    try:
        os.makedirs(os.path.join(output_dir,'all'))                                                            # create pRF output data all pRF
        os.makedirs(os.path.join(output_dir,'pos'))                                                            # create pRF output data positive pRF
        os.makedirs(os.path.join(output_dir,'neg'))                                                            # create pRF output data negative pRF
    except:
        pass

    # Get data details
    # ----------------
    aff                             =   nb.load(prf_filenames[0]).affine                                        # get image affine
    hdr                             =   nb.load(prf_filenames[0]).header                                        # get image header
    prf_data                        =   np.array([nb.load(pf).get_data() for pf in prf_filenames])              # load all data output of pRF fitting

    # Create derived measures from prfs
    # ---------------------------------
    print('\nextracting %s pRF parameters...'%typeData)

    # pRF sign
    # --------
    prf_mean                        =   prf_data.mean(axis=0)                                                   # mean across cv folds of the 
    prf_sign_all                    =   np.sign((prf_mean[..., 3]))                                             # sign of prf amplitude (beta)

    pos_mask                        =   prf_sign_all > 0.0                                                      # positive pRF mask
    neg_mask                        =   prf_sign_all < 0.0                                                      # negative pRF mask
    all_mask                        =   pos_mask | neg_mask #  pos_mask + neg_mask                              # all pRF mask
    
    # Averaging 
    # ---------
    prf_mean_all                    =   prf_mean                                                                # mean across cv folds of the 
    
    # Cross-validated r-square or r-square
    # -------------------------------------
    prf_cv_rsq_all                  =   prf_mean[..., 6]                                                    # define mean rsquared values

    # pRF eccentricity
    # ----------------
    prf_ecc_all                     =   np.nan_to_num(np.sqrt(prf_mean[..., 0]**2 + prf_mean[..., 1]**2))       # determine eccentricity from fovea

    # pRF polar angle
    # ---------------   
    complex_polar                   =   prf_mean[..., 0] + 1j * prf_mean[..., 1]                                # define complex number with coordinates of pRF
    normed_polar                    =   complex_polar / np.abs(complex_polar)                                   # normalized them to give a unity vector
    prf_polar_real_all              =   np.real(normed_polar)                                                   # real element of polar angle
    prf_polar_imag_all              =   np.imag(normed_polar)                                                   # imaginary element of polar angle
    
    # pRF size
    # --------
    prf_size_all                    =   prf_mean[..., 2].astype(np.float64)                                     # take the sigma value of the Gaussian
    prf_size_all[prf_size_all<1e-4] =   1e-4                                                                    # convert 0 in 1e-4

    # pRF amplitude
    # -------------
    prf_amp_all                     =   prf_mean[..., 3]

    # pRF baseline
    # ------------
    prf_baseline_all                =   prf_mean[..., 4]

    # pRF proportion of stimulus region
    # ---------------------------------
    deg_x, deg_y                    =   np.meshgrid(np.linspace(-30, 30, 60), np.linspace(-30, 30, 60))       # define prfs in visual space
    flat_prf_pars                   =   prf_mean.reshape((-1, prf_mean.shape[-1])).astype(np.float64)           # (?)
    
    rfs                             =   generate_og_receptive_fields(                                           # generate gaussian
                                                                      flat_prf_pars[..., 0],                    # coordinate of the center x of the Gaussian
                                                                      flat_prf_pars[..., 1],                    # coordinate of the center y of the Gaussian
                                                                      prf_size_all.ravel(),                     # dispersion of the Gaussian
                                                                      np.ones(np.prod(flat_prf_pars.shape[0])), # amplitude of the Gaussian
                                                                      deg_x,                                    # coordinate matrix along the horizontal dimension of the display (degrees)
                                                                      deg_y)                                    # coordinate matrix along the vertical dimension of the display (degrees)

    total_prf_content               =   rfs.reshape((-1, flat_prf_pars.shape[0])).sum(axis=0)                   # total vs "inside the stimulus region" volume of the prf
    stim_vignet                     =   np.sqrt(deg_x ** 2 + deg_y**2) < stim_radius                            # stim coordinates
    
    prf_stim_ratio_all              =   rfs[stim_vignet, :].sum(axis=0) / total_prf_content                     # compute in stim ratio
    prf_stim_ratio_all              =   prf_stim_ratio_all.reshape(prf_size_all.shape)                          # reshape the data to original

    # Saving
    # ------
    for output_type in ['prf_sign','prf_mean','prf_cv_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_amp','prf_baseline','prf_stim_ratio']:
        for mask_dir in ['all','pos','neg']:

            print('writing: %s'%(os.path.join(output_dir,"{mask_dir}","{typeData}_{output_type}_{mask_dir}.nii.gz").format(typeData = typeData, mask_dir = mask_dir, output_type = output_type)))
            exec('{output_type}_{mask_dir} = np.copy({output_type}_all)'.format(mask_dir = mask_dir, output_type = output_type))
            exec('{output_type}_{mask_dir}[~{mask_dir}_mask] = np.nan'.format(mask_dir = mask_dir, output_type = output_type))

            exec('out_img = nb.Nifti1Image(dataobj = {output_type}_{mask_dir}, affine = aff, header = hdr)'.format(mask_dir = mask_dir, output_type = output_type))
            exec('out_img.to_filename(os.path.join(output_dir,"{mask_dir}","{typeData}_{output_type}_{mask_dir}.nii.gz"))'.format(typeData = typeData, mask_dir = mask_dir, output_type = output_type))
             
    return None

def get_mask_from_roi(roi_filename,
                      all_fit_prf_filename,
                      cv_rsq_diff_filename,
                      rsq_threshold=-0.1,
                      rsq_diff_threshold=0.0,
                      nr_voxels=0):
    import nibabel as nb
    import numpy as np

    # the anatomical mask for both hemispheres
    roi_mask = np.array(nb.load(roi_filename).get_data()).astype(bool) 
    # the rsq diff for CV
    cv_rsq_diff = nb.load(cv_rsq_diff_filename).get_data()
    cv_rsq_mask = cv_rsq_diff > rsq_diff_threshold

    all_prf_data = nb.load(all_fit_prf_filename).get_data().astype(np.float64)
    all_prf_rsq = all_prf_data[..., -2]
    # sign of prf amplitude
    prf_sign = np.sign(np.nan_to_num(all_prf_data[..., 4]))

    # create a mask
    signed_rsq = np.nan_to_num(all_prf_rsq * prf_sign)

    # if we need a specified # of voxels instead of an rsq threshold,
    # we recompute the mask based on this
    if nr_voxels != 0:
        anat_mask = roi_mask & cv_rsq_mask
        signed_rsq[~anat_mask] = 0
        nr_voxels_in_mask = 1
        if rsq_threshold > 0:
            rsq_threshold = np.max(signed_rsq)
            while np.sum(signed_rsq > rsq_threshold) < nr_voxels:
                rsq_threshold -= 0.001
                # print('rsq_threshold %1.3f, nr of voxels: %i'%(rsq_threshold, np.sum(signed_rsq < rsq_threshold)))
        elif rsq_threshold < 0:
            rsq_threshold = np.min(signed_rsq)
            while np.sum(signed_rsq < rsq_threshold) < nr_voxels:
                rsq_threshold += 0.001
                # print('rsq_threshold %1.3f, nr of voxels: %i'%(rsq_threshold, np.sum(signed_rsq < rsq_threshold)))
        print('Estimated rsq threshold for %i voxels is %1.3f' %
              (nr_voxels, rsq_threshold))

    # do this regardless of nr_voxels
    if rsq_threshold > 0:
        rsq_mask = signed_rsq > rsq_threshold
    else:
        rsq_mask = signed_rsq < rsq_threshold

    all_mask = roi_mask & rsq_mask & cv_rsq_mask

    if nr_voxels == 0:
        print('Set rsq threshold %1.3f returns %i voxels' %
              (rsq_threshold, all_mask.sum()))

    return all_mask


def draw_cortex_volume(
            subject,                                                                                            # subject id (e.g. 'sub-001')
            xfmname,                                                                                            # xfm name used in pycortex database
            data,                                                                                               # the data you would like to plot on a flatmap
            cmap,                                                                                               # colormap that shoudl be used for plotting
            vmin,                                                                                               # minimal value iqn colormap
            vmax,                                                                                               # maximal value in colormap
            cbar                    =   'discrete',                                                             # color bar layout
            cmap_steps              =   255,                                                                    # number of colormap bins
            alpha                   =   None,                                                                   # alpha map
            depth                   =   1,                                                                      # Value between 0 and 1 for how deep to sample the surface for the flatmap (0 = gray/white matter boundary, 1 = pial surface)
            thick                   =   1,                                                                      # Number of layers through the cortical sheet to sample. Only applies for pixelwise = True
            height                  =   1024,                                                                   # Height of the image to render.
            sampler                 =   'nearest',                                                              # Name of the sampling function used
            with_curvature          =   True,                                                                   # Display the rois, labels, colorbar, annotated flatmap borders, or cross-hatch dropout?
            with_dropout            =   True,                                                                   # Display the cross-hatch dropout?
            with_labels             =   True,                                                                  # Display labels?
            with_colorbar           =   False,                                                                  # Display pycortex' colorbar?
            with_borders            =   False,                                                                  # Display borders?
            curv_brightness         =   0.95,                                                                   # Mean brightness of background. 0 = black, 1 = white, intermediate values are corresponding grayscale values.
            curv_contrast           =   0.05,                                                                   # Contrast of curvature. 1 = maximal contrast (black/white), 0 = no contrast (solid color for curvature equal to curvature_brightness).
            add_roi                 =   False,                                                                  # add roi to overlay.svg
            roi_name                =   'empty',                                                                # roi name
            col_offset              =   0):                                                                      # color offset (0 to 1)

    """
    Plot brain data onto a previously saved flatmap.

    Parameters
    ----------
    subject             : subject id (e.g. 'sub-001')
    xfmname             : xfm name used in pycortex database
    data                : the data you would like to plot on a flatmap
    cmap                : colormap that shoudl be used for plotting
    vmin                : minimal value iqn colormap
    vmax                : maximal value in colormap
    cbar                : color bar layout
    cmap_steps          : number of colormap bins
    alpha               : alpha map
    depth               : Value between 0 and 1 for how deep to sample the surface for the flatmap (0 = gray/white matter boundary, 1 = pial surface)
    thick               : Number of layers through the cortical sheet to sample. Only applies for pixelwise = True
    height              : Height of the image to render. Automatically scales the width for the aspect of the subject's flatmap
    sampler             : Name of sampling function used to sample underlying volume data. Options include 'trilinear', 'nearest', 'lanczos'
    with_curvature      : Display the rois, labels, colorbar, annotated flatmap borders, or cross-hatch dropout?
    with_dropout        : Display the cross-hatch dropout?
    with_labels         : Display labels?
    with_colorbar       : Display pycortex' colorbar?
    with_borders        : Display borders?
    curv_brightness     : Mean brightness of background. 0 = black, 1 = white, intermediate values are corresponding grayscale values.
    curv_contrast       : Contrast of curvature. 1 = maximal contrast (black/white), 0 = no contrast (solid color for curvature equal to curvature_brightness).
    add_roi             : add roi to overlay.svg
    roi_name            : roi name
    col_offset          : colormap offset between 0 and 1

    Returns
    -------
    volrgb - pycortex volume file
    """

    import cortex
    import pdb
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from matplotlib import cm
    import matplotlib as mpl

    # define colormap
    base                            =   cortex.utils.get_cmap(cmap)
    
    val                             =   np.fmod(np.linspace(0+col_offset, 1+col_offset,cmap_steps+1,endpoint=False),1.0)

    colmap                          =   colors.LinearSegmentedColormap.from_list(
                                                        'my_colmap',                                            # new colormap name
                                                        base(val),                                              # segmented data
                                                        N               =   cmap_steps)                         # number of steps

    # convert data to RGB
    vrange                          =   float(vmax) - float(vmin)                                                 
    norm_data                       =   ((data-float(vmin))/vrange)*cmap_steps
    mat                             =   colmap(norm_data.astype(int))

    volrgb                          =   cortex.VolumeRGB(                                                       # RGB colors for each voxel in a volumetric dataset
                                                    red             =   mat[...,0],                             # volume that represent the red component (1D or 3D)
                                                    green           =   mat[...,1],                             # volume that represent the green component (1D or 3D)
                                                    blue            =   mat[...,2],                             # volume that represent the blue component (1D or 3D)
                                                    alpha           =   alpha,                                  # alpha color
                                                    subject         =   subject,                                # subject identifier in database
                                                    xfmname         =   xfmname)                                # transform name in database
            
    volrgb_fig                      =   cortex.quickflat.make_figure(                                           # create a quickshow
                                                    braindata           =   volrgb,                             # the data you would like to plot on a flatmap
                                                    depth               =   depth,                              # walue between 0 and 1 for how deep to sample the surface for the flatmap (0 = gray/wm boundary, 1 = pial surface)
                                                    thick               =   thick,                              # number of layer through the cortical shet to sample
                                                    height              =   height,                             # height of the image to render
                                                    sampler             =   sampler,                            # sampling function used to sample underlying volume data
                                                    with_curvature      =   with_curvature,                     # display the curvature
                                                    with_dropout        =   with_dropout,                       # display the dropout
                                                    with_labels         =   with_labels,                        # display the labels
                                                    with_colorbar       =   with_colorbar,                      # display the colorbar
                                                    with_borders        =   with_borders,                       # show borders
                                                    curvature_brightness=   curv_brightness,                    # Mean brightness of background. 0 = black, 1 = white, intermediate values are corresponding grayscale values.
                                                    curvature_contrast  =   curv_contrast)                      # Contrast of curvature. 1 = maximal contrast (black/white), 0 = no contrast (solid color for curvature equal to curvature_brightness) 

        
    ## COLOR BARS
    ## ----------

    colorbar_location=[0.5, 0.07, 0.8, 0.2]                                                                     # specifices the location and size of the color bar on the flatmap                                                                                 # 
    n = 200                                                                                                     # the number of secants for the mesh

    if cbar == 'polar':

        ## POLAR ANGLE COLOR BAR
        ## ---------------------                                                                                
        cbar_axis = volrgb_fig.add_axes(colorbar_location, projection='polar')                                  # adds the color bar axis to the flatmap
        norm = mpl.colors.Normalize(0, 2*np.pi)                                                                 # defines colormap normalization for 0 to 2*pi

        # Plot a color mesh on the polar plot
        # with the color set by the angle
        #t = np.linspace(0,2*np.pi,n)                                                                            # theta values
        t = np.linspace(2*np.pi,0,n)                                                                            # theta values
        r = np.linspace(1,0,2)                                                                                  # raidus values - 0 for full (filled) circle
        rg, tg = np.meshgrid(r,t)                                                                               # create a radius,theta meshgrid
        c = tg                                                                                                  # define color values as theta value
        im = cbar_axis.pcolormesh(t, r, c.T,norm= norm, cmap = colmap)                                          # plot the colormesh on axis with specified colormap
        cbar_axis.set_theta_zero_location("W")                                                                         # rotates the start (theta 0)
        #cbar_axis.set_theta_offset(pi)
        cbar_axis.set_yticklabels([])                                                                           # turn off radial tick labels (yticks)
        cbar_axis.set_xticklabels([])                                                                           # turn off xticks
        #cbar_axis.tick_params(pad=15,labelsize=24)                                                             # cosmetic changes to tick labels
        cbar_axis.spines['polar'].set_visible(False)                                                            # turn off the axis spine.

        ## Doing the horizontal Meridian
    elif cbar == 'ecc':
    
        ## ECC COLOR BAR
        ## -------------
        colorbar_location=[0.5, 0.07, 0.8, 0.2]                                                                 # specifices the location and size of the color bar on the flatmap
        cbar_axis = volrgb_fig.add_axes(colorbar_location, projection='polar')                                  # adds the color bar axis to the flatmap   

        t = np.linspace(0,2*np.pi, n)                                                                           # theta values
        r = np.linspace(0,1, n)                                                                                 # radius values change 0.6 to 0 for full circle
        rg, tg = np.meshgrid(r,t)                                                                               # create a r,theta meshgrid
        c = tg                                                                                                  # define color values as theta meshgrid values
        #bounds = np.linspace(vmin, vmax, cmap_steps + 1)
        #bounds_label = np.linspace(vmin, vmax, 5)

        im = cbar_axis.pcolormesh(t, r, c, norm= mpl.colors.Normalize(0, 2*np.pi), cmap = colmap)               # plot the colormesh on axis with colormap
        cbar_axis.tick_params(pad=1,labelsize=15)                                                               # adjust padding between labels and labelsize
        cbar_axis.spines['polar'].set_visible(False)                                                            # turn off the axis spine.
        
        ## superimpose new axis for dva labeling
        box= cbar_axis.get_position()                                                                           # get position matrix of colorbar axis
        cbar_axis.set_yticklabels([])                                                                           # delete 'standard' yticks
        cbar_axis.set_xticklabels([])                                                                           # delete 'standard' xticks
        axl=volrgb_fig.add_axes([1.8*box.xmin,                                                                  # put it half way between the edge of the 1st subplot and the left edge of the figure, 1.91 for having the axis next to figure
                                 0.5*(box.ymin+box.ymax),                                                       # put the origin at the same height of the origin of the polar plots
                                 box.width/600,                                                                 # distance between spine and ticks
                                 box.height*0.5],                                                               # height of ticks relative to color bar
                                 axisbg=None)                                                                   # transparent background.
        axl.spines['top'].set_visible(False)                                                                    # invisible top spine
        axl.spines['right'].set_visible(False)                                                                  # invisible right spine
        axl.spines['bottom'].set_visible(False)                                                                 # invisible bottom spine
        axl.yaxis.set_ticks_position('right')                                                                   # specify ticks position relative to spine
        axl.xaxis.set_ticks_position('none')                                                                    # 'none' position since x-axis is not shown
        axl.set_xticklabels([])                                                                                 # empty list for xticklabels to hide it
        axl.set_yticklabels([0,3,6], size='x-large')                                                            # specify yticklabels
        axl.set_ylabel('$dva$\t\t', rotation=0, size='x-large')                                                 # set name, rotation and size of ylabel
        axl.yaxis.set_label_coords(box.xmax+30,0.4)                                                             # set location of ylabel
        axl.patch.set_alpha(0.5)                                                                                # set alpha of axis, 0.5 for transparent


    elif cbar == 'discrete':

        ## DISCRETE COLOR BARS
        ## -------------------
        colorbar_location= [0.9, 0.05, 0.03, 0.25] #[0.7, 0.07, 0.03, 0.3]                                      # specifices the location and size of the color bar on the flatmap

        cmaplist = [colmap(i) for i in range(colmap.N)]                                                         # extract all colors from colormap

        # define the bins and normalize
        bounds = np.linspace(vmin, vmax, cmap_steps + 1)                                                        # define bins of colors (default = 255 bins)
        bounds_label = np.linspace(vmin, vmax, 3)                                                               # create labels for the colors
        norm = mpl.colors.BoundaryNorm(bounds, colmap.N)                                                        # Generate a colormap index based on discrete intervals

        
        cbar_axis = volrgb_fig.add_axes(colorbar_location)                                                      # create a second axes to add colorbar                       
        cb = mpl.colorbar.ColorbarBase(cbar_axis,cmap=colmap,norm=norm,ticks=bounds_label,boundaries=bounds)    # add colorbar to axis


    if add_roi == True:
        _                               =   cortex.utils.add_roi(
                                                        data                =   volrgb,                         # data used to generate the flatmap image
                                                        name                =   roi_name,                       # name that will be assigned to the sub-layer in the rois.svg file
                                                        open_inkscape       =   False,                          # open inkscape
                                                        add_path            =   False,                          # If True also adds a sub-layer to the `rois` new SVG layer will automatically be created in the ROI group with the same `name` as the overlay
                                                        depth               =   depth,                          # walue between 0 and 1 for how deep to sample the surface for the flatmap (0 = gray/wm boundary, 1 = pial surface)
                                                        thick               =   thick,                          # number of layer through the cortical shet to sample
                                                        sampler             =   sampler,                        # sampling function used to sample underlying volume data
                                                        with_curvature      =   with_curvature,                 # display the curvature
                                                        with_dropout        =   with_dropout,                   # display the dropout
                                                        with_colorbar       =   with_colorbar,                  # display the colorbar
                                                        with_borders        =   with_borders,                   # show borders
                                                        curvature_brightness=   curv_brightness,                # Mean brightness of background. 0 = black, 1 = white, intermediate values are corresponding grayscale values.
                                                        curvature_contrast  =   curv_contrast)                  # Contrast of curvature. 1 = maximal contrast (black/white), 0 = no contrast (solid color for curvature equal to curvature_brightness) 
    
    return volrgb


def weighted_avg_std(matrix_input,matrix_weight,axisnum=0):
    import numpy as np

    wavg                        =    np.average(                                                                 # compute average weighted by cv rsq
                                                a                    =    matrix_input,                          # matrix to average
                                                axis                 =    axisnum,                               # axis to mean
                                                weights              =    matrix_weight)                         # weight of the average



    wstd                        =    np.sqrt(np.average(                                                         # compute weighted std by cv rsq
                                                a                    =    (matrix_input-wavg)**2,                # matrix to average (variance)
                                                axis                 =    axisnum,                               # axis to mean
                                                weights              =    matrix_weight))                        # weight of the average
    
    numvox                      =    matrix_input.shape[0]                                                       # number of voxel in input matrix
    return wavg, wstd, numvox
