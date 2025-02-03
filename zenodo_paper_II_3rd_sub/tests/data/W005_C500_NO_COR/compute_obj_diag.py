

import glob
import os
import numpy as np
import xarray as xr


# run_name = 'WANG1_FR_lat30'
#run_name = 'W05_C500'
# run_name = 'W05_C500_NO_COR'
# run_name = 'test'
run_name = 'W005_C500_NO_COR'
#sampling = 'Csv_m1'
# sampling = 'w_sieb_p01'
# sampling = 'Cw_m1'
sampling = 'Cw_m2'
path = '../data/'+run_name+'/'
# saving_name = '_object_diags_'+sampling+'_276h.nc'
# saving_name = '_object_diags_'+sampling+'_264h.nc'
saving_name = '_object_diags_'+sampling+'_72h.nc'


files = ['GN_01.1.OC_01.036.nc']
# files = ['GN_01.1.OC_01.022.nc']
# files = ['GN_01.1.OC_01.023.nc']

#files = ['GN_01.1.OC_01.012.nc']

saving_path = '../data/'+run_name+'/'+run_name+saving_name
n_snap = 1  # number of snapshots starting from the last

nu = 1.15e-6  # m^2.s-1 kinematic viscosity


m = 1  # standard deviation treshold in the CS sampling

percent = 1  # fraction of minimum STD in CS and Cw sampling

alpha = 2.0 * 10 ** (-4)
g = 9.81
rho0 = 1024.
cp = 3900.

Q0 = 200
B0 = g * alpha * Q0 / (rho0 * cp)

# ======================================FUNCTIONS=========================================


def subdomain_variables(mask_name):

    def MEAN(X_name):
        return variables[X_name].where(masks[mask_name]).mean(dim=('ni', 'nj')).rename(mask_name+'_'+X_name)
  # Compute MEANS
    MEAN_VARS[mask_name] = {var: MEAN(var) for var in variables}
    MEAN_VARS[mask_name]['FRAC'] = masks[mask_name].mean(
        dim=('ni', 'nj')).rename(mask_name+'_'+'FRAC')
    MEAN_VARS[mask_name]['MASS_FLUX'] = (
        MEAN_VARS[mask_name]['FRAC']*MEAN_VARS[mask_name]['WT']).rename(mask_name+'_MASS_FLUX')

    ## Add a computation of vorticity and divergence
    MEAN_VARS[mask_name]['VORT_x'] = ( variables['dyWT'] - variables['dzVT'] ).where(masks[mask_name]).mean(dim=('ni', 'nj')).rename(mask_name+'_'+'VORT_x')
    MEAN_VARS[mask_name]['VORT_y'] = (-variables['dxWT'] + variables['dzUT'] ).where(masks[mask_name]).mean(dim=('ni', 'nj')).rename(mask_name+'_'+'VORT_y')
    MEAN_VARS[mask_name]['VORT_z'] = ( variables['dxVT'] - variables['dyUT'] ).where(masks[mask_name]).mean(dim=('ni', 'nj')).rename(mask_name+'_'+'VORT_z')
    MEAN_VARS[mask_name]['DIV_h']  = ( variables['dxUT'] + variables['dyVT'] ).where(masks[mask_name]).mean(dim=('ni', 'nj')).rename(mask_name+'_'+'DIV_h')



    def intra_covariance(X_name, Y_name):
        return (MEAN_VARS[mask_name]['FRAC'] * ((variables[X_name] - MEAN_VARS[mask_name][X_name]) * (variables[Y_name] - MEAN_VARS[mask_name][Y_name])).where(masks[mask_name]).mean(dim=('ni', 'nj'))).rename(mask_name+'_intra'+'_'+X_name[:-1]+Y_name[:-1])

    def inter_covariance(X_name, Y_name):
        """Computes the inter-domain covariance between curent mask and 'TOT' mask"""
        return (MEAN_VARS[mask_name]['FRAC'] * (MEAN_VARS[mask_name][X_name] - MEAN_VARS['TOT'][X_name]) * (MEAN_VARS[mask_name][Y_name] - MEAN_VARS['TOT'][Y_name])).rename(mask_name+'_inter'+'_'+X_name[:-1]+Y_name[:-1])

    def triple_intra_covariance(X_name, Y_name, Z_name):
        return (MEAN_VARS[mask_name]['FRAC'] * ((variables[X_name] - MEAN_VARS[mask_name][X_name]) * (variables[Y_name] - MEAN_VARS[mask_name][Y_name]) * (variables[Z_name] - MEAN_VARS[mask_name][Z_name])).where(masks[mask_name]).mean(dim=('ni', 'nj'))).rename(mask_name+'_intra'+X_name[:-1]+Y_name[:-1]+Z_name[:-1])

    # Compute intra domain covariances
    intra_cov[mask_name] = {}
    intra_cov[mask_name]['WTH'] = intra_covariance('WT', 'THT')
    intra_cov[mask_name]['WSV'] = intra_covariance(
        'WT', 'SVT001').rename(mask_name+'_intra'+'_WSV')
    intra_cov[mask_name]['WU'] = intra_covariance('WT', 'UT')
    intra_cov[mask_name]['WV'] = intra_covariance('WT', 'VT')
    intra_cov[mask_name]['WPHI_over_RHO_0'] = (
        1/rho0 * intra_covariance('WT', 'PHIT')).rename(mask_name+'_intra'+'_WPHI_over_RHO_0')
    intra_cov[mask_name]['TH2'] = intra_covariance('THT', 'THT')
    intra_cov[mask_name]['SV2'] = intra_covariance(
        'SVT001', 'SVT001').rename(mask_name+'_intra'+'_SV2')
    intra_cov[mask_name]['U2'] = intra_covariance('UT', 'UT')
    intra_cov[mask_name]['V2'] = intra_covariance('VT', 'VT')
    intra_cov[mask_name]['W2'] = intra_covariance('WT', 'WT')
    intra_cov[mask_name]['TKE'] = (1/2 * (intra_cov[mask_name]['U2'] + intra_cov[mask_name]
                                   ['V2'] + intra_cov[mask_name]['W2'])).rename(mask_name+'_intra'+'_TKE')

    
    ## Add a computation of dissipation
    # intra_cov[mask_name]['diss'] = (nu*(
    #     intra_covariance('dxUT', 'dxUT') + intra_covariance('dxVT', 'dxVT') + intra_covariance('dxWT', 'dxWT') +
    #     intra_covariance('dyUT', 'dyUT') + intra_covariance('dyVT', 'dyVT') + intra_covariance('dyWT', 'dyWT') +
    #     intra_covariance('dzUT', 'dzUT') + intra_covariance('dzVT',
    #                                                         'dzVT') + intra_covariance('dzWT', 'dzWT')
    # )).rename(mask_name+'_intra'+'_diss')



    if mask_name == 'TOT':
        triple_cov[mask_name] = {}
        triple_cov[mask_name]['WKE'] = (1/2 * (triple_intra_covariance('WT', 'UT', 'UT') + triple_intra_covariance(
            'WT', 'VT', 'VT') + triple_intra_covariance('WT', 'WT', 'WT'))).rename(mask_name+'_intra'+'_WKE')

    if mask_name != 'TOT':
        # compute "mass-flux type" (ie inter domain) covariances

        inter_cov[mask_name] = {}
        inter_cov[mask_name]['WTH'] = inter_covariance('WT', 'THT')
        inter_cov[mask_name]['WSV'] = inter_covariance(
            'WT', 'SVT001').rename(mask_name+'_inter'+'_WSV')
        inter_cov[mask_name]['WU'] = inter_covariance('WT', 'UT')
        inter_cov[mask_name]['WV'] = inter_covariance('WT', 'VT')
        inter_cov[mask_name]['WPHI_over_RHO_0'] = (
            1/rho0 * intra_covariance('WT', 'PHIT')).rename(mask_name+'_inter'+'_WPHI_over_RHO_0')
        inter_cov[mask_name]['TH2'] = inter_covariance('THT', 'THT')
        inter_cov[mask_name]['SV2'] = inter_covariance(
            'SVT001', 'SVT001').rename(mask_name+'_inter'+'_SV2')
        inter_cov[mask_name]['U2'] = inter_covariance('UT', 'UT')
        inter_cov[mask_name]['V2'] = inter_covariance('VT', 'VT')
        inter_cov[mask_name]['W2'] = inter_covariance('WT', 'WT')
        inter_cov[mask_name]['TKE'] = (1/2 * (inter_cov[mask_name]['U2'] + inter_cov[mask_name]
                                       ['V2'] + inter_cov[mask_name]['W2'])).rename(mask_name+'_inter'+'_TKE')

        # compute the 4 terms "mass-flux" WKE
        triple_cov[mask_name] = {}
        triple_cov[mask_name]['WKE1'] = (1/2 * (triple_intra_covariance('WT', 'UT', 'UT') +
                                         triple_intra_covariance('WT', 'VT', 'VT') + triple_intra_covariance('WT', 'WT', 'WT')))
        triple_cov[mask_name]['WKE1'] = triple_cov[mask_name]['WKE1'].rename(
            mask_name+'_WKE1')

        # FRAC is already contained in intra cov
        triple_cov[mask_name]['WKE2'] = (
            (MEAN_VARS[mask_name]['UT'] - MEAN_VARS['TOT']['UT']) * intra_cov[mask_name]['WU'])
        triple_cov[mask_name]['WKE2'] += (
            (MEAN_VARS[mask_name]['VT'] - MEAN_VARS['TOT']['VT']) * intra_cov[mask_name]['WV'])
        triple_cov[mask_name]['WKE2'] += (
            (MEAN_VARS[mask_name]['WT'] - MEAN_VARS['TOT']['WT']) * intra_cov[mask_name]['W2'])
        triple_cov[mask_name]['WKE2'] = triple_cov[mask_name]['WKE2'].rename(
            mask_name+'_WKE2')

        triple_cov[mask_name]['WKE3'] = (
            MEAN_VARS[mask_name]['WT'] - MEAN_VARS['TOT']['WT']) * intra_cov[mask_name]['TKE']
        triple_cov[mask_name]['WKE3'] = triple_cov[mask_name]['WKE3'].rename(
            mask_name+'_WKE3')

        triple_cov[mask_name]['WKE4'] = (
            MEAN_VARS[mask_name]['WT'] - MEAN_VARS['TOT']['WT']) * inter_cov[mask_name]['TKE']
        triple_cov[mask_name]['WKE4'] = triple_cov[mask_name]['WKE4'].rename(
            mask_name+'_WKE4')
    print('Computation of '+mask_name + ' variables done')


def mixed_layer():
    # mld from min of turb flux
    # /!\ has to be called AFTER subdomain_variables('TOT')
    mld_f = z_r[(-intra_cov['TOT']['WTH']).argmax()]
    mld_f = mld_f.rename('MLD_FLUX').drop_vars(
        ['level', 'level_w']).expand_dims(dim={'time': ds.time.data})
    level_adim = (
        z_r/mld_f).rename('z_r_adimensional').drop_vars(['level_w', 'time'])

    print('Computation of mld_flux done, mld='+str(mld_f.data)+'m')

    return mld_f, level_adim


def merger(variable_group, exclude_TOT=False):
    '''merges variable_group which is a dict [mask_name] of dict [var] into a single dataset
        /!\ using .name of the variable as the NAME, NOT the var_name key of the dict '''
    dsnew = xr.Dataset({})

    for mask_name in masks:
        if (exclude_TOT == True) & (mask_name == 'TOT'):
            continue  # skip

        for var_name in variable_group[mask_name]:
            dsnew = xr.merge([dsnew, variable_group[mask_name][var_name]])
    return dsnew


# ============================================================================================================


# ==========================RUNNING the functions


#files = [os.path.basename(f) for f in sorted(glob.glob(path + 'GN*.nc'))][1:]
#files = [files[i] for i in range(len(files) - n_snap, len(files) - 1)] + [files[-1]]

# files = ['GN_01.1.OC_01.036.nc']
# # files = ['GN_01.1.OC_01.022.nc']
# # files = ['GN_01.1.OC_01.023.nc']

# #files = ['GN_01.1.OC_01.012.nc']

dsout = xr.Dataset({'initialize': 0})

for file in files:
    print('opening '+path+file)
    ds = xr.open_dataset(path + file)

    H = ds.level_w[-1]
    # removing ghost points and mapping on ocean standar coord
    z_r = ds.level[1:-1] - H
    # removing ghost points and mapping on ocean standar coord
    z_w = ds.level_w[1:] - H

    vars = ['UT', 'VT', 'WT', 'THT', 'SVT001', 'PHIT', 'TKET']
    variables = {var: ds[var] for var in vars}

    # interpolation of velocities on the grid (level, n_j, n_i), ie center of the cells. It interpolates with half sums, and so inter
    # typically U[i] = (UT[i]+UT[i+1])/2
    variables['UT'] = variables['UT'].interp(ni_u=ds.ni, nj_u=ds.nj)
    variables['VT'] = variables['VT'].interp(ni_v=ds.ni, nj_v=ds.nj)
    variables['WT'] = variables['WT'].interp(level_w=ds.level)

    # removing Nans and ghosts points, every field live on the X,Y,Z grid consiting of center of grids boxes
    variables = {var: variables[var][:, 1:-1, 3:-3, 3:-3] for var in variables}

    # ==================================
    # Compute Dissipation & vorticity
    # ==================================

    u = variables['UT']
    v = variables['VT']
    w = variables['WT']

    variables['dxUT'] = variables['UT'].differentiate('ni')
    variables['dxVT'] = variables['VT'].differentiate('ni')
    variables['dxWT'] = variables['WT'].differentiate('ni')

    variables['dyUT'] = variables['UT'].differentiate('nj')
    variables['dyVT'] = variables['VT'].differentiate('nj')
    variables['dyWT'] = variables['WT'].differentiate('nj')

    variables['dzUT'] = variables['UT'].differentiate('level')
    variables['dzVT'] = variables['VT'].differentiate('level')
    variables['dzWT'] = variables['WT'].differentiate('level')

    # intializing outputs
    MEAN_VARS = {}
    intra_cov = {}
    inter_cov = {}
    triple_cov = {}
    # initializing masks with the TOTAL domain (always True)
    masks = {'TOT': 1+xr.zeros_like(variables['THT'])}

    # compute TOT variables
    subdomain_variables('TOT')

    # compute mld
    mld_f, level_adim = mixed_layer()

    # CONDITIONAL SAMPLING

    if sampling == 'Csv_m1':
        # computation of minimum stdev from Couvreux, obtained as 5% of vertical mean stdev
        sigmaSVmin = np.abs(intra_cov['TOT']['SV2'])
        for k in range(z_r.size):
            sigmaSVmin[:, k] = percent/(- z_r[k]) * (np.sqrt(intra_cov['TOT']['SV2'].isel(
                level=[i for i in range(k, z_r.size)]))).integrate('level')

        # Compute sampling masks FROM TRACER: Couvreux 2010, Pergaud 2009, Berg,Stull 2002 on SV

        m = m
        masks['DW'] = (((variables['WT'] - MEAN_VARS['TOT']['WT']) < 0) & ((variables['SVT001'] - MEAN_VARS['TOT']
                       ['SVT001']) > m*np.maximum(sigmaSVmin, np.sqrt(intra_cov['TOT']['SV2']))))  # for down drafts
        masks['UP'] = ~masks['DW']

    if sampling == 'Cw_m1':
        # computation of minimum stdev from Couvreux, obtained as 5% (10%) of vertical mean stdev
        sigmaWmin = np.abs(intra_cov['TOT']['W2'])
        for k in range(z_r.size):
            # sigmaWmin[:, k] = 0.05/(- z_r[k]) * (np.sqrt(intra_cov['TOT']['W2'].isel(
            sigmaWmin[:, k] = percent/(- z_r[k]) * (np.sqrt(intra_cov['TOT']['W2'].isel(
                level=[i for i in range(k, z_r.size)]))).integrate('level')

        # Compute sampling masks FROM W: Couvreux 2010, Pergaud 2009, Berg,Stull 2002 on SV

        m = m
        masks['DW'] = (((variables['WT'] - MEAN_VARS['TOT']['WT']) < 0) & ((variables['WT'] - MEAN_VARS['TOT']
                                                                            ['WT']) < - m*np.maximum(sigmaWmin, np.sqrt(intra_cov['TOT']['W2']))))  # for down drafts
        masks['UP'] = ~masks['DW']

    if sampling == 'Cw_m2':
        # computation of minimum stdev from Couvreux, obtained as 5% (10%) of vertical mean stdev
        sigmaWmin = np.abs(intra_cov['TOT']['W2'])
        for k in range(z_r.size):
            # sigmaWmin[:, k] = 0.05/(- z_r[k]) * (np.sqrt(intra_cov['TOT']['W2'].isel(
            sigmaWmin[:, k] = percent/(- z_r[k]) * (np.sqrt(intra_cov['TOT']['W2'].isel(
                level=[i for i in range(k, z_r.size)]))).integrate('level')

        # Compute sampling masks FROM W: Couvreux 2010, Pergaud 2009, Berg,Stull 2002 on SV

        m = 2
        masks['DW'] = (((variables['WT'] - MEAN_VARS['TOT']['WT']) < 0) & ((variables['WT'] - MEAN_VARS['TOT']
                                                                            ['WT']) < - m*np.maximum(sigmaWmin, np.sqrt(intra_cov['TOT']['W2']))))  # for down drafts
        masks['UP'] = ~masks['DW']

    if sampling == 'w_sieb_p01':
        # /!\ THIS SETS BY DEFINITION DW_p_FRAC = percentile !!!!
        #     #p-percentile velocity (van Ulden and Siebesma 1997, Siebesma et al 2007)

        p = 0.10
        wp = variables['WT'].quantile(p, dim=('ni', 'nj'))
        masks['DW'] = ((variables['WT']-wp) < 0)
        masks['UP'] = ~masks['DW']

    # Compute subdomain variables
    subdomain_variables('DW')
    subdomain_variables('UP')

    # # ------------------------------------------------------------
    # # Computing E,D
    # TODO, à adapter car ds il faut vérifier que ds c'est le bon dataset...
    # # ------------------------------------------------------------
    # dz = ds['level'].data[1] - ds['level'].data[0]


    # # interpolate on the level_w grid
    # UP_SVT_interp = MEAN_VARS['UP']['SVT001'].interp(level=ds.level_w).data
    # DW_SVT_interp = MEAN_VARS['UP']['SVT001'].interp(level=ds.level_w).data
    # DW_FRAC_interp = ds['DW_FRAC'].interp(level=ds.level_w).data
    # DW_WT_interp = ds['DW_WT'].interp(level=ds.level_w).data

    # E_minus_D = 1/dz * (ds['DW_FRAC'][:, 1:].data * ds['DW_WT'][:, 1:].data -
    #                     ds['DW_FRAC'][:, :-1].data * ds['DW_WT'][:, :-1].data)


    # Etilde = (DW_FRAC_interp[:, 1:-1] * DW_WT_interp[:, 1:-1]) / (UP_SVT_interp[:, 1:-1] -
    #                                                             DW_SVT_interp[:, 1:-1]) * 1/dz * (ds['DW_SVT001'][:, 1:].data - ds['DW_SVT001'][:, :-1].data)

    # # Etilde est le proxy pour E calculer avec traceur ou température
    # # Astuce pour éviter d'avoir D < 0:
    # # D = max (E-D, 0)  ,  E = max(Etilde, E-D) (car quand E-D > Etilde, on impose D=0 donc ça implique que E = E-D)

    # D = np.maximum(Etilde - E_minus_D, np.zeros_like(Etilde))
    # E = np.maximum(Etilde, E_minus_D)


#===========================================================================================

    new_MEAN_VARS = merger(MEAN_VARS)
    new_intra_cov = merger(intra_cov)
    new_inter_cov = merger(inter_cov, exclude_TOT=True)
    new_triple_cov = merger(triple_cov)

    dsout = xr.merge([dsout, new_MEAN_VARS, new_inter_cov,
                     new_intra_cov, new_triple_cov, mld_f])

    ds.close()
    print('closing '+file)
# WRITTING NETCDF FILE
#dsout = dsout.assign_coords(level_adim = level_adim)
dsout['level'] = z_r
dsout['level_w'] = z_w

dsout.to_netcdf(saving_path, mode='w', format="NETCDF4")
print('objects diagnostics saved at '+saving_path)












