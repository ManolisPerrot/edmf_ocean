import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoLocator
from scipy.signal import savgol_filter

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
plt.rcParams["axes.formatter.limits"] = [-3,3]
# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
    color=['#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'])


# CHOOSE CASE NAME

case = 'W005_C500_NO_COR'


# Parameters of the case for w scaling, ATTENTION if you change the case

Q0 = 500
alpha = 2.0 * 10 ** (-4)
g = 9.81
rho0 = 1024.
cp = 3900.

deltaT0 = 1/1000
T0 = 2 + 273.16  # surface init temp
f0 = 2*(1/(24*3600))*np.sin(2*np.pi/360 * 45)

B0 = g * alpha * Q0 / (rho0 * cp)
WT0 = Q0/(rho0 * cp)


saving_path = '../figures/'

# saving_name = case+'_object_diags_Cw_m1_all'
# saving_name = case+'_object_diags_Cw_m1_72h'
saving_name = case+'_object_diags_Cw_m05_72h'

# saving_name = case+'_object_diags_Cw_m2_72h'

path = '../data/'+case+'/'


file = saving_name+'.nc'


ds = xr.open_dataset(path+file)

time = ds.time
wstar = (-B0 * ds.MLD_FLUX)**(1/3)

instants = [-1]
instant = instants[0]


# CHOOSE number of snapshots you want to average


plt.rcParams.update({'font.size': 18})

# ------------------------------------------------------------
# FC500
# ------------------------------------------------------------

zlim = ds.MLD_FLUX - 50

zadim = ds.level/(-ds.MLD_FLUX[instant])

var_names = ['THT', 'WT', 'FRAC']
var_titles = [r'$\theta$', r'$w$', r'$a_p$']
units = [r'$\mathrm{K}$', r'$\mathrm{m\,s^{-1}}$', '']

masks = ['TOT', 'DW', 'UP']
style = {'TOT': '-k', 'DW': '-', 'UP': '-'}
linewidth = {'TOT': 4, 'DW': 2.5, 'UP': 2.5}
subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                 r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}']

variables = {}
for mask in masks:
    variables[mask] = [ds[mask+'_'+var] for var in var_names]


fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(13,11), sharey=True)

for i in range(len(var_names)-1):

    for instant in instants:
        for mask in masks:
            axs.flat[i].plot(variables[mask][i][instant, :], zadim,
                             style[mask], linewidth=linewidth[mask], label=mask)
    axs.flat[i].set_title(var_titles[i])
    axs.flat[i].grid(False)
    axs.flat[i].set_ylim(-1.25, 0)

    axs.flat[i].set_xlabel(units[i])

    if var_names[i] == 'THT':
        axs.flat[i].set_xlim(274.8, 274.91)
        for tick in axs.flat[i].get_xticklabels():
            tick.set_rotation(45)
# ----------------------
i = 2  # 'FRAC'
mask = 'DW'
axs.flat[i].plot(variables[mask][i][instant, :], zadim,
                 linestyle='-', linewidth=linewidth[mask], color='#0db4c3')
axs.flat[i].set_title(var_titles[i])
axs.flat[i].grid(False)
axs.flat[i].set_xlabel(units[i])

# -----------------------
i=1 # 'WT'
ax = axs.flat[i] 
mask = 'DW'
epsilon = 0.125

ax.plot(ds['DW_FRAC'][instant]
         * ds['DW_WT'][instant], zadim,
                 linestyle='-.', linewidth=linewidth[mask], color='#0db4c3', label=r'$a_p w_p$')
# wstar line
ax.axvline(x=-wstar[instant], alpha=1,
                    color='k', linestyle='--', label=r'$w_*$')
# ------------------------------------
ax = axs.flat[3]
ax.set_title(r'$\overline{w^\prime \theta^\prime}$')
ax.set_xlabel(r'$ \mathrm{m s^{-1} K}$')
ax.grid(False)


ax.plot(ds.TOT_intra_WTH[instant], zadim, style['TOT'], linewidth=linewidth['TOT'], label = r'$ \overline{w^\prime \theta^\prime}$')
ax.plot(ds.DW_inter_WTH[instant] , zadim, linestyle='--', color = blue  , linewidth=linewidth['DW'] , label=r'$a_p w_p (\theta_p - \overline{\theta})$')
ax.plot(ds.DW_intra_WTH[instant] , zadim, linestyle=':', color = blue , linewidth=linewidth['DW'] , label=r'$a_p \overline{w_p^\prime \theta_p^\prime}$')
ax.plot(ds.UP_inter_WTH[instant] , zadim, linestyle='--', color = orange , linewidth=linewidth['DW'] , label=r'$a_e w_e (\theta_e - \overline{\theta})$')
ax.plot(ds.UP_intra_WTH[instant] , zadim, linestyle=':', color = orange , linewidth=linewidth['DW'] , label=r'$a_e \overline{w_e^\prime \theta_e^\prime}$')

#---------------------------------------------------------------------------------------

# ----------------------------------------
ax = axs.flat[4]
ax.set_title(r'$1/2\overline{w^\prime \mathbf{u}^\prime\cdot\mathbf{u}^\prime}$')
ax.plot(ds.TOT_intra_WKE[instant], zadim, color='k', linestyle='-', linewidth=4,
        label=r'$\overline{w^\prime \frac{1}{2} \mathbf{u}^\prime \cdot \mathbf{u}^\prime  }$')

# I_e + III_e + III_p + IV_p =  \overline{w_e^\prime \frac{1}{2} \mathbf{u}_e^\prime \cdot \mathbf{u}_e^\prime } +  a_p w_p (k_p - k_e) + a_p w_p \frac{1}{2} (\mathbf{u}_p - \overline{\mathbf u})^2
# II_e + II_p = $ a_p(\mathbf u_p - \overline{\mathbf u}) \cdot \overline{\mathbf u_p^\prime w_p^\prime} + a_e(\mathbf u_e - \overline{\mathbf u}) \cdot \overline{\mathbf u_e^\prime w_e^\prime} $
# I_p = \overline{w_p^\prime \frac{1}{2} \mathbf{u}_p^\prime \cdot \mathbf{u}_p^\prime }
# IV_e = a_e w_e \frac{1}{2} (\mathbf{u}_e - \overline{\mathbf u})^2

# Param:
ax.plot((ds.UP_WKE1+ds.DW_WKE3+ds.UP_WKE3+ds.DW_WKE4)
        [instant], zadim, color=green, linestyle='-', linewidth=2.5, label=r'${\rm I}_e + {\rm III}_e + {\rm III}_p + {\rm IV}_p $')

factor = 1/(1-ds['DW_FRAC'])
ap = ds['DW_FRAC']
wp = ds['DW_WT']
kp = ds['DW_intra_TKE']/ap
kp_conv = ds['DW_inter_TKE']/ap
k = ds['TOT_intra_TKE']
mf = ds['DW_FRAC'] * ds['DW_WT']
factor2 = (1-ds['DW_FRAC']**2+ds['DW_FRAC'])/(1-ds['DW_FRAC'])**2
up = ds['DW_UT']
ubar = ds['TOT_UT']
wu_p = ds['DW_intra_WU']
ww_p = ds['DW_intra_WW']

ue = ds['UP_UT']
we = ds['UP_WT']
ubar = ds['TOT_UT']
wu_e = ds['UP_intra_WU']
ww_e = ds['UP_intra_WW']
ae = ds['UP_FRAC']

# # Negkected terms:
# ax.plot( (ap*(up-ubar)*wu_p)[instant], zadim, color='r',
#         linestyle='--', linewidth=2.5, alpha=0.5, label=r'ap (up - u) wup')
# ax.plot( (ap*(wp)*ww_p)[instant], zadim, color='blue',
#         linestyle='--', linewidth=2.5, alpha=0.5, label=r'ap (wp - w) wwp')

# ax.plot( (ae*(ue-ubar)*wu_e)[instant], zadim, color='r',
#         linestyle='-', linewidth=2.5, alpha=0.5, label=r'ae (ue - u) wue')
# ax.plot( (ae*(we)*ww_e)[instant], zadim, color='blue',
#         linestyle='--', linewidth=2.5, alpha=0.5, label=r'ae (we - w) wwe')

ax.plot((ds['DW_WKE2'] + ds['UP_WKE2'])[instant], zadim, color='k',
        linestyle='--', linewidth=2.5, alpha=0.5, label=r'${\rm II}_e + {\rm II}_p $')
ax.plot((ds['DW_WKE1'])[instant], zadim, color='k',
        linestyle='-.', linewidth=2.5, alpha=0.5, label=r'${\rm I}_p $')
ax.plot((ds['UP_WKE4'])[instant], zadim, color='k',
        linestyle=':', linewidth=2.5, alpha=0.5, label=r'${\rm IV}_e$')

ax.set_ylim(-1.2, 0)
ax.set_xlabel(r'$m^3 s^{-3}$')

# ----------------------------------------
axs.flat[5].set_visible(False)

#=================================================================================

# adding subplot labels
for i in range(len(axs.flat)-1):
    ax = axs.flat[i]
    ax.text(0.15, 0.98, subplot_label[i], transform=ax.transAxes,
            fontsize=16, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0), fontweight='bold', va='top', ha='right')

axs.flat[0].set_ylabel(r'$z/h$')


handles, labels = axs.flat[0].get_legend_handles_labels()
labels = [r'$\overline{\theta}$', r'$\theta_p$', r'$\theta_e$']


axs.flat[0].legend(handles, labels, loc='lower center',
           fancybox=False, shadow=False, ncol=1, bbox_to_anchor=(0.25,0.05))

handles, labels = axs.flat[1].get_legend_handles_labels()
labels = [r'$\overline{w}$', r'$w_p$', r'$w_e$',r'$a_p w_p$',r'$w_*$']


axs.flat[1].legend(handles, labels, loc='lower center',
           fancybox=False, shadow=False, ncol=1, bbox_to_anchor=(0.4,0.1))

axs.flat[3].legend(fontsize='14',fancybox=False, shadow=False, ncol=1, loc='lower center', bbox_to_anchor=(0.7,-0.03))

axs.flat[4].legend(loc='lower left', bbox_to_anchor=(0.98,0))

for ax in axs.flat:
    ax.set_box_aspect(1)
  


plt.savefig(saving_path+'LES_diags_W005_C500_Cwm05', bbox_inches='tight', dpi=600)
print('fig saved at '+saving_path+'LES_diags_W005_C500')

plt.show()
