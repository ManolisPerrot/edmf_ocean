import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoLocator
from scipy.signal import savgol_filter
import subprocess
plt.rcParams['font.family'] = 'serif'
plt.rcParams['text.usetex'] = True
# plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
plt.rcParams["axes.formatter.limits"] = [-3,3]
# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
    color=['#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'])


# CHOOSE CASE NAME

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

cases = ['WANG1_NR_new','WANG1_FR_lat60']
# cases = ['WANG1_FR_lat60']
rerun_case=False
save_plot=True
saving_path='../figures/WANG1_NR_FR_les_scm_condsamp.pdf'
# rerun_case=False
if rerun_case:
    subprocess.run(["python", "WANG1_NR_FR306090_scm_les.py"])
#== Opening LES and SCM
les,scm = {},{}
for case in cases:
    file = case+'_object_diags_Cw_m1_276h.nc'
    path = '../data/'+case+'/'
    les[case] = xr.open_dataset(path+file)
    scm[case] = xr.open_dataset('scm_'+case+'.nc')
#----------------------------------------
time = les[cases[0]].time
wstar = (-B0 * les[cases[0]].MLD_FLUX)**(1/3)
instant = -1

MF_div = {}
for case in cases:
    MF_div[case] = (les[case]['DW_FRAC']*les[case]['DW_WT']*les[case]['DW_THT']).diff('level')

def E_D_LES(les):
# ------------------------------------------------------------
# Computing E,D on LES from passive tracer SVT001
# ------------------------------------------------------------
# $$ E\_ m\_ D = \partial_z (a_p w_p)$$
# $$\tilde{E} = \frac{a_p w_p}{\phi_e - \phi_p} \partial_z \phi_p
# $$

    dz = les['level'].data[1] - les['level'].data[0]
    regul = -0
    # interpolate on the level_w grid
    UP_SVT_interp = les['UP_SVT001'].interp(level=les.level_w).data
    DW_SVT_interp = les['DW_SVT001'].interp(level=les.level_w).data
    DW_FRAC_interp = les['DW_FRAC'].interp(level=les.level_w).data
    DW_WT_interp = les['DW_WT'].interp(level=les.level_w).data

    E_minus_D = 1/dz * (les['DW_FRAC'][:, 1:].data * les['DW_WT'][:, 1:].data -
                        les['DW_FRAC'][:, :-1].data * les['DW_WT'][:, :-1].data)

    Etilde = (DW_FRAC_interp[:, 1:-1] * DW_WT_interp[:, 1:-1]) / (UP_SVT_interp[:, 1:-1] -
                                                                DW_SVT_interp[:, 1:-1]) * 1/dz * (les['DW_SVT001'][:, 1:].data - les['DW_SVT001'][:, :-1].data)

    # Etilde est le proxy pour E calculer avec traceur ou température
    # Astuce pour éviter d'avoir D < 0:
    # D = max (Etilde - (E-D), 0)  
    # Ce qui imlique que E = max(Etilde, E-D) (car quand E-D > Etilde, on impose D=0 donc ça implique que E = E-D)

    D = np.maximum(Etilde - E_minus_D, np.zeros_like(Etilde))
    E = np.maximum(Etilde, E_minus_D)

    return E,D

E_les = { case:E_D_LES(les[case])[0] for case in cases }
D_les = { case:E_D_LES(les[case])[1] for case in cases }

zlim = les[cases[0]].MLD_FLUX - 50

zadim = les[cases[0]].level/(-les[cases[0]].MLD_FLUX[instant])
mld = (-les[cases[0]].MLD_FLUX[instant])

masks = ['TOT', 'DW']
color =     {'TOT': 'k','DW': 'k','UP': orange}
linewidth = {'TOT': 2,  'DW': 2, 'UP': 2}
linestyle = {cases[0]: ':', cases[1]: ':'}
linestyle_scm = {cases[0]: '-', cases[1]: '-'}

alpha = {cases[0]: 0.3, cases[1]: 1}


subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                 r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}', r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}']


fig, axs = plt.subplots(nrows=2, ncols=3, sharey=True,constrained_layout=True)
i_ax = -1
#-----------------------------------------------------------------------
i_ax+=1
ax=axs.flat[i_ax]
ax.set_title(r'$\theta_p$')
ax.set_xlabel(r'$\mathrm{^oC}$')
ax.set_ylabel(r'$z/h$')

for case in cases:
    mask='DW'
    ax.plot(les[case][mask+'_THT'][instant]-273.15 , zadim, color[mask], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case],label=case+'LES')
    #
    ax.plot(scm[case]['temp_p'][instant] , scm[case]['z_w']/mld, color[mask], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+'SCM')
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(loc='lower left')

# labels=[r'$\overline{\theta}$ (NR)',r'$\theta_p$ (NR)', r'$\overline{\theta}$ (FR)',r'$\theta_p$ (NR)']
# ax.legend(handles,labels,fancybox=False)
# #-----------------------------------------------------------------------
i_ax+=1
ax=axs.flat[i_ax]
ax.set_title(r'$w_p$')
ax.set_xlabel(r'$\mathrm{m s^{-1}}$')

for case in cases:
    mask='DW'
    ax.plot(les[case][mask+'_WT'][instant] , zadim, color[mask], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case],label=case)
    ax.plot(scm[case]['w_p'][instant] , scm[case]['z_w']/mld, color[mask], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
#-----------------------------------------------------------------------
i_ax+=1
ax=axs.flat[i_ax]
ax.set_title(r'$a_p$')

for case in cases:
    mask='DW'
    ax.plot(les[case][mask+'_FRAC'][instant] , zadim, color[mask], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case])
    ax.plot(scm[case]['a_p'][instant] , scm[case]['z_w']/mld, color[mask], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
# # #-----------------------------------------------------------------------
# i_ax+=1
# ax=axs.flat[i_ax]
# ax.set_title(r'$\partial_z(a_p w_p \theta_p)$')
# ax.set_xlabel(r'$\mathrm{K s^{-1}}$')

# for case in cases:
#     mask='DW'
#     ax.plot(MF_div[case][instant] , zadim[1:], color[mask], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case])

# ax.set_xlim(-0.06,0.2)
# # #-----------------------------------------------------------------------

# #-----------------------------------------------------------------------
i_ax+=1
ax=axs.flat[i_ax]
ax.set_title(r'$\omega_p$')
ax.set_xlabel(r'$\mathrm{s^{-1}}$')
ax.set_ylabel(r'$z/h$')

for case in cases:
    mask='DW'
    ax.plot(les[case][mask+'_VORT_z'][instant] , zadim, color[mask], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case])
    ax.plot(scm[cases[1]]['vort_p'][instant] , scm[case]['z_w']/mld, color[mask], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
i_ax+=1
ax=axs.flat[i_ax]
ax.set_title(r'$E$')
ax.set_xlabel(r'$\mathrm{s^{-1}}$')

for case in cases:
    mask='DW'
    # ax.plot(E_les[case][instant]-D_les[case][instant] , zadim[1:], color[mask], linewidth=linewidth[mask], linestyle='-',alpha=alpha[case])
    # ax.plot(scm[case]['Ent'][instant]-scm[case]['Det'][instant] , scm[case]['z_r']/mld, color[mask], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    ax.plot(E_les[case][instant] , zadim[1:], color[mask], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case])
    ax.plot(scm[case]['Ent'][instant] , scm[case]['z_r']/mld, color[mask], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
# ax.set_xlim(-1e-5,1e-4)
ax.set_xscale('log')
# #-----------------------------------------------------------------------
# #-----------------------------------------------------------------------
i_ax+=1
ax=axs.flat[i_ax]
ax.set_title(r'$D$')
ax.set_xlabel(r'$\mathrm{s^{-1}}$')

for case in cases:
    mask='DW'
    ax.plot(D_les[case][instant] , zadim[1:], color[mask], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case])
    ax.plot(scm[case]['Det'][instant] , scm[case]['z_r']/mld, color[mask], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
ax.set_xscale('log')

# # ax.set_xlim(-1e-5,1e-4)
# ax.set_xscale('log')
# i_ax+=1
# ax=axs.flat[i_ax]
# ax.remove()
# # # #-----------------------------------------------------------------------
# i_ax+=1
# ax=axs.flat[i_ax]
# ax.set_title(r'$D$')
# ax.set_xlabel(r'$\mathrm{s^{-1}}$')

# for case in cases:
#     mask='DW'
#     ax.plot(D_les[case][instant] , zadim[1:], color[mask], linewidth=linewidth[mask], linestyle='--',alpha=alpha[case])
# ax.set_xlim(-3e-6,5e-5)

subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                    r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}',r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}',r'\rm{(j)}',r'\rm{(k)}',r'\rm{(l)}']

for i,ax in enumerate(axs.flat):
    # ax.set_box_aspect(1)
    ax.grid(alpha=0.5)
    ax.set_ylim(-1.2,0)
    # adding subplot labels
    ax.text(0.15, 0.97, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none',alpha=0.5), fontweight='bold', va='top', ha='right')




handles, labels = axs.flat[0].get_legend_handles_labels()
labels = [r'NR, LES',r'NR, SCM',r'FR90, LES',r'FR90, SCM']

# fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.01, 0.8),
#            fancybox=False, shadow=False, ncols=2)
axs.flat[0].legend(handles, labels,fancybox=False)
fig.suptitle(r'Plume variables')
# fig.legend()
if save_plot:
    plt.savefig(saving_path, bbox_inches='tight', dpi=600)
    print('figure saved at ' + saving_path)
plt.show()
