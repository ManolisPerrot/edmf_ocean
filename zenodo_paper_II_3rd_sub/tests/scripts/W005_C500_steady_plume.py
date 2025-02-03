import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoLocator
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor':'white'})
plt.rcParams.update({'savefig.facecolor':'white'})

##colors
blue, orange, magenta, grey, green ='#0db4c3','#eea021','#ff0364','#606172','#3fb532'
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#0db4c3','#eea021','#ff0364','#606172','#3fb532']) 


### CHOOSE CASE NAME

case = 'W005_C500_NO_COR'

#### Parameters of the case for w scaling, ATTENTION if you change the case

Q0 = 500
alpha = 2.0 * 10 ** (-4)
g = 9.81
rho0 = 1024.
cp = 3900.

deltaT0 = 1/1000 # 
T0 = 2 + 273.16 #surface init temp
f0 = 2*(1/(24*3600))*np.sin(2*np.pi/360 * 45)
N0 = np.sqrt(alpha*g*(deltaT0))


B0 = g * alpha * Q0 / (rho0 * cp)
WT0 = Q0/(rho0 * cp)
L0 = (B0/N0**3)**(1/2)


saving_path = '../figures/'


path = '../data/'+case+'/'
file = 'W005_C500_NO_COR_object_diags_Cw_m1_all.nc'

ds = xr.open_dataset(path+file)

time = ds.time
wstar = (-B0 * ds.MLD_FLUX)**(1/3)




### CHOOSE which instants you want to visualize, whith a method (start-stop/manual/evenly space)

#start-end
start = 0 #ATTENTION this a backup number, has to be multiplied by backup spacing to get hours!
stop   = 36

instants =  [i for i in range(start,stop)]



###CHOOSE number of snapshots you want to average

nb_avg = len(instants)

def time_avg(variable):
    return variable[stop-nb_avg:stop,:].sum('time')/(nb_avg)


tt = (time - time[0])/np.timedelta64(1,'h')+2 #array of hours
z_enc =  np.sqrt(2) * (WT0/(deltaT0) * tt*3600)**(1/2)

## /!\ differentiate is working for 'level', but for 'time' the result is multiplied by ~ 10-9 (because time is in ns?)
#so we use handmade differentiation...

def d_z (X):
    dz = ds['level'][1].data - ds['level'][0].data
    return (X.data[1:, 1:] - X.data[1:,:-1]) / dz

def d_t (X):
    dt = (tt.data[1] - tt.data[0])*3600
    return (X.data[ 1:,1:] - X.data[:-1,1:]) / dt

###########################

tendTH = np.abs(d_t(ds['DW_FRAC']*ds['DW_THT']))
advTH  = np.abs(d_z(ds['DW_FRAC']*ds['DW_WT']*ds['DW_THT']))
fluxTH = np.abs(d_z(ds['DW_intra_WTH']))


tendW = np.abs(d_t(ds['DW_FRAC']*ds['DW_WT']))
advW  = np.abs(d_z(ds['DW_FRAC']*ds['DW_WT']*ds['DW_WT']) )
fluxW = np.abs(d_z(ds['DW_intra_WW']) )

tendU = np.abs(d_t(ds['DW_FRAC']*ds['DW_UT']))
advU  = np.abs(d_z(ds['DW_FRAC']*ds['DW_WT']*ds['DW_UT']) )
fluxU = np.abs(d_z(ds['DW_intra_WU']) )

#------------------------------------------------------------
#FC500
#------------------------------------------------------------



time_mesh, z_mesh = np.meshgrid(tt[1:], ds.level[1:])


fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15,4),sharey=True, layout='constrained')



def t_to_z_over_L0 (t):
    return 1/L0 * np.sqrt(2) * (WT0/(deltaT0) * t*3600)**(1/2)

def z_over_L0_to_t(z_over_L0):
    return 1/2 * (deltaT0)/WT0 * 1/3600 * (z_over_L0 * L0)**2

#------------------------------------------------------------

# Set the levels for contourf, using a logarithmic scale
levels = np.logspace(np.log10(1e-5), np.log10(1e-0), 16)

cmap = 'Greys'
norm = colors.LogNorm()
# contour_colors=[blue,'silver',orange,magenta]
contour_colors=['tab:green','silver','tab:pink','tab:orange']
ax = axs.flat[0]

#------------------------------------------------------------

im = ax.contourf(time_mesh, z_mesh, (tendTH/ ((T0)*N0)).T, levels=levels, norm=norm,cmap=cmap )
CS2 = ax.contour(im, levels=im.levels[::3],colors=contour_colors)


ax.set_title(r'$ |\partial_t (a_p \theta_p)| / (\theta_0 N_0) $')
ax.set_ylabel(r'$z (m)$')
ax.set_xlabel(r'\textrm{time} ($h$)')

ax.set_ylim(-400,0)

#------------------------------------------------------------

ax = axs.flat[1]

im = ax.contourf(time_mesh, z_mesh, (advTH/ ((L0*N0)*(T0)/L0) ).T, levels=levels, norm=norm,cmap=cmap  )
CS2 = ax.contour(im, levels=im.levels[::3],colors=contour_colors)

ax.set_title(r'$|\partial_z(a_p w_p \theta_p) | / ( \theta_0 N_0 ) $')
ax.set_xlabel(r'\textrm{time} ($h$)')

ax.set_ylim(-400,0)


#------------------------------------------------------------

ax = axs.flat[2]

im = ax.contourf(time_mesh, z_mesh, (tendW/ ((L0*N0)*N0)).T,levels=levels, norm=norm,cmap=cmap  )
CS2 = ax.contour(im, levels=im.levels[::3],colors=contour_colors)

ax.set_title(r'$ |\partial_t (a_p w_p)| / ( L_0 N_0^2) $')
ax.set_xlabel(r'\textrm{time} ($h$)')

ax.set_ylim(-400,0)

#------------------------------------------------------------

ax = axs.flat[3]

im = ax.contourf(time_mesh, z_mesh, (advW/ ((L0*N0)**2/L0) ).T,10, levels=levels, norm=norm,cmap=cmap  )
CS2 = ax.contour(im, levels=im.levels[::3],colors=contour_colors)

ax.set_title(r'$|\partial_z(a_p w_p^2) | / ( L_0 N_0^2 ) $')
ax.set_xlabel(r'\textrm{time} ($h$)')

ax.set_ylim(-400,0)


 

cbar=fig.colorbar(im, orientation="vertical" , ticks=levels)
cbar.add_lines(CS2)



#adding subplot labels
subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}', r'\rm{(d)}',r'\rm{(e)}',r'\rm{(f)}',r'\rm{(g)}']
for i in range(axs.size):
    ax = axs.flat[i]
    ax.set_box_aspect(1)
    ax.text(0.15, 0.98, subplot_label[i], transform=ax.transAxes,
      fontsize=16, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0, alpha=0.5), fontweight='bold', va='top', ha='right')


plt.savefig(saving_path+'W005_C500_steady_plume.png',dpi=300)


#plt.show()
