#******************************************************************************
#MODULES
import os
import netCDF4 as nc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

#******************************************************************************
#VARIABLES
# Inputs
case_name = 'AY01_test08_0p7M'
D         = 0.04
T         = 4.4
N         = 5
dt        = 0.1

# Directories
case_path = os.path.join('..',case_name)
post_path = os.path.join(case_path,'postProcessing')
exp_path = os.path.join(post_path,'expdata')
nc_path = os.path.join(post_path,'ncfiles')
fig_path = os.path.join(post_path,'figures')
os.makedirs(fig_path, exist_ok = True)


#******************************************************************************
#PROCESSING
# Read experiment data
test_no   = int(case_name.split('_')[1][-2:])
if test_no == 8:
    S1997L_name = 'SumerEtAl1997_Fig17c_left.csv'
    S1997R_name = 'SumerEtAl1997_Fig17c_right.csv'
    B2017L_name = 'BaykalEtAl2017_Fig5a_left.csv'
    B2017R_name = 'BaykalEtAl2017_Fig5a_right.csv'
    J2021L_name = 'JangEtAl2021_Fig11b_left.csv'
    J2021R_name = 'JangEtAl2021_Fig11b_right.csv'
    tau0m_exp = 1000*0.013*0.013
    zero_cross_idx = int(2.3/dt)
    idxC = int(0.70/dt)
    idxT = int(3.80/dt)
#    zero_cross_idx = int(2.4/dt)
#    idxC = int(0.80/dt)
#    idxT = int(3.85/dt)
elif test_no == 14:
    tau0m_exp = 1000*0.019*0.019
    zero_cross_idx = 20
    idxC = int(0.55/dt)
    idxT = int(3.75/dt)
else:
    print("Wrong test no.")
S1997L = np.genfromtxt(os.path.join(exp_path,S1997L_name), delimiter=',')
S1997R = np.genfromtxt(os.path.join(exp_path,S1997R_name), delimiter=',')
B2017L = np.genfromtxt(os.path.join(exp_path,B2017L_name), delimiter=',')
B2017R = np.genfromtxt(os.path.join(exp_path,B2017R_name), delimiter=',')
J2021L = np.genfromtxt(os.path.join(exp_path,J2021L_name), delimiter=',')
J2021R = np.genfromtxt(os.path.join(exp_path,J2021R_name), delimiter=',')

# Read OpenFOAM data
taubx_nc = nc.Dataset(os.path.join(nc_path,'taubx.nc'), 'r')
x = taubx_nc.variables['x'][:]
t = taubx_nc.variables['time'][:]
taubx = taubx_nc.variables['taubx'][:]

# Process data
N1 = np.floor(np.round(t[-1]/T, decimals=3))
N0 = N1-N
idx_halfx = int((len(x))/2-1)
idx = np.arange(N0*T/dt, N1*T/dt, 1).astype(int)
taub = taubx.take(idx,axis=0)
taub_phase = np.empty((int(T/dt),taub.shape[1],N))
for k in np.arange(N):
    idx = (np.arange(0, T/dt, 1)+k*T/dt).astype(int)
    tmp = taub.take(idx,axis=0)
    taub_phase[:,:,k] = tmp
taub_crestphase = taub_phase.take(np.arange(zero_cross_idx),axis=0)
taub_troughphase = taub_phase.take(np.arange(zero_cross_idx,int(T/dt)),axis=0)
taub_crestphase_max = np.nanmax(np.absolute(taub_crestphase), axis=0)
taub_troughphase_max = np.nanmax(np.absolute(taub_troughphase), axis=0)
idx1 = x<-0.25
idx2 = x>0.25
tau0m_crest = np.nanmean((np.nanmean(taub_crestphase_max[idx1,:], axis=0)+np.nanmean(taub_crestphase_max[idx2,:], axis=0))/2)
tau0m_trough = np.nanmean((np.nanmean(taub_troughphase_max[idx1,:], axis=0)+np.nanmean(taub_troughphase_max[idx2,:], axis=0))/2)
tau0m = np.max(np.array([tau0m_crest, tau0m_trough]))

# Plot
plt.figure(figsize=(9,5))
sp = plt.subplot(1, 1, 1)
plt.plot(S1997L[:,0], S1997L[:,1], 'ro', mfc='none', markersize=5, zorder=4)
plt.plot(S1997R[:,0], S1997R[:,1], 'ro', mfc='none', markersize=5, zorder=4)
plt.plot(B2017L[:,0], B2017L[:,1], 'g--', zorder=2)
plt.plot(B2017R[:,0], B2017R[:,1], 'g--', zorder=2)
plt.plot(J2021L[:,0], J2021L[:,1], 'b:', zorder=3)
plt.plot(J2021R[:,0], J2021R[:,1], 'b:', zorder=3)
for k in np.arange(N):
    plt.plot(x[0:idx_halfx]/D, taub_phase[idxC,0:idx_halfx,k]/tau0m,'-', color='lightgrey', zorder=1)
    plt.plot(x[(idx_halfx+1):]/D, -taub_phase[idxT,(idx_halfx+1):,k]/tau0m,'-', color='lightgrey', zorder=1)
plt.plot(x[0:idx_halfx]/D, np.nanmean(taub_phase[idxC,0:idx_halfx,:], axis=1)/tau0m,'k-', zorder=5)
plt.plot(x[(idx_halfx+1):]/D, -np.nanmean(taub_phase[idxT,(idx_halfx+1):,:], axis=1)/tau0m,'k-', zorder=5)
rect = plt.Rectangle((-0.5,-2.5), 1, 5, color='gray', alpha=1, zorder=6)
sp.add_patch(rect)
sp.axis([-3, 3, -2.0, 1.0])
sp.set_xlabel(r'$x/D$'+' [-]',fontsize=24)
sp.set_ylabel(r'$\tau_b \, / \, \tau_{bm}$'+' [-]',fontsize=24)
sp.set_axisbelow(True)
sp.xaxis.grid(color='whitesmoke')
sp.yaxis.grid(color='whitesmoke')

# Save
case_no = case_name.split('_')[0]
filename = case_no + '_bedShearStress.tif'
plt.savefig(os.path.join(fig_path,filename), dpi=300)

# Close netCDF file
taubx_nc.close()
