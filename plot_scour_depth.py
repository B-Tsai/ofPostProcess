#******************************************************************************
#MODULES
import os
import netCDF4 as nc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


#******************************************************************************
#VARIABLES
# Inputs
case_name = 'AL01_test09_9M'
D         = 0.04
# Test 6
#SD        = 0.45
#T         = 400
# Test 9
SD        = 0.28
T         = 90

# Directories
case_path = os.path.join('..',case_name)
post_path = os.path.join(case_path,'postProcessing')
os.makedirs(post_path, exist_ok = True)
fig_path = os.path.join(post_path,'figures')
os.makedirs(fig_path, exist_ok = True)
nc_path = os.path.join(post_path,'ncfiles')
os.makedirs(nc_path, exist_ok = True)


#******************************************************************************
#PROCESSING
# Read data
Z_nc = nc.Dataset(os.path.join(nc_path,'Z.nc'), 'r')
t = Z_nc.variables['time'][:]
Z = Z_nc.variables['Z'][:]

# Plot
fig = plt.figure(figsize=(20,5))
ax = fig.add_subplot(1,1,1)
cmap = plt.get_cmap('hsv', Z.shape[1])
for i in np.arange(Z.shape[1]):
    plt.plot(t, (Z[0,i]-Z[:,i])/D, linestyle='dashed', linewidth=1.0, color=cmap(i))
t_exp = np.arange(0.0, 3*T, 0.05)
plt.plot(t_exp, SD*(1-np.exp(-t_exp/T)), linestyle='solid', linewidth=1.5, color='grey')
plt.plot(t_exp, SD*np.ones(t_exp.shape), linestyle='dashed', linewidth=1.5, color='grey')
plt.plot(t, (np.nanmean(Z[0,:])-np.nanmean(Z, axis=1))/D, linestyle='solid', linewidth=1.5, color='black')
ax.set_xlabel('$\it{t}$ [sec]',fontsize=24)
ax.set_ylabel('$\it{S}$/$\it{D}$',fontsize=24)
ax.axis([0, 200, 0, 0.4])

# Save
case_no = case_name.split('_')[0]
filename = case_no + '_t_vs_SD.tif'
plt.savefig(os.path.join(fig_path,filename), dpi=300)

# Close netCDF file
Z_nc.close()
