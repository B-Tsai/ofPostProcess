#******************************************************************************
#MODULES
import math
import os
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from scipy.interpolate import LinearNDInterpolator


#******************************************************************************
#FUNCTIONS
# Give the index of the interpolation area
def interpIdx(x, y, xi, yi, dx, dy, ngrid):
    xidx = ((np.min(xi)-ngrid*dx) < x) & ((np.max(xi)+ngrid*dx) > x)
    yidx = ((np.min(yi)-ngrid*dy) < y) & ((np.max(yi)+ngrid*dy) > y)
    idx = xidx & yidx
    return idx

# Extract the interpolation area
def extractInterp(x, y, v, idx):
    xi = x[idx]
    yi = y[idx]
    vi = v[idx]
    return xi, yi, vi


#******************************************************************************
#VARIABLES
# Inputs
case_name = 'AY01_test08_0p7M'
bnd       = 'bottom'
rhof      = 1000
D         = 0.04
dx        = 0.0003*2
dy        = dx
x1        = np.arange(-0.3+dx/2, -D/2-dx/2, dx)
x2        = np.arange(D/2+dx/2, 0.3-dx/2, dx)
xi        = np.concatenate((x1, x2), axis=0)
yi        = 0
ngrid     = 3

# Directories
case_path = os.path.join('..', case_name)
nc_path = os.path.join(case_path, 'postProcessing', 'ncfiles')
fig_path = os.path.join(case_path, 'postProcessing','figures')
os.makedirs(nc_path, exist_ok = True)
os.makedirs(fig_path, exist_ok = True)


#******************************************************************************
#PROCESSING
# Read data
F_path = os.path.join(nc_path,'F_' + bnd + '.nc')
wallShearStress_path = os.path.join(nc_path,'wallShearStress_' + bnd + '.nc')
F = nc.Dataset(F_path, 'r')
wallShearStress_nc = nc.Dataset(wallShearStress_path, 'r')
Fx = F.variables['Fx'][:]
Fy = F.variables['Fy'][:]
Fz = F.variables['Fz'][:]
t_new = wallShearStress_nc.variables['time'][:]
wallShearStress = wallShearStress_nc.variables['wallShearStress'][:]

# Start processing data
print('')
print('Case: ' + case_name)
ncfile_path = os.path.join(nc_path, 'taubx.nc')

# Create netCDF files if not exist
if not os.path.exists(ncfile_path):
    print('    ' + ncfile_path + ' does not exist, creating... ')
    ncfile = nc.Dataset(ncfile_path, mode='w', format='NETCDF4')
    t_dim = ncfile.createDimension('t', None)
    x_dim = ncfile.createDimension('x', len(xi))
    time = ncfile.createVariable('time', np.float64, ('t',))
    time.units = 's'
    time.long_name = 'time'
    x = ncfile.createVariable('x', np.float64, ('x',))
    x.units = 'm'
    x.long_name = 'x_coordinate'
    x[:] = xi
    taubx = ncfile.createVariable('taubx', np.float64, ('t', 'x'))
    taubx.units = 'kg m-1 s-2'
    taubx.long_name = 'bed_shear_stress_along_x_coordinate'
    ncfile.close()
    print('    Done')

# Attach new OpenFOAM data to netCDF files
ncfile = nc.Dataset(ncfile_path, mode='a')
t_old = ncfile.variables['time'][:]
count = 0
for it in t_new:
    prog_perc = "%0.2f" % (count/len(t_new)*100)
    if (count >= len(t_old) and count <= len(t_new)):
        print('[' + prog_perc + '%] Processing taubx at t = ' + str(it) + ' sec')
        ncfile.variables['time'][count] = it
        v = wallShearStress[count,:,0]
        idx = interpIdx(Fx, Fy, xi, yi, dx, dy, ngrid)
        xv, yv, vv = extractInterp(Fx, Fy, v, idx)
        xp = xi
        yp = yi*np.ones(xi.shape)
        xyp = np.vstack((xp.flatten(), yp.flatten())).T
        interp = LinearNDInterpolator((xv, yv), vv)
        vp = interp(xyp)
        print(vp)
        ncfile.variables['taubx'][count,:] = -rhof*vp
    else:
        print('[' + prog_perc + '%] taubx at t = ' + str(it) + ' sec exists, skip')
    count += 1
print('[100.00%] Done with writing taubx.nc for ' + case_name)
print('')

# Close netCDF file
F.close()
wallShearStress_nc.close()
ncfile.close()
