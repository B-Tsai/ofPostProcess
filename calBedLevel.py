#******************************************************************************
#MODULES
import fluidfoam
import math
import os
import netCDF4 as nc
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d


#******************************************************************************
#FUNCTIONS
# Give the index of the interpolation area
def interpIdx(x, y, z, xi, yi, zi, dx, dy, dz, ngrid):
    xidx = ((np.min(xi)-ngrid*dx) < x) & ((np.max(xi)+ngrid*dx) > x)
    yidx = ((np.min(yi)-ngrid*dy) < y) & ((np.max(yi)+ngrid*dy) > y)
    zidx = ((np.min(zi)-ngrid*dz) < z) & ((np.max(zi)+ngrid*dz) > z)
    idx = xidx & yidx & zidx
    return idx

# Extract the interpolation area
def extractInterp(x, y, z, v, idx):
    xi = x[idx]
    yi = y[idx]
    zi = z[idx]
    vi = v[idx]
    return xi, yi, zi, vi


#******************************************************************************
#VARIABLES
# Inputs
case_name = 'AL01_test09_9M'
dx        = 0.0003*1
dy        = dx
dz        = dx
alpha_min = 0.57
D         = 0.04
xi        = (D/2 + dx)*np.cos(np.linspace(0,2*math.pi,17))
yi        = (D/2 + dx)*np.sin(np.linspace(0,2*math.pi,17))
zi        = np.arange(-D/2, D/2, dz)
ngrid     = 3

# Directories
case_path = os.path.join('..',case_name)
nc_path = os.path.join(case_path, 'postProcessing', 'ncfiles')
os.makedirs(nc_path, exist_ok = True)


#******************************************************************************
#PROCESSING
# Read data
C_path = os.path.join(nc_path,'C.nc')
alphaa_path = os.path.join(nc_path,'alphaa.nc')
C = nc.Dataset(C_path,'r')
alphaa_nc = nc.Dataset(alphaa_path,'r')
x = C.variables['Cx']
y = C.variables['Cy']
z = C.variables['Cz']
t_new = alphaa_nc.variables['time'][:]
alphaa = alphaa_nc.variables['alphaa'][:]

# Start processing data
print('')
print('Case: ' + case_name)
ncfile_path = os.path.join(nc_path, 'Z.nc')

# Create netCDF files if not exist
if not os.path.exists(ncfile_path):
    print('    ' + ncfile_path + ' does not exist, creating... ')
    ncfile = nc.Dataset(ncfile_path, mode='w', format='NETCDF4')
    t_dim = ncfile.createDimension('t', None)
    n_dim = ncfile.createDimension('n', len(xi)-1)
    time = ncfile.createVariable('time', np.float64, ('t',))
    time.units = 's'
    time.long_name = 'time'
    Z = ncfile.createVariable('Z', np.float64, ('t', 'n'))
    Z.units = 'm'
    Z.long_name = 'Bed_level'
    ncfile.close()
    print('    Done')

# Attach new OpenFOAM data to netCDF files
ncfile = nc.Dataset(ncfile_path, mode='a')
t_old = ncfile.variables['time'][:]
count = 0
for it in t_new:
    prog_perc = "%0.2f" % (count/len(t_new)*100)
    if (count >= len(t_old) and count <= len(t_new)):
        print('[' + prog_perc + '%] Processing Z at t = ' + str(it) + ' sec')
        ncfile.variables['time'][count] = it
        v = alphaa[count,:]
        Zj = np.empty(len(xi)-1)
        Zj[:] = np.NaN
        for j in np.arange(len(xi)-1):
            idx = interpIdx(x, y, z, xi[j], yi[j], zi, dx, dy, dz, ngrid)
            xv, yv, zv, vv = extractInterp(x, y, z, v, idx)
            xp = xi[j]*np.ones(zi.shape)
            yp = yi[j]*np.ones(zi.shape)
            zp = zi
            xyzp = np.vstack((xp.flatten(), yp.flatten(), zp.flatten())).T
            interp3 = LinearNDInterpolator((xv, yv, zv), vv)
            vp = interp3(xyzp)
            idx_z = np.argmin(np.abs(vp-alpha_min))
            Zj[j] = zp[idx_z]
        ncfile.variables['Z'][count,:] = Zj
    else:
        print('[' + prog_perc + '%] Z at t = ' + str(it) + ' sec exists, skip')
    count += 1
print('[100.00%] Done with writing Z.nc for ' + case_name)

# Close netCDF file
C.close()
alphaa_nc.close()
ncfile.close()
