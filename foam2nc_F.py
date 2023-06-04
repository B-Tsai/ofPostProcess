#******************************************************************************
#MODULES
import fluidfoam
import os
import netCDF4 as nc
import numpy as np


#******************************************************************************
#VARIABLES
# Inputs
case_name = 'AK03_test08_88k'
bnd_name  = 'bottom'

# Directories
case_path = os.path.join('..', case_name)
nc_path = os.path.join(case_path, 'postProcessing', 'ncfiles')
ncfile_path = os.path.join(nc_path,'F_' + bnd_name + '.nc')
os.makedirs(nc_path, exist_ok = True)


#******************************************************************************
#PROCESSING
# Close netCDF file if any
try:
    ncfile.close()
except:
    pass

# Read OpenFOAM mesh data
x, y, z = fluidfoam.readmesh(case_path, boundary=bnd_name)
n = len(x)

# Create netCDF files
print('')
print('Writing ' + ncfile_path)
ncfile = nc.Dataset(ncfile_path, mode='w', format='NETCDF4')
ncfile.title='Face centres of ' + bnd_name

# Create dimensions
n_dim = ncfile.createDimension('n', n)

# Create variables
Fx = ncfile.createVariable('Fx', np.float64, ('n',))
Fx.units = 'm'
Fx.long_name = 'face_centres_x_coordinate'
Fx[:] = x
Fy = ncfile.createVariable('Fy', np.float64, ('n',))
Fy.units = 'm'
Fy.long_name = 'face_centres_y_coordinate'
Fy[:] = y
Fz = ncfile.createVariable('Fz', np.float64, ('n',))
Fz.units = 'm'
Fz.long_name = 'face_centres_z_coordinate'
Fz[:] = z

# Close ncfile
ncfile.close()
print('Done')
print('')