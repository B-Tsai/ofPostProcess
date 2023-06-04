#******************************************************************************
#MODULES
import fluidfoam
import os
import netCDF4 as nc
import numpy as np


#******************************************************************************
#VARIABLES
# Inputs
case_name = 'AT03_test09_134k'

# Directories
case_path = os.path.join('..', case_name)
nc_path = os.path.join(case_path, 'postProcessing', 'ncfiles')
ncfile_path = os.path.join(nc_path, 'C.nc') 
os.makedirs(nc_path, exist_ok = True)


#******************************************************************************
#PROCESSING
# Close netCDF file if any
try: 
    ncfile.close()
except: 
    pass

# Read OpenFOAM mesh data
x, y, z = fluidfoam.readmesh(case_path)
n = len(x)

# Create netCDF files
print('')
print('Writing ' + ncfile_path)
ncfile = nc.Dataset(ncfile_path, mode='w', format='NETCDF4')
ncfile.title='Cell centres'

# Create dimensions
n_dim = ncfile.createDimension('n', n)

# Create variables
Cx = ncfile.createVariable('Cx', np.float64, ('n',))
Cx.units = 'm'
Cx.long_name = 'cell_centres_x_coordinate'
Cx[:] = x
Cy = ncfile.createVariable('Cy', np.float64, ('n',))
Cy.units = 'm'
Cy.long_name = 'cell_centres_y_coordinate'
Cy[:] = y
Cz = ncfile.createVariable('Cz', np.float64, ('n',))
Cz.units = 'm'
Cz.long_name = 'cell_centres_z_coordinate'
Cz[:] = z

# Close ncfile
ncfile.close()
print('Done')
print('')