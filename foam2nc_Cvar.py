#******************************************************************************
#MODULES
import csv
import fluidfoam
import os
import netCDF4 as nc
import numpy as np


#******************************************************************************
#FUNCTIONS
# Reading time list
def readOfTime(case_path):
    stream = os.popen('foamListTimes -case ' + case_path)
    time_list = stream.readlines()
    t = np.empty(len(time_list))
    t[:] = np.NaN
    for it in range(len(t)):
        t[it] = float(time_list[it])
    return t

# Obtaining time name
def t2str(t):
    if np.mod(t, 1) == 0:
        t_name = str(int(t))
    else:
        t_name = str(t)
    return t_name

# Reading netCDF variable table
def readNcTable():
    file_name = 'nc_variable_table.csv'
    fields = []
    rows = []
    with open(file_name, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
        for row in csvreader:
            rows.append(row)
    return fields, rows

# Finding OpenFOAM variable attributes
def findNcVar(of_name):
    fields, rows = readNcTable()
    col = fields.index('of_name') 
    of_name_list = []
    for row in range(len(rows)):
        of_name_list.append(rows[row][col])
    idx = of_name_list.index(of_name) 
    return rows[idx]

# Reading OpenFOAM data
def readOf(case_path, t, of_name):
    fields = readNcTable()[0]
    offields = findNcVar(of_name)
    nc_name = offields[fields.index('nc_name')]
    of_variable_type = offields[fields.index('of_variable_type')]
    readtype = 'read' + of_variable_type
    ofdata = eval("fluidfoam." + readtype + "(case_path, t2str(t), '" + of_name + "')")
    return ofdata


#******************************************************************************
#VARIABLES
# Inputs
case_name     = 'AL01_test09_9M'
of_name_list  = ['alpha.a']

# Directories
case_path = os.path.join('..', case_name)
nc_path = os.path.join(case_path, 'postProcessing', 'ncfiles')
os.makedirs(nc_path, exist_ok = True)


#******************************************************************************
#PROCESSING
# Close netCDF file if any
try:
    ncfile.close()
except:
    pass

# Start processing data
print('')
print('Case: ' + case_name)
t_new = readOfTime(case_path)    
fields = readNcTable()[0]
for of_name in of_name_list: # Loop over all assigned OpenFOAM variables
    print('')
    nc_attrib = findNcVar(of_name)
    of_variable_type = nc_attrib[fields.index('of_variable_type')]
    nc_name = nc_attrib[fields.index('nc_name')]
    ncfile_path = os.path.join(nc_path, nc_name + '.nc')
    print('Writing ' + ncfile_path + ':')
    
    # Create netCDF files if not exist
    if not os.path.exists(ncfile_path): 
        print('    ' + ncfile_path + ' does not exists, creating... ')
        C = nc.Dataset(os.path.join(nc_path, 'C.nc'), 'r')
        n = len(C.variables['Cx'])
        C.close()
        ncfile = nc.Dataset(ncfile_path, mode='w', format='NETCDF4')
        t_dim = ncfile.createDimension('t', None)
        n_dim = ncfile.createDimension('n', n)
        match of_variable_type:
            case "vector":
                idx_dim = ncfile.createDimension('idx', 3)
            case "symmTensor":
                idx_dim = ncfile.createDimension('idx', 6)
            case "tensor":
                idx_dim = ncfile.createDimension('idx', 9)
        dim_list = nc_attrib[fields.index('dimensions')].split( )        
        dim_str = ''
        for dim in dim_list:
            dim_str += ("'" + dim + "',")
        time = ncfile.createVariable('time', np.float64, ('t',))
        time.units = 's'
        time.long_name = 'time'
        exec(nc_name + " = ncfile.createVariable('" + nc_name + "', np.float64, (" + dim_str[0:-1] + "))")
        exec(nc_name + ".units = nc_attrib[fields.index('units')]")
        exec(nc_name + ".long_name = nc_attrib[fields.index('long_name')]")
        ncfile.close()
        print('    Done')

    # Attach new OpenFOAM data to netCDF files
    ncfile = nc.Dataset(ncfile_path, mode='a')
    t_old = ncfile.variables['time'][:]
    idx = 0
    for it in t_new:
        prog_perc = "%0.2f" % ((of_name_list.index(of_name)*len(t_new)+idx)/(len(of_name_list)*len(t_new))*100)
        if (idx >= len(t_old) and idx <= len(t_new)):
            print('[' + prog_perc + '%] Processing ' + nc_name + ' at t = ' + str(it) + ' sec')
            ncfile.variables['time'][idx] = it
            ofdata = readOf(case_path, it, of_name)
            if of_variable_type == 'scalar':
                exec("ncfile.variables['" + nc_name + "'][idx,:] = ofdata")
            else:
                exec("ncfile.variables['" + nc_name + "'][idx,:,:] = np.transpose(ofdata)")
        else:
            print('[' + prog_perc + '%] ' + nc_name + ' at t = ' + str(it) + ' sec exists, skip')
        idx += 1
    ncfile.close()
    
print('')
print('[100.00%] Done with writing ncfiles for ' + case_name)
