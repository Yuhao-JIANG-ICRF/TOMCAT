import numpy as np

'''
  'Device name': (minor radius; major radius; magnetic axis)
'''
parameter_dict = {    
    'CFEDR': (2.5, 7.8, 8.2),
    'CFETR': (2.5, 7.8, 8.2),
    'WEST': (0.43, 2.5, 2.5+0.05)
    # add new device here
}

species_dict = {
    'Li' : '^7Li',
    'He3': '^3He',
    'He4': '^4He'
}

#%%
def get_device_input(name):
    global N0, N1, ap, R0, B0, xkap, del0, delt, T0e, T0i, T1, exponn, expont
    global  RRR0, ZZZ0, Rsep0, Zsep0, rho0, Nprf0, TprfE0, TprfI0, FWR, FWZ
    FWR = []
    FWZ = []
    name = name.lower()
    if name == 'west':
        N0 = 6.02e19
        N1 = 2.25e19
        ap = 0.43
        R0 = 2.5
        B0 = 3.657
        xkap = 1.4
        del0 = 0.05
        delt = 0.3 
        T0e = 1.609e3
        T0i = 1.37e3
        T1 = 0.25e3
        exponn = 0.3
        expont = 1.0
        
        pre = 'input/WEST/input_file/'
        import netCDF4    
        global  RRR1, ZZZ1, Rsep1, Zsep1, rho1, Nprf1, TprfE1, TprfI1
        # - Reference - #
        d0 = netCDF4.Dataset(pre+'west_59141_t6000ms_exp.inp.nc')
        # - Analytical from eve- #
        d1 = netCDF4.Dataset(pre+'west_55605_analytical.inp.nc')

        nsgpls0 = d0.dimensions['nsgpls'].size
        RR0 = d0.variables['sRR'][:, :, 0]
        ZZ0 = d0.variables['sZZ'][:, :, 0]
        sg0 = d0.variables['sg'][:nsgpls0]
        FT0 = d0.variables['sFT'][:nsgpls0, 0]
        psi0 = d0.variables['spsi'][:nsgpls0, 0]
        ne0 = d0.variables['sns'][0, :nsgpls0, 0]
        Te0 = d0.variables['sTs'][0, :nsgpls0, 0]
        Ti0 = d0.variables['sTs'][1, :nsgpls0, 0]

        nsgpls1 = d1.dimensions['nsgpls'].size
        RR1 = d1.variables['sRR'][:, :, 0]
        ZZ1 = d1.variables['sZZ'][:, :, 0]
        sg1 = d1.variables['sg'][:nsgpls1]
        FT1 = d1.variables['sFT'][:nsgpls1, 0]
        psi1 = d1.variables['spsi'][:nsgpls1, 0]
        ne1 = d1.variables['sns'][0, :nsgpls1, 0]
        Te1 = d1.variables['sTs'][0, :nsgpls1, 0]
        Ti1 = d1.variables['sTs'][1, :nsgpls1, 0]
        
        RRR0 = RR0[nsgpls0, :]
        ZZZ0 = ZZ0[nsgpls0, :]
        RRR1 = RR1[nsgpls1, :]
        ZZZ1 = ZZ1[nsgpls1, :]  
        
        Rsep0 = RR0[0,0] 
        Zsep0 = ZZ0[0,0] 
        Rsep1 = RR1[0,0] 
        Zsep1 = ZZ1[0,0] 
        
        
        rho0 = sg0
        Nprf0 = ne0
        rho1 = sg1
        Nprf1 = ne1
        TprfE0 = Te0
        TprfE1 = Te1
        TprfI0 = Ti0
        TprfI1 = Ti1
      
        
    elif name == 'cfedr':
        N0 = 14e19
        N1 = 4e19
        ap = 2.5
        R0 = 7.8
        B0 = 6.3
        xkap = 1.8
        del0 = 0.40
        delt = 0.25
        T0e = 30e3
        T0i = 24e3
        T1 = 0.8e3
        exponn = 0.7
        expont = 1.3
        
        pre = 'input/CFEDR/CFEDR_Conventional_H_mode_V1_20240522/'
        from scipy.io import loadmat
        from read_gfile_func import read_gfile_func

        # Load MAT file with better handling of MATLAB structs
        pre = 'input/CFEDR/CFEDR_Conventional_H_mode_V1_20240522/'
        data = loadmat(pre+'CFETR_PROFILES.mat', squeeze_me=True, struct_as_record=False)  
        # Access MATLAB struct fields using dot notation
        ELECTRON = data['ELECTRON']
        IONS_1 = data['IONS_1']
        rho0 = data['rho']
        Nprf0 = ELECTRON.density
        TprfE0 = ELECTRON.temperature*1000
        TprfI0 = IONS_1.temperature*1000
        # read data from gfile
        sfile = 'input/CFEDR/CFEDR_Conventional_H_mode_V1_20240522/gfile_efit'
        gdata, ireadok = read_gfile_func(sfile, 12, 0, 0)
        gvar = gdata
        
        RRR0 = gvar['rbbbs']
        ZZZ0 = gvar['zbbbs']        
        Rsep0 = gvar['rmaxis']
        Zsep0 = gvar['zmaxis']
        FWR = gvar['xlim']
        FWZ = gvar['ylim']
        
    else:
        raise ValueError("unknow device: {}".format(name))


# %%

def get_model_data(pre: str):
    """
    Reads the file '{pre}fort.70' and returns the device name, number of species,
    species list (with conversions applied), concentrations (as a NumPy array),
    and the tokamak geometric parameters corresponding to the device.
    
    Parameters:
        pre (str): The file path prefix (e.g., "/path/to/"). Ensure the path ends with a slash.
    
    Returns:
        tuple: (device_name, num_species, name_species, concentrations, minor_radius, major_radius)
    """
 
    filename = pre + 'fort.70'
    with open(filename, "r") as f:
        lines = f.readlines()

    device_name = None
    num_species = None
    species = None
    concentrations = None

    # Loop over the lines to find the keywords and extract the data on the next line.
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("Name of the case:"):
            device_name = lines[i+1].strip()
        elif line.startswith("Number of species:"):
            num_species = int(lines[i+1].strip())
        elif line.startswith("Names of the species:"):
            species = lines[i+1].split()
        elif line.startswith("Concentrations:"):
            conc_line = lines[i+1].strip()
            concentrations = np.array([float(val) for val in conc_line.split()]) * 100

    # Replace species names according to species_dict without changing their case.
    name_species = [species_dict.get(item, item) for item in species]

    # Retrieve the tokamak parameters for the device from the parameter database.
    if device_name in parameter_dict:
        minor_radius, major_radius, position_axis = parameter_dict[device_name]
    else:
        available_devices = ", ".join(parameter_dict.keys())
        raise ValueError(
            f"Unsupported device: '{device_name}'. Valid options: {available_devices}"
        )

    return device_name, num_species, name_species, concentrations, minor_radius, major_radius, position_axis



