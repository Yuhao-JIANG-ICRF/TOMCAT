import numpy as np

'''
  'Device name': (minor radius; major radius; magnetic axis)
'''
parameter_dict = {    
    'CFEDR': (2.5, 7.8, 8.2),
    'CFETR': (2.5, 7.8, 8.2),
    'WEST': (0.43, 2.5, 2.5)
    # add new device here
}

species_dict = {
    'Li' : '^7Li',
    'He3': '^3He',
    'He4': '^4He'
}


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




# %%

# filename = pre + 'fort.70'
# with open(filename, "r") as f:
#     lines = f.readlines()

# # case = None
# # num_species = None
# # species = None
# # concentrations = None

# for i, line in enumerate(lines):
#     line = line.strip()
#     if line.startswith("Name of the case:"):
#         device_name = lines[i+1].strip()
#     elif line.startswith("Number of species:"):
#         num_species = int(lines[i+1].strip())
#     elif line.startswith("Names of the species:"):
#         species = lines[i+1].split()
#     elif line.startswith("Concentrations:"):
#         conc_line = lines[i+1].strip()
#         concentrations = np.array([float(val) for val in conc_line.split()])*100

# name_species = [species_dict.get(item, item) for item in species]

# print("Device name:", device_name)
# print("Number of species:", num_species)
# print("Species:", name_species)
# print("Concentrations:", concentrations,'%')


# def get_tokamak_parameters(device_name: str) -> tuple[float, float]:
#     """
#     Retrieve tokamak geometric parameters
    
#     Returns:
#         tuple: (minor_radius[m], major_radius[m])
    
#     Raises:
#         ValueError: For unsupported device names
    
#     Example:
#         >>> a, R0 = get_tokamak_parameters('CFEDR')
#         >>> print(f"Minor radius: {a}, Major radius: {R0}")
#         Minor radius: 2.5, Major radius: 7.8
#     """
    

    
#     if device_name in parameter_db:
#         return parameter_db[device_name]
#     else:
#         available_devices = ", ".join(parameter_db.keys())
#         raise ValueError(
#             f"Unsupported device: '{device_name}'. "
#             f"Valid options: {valid_devices}"
#         )



