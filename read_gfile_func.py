import numpy as np
import re

def read_gfile_func(filename, nfcoil=18, nesum=6, old_efit=0):

    if not isinstance(filename, str):
        print(f"error: {filename}")
        return {}, 0
    
    try:
        with open(filename, 'r') as f:
            content = f.read()
            lines = content.strip().split('\n')
    except Exception as e:
        print(f"error {filename}: {e}")
        return {}, 0
    
    line_idx = 0
    
    ecase = lines[line_idx][:48]
    parts = lines[line_idx][48:].strip().split()
    imfit = int(parts[0])
    nw = int(parts[1])
    nh = int(parts[2])
    line_idx += 1
    

    vals = parse_floats(lines[line_idx])
    xdim, zdim, rzero, rgrid1, zmid = vals[:5]
    line_idx += 1
    
    vals = parse_floats(lines[line_idx])
    rmaxis, zmaxis, ssimag, ssibry, bzero = vals[:5]
    line_idx += 1
    
    vals = parse_floats(lines[line_idx])
    cpasma, ssimag, xdum, rmaxis, xdum = vals[:5]
    line_idx += 1
    
    vals = parse_floats(lines[line_idx])
    zmaxis, xdum, ssibry, xdum, xdum = vals[:5]
    line_idx += 1
    
    fpol = read_array_from_lines(lines, line_idx, nw)
    line_idx += count_lines_for_n_values(nw)
    
    pres = read_array_from_lines(lines, line_idx, nw)
    line_idx += count_lines_for_n_values(nw)
    
    ffprim_raw = read_array_from_lines(lines, line_idx, nw)
    ffprim = -np.sign(cpasma) * ffprim_raw
    line_idx += count_lines_for_n_values(nw)
    
    pprime_raw = read_array_from_lines(lines, line_idx, nw)
    pprime = -np.sign(cpasma) * pprime_raw
    line_idx += count_lines_for_n_values(nw)
    
    psirz_flat = read_array_from_lines(lines, line_idx, nw * nh)
    psirz = np.reshape(psirz_flat, (nw, nh), order='F').T
    line_idx += count_lines_for_n_values(nw * nh)
    
    qpsi = read_array_from_lines(lines, line_idx, nw)
    line_idx += count_lines_for_n_values(nw)
    
    vals = parse_floats(lines[line_idx])
    nbbbs = int(vals[0])
    limitr = int(vals[1])
    line_idx += 1
    
    boundary_vals = read_array_from_lines(lines, line_idx, nbbbs * 2)
    rbbbs = boundary_vals[0::2]  
    zbbbs = boundary_vals[1::2]  
    line_idx += count_lines_for_n_values(nbbbs * 2)
    
    limiter_vals = read_array_from_lines(lines, line_idx, limitr * 2)
    xlim = limiter_vals[0::2]  
    ylim = limiter_vals[1::2]  
    
    sn_ip = np.sign(cpasma)
    psirz = -sn_ip * psirz * 2 * np.pi
    ssimag = -sn_ip * ssimag * 2 * np.pi
    ssibry = -sn_ip * ssibry * 2 * np.pi
    
    rg = np.linspace(rgrid1, rgrid1 + xdim, nw)
    zg = np.linspace(zmid - zdim/2, zmid + zdim/2, nh)
    
    gdata = {
        'bzero': bzero,
        'cpasma': cpasma,
        'ecase': ecase,
        'ffprim': ffprim,
        'fpol': fpol,
        'limitr': limitr,
        'nbbbs': nbbbs,
        'nh': nh,
        'nw': nw,
        'pprime': pprime,
        'pres': pres,
        'psirz': psirz,
        'qpsi': qpsi,
        'rbbbs': rbbbs,
        'rgrid1': rgrid1,
        'rmaxis': rmaxis,
        'rzero': rzero,
        'ssibry': ssibry,
        'ssimag': ssimag,
        'xdim': xdim,
        'xlim': xlim,
        'ylim': ylim,
        'zbbbs': zbbbs,
        'zdim': zdim,
        'zmaxis': zmaxis,
        'zmid': zmid,
        'rg': rg,
        'zg': zg
    }
    
    return gdata, 1

def parse_floats(line):

    try:
        return [float(x) for x in line.split()]
    except ValueError:

        pattern = r'[+-]?\d+\.\d+E[+-]\d+'
        matches = re.findall(pattern, line)
        if matches:
            return [float(match) for match in matches]
        else:
            result = []
            chars_per_value = 16
            for i in range(0, len(line), chars_per_value):
                if i + chars_per_value <= len(line):
                    val_str = line[i:i+chars_per_value]
                    try:
                        result.append(float(val_str))
                    except ValueError:
                        pass
            return result

def count_lines_for_n_values(n, values_per_line=5):
    return (n + values_per_line - 1) // values_per_line

def read_array_from_lines(lines, start_line, n):
    values = []
    line_idx = start_line
    
    while len(values) < n and line_idx < len(lines):
        new_values = parse_floats(lines[line_idx])
        values.extend(new_values)
        line_idx += 1
    
    return np.array(values[:n])

# =============================================================================
# 
# sfile = 'input/CFEDR/CFEDR_Conventional_H_mode_V1_20240522/gfile_efit'
# gdata, ireadok = read_gfile_func(sfile, 12, 0, 0)
# gvar = gdata
# 
# # 
# dr = gvar['xdim'] / (gvar['nw'] - 1)
# dz = gvar['zdim'] / (gvar['nh'] - 1)
# 
# # 
# r1d = np.arange(gvar['rgrid1'], gvar['rgrid1'] + gvar['xdim'] + dr/2, dr)  
# z1d = np.arange(gvar['zmid'] - gvar['zdim']/2, gvar['zmid'] + gvar['zdim']/2 + dz/2, dz)
# psirz = gvar['psirz']
# 
# plt.figure(figsize=(8, 8))
# nlevel = 80
# contour = plt.contourf(r1d, z1d, -psirz, nlevel)
# plt.colorbar()
# plt.set_cmap('jet')
# 
# # Add additional plots
# plt.plot(gvar['xlim'], gvar['ylim'], 'k')
# plt.plot(gvar['rbbbs'], gvar['zbbbs'], 'r--', linewidth=2)
# plt.scatter(gvar['rmaxis'], gvar['zmaxis'], marker='*', c='r', linewidth=2)
# 
# plt.clim(-150, 150)  # equivalent to caxis in MATLAB
# 
# =============================================================================
