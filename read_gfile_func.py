
globals().clear()

import matplotlib.pyplot as plt
import numpy as np
import re

def read_gfile_func(filename, nfcoil=18, nesum=6, old_efit=0):
    """
    读取 EFIT ASCII g-eqdsk 文件(使用 iecurr=2 标志)
    
    参数:
    -----------
    filename : str
        g文件文件名
    nfcoil : int, optional
        F线圈数量(DIII-D为18)
    nesum : int, optional
        E线圈数量(DIII-D为6)
    old_efit : int, optional
        对于旧版本的EFIT(例如NSTX)，使用old_efit=1
        
    返回:
    --------
    gdata : dict
        包含EFIT G文件变量的数据结构
    ireadok : int
        指示g文件是否成功读取的标志(0=失败, 1=成功)
    """
    # 输入验证
    if not isinstance(filename, str):
        print(f"错误: 文件名必须是字符串: {filename}")
        return {}, 0
    
    # 打开文件
    try:
        with open(filename, 'r') as f:
            content = f.read()
            lines = content.strip().split('\n')
    except Exception as e:
        print(f"错误: 无法打开文件 {filename}: {e}")
        return {}, 0
    
    line_idx = 0
    
    # 第一行包含头部和维度
    ecase = lines[line_idx][:48]
    parts = lines[line_idx][48:].strip().split()
    imfit = int(parts[0])
    nw = int(parts[1])
    nh = int(parts[2])
    line_idx += 1
    
    # 第二行: 维度
    vals = parse_floats(lines[line_idx])
    xdim, zdim, rzero, rgrid1, zmid = vals[:5]
    line_idx += 1
    
    # 第三行
    vals = parse_floats(lines[line_idx])
    rmaxis, zmaxis, ssimag, ssibry, bzero = vals[:5]
    line_idx += 1
    
    # 第四行
    vals = parse_floats(lines[line_idx])
    cpasma, ssimag, xdum, rmaxis, xdum = vals[:5]
    line_idx += 1
    
    # 第五行
    vals = parse_floats(lines[line_idx])
    zmaxis, xdum, ssibry, xdum, xdum = vals[:5]
    line_idx += 1
    
    # 读取数组
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
    
    # 读取psirz矩阵 - 修复方法
    psirz_flat = read_array_from_lines(lines, line_idx, nw * nh)
    # 重要修复: MATLAB中psirz的存储是按列优先排列的
    # 首先将数据重新排列为nw×nh数组，然后进行转置
    psirz = np.reshape(psirz_flat, (nw, nh), order='F').T
    line_idx += count_lines_for_n_values(nw * nh)
    
    # 读取qpsi
    qpsi = read_array_from_lines(lines, line_idx, nw)
    line_idx += count_lines_for_n_values(nw)
    
    # 读取边界和限制器尺寸
    vals = parse_floats(lines[line_idx])
    nbbbs = int(vals[0])
    limitr = int(vals[1])
    line_idx += 1
    
    # 读取边界和限制器点
    boundary_vals = read_array_from_lines(lines, line_idx, nbbbs * 2)
    rbbbs = boundary_vals[0::2]  # 偶数索引
    zbbbs = boundary_vals[1::2]  # 奇数索引
    line_idx += count_lines_for_n_values(nbbbs * 2)
    
    limiter_vals = read_array_from_lines(lines, line_idx, limitr * 2)
    xlim = limiter_vals[0::2]  # 偶数索引
    ylim = limiter_vals[1::2]  # 奇数索引
    
    # 创建有用的通量
    sn_ip = np.sign(cpasma)
    psirz = -sn_ip * psirz * 2 * np.pi
    ssimag = -sn_ip * ssimag * 2 * np.pi
    ssibry = -sn_ip * ssibry * 2 * np.pi
    
    # 定义网格数组
    rg = np.linspace(rgrid1, rgrid1 + xdim, nw)
    zg = np.linspace(zmid - zdim/2, zmid + zdim/2, nh)
    
    # 创建输出结构
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
    """解析可能没有空格分隔的科学计数法浮点数行"""
    # 首先尝试简单的空格分隔
    try:
        return [float(x) for x in line.split()]
    except ValueError:
        # 尝试科学计数法模式
        pattern = r'[+-]?\d+\.\d+E[+-]\d+'
        matches = re.findall(pattern, line)
        if matches:
            return [float(match) for match in matches]
        else:
            # 最后尝试固定宽度解析
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
    """估计n个值需要的行数"""
    return (n + values_per_line - 1) // values_per_line

def read_array_from_lines(lines, start_line, n):
    """从start_line开始的行中读取n个值"""
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
# # 计算网格间距
# dr = gvar['xdim'] / (gvar['nw'] - 1)
# dz = gvar['zdim'] / (gvar['nh'] - 1)
# 
# # 创建坐标网格
# r1d = np.arange(gvar['rgrid1'], gvar['rgrid1'] + gvar['xdim'] + dr/2, dr)  # 添加 dr/2 确保包含末端点
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
