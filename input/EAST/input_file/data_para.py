#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:23:50 2025

@author: YJ281217
"""
import re
import numpy as np
def parse_data_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    n_y = None
    data_follow_index = None
    for i, line in enumerate(lines):
        if "# OF Y PTS" in line:
            match = re.search(r'(\d+)', line)
            if match:
                n_y = int(match.group(1))
        if "DATA FOLLOW:" in line:
            data_follow_index = i
            break

    if n_y is None:
        raise ValueError("No n_y")
    if data_follow_index is None:
        raise ValueError("No DATA FOLLOW")

    data_str = " ".join(lines[data_follow_index+1:])
    
    nums_str = re.findall(r'[-+]?\d*\.\d+e[+-]\d+', data_str)
    nums = [float(num) for num in nums_str]
    
    if len(nums) < 1 + 2 * n_y:
        raise ValueError("No enough number")
    
    y_array = nums[1:1+n_y]
    f_array = nums[1+n_y:1+n_y+n_y]
    rho = np.array(y_array)
    NT = np.array(f_array)
    
    return rho, NT