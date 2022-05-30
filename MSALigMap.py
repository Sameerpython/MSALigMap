#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 21:34:40 2022

@author: xhasam
"""

import sys
import subprocess


feature= sys.argv[1]
sequence= sys.argv[2]

if feature == 'Lig':
    lig= sys.argv[3]
    
    subprocess.call("python LigandBinding.py %s %s" %(sequence, lig), shell=True)

elif feature == 'DNA':
    subprocess.call("python DNABinding.py %s" %(sequence), shell=True)
    
    
    




