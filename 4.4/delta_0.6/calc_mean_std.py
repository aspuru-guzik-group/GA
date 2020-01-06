#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 11:20:08 2019

@author: akshat
"""
import numpy as np

filename = './IMPROVEMENT.txt'

with open(filename) as f:
    content = f.readlines()
content = [x.strip() for x in content] 
content = [x[1:len(x)-1] for x in content]
content = [float(x) for x in content]



fail_cases = [x for x in content if x == 0]


# Success percentage:
success_P = 1 - (len(fail_cases) / len(content))
print('Success percentage: ', success_P * 100)


# Average Improvement and standard deviation:
A = np.array(content)
print('Mean: ', A.mean())
print('Std : ', A.std())


# Get the indices of the 3 best elements: 
print(A.argsort()[-3:][::-1])
# Output: [649 350 762]