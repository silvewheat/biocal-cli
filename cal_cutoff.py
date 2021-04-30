# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:07:23 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import typer
import numpy as np
import pandas as pd


def cal_outlier(counts, x=3):
    Q1 = np.percentile(counts, 25)
    Q3 = np.percentile(counts, 75)
    IQR = Q3 - Q1
    upper_boundary = Q3 + IQR * x
    lower_boundary = Q1 - IQR * x
    return upper_boundary, lower_boundary

