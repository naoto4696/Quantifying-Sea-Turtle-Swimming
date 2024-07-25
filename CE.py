import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import glob
import cv2
import csv
import random
import warnings
import time
from progress.bar import Bar

# ignore warnings
warnings.filterwarnings("ignore")

# variables
Maxfilenum = 18000
sample_interval = 4
lumi_min = 200
lumi_max = 255

# set pass to csv files
output_csv_file = "output.csv"
with open(output_csv_file, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(['coords', 'chain', 'x', 'y'])

# contour extraction
bar = Bar('Extracting contours...', max=Maxfilenum)
for i in range(1, Maxfilenum):
    bar.next()
    filename = f"{i:05d}.jpeg"
    Rowname = f"{i-1:05d}"
    image = cv2.imread(filename)
    if image is not None:
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        gray = 255 - gray
        ret, thresh = cv2.threshold(gray, lumi_min, lumi_max, cv2.THRESH_BINARY)
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        largest_contour = None
        max_area = 0
        
        for contour in contours:
            area = cv2.contourArea(contour)
            if area > max_area:
                max_area = area
                largest_contour = contour

        with open(output_csv_file, mode='a', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            J = -1
            for p in range(0,len(largest_contour), sample_interval):
                x, y = largest_contour[p][0]
                J = J + 1
                csv_writer.writerow([Rowname,J, x, y])
                
bar.finish()
                            
# normalization and rotation
df = pd.read_csv('output.csv')
df.set_index(['coords', 'chain'], inplace=True)
datalen = len(df.index.get_level_values('coords').unique())
ct = pd.read_csv('CT.csv')
ct.set_index(['coords'], inplace=True)
ctlen = len(ct)
minlen = min(datalen,ctlen)

df = df.loc[list(range(0,minlen))]
ct = ct.loc[list(range(0,minlen))]

bar = Bar('Centering and Rotating...', max=datalen)
for i in range(0,datalen):
    bar.next()
    df_selected = df.loc[i]
    ct_selected = ct.loc[i]
    center_x = (ct_selected['x_head'] + ct_selected['x_tail'])/2
    center_y = (ct_selected['y_head'] + ct_selected['y_tail'])/2
    df_centered_x = df_selected['x'] - center_x
    df_centered_y = df_selected['y'] - center_y
    tail_centered_x = ct_selected['x_tail'] - center_x
    tail_centered_y = ct_selected['y_tail'] - center_y
    r = math.sqrt(tail_centered_x**2 + tail_centered_y**2)
    cos = tail_centered_x/r
    sin = tail_centered_y/r
    df_centered_rotated_x = (df_centered_x * cos + df_centered_y * sin)/r
    df_centered_rotated_y = (df_centered_y * cos - df_centered_x * sin)/r
    if i == 0:
        norm_rotated_x = df_centered_rotated_x.tolist()
        norm_rotated_y = df_centered_rotated_y.tolist()
    else:
        norm_rotated_x = norm_rotated_x + df_centered_rotated_x.tolist()
        norm_rotated_y = norm_rotated_y + df_centered_rotated_y.tolist()
        
df['norm_rotated_x'] = norm_rotated_x
df['norm_rotated_y'] = norm_rotated_y

df.drop('x', axis=1, inplace=True)
df.drop('y', axis=1, inplace=True)
df.to_csv('output_norm.csv', index=True)

bar.finish()

# rechain of contours
df = pd.read_csv('output_norm.csv')
df.set_index(['coords', 'chain'], inplace=True)
df_rechain = pd.DataFrame(columns=["coords","chain","norm_rotated_x","norm_rotated_y"])
datalen = len(df.index.get_level_values('coords').unique())

bar = Bar('Rechaining...', max=datalen)
for i in range(0,datalen):
    bar.next()
    dist = np.sqrt((df.loc[i]["norm_rotated_x"]-1)**2+(df.loc[i]["norm_rotated_y"])**2)
    startpoint = dist.idxmin(axis=0)
    point_num = len(dist)
    newchain = list(range(startpoint,point_num)) + list(range(0,startpoint))
    df_part = df.loc[i].reindex(newchain).reset_index(drop=True)
    new_col1 = list(range(point_num))
    new_col2 = [i] * point_num
    df_part.insert(0,"chain",new_col1)
    df_part.insert(0,"coords",new_col2)
    df_rechain = pd.concat([df_rechain,df_part], ignore_index=True)
    
df_rechain.set_index(['coords', 'chain'], inplace=True)
df_rechain.to_csv('output_rechain.csv', index=True)

bar.finish()
            



