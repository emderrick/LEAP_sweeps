#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import os as os
import numpy as np
import csv
import pickle


# In[40]:


depth_files = os.listdir('Documents/instrain_mag_depth/')


# In[41]:


for file in depth_files:
    input_file = "Documents/instrain_mag_depth/" + file
    pond = input_file[29:41]
    print(pond)
    with open(input_file, 'rb') as handle:
        depth_dict = pickle.load(handle)
        
    pond_df_list = []
    for scaffold in depth_dict.keys():
        scaf_dict = depth_dict[scaffold]
        scaf_df = pd.DataFrame.from_dict(scaf_dict)
        scaf_df = scaf_df.transpose()
        scaf_df['position'] = scaf_df.index
        scaf_df.reset_index(drop=True, inplace=True)
        scaf_df.columns =['coverage', 'position']
        scaf_df['scaffold'] = scaffold
        pond_df_list.append(scaf_df)

    pond_df = pd.concat(pond_df_list)
    output_file = pond + '_depth.csv'
    pond_df.to_csv(output_file, index = False)

