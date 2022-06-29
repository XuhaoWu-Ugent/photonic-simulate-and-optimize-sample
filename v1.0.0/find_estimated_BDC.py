import pandas as pd
import numpy as np
import sys, os
#from utils import * #necessary if want to extract data from .csv to array

# find estimated points for 50/50
df = pd.read_csv("max_min_coupling_v1.txt",delimiter=" ")
df.to_csv( "max_min_coupling_v1.csv", encoding='utf-8', index=False) 
df_50_target = df[(df["average_power"] > 0.45) & (df["average_power"] < 0.55) & (df["min_max_coupling"] < 2.0)]

# find estimated points for 60/40
df_60_target = df[(df["average_power"] > 0.55) & (df["average_power"] < 0.65) & (df["min_max_coupling"] < 2.0)]
print df_50_target, df_60_target

df_50_target.to_csv("estimated_target_points_0.5.csv", encoding='utf-8', index=False)
df_60_target.to_csv("estimated_target_points_0.6.csv", encoding='utf-8', index=False)

