import numpy as np
import scipy as sp
import pandas as pd
import math
import matplotlib.pyplot as plt
import statistics as sta
import sys

data = []
for arg in sys.argv:
	if arg!="errori.py":
		df=pd.read_csv(arg, sep=' ')
		df.name = arg
		data.append(df)
		print(df)

"""
for df in data:
	for i in df.index:
		if i < df.index.max():
			df.iloc[i]=df.iloc[i+1]-df.iloc[i]
	df=df.drop(df.index.max())
	print(df.iloc[df.index*2//3])
"""

if len(data) ==2:	
	for col in data[0]:
		data[0][col]=data[0][col]-data[1][col]

		plt.figure(figsize=(18,8))
		plt.title(str(col)+" - differenza tra vari metodi di approssimazione")
		plt.xlabel("tempo [step]")
		plt.ylabel("distanza [m]")
		plt.plot(data[0].index, data[0].loc[:,col])	
		plt.show()

else:
	for col in data[0]:
		plt.figure(figsize=(18,8))
		plt.title(str(col)+" - Confronto tra vari metodi di approssimazione")
		plt.xlabel("tempo [step]")
		plt.ylabel("distanza [m]")
		for dt in data:
			#print(dt)
			if dt.name.find("1800")!=-1:
				lol="1800s - mode "+dt.name[-5]
				plt.plot(dt.index, dt.loc[:,col], label=lol, )	
			if dt.name.find("3600")!=-1:
				lol="3600s - mode "+dt.name[-5]
				plt.plot(dt.index*2, dt.loc[:,col], label=lol)
		plt.legend()
		plt.show()
