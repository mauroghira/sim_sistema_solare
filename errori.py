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

fig = plt.figure(figsize=(19,9))
for dt in data:
	#print(dt)
	if dt.name.find("1800")!=-1:
		lol="1800s - mode "+dt.name[-5]
		plt.plot(dt.index, dt.loc[:,'Terra'], label=lol)	
	elif dt.name.find("3600")!=-1:
		lol="3600s - mode "+dt.name[-5]
		plt.plot(dt.index*2, dt.loc[:,'Terra'], label=lol)

plt.title("Confronto tra vari metodi di approssimazione")
plt.xlabel("tempo [step]")
plt.ylabel("distanza [m]")
plt.legend()
plt.show()
