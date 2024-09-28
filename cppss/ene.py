import numpy as np
import scipy as sp
import pandas as pd
import math
import matplotlib.pyplot as plt
import statistics as sta
import sys

def cos_func(x, C, D, E, F, a, b,g, h):
    y =  C + a*np.sin(x*2*np.pi/b-g) + D*np.sin(x*2*np.pi/E-F) + h*x
    return y

def f2(x, C, b,g,a, h):
    y =  C + a*np.cos(x*2*np.pi/b-g) + h*x
    return y

def R_squared(x, y, func, *params):
    res = y - func(x, *params)
    ss_res = np.sum(res**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res/ss_tot)
    return r_squared

SMALL_SIZE = 14
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

df = []
for arg in sys.argv:
	if arg!="ene.py":
		d=pd.read_csv(arg, sep=' ')
		d.name = arg
		df.append(d)
		print(df)

#df=df.drop('0', axis=1)
ene=df[0];
dist=df[1]
dist.set_index('0', inplace=True)
ene.set_index('0', inplace=True)
ene.columns = ['Epot/2', '-Ecin', 'Emec']
dist.columns = ['ddd']
#print(df)

"""
#ene=ene.head(int(len(ene.index)/15))
mec=ene['Emec']
ndata = len(ene.index)
guess=[-1.958621e35, 5.8e29, 12, -np.pi, 2.5e29, 0.2405, 0, -5e26]
inf=[min(mec), 3e29, 10, -np.pi, 1.5e29, 0.235, -np.pi, -1e27]
sup=[max(mec), 10e29, 15, np.pi, 10e29, 0.25, np.pi, -1e26]
popt, pcov = sp.optimize.curve_fit(cos_func, ene.index, ene['Emec'], p0=guess, bounds=[inf,sup])
sigma = np.sqrt(ndata * np.diag(pcov))
R2 = R_squared(ene.index, mec, cos_func, *popt)

print(f"Valore di partenza dell'oscillazione = {popt[0]} +/- {sigma[1]}")
print(f"Ampiezza portante = {popt[1]} +/- {sigma[1]}")
print(f"Periodo portante = {popt[2]} +/- {sigma[2]}")
print(f"fase portante = {popt[3]*180/np.pi} +/- {sigma[3]*180/np.pi}")
print(f"Ampiezza modulante = {popt[4]} +/- {sigma[4]}")
print(f"Periodo modulante = {popt[5]*365.26} +/- {sigma[5]*365.26}")
print(f"fase modulante= {popt[6]*180/np.pi} +/- {sigma[6]*180/np.pi}")
print(f"pendenza decrescita= {popt[7]} +/- {sigma[7]}")
print(f"R^2 = {R2}")

plt.figure(figsize=(18,8))
plt.title("Energia meccanica nel tempo")
plt.xlabel("tempo [step]")
plt.ylabel("Energia [J]")
plt.plot(ene.index, mec, label='Emec')
plt.plot(ene.index, cos_func(ene.index, *popt), '-', label='Fit')
plt.legend()
plt.show()

fig, axs = plt.subplots(2, figsize=(19,10))
fig.suptitle('Energia meccanica nel tempo')
axs[1].plot(ene.index, cos_func(ene.index, *popt), '-', label='Fit')
axs[0].plot(ene.index, mec, label='Emec')
axs[1].set_title('f(t)=A+Bcos($2 \pi t/C$-D)+Ecos($2\pi t/F$-G)+Ht')
for ax in axs.flat:
    ax.set(xlabel='Tempo [anni]', ylabel='Energia [J]')
for ax in axs.flat:
    ax.label_outer()
plt.legend()
plt.show()
#"""

"""
#ene=ene.head(int(len(ene.index)/10))
mec=ene['-Ecin']
ndata = len(ene.index)
guess=[mec.mean(), 15e33, 12.03, np.pi/2, 2.5e27	, 0.2405, 0, -0]
inf=[min(mec), 5e33, 10, -np.pi, 1.5e27, 0.235, -np.pi, -np.inf]
sup=[max(mec), 5e35, 15, np.pi, 10e27, 0.25, np.pi, 0]
popt, pcov = sp.optimize.curve_fit(cos_func, ene.index, ene['Emec'], p0=guess, bounds=[inf,sup])
sigma = np.sqrt(ndata * np.diag(pcov))
R2 = R_squared(ene.index, mec, cos_func, *popt)

print(f"Valore di partenza dell'oscillazione = {popt[0]} +/- {sigma[1]}")
print(f"Ampiezza portante = {popt[1]} +/- {sigma[1]}")
print(f"Periodo portante = {popt[2]} +/- {sigma[2]}")
print(f"fase portante = {popt[3]*180/np.pi} +/- {sigma[3]*180/np.pi}")
print(f"Ampiezza modulante = {popt[4]} +/- {sigma[4]}")
print(f"Periodo modulante = {popt[5]*365.26} +/- {sigma[5]*365.26}")
print(f"fase modulante= {popt[6]*180/np.pi} +/- {sigma[6]*180/np.pi}")
print(f"pendenza decrescita= {popt[7]} +/- {sigma[7]}")
print(f"R^2 = {R2}")

plt.figure(figsize=(18,8))
plt.title("Energia meccanica nel tempo")
plt.xlabel("tempo [step]")
plt.ylabel("Energia [J]")
plt.plot(ene.index, mec, label='Emec')
plt.plot(ene.index, cos_func(ene.index, *popt), '-', label='Fit')
plt.legend()
plt.show()

fig, axs = plt.subplots(2, figsize=(19,10))
fig.suptitle('Energia meccanica nel tempo')
axs[1].plot(ene.index, cos_func(ene.index, *popt), '-', label='Fit')
axs[0].plot(ene.index, mec, label='Emec')
axs[1].set_title('f(t)=A+Bcos($2 \pi t/C$-D)+Ecos($2\pi t/F$-G)+Ht')
for ax in axs.flat:
    ax.set(xlabel='Tempo [anni]', ylabel='Energia [J]')
for ax in axs.flat:
    ax.label_outer()
plt.legend()
plt.show()
#"""

"""
#dist=dist.head(int(len(ene.index)/10))
d=dist['ddd']
ndata = len(dist.index)
guess=[d.mean(), 4.2e10, 12.03, np.pi/2, 0	, 0.2405, 0, -0]
inf=[min(d), 1e10, 10, -np.pi, 0, 0.235, -np.pi, -np.inf]
sup=[max(d), 5e10, 15, np.pi, 1, 0.25, np.pi, np.inf]
popt, pcov = sp.optimize.curve_fit(cos_func, dist.index, d, p0=guess, bounds=[inf,sup])
sigma = np.sqrt(ndata * np.diag(pcov))
R2 = R_squared(dist.index, d, cos_func, *popt)

print(f"Valore di partenza dell'oscillazione = {popt[0]} +/- {sigma[1]}")
print(f"Ampiezza portante = {popt[1]} +/- {sigma[1]}")
print(f"Periodo portante = {popt[2]} +/- {sigma[2]}")
print(f"fase portante = {popt[3]*180/np.pi} +/- {sigma[3]*180/np.pi}")
print(f"Ampiezza modulante = {popt[4]} +/- {sigma[4]}")
print(f"Periodo modulante = {popt[5]*365.26} +/- {sigma[5]*365.26}")
print(f"fase modulante= {popt[6]*180/np.pi} +/- {sigma[6]*180/np.pi}")
print(f"pendenza decrescita= {popt[7]} +/- {sigma[7]}")
print(f"R^2 = {R2}")

plt.figure(figsize=(18,8))
plt.title("Distanza di Giove dal Sole nel tempo")
plt.xlabel("tempo [step]")
plt.ylabel("Distanza [m]")
plt.plot(dist.index, d, label='Emec')
plt.plot(dist.index, cos_func(dist.index, *popt), '-', label='Fit')
plt.legend()
plt.show()

fig, axs = plt.subplots(2, figsize=(19,10))
fig.suptitle('Distanza Giove-SOle nel tempo')
axs[1].plot(dist.index, cos_func(dist.index, *popt), '-', label='Fit')
axs[0].plot(dist.index, d, label='dist')
axs[1].set_title('f(t)=A+Bcos($2 \pi t/C$-D)+Ecos($2\pi t/F$-G)+Ht')
for ax in axs.flat:
    ax.set(xlabel='Tempo [anni]', ylabel='Distanza [m]')
for ax in axs.flat:
    ax.label_outer()
plt.legend()
plt.show()
#"""

#parametri fit su brevbe periodo a sole fisso
#"""
ene=ene.head(int(len(ene.index)/25))
mec=ene['Emec']
ndata = len(ene.index)
guess=[-1.9586215e35, 0.7e28, 0.71, -np.pi, 0.7e28, 0.07, 0, -4e26]
inf=[min(mec), 1e27, 0.235, -np.pi, 1.5e27, 0, -np.pi, -1e27]
sup=[max(mec), 10e28, 1, np.pi, 10e28, 0.1, np.pi, -1e26]
popt, pcov = sp.optimize.curve_fit(cos_func, ene.index, ene['Emec'], p0=guess, bounds=[inf,sup])
sigma = np.sqrt(ndata * np.diag(pcov))
R2 = R_squared(ene.index, mec, cos_func, *popt)

print(f"Valore di partenza dell'oscillazione = {popt[0]} +/- {sigma[1]}")
print(f"Ampiezza portante = {popt[1]} +/- {sigma[1]}")
print(f"Periodo portante = {popt[2]} +/- {sigma[2]}")
print(f"fase portante = {popt[3]*180/np.pi} +/- {sigma[3]*180/np.pi}")
print(f"Ampiezza modulante = {popt[4]} +/- {sigma[4]}")
print(f"Periodo modulante = {popt[5]*365.26} +/- {sigma[5]*365.26}")
print(f"fase modulante= {popt[6]*180/np.pi} +/- {sigma[6]*180/np.pi}")
print(f"pendenza decrescita= {popt[7]} +/- {sigma[7]}")
print(f"R^2 = {R2}")

plt.figure(figsize=(18,8))
plt.title("Energia meccanica nel tempo")
plt.xlabel("tempo [step]")
plt.ylabel("Energia [J]")
plt.plot(ene.index, mec, label='Emec')
plt.plot(ene.index, cos_func(ene.index, *popt), '-', label='Fit')
plt.legend()
plt.show()

fig, axs = plt.subplots(2, figsize=(19,10))
fig.suptitle('Energia meccanica nel tempo')
axs[1].plot(ene.index, cos_func(ene.index, *popt), '-', label='Fit')
axs[0].plot(ene.index, mec, label='Emec')
axs[1].set_title('f(t)=A+Bcos($2 \pi t/C$-D)+Ecos($2\pi t/F$-G)+Ht')
for ax in axs.flat:
    ax.set(xlabel='Tempo [anni]', ylabel='Energia [J]')
for ax in axs.flat:
    ax.label_outer()
plt.legend()
plt.show()
#"""



#sr = df[0][df[0]['Emec'] == max(df[0]['Emec'])]
#print(sr)
