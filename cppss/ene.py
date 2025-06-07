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

######

def fit_mobile(mec):
    ndata = len(mec.index)
    serie=mec.name
    if serie=="Emec":
        guess=[-1.958621e35, 5.8e29, 12, -np.pi, 2.5e29, 0.2405, 0, -5e26]
        inf=[min(mec), 3e29, 10, -np.pi, 1.5e29, 0.235, -np.pi, -1e27]
        sup=[max(mec), 10e29, 15, np.pi, 10e29, 0.25, np.pi, -1e26]
    elif serie=="-Ecin" or serie=="Epot/2":
        guess=[mec.mean(), 15e33, 12.03, np.pi/2, 2.5e27	, 0.2405, 0, -0]
        inf=[min(mec), 5e33, 10, -np.pi, 1.5e27, 0.235, -np.pi, -np.inf]
        sup=[max(mec), 5e35, 15, np.pi, 10e27, 0.25, np.pi, 0]
    elif serie=="ddd":
        guess=[mec.mean(), 4.2e10, 12.03, np.pi/2, 0	, 0.2405, 0, -0]
        inf=[min(mec), 1e10, 10, -np.pi, 0, 0.235, -np.pi, -np.inf]
        sup=[max(mec), 5e10, 15, np.pi, 1, 0.25, np.pi, np.inf]

    popt, pcov = sp.optimize.curve_fit(cos_func, mec.index, mec, p0=guess, bounds=[inf,sup])
    R2 = R_squared(mec.index, mec, cos_func, *popt)
    return popt, pcov, R2

#######

def fit_fisso(mec):
    ndata = len(mec.index)
    serie=mec.name
    guess=[-1.9586215e35, 0.7e28, 0.71, -np.pi, 0.7e28, 0.07, 0, -4e26]
    inf=[min(mec), 1e27, 0.235, -np.pi, 1.5e27, 0, -np.pi, -1e27]
    sup=[max(mec), 10e28, 1, np.pi, 10e28, 0.1, np.pi, -1e26]

    popt, pcov = sp.optimize.curve_fit(cos_func, mec.index, mec, p0=guess, bounds=[inf,sup])
    sigma = np.sqrt(ndata * np.diag(pcov))
    R2 = R_squared(mec.index, mec, cos_func, *popt)
    return popt, pcov, R2

#########

def stampa(popt, pcov, R2, ndata):
    sigma = np.sqrt(ndata * np.diag(pcov))
    print(f"Valore di partenza dell'oscillazione = {popt[0]} +/- {sigma[1]}")
    print(f"Ampiezza portante = {popt[1]} +/- {sigma[1]}")
    print(f"Periodo portante = {popt[2]} +/- {sigma[2]}")
    print(f"fase portante = {popt[3]*180/np.pi} +/- {sigma[3]*180/np.pi}")
    print(f"Ampiezza modulante = {popt[4]} +/- {sigma[4]}")
    print(f"Periodo modulante = {popt[5]*365.26} +/- {sigma[5]*365.26}")
    print(f"fase modulante= {popt[6]*180/np.pi} +/- {sigma[6]*180/np.pi}")
    print(f"pendenza decrescita= {popt[7]} +/- {sigma[7]}")
    print(f"R^2 = {R2}")

#######

def Evst(mec, popt):
    serie=mec.name
    ndata = len(mec.index)
    #"""confronto sullo stesso grafico
    plt.figure(figsize=(10,6))
    plt.title("Energia nel tempo")
    plt.xlabel("tempo [anni]")
    plt.ylabel("Energia [J]")
    plt.plot(mec.index, mec, label=serie)
    #plt.plot(mec.index, cos_func(mec.index, *popt), '-', label='Fit')
    plt.legend()
    plt.grid(True)
    plt.show()
    """
    fig, axs = plt.subplots(2, figsize=(10,6))
    fig.suptitle('Energia nel tempo')
    axs[1].plot(mec.index, cos_func(mec.index, *popt), '-', label='Fit')
    axs[0].plot(mec.index, mec, label=serie)
    axs[1].set_title('f(t)=A+Bcos($2 \pi t/C$-D)+Ecos($2\pi t/F$-G)+Ht')
    for ax in axs.flat:
        ax.set(xlabel='Tempo [anni]', ylabel='Energia [J]')
        ax.grid(True)
        ax.set_xlim(0,mec.index.max())
        ax.label_outer()
    plt.legend()
    plt.show()
    #"""

#####

def delta(mec):
    ndata = len(mec.index)
    plt.figure(figsize=(10,6))
    plt.title("variazione Energia meccanica nel tempo")
    plt.xlabel("tempo [anni]")
    plt.ylabel("$\Delta$E [J]")
    plt.plot(mec.index, mec-mec.iloc[0])
    plt.show()

#####

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

#passo i file coi dati dq linea di comando delle energie e distanze dal sole di un pianeta
df = []
for arg in sys.argv:
	if arg!="ene.py":
		d=pd.read_csv(arg, sep=' ')
		d.name = arg
		df.append(d)
		print(df)

#df=df.drop('0', axis=1)
ene=df[0];
#dist=df[1]
#dist.set_index('0', inplace=True)
ene.set_index('0', inplace=True)
ene.columns = ['Epot/2', '-Ecin', 'Emec']
#dist.columns = ['ddd']
#print(ene["Emec"].name)
#ene=ene.head(int(len(ene.index)/15))

#"""fit con sole in moto, tiene periodo giove ecc
mec=ene["Emec"]
popt, pcov, R2=fit_mobile(mec)
stampa(popt, pcov, R2, len(mec.index))
Evst(mec, popt)
#plot dello scostamento dell'energia nel tempo dal valore iniziale
delta(mec)
#"""
"""fit e grafici energia cinetica e potenziale
mec=ene['-Ecin']
popt, pcov, R2=fit_mobile(mec)
stampa(popt, pcov, R2, len(mec.index))
Evst(mec, popt)

mec=ene['Epot/2']
popt, pcov, R2=fit_mobile(mec)
stampa(popt, pcov, R2, len(mec.index))
Evst(mec, popt)
#"""

"""parametri fit su brevbe periodo a sole fisso
popt, pcov, R2=fit_(mec)
stampa(popt, pcov, R2, len(mec.index))
Evst(mec, popt)
"""