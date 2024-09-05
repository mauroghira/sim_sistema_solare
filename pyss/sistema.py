import numpy as np
import math
import matplotlib.pyplot as plt
from funzioni import *
#from corpo import CorpoCeleste as cc

# Classe Sistema
class sistema:
	"""Classe per la descrizione del sistema solare"""

	def __init__(self, pianeti, deltaT, steps, outFile):
		""" Costruttore: lista pianeti, delta T e numero di step """
		self.pianeti = pianeti   # List di pianeti
		self.deltaT  = deltaT    # Granularità temporale
		self.steps   = steps     # Numero di step
		self.n       = 0
		print("Durata:", steps, "step da ",  self.deltaT, "s, ",
		      steps*self.deltaT//(3600*24*365), "anni" )
		for p in self.pianeti:   # Preallocazione della memoria
			p.allocate(steps-1)
			
		self.ES = np.zeros(steps-1, dtype=float) #per i plot di E ed L totali
		self.LS = np.zeros(steps-1, dtype=float)
		self.file = outFile

	####à
	
	def add(self, corpo):
		self.pianeti.append(corpo)
		self.pianeti[-1].allocate(self.steps-1)
		#print(angolo(corpo.pos, corpo.vel))

	#####
	
	def modD(self, anni):
		self.steps = anni*365*24*3600//self.deltaT
		self.extend(self.steps)
		for p in self.pianeti:
			p.ext(self.steps)

	#####
	
	def extend(self, dim):
		self.ES = np.append(self.ES, np.zeros(dim-1))
		self.LS = np.append(self.LS, np.zeros(dim-1))

	#####

	def simulazione(self, mode):
		""" Simulazione per il numero di step impostati """
		cnt = 1
		for i in range(1, self.steps):
			if i%(self.steps/20) == 0:						# Stampa avanzamento
				print("%3d%%" % (cnt*5), end=" "); cnt+=1
		    	
				self.print()
			
			self.step(mode)
			#self.n += 1
			
	#####
	
	def step(self, mode):
		""" Evoluzione di uno step """
		EE=0
		LL=0
		for planet in self.pianeti:
			planet.update(self.pianeti, self.deltaT, self.n, mode)
	
		for planet in self.pianeti:
		#calcolo L ed E tot
			LL+= planet.L[self.n]
			EE+= planet.Ek[self.n] + planet.Ep[self.n]/2
			
		self.LS[self.n]=LL
		self.ES[self.n]=EE
		
		#implementa l'output nel file per i grafici
		
		self.n += 1
	
	#####
	
	def selectplanet(self, cmd, val):
		i=-1
		for p in self.pianeti:
			if p.nome.find(cmd) == 0:
				p.selplot(val, self.deltaT)
				i=0
				break
		if i == -1:
			print("pianeta non riconosciuto")
	
	#####

	def selectplot(self, num):
		match num:
			case 0:
				self.plotLL()			
			case 1:
				self.plotEE()
			case 2:
				self.plotP()
			case _:
				print("indice non valido")

	#####

	def plotP(self, outfile = None):
		"""Disegno xy"""
		fig = plt.figure(figsize=(9,9))
		for p in self.pianeti:
			if p.nome=="Sole":
				plt.plot(p.Xs, p.Ys, "*", markersize=10, label = p.nome)
			else:
				plt.plot(p.Xs, p.Ys, label = p.nome)
		plt.xlabel("x [m]")
		plt.ylabel("y [m]")
		plt.legend()
		plt.show()
	
	#####

	def plotLL(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		st = "Media: "+str(self.LS.mean())+"\n Dev Std: "+str(np.std(self.LS))
		plt.hist(self.LS, bins=200, range=(self.LS.mean()-np.std(self.LS)*4, self.LS.mean()+np.std(self.LS)*4), label = st)
		plt.axvline(x = self.LS.mean(), color = 'red', linestyle = '--', alpha = 0.5, label="media")
		plt.axvline(x = self.LS.mean()-np.std(self.LS), color = 'green', linestyle = 'dotted', alpha = 1)
		plt.axvline(x = self.LS.mean()+np.std(self.LS), color = 'green', linestyle = 'dotted', alpha = 1)
		plt.axvline(x = self.LS.mean()-np.std(self.LS)*3, color = 'black', linestyle = 'dashdot', alpha = 1)
		plt.axvline(x = self.LS.mean()+np.std(self.LS)*3, color = 'black', linestyle = 'dashdot', alpha = 1)
		plt.title("Sistema: momento angolare totale [kg*m^2/s]")
		plt.xlabel("L [kg*m^2/s]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()
		
	#####

	def plotEE(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		st = "Media: "+str(self.ES.mean())+"\n Dev Std: "+str(np.std(self.ES))
		plt.hist(self.ES, bins=200, range=(self.ES.mean()-np.std(self.ES)*5, self.ES.mean()+np.std(self.ES)*5), label = st)
		plt.axvline(x = self.ES.mean(), color = 'red', linestyle = '--', alpha = 0.5, label="Media")
		plt.axvline(x = self.ES.mean()-np.std(self.ES), color = 'green', linestyle = 'dotted', alpha = 1)
		plt.axvline(x = self.ES.mean()+np.std(self.ES), color = 'green', linestyle = 'dotted', alpha = 1)
		plt.axvline(x = self.ES.mean()-np.std(self.ES)*3, color = 'black', linestyle = 'dashdot', alpha = 1)
		plt.axvline(x = self.ES.mean()+np.std(self.ES)*3, color = 'black', linestyle = 'dashdot', alpha = 1)
		plt.title("Sistema: energia meccanica totale [J]")
		plt.xlabel("E [J]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()
	
	#####
	
	def save(self, out):
		for p in self.pianeti:
			o = open("res/"+p.nome+"_"+out, 'w')  # 'w' means writing
			for i in range(p.L.size):
				o.write(str(p.L[i])+" "+str(p.Ek[i]+p.Ep[i])+" "+str(p.ecc[i])+" "+str(p.teta[i])+" "+str(p.DS[i])+" "+str(p.V[i])+str(p.Xs[i])+" "+str(p.Ys[i])+"\n")
			o.close()
		o = open("res/sistema_"+out, 'w')  # 'w' means writing	
		for i in range(p.L.size):
			o.write(str(self.LS[i])+" "+str(self.ES[i])+"\n")
		o.close()
	#####

	def print(self):
		print("====================")
		for pla in self.pianeti:
			print(pla)
		
		
		if self.n==0:
			print("L0= ", self.LS[0], "kg*m^2/s")
			print("E0= ", self.ES[0], "J")
		else:
			print("Ltot: ", self.LS[self.n-1],"kg*m^2/s")
			print("Etot: ", self.ES[self.n-1], "J")			
		
