import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from copy import copy
from funzioni import *
#from corpo import CorpoCeleste as cc

G = 6.67408e-11    # m3 kg-1 s-2 (Costante di gravitazione)

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
			p.allocate(1)
			
		self.ES = np.zeros(1, dtype=float) #per i plot di E ed L totali
		self.LS = np.zeros(1, dtype=float)
		self.file = outFile

	####à
	
	def add(self, corpo):
		self.pianeti.append(corpo)
		self.pianeti[-1].allocate(1)
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
	
	def cfr(self):
		m=[]
		for p in self.pianeti:
			m.append(p.massa)
			p.modE(0, self.pianeti)
			self.LS[0]+= p.L[0]
			self.ES[0]+= p.Ek[0] + p.Ep[0]/2
		mu=m[0]*m[1]/(m[0]+m[1])
		K=G*m[0]*m[1]
		
		E0=self.ES[0]
		L0=np.linalg.norm(self.LS[0])
		
		a=-K/(2*E0) #semiasse magiore
		T=2*np.pi*math.sqrt(mu*a**3/K) #eriodo
		e=math.sqrt(1+2*E0*L0*L0/(mu*K*K)); #eccentircita
		
		#calcolo le velocità orbitali associate a afelio e perielio
		rp=a*(1-e)
		ra=a*(1+e)
		vp=math.sqrt(2*(E0+K/rp)/mu)
		va=math.sqrt(2*(E0+K/ra)/mu)
		print(e, rp, vp, ra, va)
		
		dim=1440
		teta=np.arange(1,2*dim+1)*np.pi/dim #angoli di un periodo orbitale per r
		teta1=teta[:dim]
		teta2=teta[-dim:]   #serparo gli array degli angoli per non avere tempi negativi
		#print(teta, teta1, teta2)
			
		tempi, teta = self.monto(teta1, teta2, e, T)	#scelgo come montare in base a afe o peri

		r=a*(1-e**2)/(1+e*np.cos(teta)) #distanze relative al primo periodo
		r, tempi = self.taglio(r, tempi, T, dim)				
		#self.stmp(r,tempi, "sol analitica")
		
		#seleziono dai dati raccolti il campione corrispondente ai tempi che ho generato
		r_num=self.dati(tempi)
		#self.stmp(r_num,tempi, "sol numerica")
		
		self.stmp(r_num,tempi, "differenza", r)
		
	#####
	
	def monto(self, teta1, teta2, e, T):
		BB=math.sqrt((1-e)/(1+e))	
		psi1=2*np.arctan(BB*np.tan(teta1/2))
		psi2=2*np.arctan(BB*np.tan(teta2/2))
		
		if self.file.find("afe") != -1: #se passo dati afelio devo montare i trempi e le distanze in modo inverso
			t2=(psi1-e*np.sin(psi1))*T/(2*np.pi)+T/2
			t1=T/2+(psi2-e*np.sin(psi2))*T/(2*np.pi)
			tempi=np.append(t1,t2) #array che raccoglie tempi relativi ad 1 periodo orbitale
			#print(tempi/(3600*24*365.24), tempi.size)		
			teta=np.append(teta2,teta1) #scambio gli angoli per comodità
		
		else:
			t1=(psi1-e*np.sin(psi1))*T/(2*np.pi)
			t2=T+(psi2-e*np.sin(psi2))*T/(2*np.pi)
			tempi=np.append(t1,t2) #array che raccoglie tempi relativi ad 1 periodo orbitale
			#print(tempi/(3600*24*365.24), tempi.size)
			teta=np.append(teta1,teta2) #scambio gli angoli per comodità
		
		return tempi, teta
	
	#####
	
	def dati(self, tempi):
		df = pd.read_csv(self.file, sep=' ', names=['sol', 'dist'], skiprows=0)
		#rimuovo i duplicati negli indici perché dan fastidio
		df= df.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')
		#print(df)
		#print(df.index)
		closest_indices = df.index.get_indexer(tempi/(365.26*24*3600), method='nearest')
		cr = df.iloc[closest_indices]
		#print(cr)
		
		d=cr['dist']
		return d
	
	#####
	
	def taglio(self, r, tempi, T, dim):
		r1=r[:dim]
		r2=r[-dim:]
		t1=tempi[:dim]
		t2=tempi[-dim:]
		
		q=365.26*3600*24*self.steps #conto quante orbite faccio nella simulazione cxonsiderata (ha self.steps anni terestri)
		n=q/T
		
		for i in range(2, int(n+1)*2): #così mi prende i valori fino alla fine del ciclo n+1, che poi taglierè
			if i%2==0:
				t0=tempi[-1]
				tempi=np.append(tempi, t0+t1)
				#print(tempi.size)
				r=np.append(r,r1)

			else:
				t0=tempi[-1]-T/2
				tempi=np.append(tempi, t0+t2)
				#print(tempi.size+1)
				r=np.append(r,r2)

		if q%T!=0: #taglio i pezzi extra se non ho multiplo intero del periodo
			tempi=tempi[:int(tempi.size*n/int(n+1))] #devo prendere la frazione del vettore contando che ho messo n+1 cicli
			r=r[:int(n*r.size/int(n+1))]
			print("taglio", tempi.size, r.size)
		
		return r, tempi
	
	#####
	
	def stmp(self, r, tempi, lab, d=np.zeros(2)):
		tt=tempi/(3600*365.26*24)
		if np.linalg.norm(d)!=0:
			fig = plt.figure(figsize=(19,9))
			plt.plot(tt, r, label = "sol numerica")
			plt.plot(tt, d, label = "sol analitica")
			plt.title("Distanza dal sole in fz del tempo")
			plt.xlabel("tempo [anni]")
			plt.ylabel("distanza [m]")
			plt.legend()
			plt.show()
			r=r-d
	
		fig = plt.figure(figsize=(19,9))
		plt.plot(tt, r, label = lab)
		plt.title("Distanza dal sole in fz del tempo - "+lab)
		plt.xlabel("tempo [anni]")
		plt.ylabel("distanza [m]")
		plt.legend()
		plt.show()
	
	#####
	
	def tempi(self):
		return 0
	
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
		plt.title("Sistema: energia meccanica totale [J]")
		plt.xlabel("E [J]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()
	
	#####
	
	def save(self, out):
		print("g")
	
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
		
