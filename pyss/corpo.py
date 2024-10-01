import numpy as np
import math
import matplotlib.pyplot as plt
import statistics as sta
import copy

G = 6.67408e-11    # m3 kg-1 s-2 (Costante di gravitazione)

def angolo(a,b):
	nn = np.dot(a,b)
	dd= np.linalg.norm(a) * np.linalg.norm(b)
	#print(nn/dd)
	return math.acos(nn/dd)*180/math.pi

# Classe corpo celeste
class CorpoCeleste:
	"""Classe che descrive un generico corpo celeste
     caratterizzato da nome, massa, posizione e velocita'"""

	def __init__(self, nome, massa, pos, vel):
		""" Costruttore: assegna i parametri ai dati membro """
		self.nome   = nome
		self.massa  = massa
		self.pos    = np.array(pos)
		self.vel    = np.array(vel)
		self.Xs     = np.array([])   # Array delle posizioni in X
		self.Ys     = np.array([])   # Array delle posizioni in Y

		self.pos0 = np.array(pos)
		""" dati che aandrebbero in istogrammi, ma che per efficienza poi verranno messi in file esterno (o con pyroot)"""
		self.Ek = np.array([])
		self.Ep = np.array([])
		self.L = np.array([])
		self.teta = np.array([])
		self.ecc = np.array([])
		self.DS = np.array([])
		self.V = np.array([])		
    
	######

	def allocate(self, size):  #size=numero di step
		"""Pre-allocazione memoria per ottimizzare le prestazioni"""
		self.Xs = np.zeros(size)
		self.Ys = np.zeros(size)
		
		self.Ek = np.zeros(size, dtype=float)
		self.Ep = np.zeros(size, dtype=float)
		self.L = np.zeros(size, dtype=float)
		self.teta = np.zeros(size, dtype=float)
		self.ecc = np.zeros(size, dtype=float)
		self.DS = np.zeros(size, dtype=float)	
		self.V = np.zeros(size, dtype=float)
    
	#####
	
	def ext(self, dim):
		self.Xs = np.append(self.Xs, np.zeros(dim-1))
		self.Ys = np.append(self.Ys, np.zeros(dim-1))
		self.Ek = np.append(self.Ek, np.zeros(dim-1))
		self.Ep = np.append(self.Ep, np.zeros(dim-1))
		self.L = np.append(self.L, np.zeros(dim-1))
		self.teta = np.append(self.teta, np.zeros(dim-1))
		self.ecc = np.append(self.ecc, np.zeros(dim-1))
		self.DS = np.append(self.DS, np.zeros(dim-1))
		self.V = np.append(self.V, np.zeros(dim-1))	
	
	#####

	def accel(self, planets):
		"""Calcolo accelerazione in base a posizione e massa degli altri Corpi Celesti"""
		acc = np.array([0.,0., 0.])    #da modif                         # Array 2D
		for other_planet in planets:                        # Loop sui pianeti
			if self.nome != other_planet.nome:                # Evitando se stesso: i==j
				r = np.linalg.norm(self.pos - other_planet.pos) # Modulo distanza tra i e j
				scalar = G * other_planet.massa / r**(3)   # Formula di newton
				acc += scalar * (other_planet.pos - self.pos)   # Vettore accelerazione
		return acc
    
	#####
	
	def modE(self, step, planets):
		self.Ek[step] = 0.5 * self.massa * np.linalg.norm(self.vel)**2
		self.Ep[step] = 0
		for p in planets:
			if self.nome != p.nome:                # Evitando se stesso: i==j
				r = np.linalg.norm(self.pos - p.pos)
				self.Ep[step] += p.massa/r
		self.Ep[step] *= -G*self.massa
		
		self.L[step] = self.massa * np.linalg.norm(np.cross(self.pos, self.vel))
	
	#####

	def update(self, planets, deltaT, step, mode):
		"""Aggiorna la posizione del Corpo Celeste"""
		match mode:
			case 0: #eulero
				acc = self.accel(planets)
				self.vel += acc * deltaT
				self.pos += self.vel * deltaT
			
			case 1: #runge kutta 1  da moficare
				v0 = copy.copy(self.vel)
				acc = self.accel(planets)
				p0 = copy.copy(self.pos)
				dp = v0*deltaT
				self.pos=self.pos+dp
				a2 = self.accel(planets)
				vv=((acc + a2)*0.5)*deltaT;
				self.vel += vv
				pp=((self.vel + v0)*0.5)*deltaT
				self.pos = p0 + pp
				
			case _:
				print("Specify the mode for calculus (0/1) \n")
				
		
		self.Xs[step] = self.pos[0]
		self.Ys[step] = self.pos[1]

		#ora devi copiare tutto il resto da evolvidT per raccogliere i dati
		self.L[step] = self.massa * np.linalg.norm(np.cross(self.pos, self.vel))
		self.modE(step, planets)
		self.V[step] = np.linalg.norm(self.vel)	
		
		sole = planets[0]
		terra = planets[3]
		sp = sole.pos
		ds=self.pos-sp
		dSole=np.linalg.norm(ds)
		
		E=0
		if self.nome=="Sole":
			E=self.Ek[step]
		else:
			Epot=-G*sole.massa*self.massa/dSole
			E=self.Ek[step]+Epot
		
		alfa = G * sole.massa * self.massa
		mr=sole.massa*self.massa/(self.massa+sole.massa)
		den = alfa * alfa * mr
		Ls=np.cross(ds,self.vel) * self.massa #per calcolare l'eccentricità devo usare il momento angolare rispetto al sole
		h2  = np.linalg.norm(Ls) * np.linalg.norm(Ls)
		num2 = 2 * h2 * E
		
		self.ecc[step] = math.sqrt(1+num2/den)  #eccentricità
		self.DS[step] = dSole
		
		tp=terra.pos
		dtt=tp-sp
		vt=terra.vel
		nt=np.cross(dtt,vt)
		if self.nome!="Luna" and self.nome!= "Sole":
			self.teta[step] = 90- angolo(ds, nt)
		elif self.nome == "Luna":
			lt=self.pos-tp;
			self.teta[step] = 90- angolo(lt, nt)
		
	####
	
	def selplot(self, num, dt):
		match num:
			case 0:
				self.plotL()			
			case 1:
				self.plotE()
			case 2:
				self.plotecc()
			case 3:
				self.plotteta()
			case 4:
				self.plotst(dt)
			case 5:
				self.plotxy()
			case 6:
				self.plotdv()
			case 7:
				self.plotD()
			case 8:
				self.plotsd(dt)
			case 9:
				self.plotV()
			case _:
				print("indice non valido")
	
	#####

	def plotL(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		st = "Media: "+str(self.L.mean())+"\n Dev Std: "+str(np.std(self.L))
		plt.hist(self.L, bins=400, range=(self.L.mean()-np.std(self.L)*4, self.L.mean()+np.std(self.L)*4), label = st)
		plt.title(self.nome+": momento angolare [kg*m^2/s]")
		plt.xlabel("L [kg*m^2/s]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()
		
	#####

	def plotE(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		E=self.Ek+self.Ep
		st = "Media: "+str(E.mean())+"\n Dev Std: "+str(np.std(E))
		plt.hist(E, bins=400, range=(E.mean()-np.std(E)*5, E.mean()+np.std(E)*5), label = st)
		plt.title(self.nome+": energia meccanica [J]")
		plt.xlabel("E [J]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()
	
	#####

	def plotecc(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		st = "Media: "+str(self.ecc.mean())+"\n Dev Std: "+str(np.std(self.ecc))
		#if self.nome != "Sole":
		plt.hist(self.ecc, bins=400, range=(self.ecc.mean()-np.std(self.ecc)*5, self.ecc.mean()+np.std(self.ecc)*5), label = st)
		"""else: 
			plt.hist(self.ecc, bins=100, label = st)	"""
		plt.title(self.nome+": eccentricita'")
		plt.xlabel("e")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()
		
	#####

	def plotteta(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		st = "Moda: " + str(sta.mode(self.teta))
		plt.hist(self.teta, bins=100, label = st)		
		plt.title(self.nome+": inclinazione orbita'")
		plt.xlabel("angolo [°]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()

	#####
	
	def plotst(self, dt, outfile = None):
		fig = plt.figure(figsize=(9,9))
		xx = np.arange(0, self.teta.size*dt, dt)
		xx = np.divide(xx, 365*24*3600)
		plt.plot(xx, self.teta, label = self.nome)
		plt.title(self.nome+": inclinazione in fz del tempo")
		plt.xlabel("tempo [anni]")
		plt.ylabel("angolo [°]")
		plt.legend()
		plt.show()
	
	#####
	
	def plotxy(self, outfile = None):
		"""Disegno xy"""
		fig = plt.figure(figsize=(9,9))
		plt.plot(self.Xs, self.Ys, label = self.nome)
		plt.title(self.nome+": traiettoria")
		plt.xlabel("x [m]")
		plt.ylabel("y [m]")
		plt.legend()
		plt.show()
		
	#####
	
	def plotdv(self, outfile = None):
		"""Disegno pv"""
		fig = plt.figure(figsize=(9,9))
		plt.plot(self.DS, self.V, label = self.nome)
		plt.title(self.nome+": |V| vs |distanza dal sole|")
		plt.xlabel("|d| [m]")
		plt.ylabel("|v| [m/s]")
		plt.legend()
		plt.show()
		
	#####

	def plotD(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		st = "Moda: "+str(sta.mode(self.DS))
		plt.hist(self.DS, bins=100, label = st)
		plt.title(self.nome+": distanza dal sole [m]")
		plt.xlabel("distanza [m]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()	
		
	#####
	
	def plotsd(self, dt, outfile = None):
		fig = plt.figure(figsize=(9,9))
		xx = np.arange(0, self.DS.size*dt, dt) / (365*24*3600)
		plt.plot(xx, self.DS, label = self.nome)
		plt.title(self.nome+": distanza dal sole in fz del tempo")
		plt.xlabel("tempo [anni]")
		plt.ylabel("distanza [m]")
		plt.legend()
		plt.show()
	
	#####

	def plotV(self, outfile = None):
		fig = plt.figure(figsize=(8,8))
		st = "Moda: "+str(sta.mode(self.V))
		plt.hist(self.V, bins=100, label = st)
		plt.title(self.nome+": |v| [m/s]")
		plt.xlabel("|v| [m/s]")
		plt.ylabel("Conteggi")
		plt.legend()
		plt.show()
	
	######
	#per stampare i copri
	def __str__(self):
		return "Posizione di "+self.nome+": "+str(self.pos)+"m"
