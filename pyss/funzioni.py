import numpy as np
import math
import matplotlib.pyplot as plt
from corpo import CorpoCeleste as cc

def leggi(inFile):
	#funzione per leggere input da file
	lista = []
	f = open(inFile)
	for line in f:
		lista.append(line.strip().split(' ')) 

	corpi = []
	for i in lista:
		nome = i[0]
		massa = float(i[1])
		if inFile.find("sf")!=-1:
			d = float(i[2])
			phi = float(i[3])
			alfa = float(i[4])
			v = float(i[5])
		
			if nome[0]!="#":
				if nome != "Luna":
					pos, vel = sf(d, phi, alfa, v)
				else:
					pos, vel = mkluna(d, phi, alfa, v, corpi[1])  ##mod indice in base ai pianeti
		
				corpi.append(cc(nome, massa, pos, vel))
		
		else:
			x = float(i[2])
			y = float(i[3])
			z = float(i[4])
			vx = float(i[5])
			vy = float(i[6])
			vz = float(i[7])
			a = float(i[8])
		
			if nome[0]!="#":
				xx=0;
				if nome != "Luna":
					xx=x*math.cos(a*math.pi/180);
					z=x*math.sin(a*math.pi/180);
				else:
					xt=corpi[3].pos[0];
					x1=abs(x-xt);
					xx=xt+x1*math.cos(a*math.pi/180);
					z=x1*math.sin(a*math.pi/180);
				
				pos = [xx, y, z]
				vel = [vx, vy, vz]
		
				corpi.append(cc(nome, massa, pos, vel))			
		
	return corpi
	
####

def sf(d, f, a, v):
	f= f*math.pi/180
	t = (90-a)*math.pi/180
	pos = np.array([d*math.cos(f)*math.sin(t), d*math.sin(f)*math.sin(t), d*math.cos(t)])
	vel = [-v*math.sin(f), v*math.cos(f), 0]
	return pos, vel
	
####

def mkluna(d, f, a, v, terra):
	rt=np.linalg.norm(terra.pos);
	rl=abs(d-rt);
	pp, vv = sf(rl, f, a, v);
	pp += terra.pos
	return pp, vv

####

def modulo(a):
	s=0
	for i in a:
		s+=i**2
	return math.sqrt(s)
