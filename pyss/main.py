from sistema import sistema
from corpo import CorpoCeleste
from funzioni import *
import sys

if len(sys.argv) != 6:
	print("Usage: configFile numeroAnni granularità mode outputfile")
	print("   configFile: file di configurazione delle condizioni iniziali")
	print("   numeroAnni: anni di evoluzione")
	print("   granularità: tempo dT di evoluzione (secondi)")
	print("   mode:       modalità di evoluzione")
	print("               0 = standard")
	print("               1 = con media dell'accelerazione ")
	exit()
		
confFile = sys.argv[1]
nA = int(sys.argv[2])
ddt = int(sys.argv[3])
mode  = int(sys.argv[4])
outFile = sys.argv[5]

listaCorpi = leggi(confFile) #####

#deltaSecondi = 3600*1               # Granularità in s
durata       = 3600*24*365*nA       # Durata totale in s
steps        = durata//ddt
raggioTerra  = 6e6 # metri
spinta       = 3397

ss = sistema(listaCorpi, ddt, steps, outFile)

#ss.add(CorpoCeleste("saturno",    1.90e27, [0. ,1352.5e9], [-10180., 0.]))
ss.print()

ss.simulazione(mode)
ss.save(outFile)

#parse
while True:
	ss.print();
	print("==============================")
	print("Inserisci: stringa intero")
	print("  Esempi:")
	print("    pianeta n   // cioè Pianeta indiceIstogramma (sis per info globali)")
	print("    list 0      // lista istogrammi disponibili")
	print("    out  0      // stampa dati istogrammi")
	print("    evo  n      // evolvi per altri n anni")
	print("    quit 0      // Uscire dal programma")
	print("______")
	print("prompt>")
	
	#Legge i valori
	cmd = input()
	val = int(input())
	
	if cmd == "quit":
		exit()
	
	elif cmd == "evo":
		ss.modD(val)
		ss.simulazione(mode)
	
	elif cmd == "list":
		print("boh") ########
		
	elif cmd == "out":
		print("boh") ######
		
	elif cmd == "sis":
		ss.selectplot(val)
		
	else:
		ss.selectplanet(cmd, val)

"""
listaCorpi = [
    CorpoCeleste("sole",     1.98e30, [0.,      0.], [     0., 0.]),
    CorpoCeleste("mercurio", 3.30e23, [0.,  4.6e10], [-58309., 0.]),
    CorpoCeleste("venere"  , 4.87e24, [0., 1.07e11], [-35267., 0.]),
    CorpoCeleste("terra",    5.97e24, [0., 1.47e11], [-30293., 0.]),
    CorpoCeleste("luna",     7.34e22, [0., 146.7E9], [-29290., 0.]),   
    CorpoCeleste("marte",    6.42e23, [0., 2.06e11], [-26457., 0.]),
    CorpoCeleste("giove",    1.90e7, [0. ,7.41e11], [-13712., 0.]),
    CorpoCeleste("saturno",    1.90e27, [0. ,1352.5e9], [-10180., 0.]),
    CorpoCeleste("urano",    8.68e25, [0. ,2.74e12], [ -7119., 0.]),
    CorpoCeleste("nettuno",  1.02e26, [0. ,4.44e12], [ -5495., 0.]),
    #CorpoCeleste("sonda", 1e6, [0. ,1.47e11+raggioTerra], [-30293. - spinta, 0.]),
]
"""
