from sistema import sistema
from corpo import CorpoCeleste
from funzioni import *
import sys

if len(sys.argv) != 4:
	print("Usage: configFile numeroAnni outputfile")
	print("   configFile: file di configurazione delle condizioni iniziali")
	print("   numeroAnni: anni di evoluzione")
	exit()
		
confFile = sys.argv[1]
nA = int(sys.argv[2])
outFile = sys.argv[3]

listaCorpi = leggi(confFile) #####

ss = sistema(listaCorpi, 1, nA, outFile)

ss.print()

ss.cfr()
