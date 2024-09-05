#include "corpo.h"
#include<vector>
#include<cstring>
#include<iostream>
#include<cmath>
#include "sistema.h"

#include <TApplication.h>
#include<TCanvas.h>
#include <TFile.h>

#define CLRSCREEN printf("\e[1;1H\e[2J");

#define CAT "\
#/bin/bash \n\
cat $(ls -Art *.txt | tail -n 1) \n\
"

#define RM "\
#/bin/bash \n\
rm *_dist_sole.txt *incl.txt \n\ 
"
//così non va bene perché sputtana tutto se faccio più sim insieme

int main(int argc, char** argv){

	if(argc!=6){
		std::cerr << "Usage: configFile numeroAnni granularità mode outputfile\n";
		std::cerr << "   configFile: file di configurazione delle condizioni iniziali\n";
		std::cerr << "   numeroAnni: anni di evoluzione\n";
		std::cerr << "   granularità: tempo dT di evoluzione (secondi)\n";		
		std::cerr << "   mode:       modalità di evoluzione\n";
		std::cerr << "               0 = standard \n";
		std::cerr << "               1 = con media dell'accelerazione \n";
		return 2;
	}
	std::string confFile = argv[1];
	float nAnni = atof(argv[2]);
	uint32_t ddt = atoi(argv[3]);
	uint32_t mode  = atoi(argv[4]);
	std::string outFile = argv[5];
	int step=0;
	
	TApplication myApp("App", &argc, argv);
	TCanvas *screen2 = new TCanvas("c1", "My Solar System", 0, 0, 900, 900);

	sistema s(nAnni, ddt, confFile, outFile);
	
	/*
	uint32_t mode;
	std::cout<<"modalità calcolo accelerazione? ";
	std::cin>>mode;
	std::string config;
	std::cout<<"file di input? ";
	std::cin>>config;
	sistema s;
	s.input(config);
	*/
	
	/*corpo X;
	X.leggi();
	corpo Y;
	Y.leggi();
	s.add(&X);
	s.add(&Y);*/
	
	s.clean();
	s.evo(mode);
	//s.print();
	//system(CAT);
	s.savehist(outFile);
	s.output(outFile);
	
	// Loop infinito per poter scegliere gli istogrammi
	std::string cmd;  // Comando
	int val;          // Valore associato al comando
	CLRSCREEN;
	while(true)
	  {
	    s.print();
	    std::cout << "============================== \n";
	    std::cout << "Inserisci: stringa intero \n";
	    std::cout << "  Esempi: \n";
	    std::cout << "    Pianeta n   // cioè Pianeta indiceIstogramma (sis per info globali) \n";
	    std::cout << "    list 0      // lista istogrammi disponibili  \n";
	    std::cout << "    graf 0      // crea i grafici fz/tempo  \n";
	    std::cout << "    out  0      // stampa dati istogrammi      \n";
	    std::cout << "    log  0/1    // scala logaritmica si/no       \n";
	    std::cout << "    evo  n      // evolvi per altri n anni       \n";
	    std::cout << "    root 0      // per passare a root            \n";
	    std::cout << "    quit 0      // Uscire dal programma          \n";
	    std::cout << "______\n";
	    std::cout << "prompt>";
	   
	    // Legge i valori
	    std::cin >> cmd >> val;

			// Parsing
	    if(cmd=="quit"){ break;}
	    else if(cmd=="log") {
	    	screen2->SetLogy(val);
	    	CLRSCREEN;
	    	continue;
	    }
	    else if(cmd=="list"){
	    	CLRSCREEN;
	    	s.PrintHistos();
	    	continue;
	    }
	    else if(cmd=="out"){
	    	CLRSCREEN;
	    	system(CAT);
	    	continue;
	    }
	    else if(cmd=="graf"){
	    	CLRSCREEN;
	    	s.mkGraf(); //creo grafici d/t
	    	//system(RM);
	    	continue;
	    }
	    else if(cmd=="evo") {
	    	step+= nAnni*365*24*3600/ddt;
	    	s.modT((float)val);
	    	s.evo(mode, step);  
			s.output(outFile);
			nAnni=val;
	    	CLRSCREEN;
	    	continue;
	    }
	    else if(cmd=="root"){
	      std::cout << "\n\n\t\t --> Comando alla finestra di root\n";
	      std::cout << "    \t\t     File->Quit per tornare al terminale\n";       
	      myApp.Run(kTRUE);
	      CLRSCREEN;
	      std::cout << "Uscito dalla finestra di root\n";
	      continue;
	    }
	    else if(cmd=="sis"){
			TH1I *h = s.getist(val);
			if(h==NULL)
			  		{ std::cerr << "Istogramma non trovato\n"; continue; }
			screen2->Clear();
			screen2->SetWindowSize(900, 900);
			h->SetFillColor(41);
			h->Draw();
			screen2->Modified();   
			screen2->Update();
			CLRSCREEN;
			continue;
	    }
	   
	//altrimenti di default CERCA ISTOGRAMMA RELATIVO AD UN PIANETA
		else{
			//cmd = cmd.substr(0,3);
			int a=s.select(cmd, val);
			
			if(a==-1){
				TH1I *h = s.getThisHisto(cmd, val);
				if(h==NULL){			
					std::cerr << "Pianeta non riconosciuto\n";
					continue;				
				}
				screen2->Clear();
				screen2->SetWindowSize(900, 900);
				h->SetFillColor(41); //41
				h->Draw();
				screen2->Modified();   
				screen2->Update();
				CLRSCREEN;
				continue;
			}
			else if(a==-2){			
				std::cerr << "Pianeta non riconosciuto\n";
				continue;
			}
			else{			
				//std::cout<<a<<std::endl;
				TGraph *h2 = s.getThisGraph(cmd, val-a);
				if(h2==NULL){
					std::cerr << "Pianeta non riconosciuto\n";
					continue;				
				}
				screen2->Clear();
				screen2->SetWindowSize(1600, 900);
				//h2->GetXaxis()->SetTitle("time");
				//h2->GetYaxis()->SetTitle("dist");
				//h->SetColor(1); //41*/
				h2->Draw();
				screen2->Modified();   
				screen2->Update();
				CLRSCREEN;		
				continue;	
			}
		
		/*	TH1I *h = s.getThisHisto(cmd, val);
			if(h==NULL)
			  { std::cerr << "Pianeta non riconosciuto\n"; continue; }

			h->SetFillColor(1); //41
			h->Draw();
			screen2->Modified();   
			screen2->Update();
			CLRSCREEN;
		*/
		
	    }
	   
	  }
	
	return 0;
}
