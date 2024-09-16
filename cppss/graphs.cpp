#include "corpo.h"
#include<vector>
#include<cstring>
#include<iostream>
#include<cmath>
#include "sistema.h"

#include <TApplication.h>
#include<TCanvas.h>
#include <TFile.h>
#include<TStyle.h>
#include<TLine.h>

#define CLRSCREEN printf("\e[1;1H\e[2J");

void linee(TH1I *h);

int main(int argc, char** argv){

	if(argc!=3){
		std::cerr << "Usage: configFile datafile\n"; 
		return 1;
	}
	std::string confFile = argv[1];
	std::string outFile = argv[2];
	
	TApplication myApp("App", &argc, argv);
	TCanvas *screen2 = new TCanvas("c1", "My Solar System", 0, 0, 900, 900);
	
	sistema ss;
	ss.modSTR(outFile);
	ss.leggi(confFile);
	
	std::string cmd;  // Comando
	int val;          // Valore associato al comando
	CLRSCREEN;
	while(true)
	  {
	    ss.print();
	    std::cout << "============================== \n";
	    std::cout << "Inserisci: stringa intero \n";
	    std::cout << "  Esempi: \n";
	    std::cout << "    Pianeta n   // cioè Pianeta indiceGrafico (max+1 per creare i grafici)\n";
	    std::cout << "    sis 0/1  	  // cioè sistema indiceGrafico \n";
	    std::cout << "    list 0      // lista grafici disponibili  \n";
	    std::cout << "    graf 0      // crea tutti i grafici  \n";	    
	    std::cout << "    root 0      // per passare a root            \n";
	    std::cout << "    quit 0      // Uscire dal programma          \n";
	    std::cout << "______\n";
	    std::cout << "prompt>";
	    
	    // Legge i valori
	    std::cin >> cmd >> val;

			// Parsing
	    if(cmd=="quit"){ break;}
	    else if(cmd=="list"){ 
	    	CLRSCREEN;
	    	ss.PrintHistos();
	    	continue;
	    }
	    else if(cmd=="graf"){
	    	CLRSCREEN;
	    	ss.mkGraf(); //creo grafici d/t
	    	//system(RM);
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
  			gStyle->SetOptStat(111111);
			std::string fx="histo"+outFile;
			TFile f(fx.c_str(), "read");
			std::string title = cmd + std::to_string(val);
			TH1I *h = f.Get<TH1I>(title.c_str());
			if(h==NULL)
			  		{ std::cerr << "Istogramma non trovato\n"; continue; }
			screen2->Clear();
			screen2->SetWindowSize(900, 900);
			h->GetXaxis()->SetNdivisions(4, 2, 0, kFALSE);
			h->SetFillColor(41);
			h->Draw();
			//linee(h);
			screen2->Modified();    
			screen2->Update();
			CLRSCREEN;
			continue;
	    }
	    
	//altrimenti di default CERCA ISTOGRAMMA RELATIVO AD UN PIANETA
		else{
			//cmd = cmd.substr(0,3);
			int a=ss.select(cmd, val);
			if(a==-1){
  				gStyle->SetOptStat(111111);
				std::string fx="histo"+outFile;
				TFile f(fx.c_str(), "read");
				std::string title = cmd + std::to_string(val);
				TH1I *h;
				if(val == 1 || val ==2) h = reinterpret_cast<TH1I*>(f.Get<TH2I>(title.c_str()));
				else h = f.Get<TH1I>(title.c_str());
				if(h==NULL){			
					std::cerr << "Pianeta non riconosciuto\n"; 
					continue;	
				}
				screen2->Clear();
				screen2->SetWindowSize(900, 900);
				h->SetFillColor(41); //41
				h->Draw();
				//if(val==4 || val==8 || val==9) linee(h);
				screen2->Modified();    
				screen2->Update();
				CLRSCREEN;
				continue;
			}
			else if(a==-2){
				CLRSCREEN;
				ss.mkGraf(cmd);
				continue;
			}
			else{	
				//std::cout<<a<<std::endl;
				TGraph *h2 = ss.getThisGraph(cmd, val-a);
				if(h2==NULL){			
					std::cerr << "Pianeta non riconosciuto\n"; 
					continue;				
				}
				screen2->Clear();
				screen2->SetWindowSize(1600, 900);
				h2->Draw();
				screen2->Modified();    
				screen2->Update();
				CLRSCREEN;		
				continue;	
			}

	    }
	    
	  }
	
	return 0;
}
void linee(TH1I *h){
	TLine *l= new TLine(h->GetMean(),0,h->GetMean(),h->GetMaximum());
	l->SetLineColor(kRed);
	l->SetLineWidth(3);
	l->Draw();
	TLine *l4= new TLine(h->GetMean()+h->GetRMS(),0,h->GetMean()+h->GetRMS(),h->GetMaximum());			
	l4->SetLineColor(kBlue);
	l4->SetLineStyle(kDashed);
	l4->SetLineWidth(4);
	l4->Draw();
	TLine *l5= new TLine(h->GetMean()-h->GetRMS(),0,h->GetMean()-h->GetRMS(),h->GetMaximum());			
	l5->SetLineColor(kBlue);
	l5->SetLineStyle(kDashed);
	l5->SetLineWidth(4);
	l5->Draw();
	TLine *l2= new TLine(h->GetMean()+3*h->GetRMS(),0,h->GetMean()+3*h->GetRMS(),h->GetMaximum());			
	l2->SetLineStyle(kDashDotted);
	l2->SetLineWidth(3);
	l2->Draw();
	TLine *l3= new TLine(h->GetMean()-3*h->GetRMS(),0,h->GetMean()-3*h->GetRMS(),h->GetMaximum());			
	l3->SetLineStyle(kDashDotted);
	l3->SetLineWidth(3);
	l3->Draw();
}
