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
#include<TMultiGraph.h>

#include <fstream>
#include <algorithm>

#define CLRSCREEN printf("\e[1;1H\e[2J");

void linee(TH1I *h);
void lin2D(TH1I *h);
void Egraf(std::string outFile);
TH1I* histo(std::string outFile);

int main(int argc, char** argv){

	if(argc!=4){
		std::cerr << "Usage: configFile planet datafile \n"; 
		return 1;
	}
	std::string confFile = argv[1];
	std::string outFile = argv[3];
	int xx = atoi(argv[2]);
	
	TApplication myApp("App", &argc, argv);
	TCanvas *screen2 = new TCanvas("c1", "My Solar System", 0, 0, 900, 900);
	
	sistema ss;
	ss.modSTR(outFile);
	ss.leggi(confFile);
	ss.ins(xx);
	
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
	    std::cout << "    sis  n  	  // cioè sistema indiceGrafico \n";
	    std::cout << "    list 0      // lista grafici disponibili  \n";
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
	    	ss.PrintHistos(xx);
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
	    	if(val==ss.numist()){
				screen2->Clear();
				screen2->SetWindowSize(1600, 900);
	      		Egraf(outFile);
	      		//screen2->BuildLegend(0.12, 0.1, 0.12, 0.1);
				screen2->Modified();    
				screen2->Update();
				CLRSCREEN;
	      		continue;
	    	}
	    	else{
	  			gStyle->SetOptStat(111111);
				std::string fx="histo"+outFile;
				TFile f(fx.c_str(), "read");
				std::string title = cmd + std::to_string(val);
				TH1I *h;
				if(val>ss.numist()) h=histo("energy_"+outFile);
				else if(val>3) h = reinterpret_cast<TH1I*>(f.Get<TH2I>(title.c_str()));
				else h = f.Get<TH1I>(title.c_str());
				if(h==NULL)
				  		{ std::cerr << "Istogramma non trovato\n"; continue; }
				screen2->Clear();
				screen2->SetWindowSize(900, 900);
				h->GetXaxis()->SetNdivisions(4, 2, 0, kFALSE);
				h->SetFillColor(41);
				h->Draw();
				if(val<4) linee(h);
				else lin2D(h);
				screen2->Modified();    
				screen2->Update();
				CLRSCREEN;
				continue;
			}
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
				if(val==4 || val==8 || val==9) linee(h);
				screen2->Modified();    
				screen2->Update();
				CLRSCREEN;
				continue;
			}
			else if(a==-2){
				CLRSCREEN;
				ss.mkGraf(xx);
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

TH1I* histo(std::string outFile){
	std::vector<double> data;
	std::ifstream in(outFile);
	if( in.fail() )
	  {
		std::cerr << "File " << outFile << " non trovato\n";
		exit(3);
	  }
	int i=0;
	double p,v,m,t;
	while(true){
		in>>t>>p>>v>>m;
		if(in.eof())break;
		data.push_back(m);	
		i++;
	}
	in.close();

	int nBins = 400; // Adjust the number of bins as needed
	double maxBin = *std::max_element(data.begin(), data.end())*199999/200000;
    double minBin = *std::min_element(data.begin(), data.end())*200001/200000;

	TH1I* histogram = new TH1I("Energia", "Energia", nBins, minBin, maxBin);

	for (auto value : data) {
		histogram->Fill(value);
	}
	return histogram;
}

void Egraf(std::string outFile){
	TMultiGraph *mg = new TMultiGraph();
	std::string s = "Energia meccanica in fz del tempo(anni)";
    mg->SetName(s.c_str());
    s+=";Anni;Energia (J)";
    mg->SetTitle(s.c_str());
	
	std::string file="energy_"+outFile;
  	std::string a="%lg";
	std::string c=" %lg";
	std::string b="";
/*  	
  	TGraph *gr1 = new TGraph(file.c_str(), (a+b+c).c_str());
  	gr1->SetTitle("E Pot/2");
   	gr1->SetLineColor(kRed);
   	gr1->SetFillStyle(kDashed);
*/	
	b+=" %*lg";
/*  	TGraph *gr2 = new TGraph(file.c_str(), (a+b+c).c_str());
  	gr2->SetTitle("-E Cin");
   	gr2->SetLineColor(kBlue);
   	gr2->SetFillStyle(kDashDotted);
  	mg->Add(gr2);
  	mg->Add(gr1);
*/	
	b+=" %*lg";
  	TGraph *gr3 = new TGraph(file.c_str(), (a+b+c).c_str());
  	gr3->SetTitle("E Mec");
  	mg->Add(gr3);
	
	mg->Draw("AL");
	
	TLine *l= new TLine(0,gr3->GetMean(2),500,gr3->GetMean(2));
	l->SetLineColor(kRed);
	l->SetLineWidth(3);
	l->Draw();
}

void lin2D(TH1I *h){
	double mx=h->GetMean(1);
	double my=h->GetMean(2);
	double lbx=h->GetXaxis()->GetBinLowEdge(1);
	double lby=h->GetYaxis()->GetBinLowEdge(1);
	double ubx=h->GetXaxis()->GetBinUpEdge(399);
	double uby=h->GetYaxis()->GetBinUpEdge(399);
	
	TLine *l= new TLine(mx,lby,mx,uby);
	l->SetLineColor(kRed);
	l->SetLineWidth(3);
	l->Draw();
	TLine *l6= new TLine(lbx,my,ubx,my);
	l6->SetLineColor(kRed);
	l6->SetLineWidth(3);
	l6->Draw();
	
	double rx=h->GetRMS(1);
	double ry=h->GetRMS(2);
	
	TLine *l5= new TLine(mx+rx,lby,mx+rx,uby);
	l5->SetLineColor(kBlue);
	l5->SetLineStyle(kDashed);
	l5->SetLineWidth(4);
	l5->Draw();
	TLine *l4= new TLine(mx-rx,lby,mx-rx,uby);
	l4->SetLineColor(kBlue);
	l4->SetLineStyle(kDashed);
	l4->SetLineWidth(4);
	l4->Draw();
	
	TLine *l3= new TLine(lbx,my+ry,ubx,my+ry);
	l3->SetLineStyle(kDashDotted);
	l3->SetLineWidth(3);
	l3->SetLineColor(kBlue);
	l3->Draw();
	TLine *l2= new TLine(lbx,my-ry,ubx,my-ry);
	l2->SetLineStyle(kDashDotted);
	l2->SetLineColor(kBlue);	
	l2->SetLineWidth(3);
	l2->Draw();
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
