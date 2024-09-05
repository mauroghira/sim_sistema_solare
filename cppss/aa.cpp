#include<cstring>
#include<vector>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<TH1I.h>
#include <TH2I.h>
#include<TGraph.h>
#include<stdlib.h>
#include <TApplication.h>
#include<TCanvas.h>

#include<TPaveText.h>

#define SHELLSCRIPT "\
#/bin/bash \n\
echo -e \"\" \n\
echo -e \"This is a test shell script inside C code!!\" \n\
read -p \"press <enter> to continue\" \n\
"

void hlabels2();

int main(int argc, char** argv){
	TApplication myApp("App", &argc, argv);
	TCanvas *screen2 = new TCanvas("c1", "My Solar System", 0, 0, 800, 800);
	std::string s =": distanza dal sole in fz del tempo(step)";
    TH1I* h = new TH1I("lol", "lol", 10, -10, 10);
    //*
    for(int i=-9;  i<10; i++){
    	h->Fill(i, i+1);
        	std::cout<<i;
    }
    /*
	TAxis* a = h->GetXaxis();
	a->SetNdivisions(-502);
	a->ChangeLabel(1,-1,-1,-1,-1,-1,"-#pi");
	a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
   //*/
   
h->GetXaxis()->SetNdivisions(9, 5, 0, kFALSE);

   
	while(true){
		std::string cmd;  // Comando
		std::cin >> cmd;

		if(cmd=="h"){
			screen2->Clear();
			h->SetFillColor(41);
			h->Draw();
			screen2->Modified();   
			screen2->Update();
		}
		
		else if(cmd=="a"){
			hlabels2();
			
		}
		
		else if(cmd=="q") break;
		
    }

	system(SHELLSCRIPT);

	return 0;
}

 
void hlabels2()
{
   const Int_t nx = 12;
   const Int_t ny = 20;
   const char *month[nx]  = {"January","February","March","April","May",
      "June","July","August","September","October","November",
      "December"};
   const char *people[ny] = {"Jean","Pierre","Marie","Odile","Sebastien",
      "Fons","Rene","Nicolas","Xavier","Greg","Bjarne","Anton",
      "Otto","Eddy","Peter","Pasha","Philippe","Suzanne","Jeff",
      "Valery"};
   TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,600,600);
   c1->SetGrid();
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);
   TH2F *h = new TH2F("h","test",3,0,3,2,0,2);
   h->SetCanExtend(TH1::kAllAxes);
   h->SetStats(0);
   for (Int_t i=0;i<15000;i++) {
      Int_t rx = random()*nx;
      Int_t ry = random()*ny;
      h->Fill(people[ry],month[rx],1);
   }
   h->LabelsDeflate("X");
   h->LabelsDeflate("Y");
   h->LabelsOption("v");
   h->Draw("text");
 
   TPaveText *pt = new TPaveText(0.6,0.85,0.98,0.98,"brNDC");
   pt->SetFillColor(18);
   pt->SetTextAlign(12);
   pt->AddText("Use the axis Context Menu LabelsOption");
   pt->AddText(" \"a\"   to sort by alphabetic order");
   pt->AddText(" \">\"   to sort by decreasing values");
   pt->AddText(" \"<\"   to sort by increasing values");
   pt->Draw();
  	c1->Modified();   
	c1->Update();
	
}
