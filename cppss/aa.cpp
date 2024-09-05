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

#define SHELLSCRIPT "\
#/bin/bash \n\
echo -e \"\" \n\
echo -e \"This is a test shell script inside C code!!\" \n\
read -p \"press <enter> to continue\" \n\
"

int main(int argc, char** argv){
	TApplication myApp("App", &argc, argv);
	TCanvas *screen2 = new TCanvas("c1", "My Solar System", 0, 0, 800, 800);
	std::string s =": distanza dal sole in fz del tempo(step)"; 
	TGraph* h = new TGraph();
    h->SetName(s.c_str());
    h->SetTitle(s.c_str());
    for(int i=0;  i<10; i++){
    	h->AddPoint(i, i+1);
    }
    
    TH1I* t = new TH1I("lol", "lol", 10, 0, 10);
    for(int i=0;  i<10; i++){
    	t->Fill(i, i+1);
    }    
    
	while(true){
		std::string cmd;  // Comando
		std::cin >> cmd;

		if(cmd=="h"){
			screen2->Clear();
			t->SetFillColor(41);
			t->Draw();
			screen2->Modified();    
			screen2->Update();
		}
		else if(cmd=="g"){
			screen2->Clear();
			h->Draw();
			screen2->Modified();    
			screen2->Update();
		}
		else if(cmd=="q") break;
		
    }

	system(SHELLSCRIPT);
 
	return 0;
}
