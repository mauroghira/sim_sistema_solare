#include "vec.h"
#include "corpo.h"
#include<vector>
#include<cstring>
#include<iostream>
#include<cmath>
#include "sistema.h"
#include<fstream>
#include <TCanvas.h>
#include <TView.h>
#include <TMarker.h>
#include <time.h>
#include <typeinfo>
#include <TFile.h>

void sistema::input(std::string config, std::string of){
	do{
	std::cout<<"tempo simulazione? (anni) ";
	std::cin>>m_T;
	}while(m_T<=0);
	std::cout<<"dT evoluzione? (secondi) ";
	std::cin>>m_dT;
	m_inc=of;
	leggi(config);
	ins();
}
sistema::sistema(float T, unsigned int dt, std::string config, std::string of){
	m_T=T;
	m_dT=dt;
	m_inc=of;
	leggi(config);
	ins();
}
sistema::sistema(){
	m_T=0;
	m_dT=0;
	m_inc=std::to_string(m_T)+"_"+std::to_string(m_dT)+"_";	
}

void sistema::leggi(std::string config){
	std::ifstream in(config);
	if( in.fail() )
	  {
		std::cerr << "File " << config << " non trovato\n";
		exit(3);
	  }

	std::string nome;
	double mas;
	int i=0;
	std::size_t found= config.find("sfe");
	
	if(found!=std::string::npos){
		float r, p, t, v, T;
		while(true){
			in>>nome>>mas>>r>>p>>t>>v>>T;
			if(in.eof())break;
			if(nome[0]!='#'){
				vettore pp(r, p, t);
				vettore vv(0-v*sin(p*M_PI/180), v*cos(M_PI*p/180), 0);
				if(nome!="Luna"){
					pp.sfToCar();
				}
				else{
					double rt=m_corpi[3]->P().modulo();
					double rl=abs(r-rt);
					vettore dp(rl, p, t);
					dp.sfToCar();
					pp = m_corpi[3]->P() + dp;
				}
				corpo *e = new corpo(nome,mas,pp,vv,T,90-t);
				add(e);
				i++;
			}
		}
	}	
	else{
		float px, py, pz, vx, vy,vz, a, T;
		while(true){
			in>>nome>>mas>>px>>py>>pz>>vx>>vy>>vz>>a>>T;
			if(in.eof())break;
			if(nome[0]!='#'){
				double xx=0;
				if(nome!="Luna"){
					xx=px*cos(a*M_PI/180);
					pz=px*sin(a*M_PI/180);
				}
				else{
					double xt=m_corpi[3]->P().x();
					double x1=abs(px-xt);
					xx=xt+x1*cos(a*M_PI/180);
					pz=x1*sin(a*M_PI/180);
				}
				vettore v(vx,vy,vz), p(xx,py,pz);
				corpo *e = new corpo(nome,mas,p,v,T,a);
				add(e);
				i++;
			}
		}
	}
	in.close();
	for(int i=0; i<m_corpi.size(); i++){
		m_corpi[i]->modE(m_corpi);
		m_corpi[i]->istEmec(m_corpi);
		//m_corpi[i]->getisto(12)->SetBinsLength(m_T*365*24*3600/m_dT+2);
	}
	print();
}

void sistema::add(corpo* c){
	m_corpi.push_back(c);
}

void sistema::evo(uint32_t mode, int st){

	unsigned long int nn=365*24*3600/m_dT;
	unsigned long int n=nn*m_T;
	for(uint64_t i=0; i<n; i++){
		evodt(mode, i+st);
		if((i+st)%(nn*5)==0){
			print();	//stampa ogni 5 ani
			std::cout<<"anno "<<(i+st)/nn<<std::endl;
		}
	}
			
	//per ora lo metto quì ma sarà poco efficiente
	for(auto p: m_corpi) p->precessione(m_corpi[3]->period());
	
}
void sistema::evodt(uint32_t mode, uint64_t j){
	for(int i=0; i<m_corpi.size(); i++){
		m_corpi[i]->evolvidT(m_corpi, m_dT, mode, j);
	}

	//*dati in file v2, così più veloce perché apre solo 1 file, inoltre velocizza campionando ogni 1000
	if(j%10000==0){
		vettore dS=m_corpi[0]->P();
		std::ofstream out;
		out.open("dist_sole_"+m_inc, std::ofstream::app);
		out<<float(j)*m_dT/(365*24*3600);
		for(auto c: m_corpi){
			vettore dd=c->P()-dS;
			float d=(float)dd.modulo();
			out<<" "<<d;
		}
		out<<std::endl;
		out.close();
	}
	
	if(j%10000==0){
		std::ofstream out;
		out.open("incl_"+m_inc, std::ofstream::app);
		out<<float(j)*m_dT/(365*24*3600);
		for(auto c: m_corpi){
			out<<" "<<c->incl();
		}
		out<<std::endl;
		out.close();
	}
	
	/*
	if(float(j)*m_dT/(365*24*3600) >= 499){
	vettore dS=m_corpi[0]->P();
	std::ofstream out;
	out.open("val_err_"+m_inc, std::ofstream::app);
	for(auto c: m_corpi){
		vettore dd=c->P()-dS;
		float d=(float)dd.modulo();
		out<<d<<" ";
	}
	out<<std::endl;
	out.close();
	}//*/
	
  	vettore L;
  	for(int i=0; i<m_corpi.size(); i++){
  		vettore dL=m_corpi[i]->LA();
  		L=L+dL;
  	}
  	double L_att=(double)L.modulo();
  	m_ist[0]->Fill(L_att);
	
	double Emec=0;
  	for(int i=0; i<m_corpi.size(); i++)
  		Emec+=m_corpi[i]->ECIN()+m_corpi[i]->EPOT()/2;
  	//std::cout<<Emec<<std::endl;
  	m_ist[1]->Fill(Emec);
}

void sistema::print(){
	std::cout<<"=============================="<<std::endl;
	for(int i=0; i<m_corpi.size(); i++){
		std::cout<<"Posizione di "<<m_corpi[i]->NAME()<<"\t"<<m_corpi[i]->P()<<"\t";
		//std::cout<<m_corpi[i]->ECIN()<<"\t"<<m_corpi[i]->EPOT()<<"\t"<<m_corpi[i]->EMEC();
		//std::cout<<m_corpi[i]->incl();
		std::cout<<std::endl;
	}
}

TH1I* sistema::getThisHisto(std::string nomeCorpo, uint32_t indice){
	nomeCorpo = nomeCorpo.substr(0,3); // Cerchiamo solo le prime 3 lettere
  	for(auto p : m_corpi)
  		if (p->NAME().find(nomeCorpo)==0)
  			return p->getisto(indice);
  	return NULL;
}

TGraph* sistema::getThisGraph(std::string nomeCorpo, uint32_t indice){
	nomeCorpo = nomeCorpo.substr(0,3); // Cerchiamo solo le prime 3 lettere
  	for(auto p : m_corpi)
  		if (p->NAME().find(nomeCorpo)==0)
  			return p->getgraf(indice);
  	return NULL;
}

void sistema::PrintHistos(){
	corpo *p = m_corpi[1];
  	std::cout << std::noshowpos;
  	for(int i=0; i<p->numHistos(); i++)
  		std::cout << "\t" << i << "\t-> " << p->getisto(i)->GetName() << std::endl;  //istogrammi di ogni corpo
  	for(int i=0; i<p->numgraf(); i++)
  		std::cout << "\t" << i+p->numHistos() << "\t-> " << p->getgraf(i)->GetName() << std::endl;  //istogrammi di ogni corpo
  	for(int i=0; i<m_ist.size(); i++){
  		std::cout<<"\t" << i << "\t-> " << getist(i)->GetName() << std::endl;    //istogrammi del sistema
  	}
}

TH1I* sistema::getist(uint32_t i){
	if( i >= m_ist.size() ) 
	    	return NULL;            // Ritorna NULL se l'indice è fuori range
  	return m_ist[i];   // 
}

//creare gli istogrammi dell'energia e del momento angolare totale
void sistema::ins(){
  int numBins = 400;
  std::string s; 

  // Histo 0
  s = "sistema: momento angolare (r x mv)";
  vettore L;
  for(int i=0; i<m_corpi.size(); i++){
  	vettore dL=m_corpi[i]->LA();
  	L=L+dL;
  }
  double L_att=L.modulo();
  m_ist.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";|L| [kg*m^2/s];Conteggi").c_str(), numBins, L_att*49999/50000, L_att*50001/50000) )   //va bene, devi aumentare pèerché strettissimo, ma funzioona
  );

  // Histo 1
  s ="sistema: energia meccanica"; //NB e' negativa!!!
  double Emec=0;
  for(int i=0; i<m_corpi.size(); i++)
  	Emec+=m_corpi[i]->ECIN()+m_corpi[i]->EPOT()/2;  //NON DEVO CONTARE L'EN POT DUE VOLTE
  std::cout<<"E0 "<<Emec<<std::endl;	
  m_ist.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, Emec*100001/100000, Emec*99999/100000) )
  );
}

int sistema::select(std::string p, int v){
	p=p.substr(0,3);
	for(int i=0; i<m_corpi.size(); i++){
		if(m_corpi[i]->NAME().find(p) == 0){
			if(v >= m_corpi[i]->numHistos() + m_corpi[i]->numgraf()) return -2;
			else if (v >= m_corpi[i]->numHistos()) return m_corpi[i]->numHistos();
		}
	}
	return -1;
}

void sistema::output(std::string file){
	std::ofstream out(file);
	
 	out<<"dati || media | dev std"<<std::endl;
	for(auto p: m_corpi){
		out<<"------------------------------"<<std::endl;
		for(int i=0; i<p->numHistos(); i++){
	 		out<<p->getisto(i)->GetName()<<" || "<<p->getisto(i)->GetMean(1)<<" | "<<p->getisto(i)->GetRMS(1)<<"\t"<<std::endl;
	 		if(i==1||i==2) out<<" || "<<p->getisto(i)->GetMean(2)<<" | "<<p->getisto(i)->GetRMS(1)<<std::endl;
	 	}
		for(int i=0; i<p->numgraf(); i++){
	 		out<<p->getgraf(i)->GetName()<<" || "<<p->getgraf(i)->GetMean(2)<<" | "<<p->getgraf(i)->GetRMS(2)<<std::endl;
	 	}
	 	out<<std::endl;
	}
	out<<"=============================="<<std::endl;
	out<<"dati generali"<<std::endl;
	for(auto h: m_ist){
	 	out<<h->GetName()<<" || "<<h->GetMean(1)<<" | "<<h->GetRMS(1)<<std::endl;		
	}	
	 
	out.close();
}

void sistema::savehist(std::string out){
	std::string fx="histo"+out;
	TFile* f(TFile::Open(fx.c_str(), "recreate"));
	for(auto P: m_corpi){
		for(int i=0; i<P->numHistos(); i++){
			std::string title = P->NAME().substr(0,3) + std::to_string(i);;
			f->WriteObject(P->getisto(i), title.c_str());
		}
	}
	for(int i=0; i<m_ist.size(); i++){
		std::string title = "sis" + std::to_string(i);;
		f->WriteObject(m_ist[i], title.c_str());
	}
	delete f;
}

void sistema::mkGraf(std::string pla){
	std::string a="%lg";
	std::string c=" %lg";
	std::string b="";
	if (pla==""){
		for(auto p: m_corpi){
			if(p->numgraf()==0){
				p->fillgraf(m_inc, a+b+c);
			}
			else std::cout<<"grafici di "<<p->NAME()<<" gia creati"<<std::endl;
			b+=" %*lg";
		}
		std::cout<<"creo tutti i grafici"<<std::endl;
	}
	else{
		pla = pla.substr(0,3); // Cerchiamo solo le prime 3 lettere
  		for(int p=0; p<m_corpi.size(); p++)
  			if (m_corpi[p]->NAME().find(pla)==0){
  				if(m_corpi[p]->numgraf()==0){
	  				for(int i=0; i<p; i++)
	  					b+=" %*lg";
	  				m_corpi[p]->fillgraf(m_inc, a+b+c);
	  				std::cout << "creo grafici di "<<m_corpi[p]->NAME()<<"\n"; 
	  			}
	  			else std::cout<<"grafici di "<<m_corpi[p]->NAME()<<" gia creati"<<std::endl;		
	  			pla="";
	  			break;
  			}
		if(pla!=""){
			std::cerr<<"Pianeta non riccnosciuto"<<std::endl;
		}
	}
}

void sistema::clean(){	
	std::ofstream out("dist_sole_"+m_inc);
	out.close();
	out.open("incl_"+m_inc);
	out.close();
	out.open("val_err_"+m_inc);
	for(auto p: m_corpi) out<<p->NAME()<<" ";
	out<<std::endl;
	out.close();	
}
