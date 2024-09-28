#include "vec.h"
#include "corpo.h"
#include<vector>
#include<cstring>
#include<iostream>
#include<iomanip>
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
sistema::sistema(float T, unsigned int dt, std::string config, std::string of, int xx){
	m_T=T;
	m_dT=dt;
	m_inc=of;
	leggi(config);
	ins(xx);
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

void sistema::evo(uint32_t mode, int xx, int st){

	unsigned long int nn=365*24*3600/m_dT;
	unsigned long int n=nn*m_T;
	for(uint64_t i=0; i<n; i++){
		evodt(mode, i+st, xx);
		if((i+st)%(nn*5)==0){
			print();	//stampa ogni 5 ani
			std::cout<<"anno "<<(i+st)/nn<<std::endl;
		}
	}
			
	//per ora lo metto quì ma sarà poco efficiente
	for(auto p: m_corpi) p->precessione(m_corpi[3]->period());
	
}
void sistema::evodt(uint32_t mode, uint64_t j, int xx){
	for(int i=1; i<m_corpi.size(); i++){
		m_corpi[i]->evolvidT(m_corpi, m_dT, mode, j);
	}
	m_corpi[0]->modE(m_corpi);
  	m_corpi[0]->getisto(9)->Fill(m_corpi[0]->EMEC());      
	
  	vettore L;
  	for(int i=0; i<m_corpi.size(); i++){
  		vettore dL=m_corpi[i]->LA();
  		L=L+dL;
  	}
  	double L_att=(double)L.modulo();
  	m_ist[0]->Fill(L_att);
	
	double Ek=0;
  	double Ep=0;
  	for(int i=0; i<m_corpi.size(); i++){
  		Ek+=m_corpi[i]->ECIN();  		
  		Ep+=m_corpi[i]->EPOT()/2;
  	}
  	//std::cout<<Emec<<std::endl;
  	m_ist[1]->Fill(Ek+Ep);
  	m_ist[4]->Fill(Ep, Ek);
  	m_ist[2]->Fill(Ek);
   	m_ist[3]->Fill(Ep);
  	m_ist[6]->Fill(Ep, Ep+Ek);
  	m_ist[5]->Fill(Ek, Ep+Ek);
  	
	vettore dS=m_corpi[0]->P()-m_corpi[xx]->P();
  	m_ist[7]->Fill(dS.modulo(), Ep+Ek);
  	vettore v=m_corpi[xx]->V();
  	m_ist[8]->Fill(v.modulo(), Ep+Ek);
  	
	std::ofstream out;
	out.open("energy_"+m_inc, std::ofstream::app);
	out<<float(j)*m_dT/(365*24*3600)<<" "<<std::setprecision(10)<<Ep/2<<" "<<-Ek<<" "<<Ep+Ek;
	out<<std::endl;
	out.close();
	
	out.open("dist_sole_"+m_inc, std::ofstream::app);
	out<<float(j)*m_dT/(365*24*3600)<<" "<<dS.modulo()<<std::endl;
	out.close();
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

void sistema::PrintHistos(int xx){
	corpo *p = m_corpi[xx];
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
	if( i > m_ist.size() ) 
	   	return NULL;            // Ritorna NULL se l'indice è fuori range
  	return m_ist[i];   // 
}

//creare gli istogrammi dell'energia e del momento angolare totale
void sistema::ins(int xx){
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
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, Emec*90001/90000, Emec*189999/190000) )
  );

    double Ek=0;
    double Ep=0;
    for(auto p: m_corpi){
    	Ek+=p->ECIN();
		Ep+=p->EPOT()/2;
    }
  
  // Histo 2
  s ="sistema: energia cinetica"; //NB e' negativa!!!
  m_ist.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, Ek*4/5, Ek*26/25) )
  );

  // Histo 3
  s ="sistema: energia potenziale"; //NB e' negativa!!!
  m_ist.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, Ep*41/40, Ep*9/10) )
  );  
 
    // Histo 4
  	s = "sistema: E_cinetica vs E_potenziale";
  	m_ist.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";E_pot [J];E_cin [J]").c_str(), numBins, Ep*41/40, Ep*9/10,
                                                         numBins, Ek*4/5, Ek*26/25) ) );
  
    // Histo 6
  	s = "sistema: E_meccanica vs E_cinetica";
  	m_ist.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";E_cin [J];E_mec [J]").c_str(), numBins, Ek*4/5, Ek*26/25,
                        								numBins, Emec*190001/190000, Emec*189999/190000) ) );
    // Histo 5
  	s = "sistema: E_meccanica vs E_potenziale";
  	m_ist.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";E_pot [J];E_mec [J]").c_str(), numBins, Ep*41/40, Ep*9/10,
                                                         numBins, Emec*190001/190000, Emec*189999/190000) ) );

	vettore dS=m_corpi[0]->P()-m_corpi[xx]->P();
	double d=dS.modulo();
    // Histo 7
  	s = "E meccanica VS distanza di "+m_corpi[xx]->NAME()+" dal Sole";
  	m_ist.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";distanza [m];E [J]").c_str(), numBins, d*99/100, d*115/100,
                                                         numBins, Emec*190001/190000, Emec*189999/190000) ) );

	vettore v=m_corpi[xx]->V();
	d=v.modulo();
    // Histo 7
  	s = "E meccanica VS velocita' di "+m_corpi[xx]->NAME();
  	m_ist.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";velocita' [m/s];E [J]").c_str(), numBins, d*85/100, d*101/100,
                                                         numBins, Emec*190001/190000, Emec*189999/190000) ) );

  for(auto h: m_ist) h->GetXaxis()->SetNdivisions(4, 2, 0, kFALSE);
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

void sistema::mkGraf(int xx){
	std::string a="%lg";
	std::string c=" %lg";
	if(m_corpi[xx]->numgraf()==0)
		m_corpi[xx]->fillgraf(m_inc, a+c);
	else std::cout<<"grafici di "<<m_corpi[xx]->NAME()<<" gia creati"<<std::endl;
}

void sistema::clean(){	
	std::ofstream out("energy_"+m_inc);
	out.close();
	out.open("dist_sole_"+m_inc);
	out.close();
}
