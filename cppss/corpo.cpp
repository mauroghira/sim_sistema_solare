#include "vec.h"
#include<cstring>
#include<vector>
#include<cmath>
#include "corpo.h"
#include<iostream>
#include<iomanip>
#include<TH1I.h>
#include <TH2I.h>
#include<TGraph.h>
#include<stdlib.h>
#include<fstream>

std::ostream& operator<<(std::ostream& o, corpo& c)
{
  vettore p = c.P()/1.E3; 
  vettore v = c.V()/1.E3; 
  o << std::setw(10) << c.NAME() << "\t" << p << " km\t" << v << " km/s";
  return o;
}

void corpo::leggi(){
	std::cout<<"nome corpo: ";
	std::cin>>m_nome;
	std::cout<<"massa copro: ";
	do{
		std::cin>>m_massa;
		if(m_massa<=0)std::cout<<"massa positiva!"<<std::endl;
	}while(m_massa<=0);
	std::cout<<"inserire la posizione iniziale del corpo: ";
	m_pos.leggi();
	m_pos0=m_pos;
	m_app=m_pos;
	m_sap=vettore(0,0,0);
	m_s0=vettore(0,0,0);
	std::cout<<"inserire la velocita' iniziale del corpo: ";
	m_vel.leggi();
	std::cout<<"inserire l'inclinazione del corpo: ";
	std::cin>>m_teta;
	m_L=m_pos*m_vel*m_massa;
	m_Ek= 0.5*m_massa*m_vel.modulo()*m_vel.modulo(); 
	m_Ep=0;	  //all'nizio la metto a zero, poi la modifico nel main quando ho aggiunto tutti i corpi del sistema
	inizia();
}
corpo::corpo(std::string n, double m, vettore r, vettore v, float Torb, float t){
	m_nome=n;
	m_massa=m;
	m_vel=v;
	m_pos=r;
	m_teta=t;
	m_pos0=r;
	m_app=r;
	m_sap=vettore(0,0,0);
	m_s0=vettore(0,0,0);
	m_TT=Torb;
	m_L=m_pos*m_vel*m_massa;
	m_Ek= 0.5*m_massa*m_vel.modulo()*m_vel.modulo(); 
	m_Ep=0;
  	inizia();	
}
corpo::corpo(){
	m_nome="";
	m_massa=0;
	m_teta=0;
	m_Ek=0;
	m_Ep=0;
	m_TT=0;
}

void corpo::ass(std::string n, double m, vettore r, vettore v, float Torb, float t){
	m_nome=n;
	m_massa=m;
	m_vel=v;
	m_pos=r;
	m_teta=t;
	m_pos0=r;
	m_app=r;
	m_s0=vettore(0,0,0);
	m_sap=vettore(0,0,0);
	m_TT=Torb;
	m_L=m_pos*m_vel*m_massa;
	m_Ek= 0.5*m_massa*m_vel.modulo()*m_vel.modulo();
	m_Ep=0;
	inizia();
}

//metodo pr computare l'energia nella configurazione iniaizle una volta aggiunti uttti i corpi e per fare l'evoluzione
void corpo::modE(std::vector<corpo*> cc){
  	m_Ep=0; //devo riazzerare l'energia ogni ovolta che la ricalcolo, altrimenti continuo a sommare tutto oidocrop
  	for(int i=0; i<cc.size(); i++){
  		if(cc[i]->NAME()!=m_nome){
  			vettore dd=cc[i]->P();
  			vettore d=m_pos-dd;
  			double D=(double)d.modulo();
  			m_Ep+= -G * cc[i]->MASS() * m_massa /D;
  		}
  	}
  	m_Ek=0.5*m_massa*m_vel.modulo()*m_vel.modulo();
  	//std::cout<<m_Ek<<", "<<m_Ep<<std::endl;
}

vettore corpo::acc(std::vector<corpo*> &cc){
	vettore a;
	float k=0;
	for(int i=0; i<cc.size(); i++){
		if(cc[i]->m_nome!=m_nome){
			vettore d=cc[i]->m_pos-m_pos;
			if(d.modulo()!=0) k=G*cc[i]->m_massa/pow(d.modulo(), 3);
			vettore A=d*k;
			a=a+A;
		}
	}
	return a;
}

void corpo::muovi(std::vector<corpo*> cc, unsigned int dt, uint32_t mode){
	switch(mode){
		case 0:
		{   //eulero 0	
			m_vel=m_vel + acc(cc)*dt;
			m_pos= m_pos + m_vel*dt;
			break;
		}
		case 1:
		{  //eulero modificato
			vettore v0=m_vel;
			vettore a0=acc(cc);
			vettore p0=m_pos;
			m_pos=m_pos+v0*dt;
			vettore a2=acc(cc);
			m_vel = m_vel + (a0 + a2)*0.5*dt;
			m_pos= p0 + (m_vel + v0)*0.5*dt;
			break;
		}
		case 2:
		{  //velocity verlet
			vettore a0=acc(cc);
			m_pos= m_pos + m_vel*dt + a0*dt*dt/2;
			vettore a2=acc(cc);
			m_vel = m_vel + (a0 + a2)*0.5*dt;
			break;
		}
		case 3:
		{  //runge-kutta ordine 4
			//step 1
			vettore k1v=acc(cc)*dt;
			vettore k1x=m_vel*dt;
			vettore p0=m_pos;
			vettore v0=m_vel;
			//step 2
			m_pos=p0+k1x/2;
			vettore k2v=acc(cc)*dt;
			m_vel=v0+k1v/2;
			vettore k2x=m_vel*dt;
			//step 3
			m_pos=p0+k2x/2;
			vettore k3v=acc(cc)*dt;
			m_vel=v0+k2v/2;
			vettore k3x=m_vel*dt;
			//step 4
			m_pos=p0+k3x;
			vettore k4v=acc(cc)*dt;
			m_vel=v0+k3v;
			vettore k4x=m_vel*dt;			
			//finale
			m_vel=v0+(k1v+k2v*2+k3v*2+k4v)/6;
			m_pos=p0+(k1x+k2x*2+k3x*2+k4x)/6;
		}
		case 4:
		{  //runge-kutta 2
			//step 1
			vettore k1v=acc(cc)*dt;
			vettore k1x=m_vel*dt;
			vettore p0=m_pos;
			vettore v0=m_vel;
			//step 2
			m_pos=p0+k1x;
			vettore k2v=acc(cc)*dt;
			m_vel=v0+k1v;
			vettore k2x=m_vel*dt;
			//finale
			m_vel=v0+k1v/2+k2v/2;
			m_pos=p0+k1x/2+k2x/2;						
			break;
		}
		case 5:
		{  //Yoshida 4th order
			//step1
			m_pos= m_pos+ m_vel*dt*W1/2;
			m_vel= m_vel+ acc(cc)*dt*W1;
			//step2
			m_pos= m_pos+ m_vel*dt*(W0+W1)/2;
			m_vel= m_vel+ acc(cc)*dt*W0;
			//step 3
			m_pos= m_pos+ m_vel*dt*(W0+W1)/2;
			m_vel= m_vel+ acc(cc)*dt*W1;
			//step 4
			m_pos= m_pos+ m_vel*dt*W1/2;	
			break;
		}
		default: 
			std::cout<<"invaild mode";
			break;
	}
}

void corpo::evolvidT(std::vector<corpo*> cc, unsigned int dt, uint32_t mode, uint64_t j){
	corpo *sole=cc[0];
	corpo *terra=cc[3];
		
	/*
	if(m_nome=="Sole" && j<2){
		//per sole m_sap non dovrebbe importare
		m_app=vettore(0,0,0);
	}
	else{
		m_app=m_pos0;
		m_sap=m_s0;
	}	
	m_s0=sole->P0();
	m_pos0=m_pos; //devo salvare la posizione che il corpo aveva prima dell'evoluzione, soprattutto per il sole, perché serve per il calcolo della precessione
	//*/
	
	//sposto corpo
	muovi(cc, dt, mode);
	
	//seleziono sole per raccoliere dati rispetto a lui
	vettore sp=sole->P();
	vettore ds=m_pos-sp;
	double dSole=ds.modulo();
	double d0=(m_pos-m_pos0).modulo();
	
	//aggiorno l'energia e il mojenot
	modE(cc);
  	double Emec = m_Ek + m_Ep;
  	//std::cout<<Ecin<<" "<<Epot<<" "<<Emec<<" ";
	m_L=m_pos*m_vel*m_massa;
	
	//calcolo solo 'energia rispetto al sole, senza contare altri corpi, per valutare meglio l'eccentricità
	double E=0;
	if(m_nome=="Sole")E=m_Ek;
	else{
		double Epot=-G*sole->MASS()*m_massa/dSole;
		E=m_Ek+Epot;
	}
	
  	double alfa = G * sole->MASS() * m_massa;
  	double mr=sole->MASS()*m_massa/(m_massa+sole->MASS());
  	double den = alfa * alfa * mr;
  	
  	//vettore Ls=ds*m_vel*m_massa; //per calcolare l'eccentricità devo usare il momento angolare rispetto al sole, non rispetto all'origine
	double h2  = m_L.modulo()*m_L.modulo();

	double num2 = 2 * h2 * E;
	double e2 = sqrt(1+num2/den);  //eccentricità2

	//valuto inclinazione
	vettore tp=terra->P();
	vettore dtt=tp-sp;
	vettore vt=terra->V();
	vettore nt=dtt*vt;
	if(m_nome!="Luna") m_teta=90-ds.angolo(nt);
	else{
		vettore lt=m_pos-tp;
		m_teta=90-lt.angolo(nt);
	}

	// Riempimento degli istogrammi
	m_histos[0]->Fill( dSole );                    // Dist dal Sole
	m_histos[1]->Fill( m_pos.x(), m_pos.y() );           // Traiettoria   
	//m_histos[2]->Fill( m_pos.modulo(), m_vel.modulo() ); // Vel vs dist
	m_histos[2]->Fill( dSole , m_vel.modulo() ); 		// Vel vs dist dal sole
	m_histos[3]->Fill(m_vel.modulo());                  // Modulo velocità
	m_histos[4]->Fill( m_L.modulo() );                           // Momento angolare
	m_histos[5]->Fill( e2 );                              // Eccentricità modo 2
	m_histos[6]->Fill( m_teta );                              // inclinazione orbita
	m_histos[8]->Fill( E );                           // Enercia solo sole  
	m_histos[9]->Fill( Emec );                           // Enercia meccanica

	//raccoglo dati perielio
	if(m_nome!="Sole"){
		//*
		float d_cfr=(m_app-m_sap).modulo();
		if(dSole<=d_cfr){
			m_app=m_pos;
			m_sap=sp;
		}
		uint64_t Tstep = m_TT *24*3600 / dt ;
		if((j+1)%Tstep == 0){
			//if(m_nome=="Mercurio") std::cout<<m_app<<m_sap<<m_app-m_sap<<std::endl;
			m_peri.push_back(m_app-m_sap);
			m_app=m_pos; //riinizializzo il vettore di confronto
			m_sap=sp;
		}
		//*/
		//van bene ambo i modi - vantaggio di questo è che posso vedere le step a cui lo becco
		/*
		float d_media=(m_pos0-m_s0).modulo();
		float d_pre=(m_app-m_sap).modulo();
		if(d_media<d_pre && d_media<dSole){
			m_peri.push_back(m_pos0-m_s0);
			//if(m_nome=="Mercurio"){
				//std::cout<<d_pre<<" "<<d_media<<" "<<dSole<<std::endl;
				//std::cout<<j<<" "<<m_pos0-m_s0<<std::endl;
			//}
		}
		//*/
	}
}

void corpo::precessione(float Tterra){
	if(m_nome=="Sole")m_histos[7]->Fill(0);
	else{
		for(int i=1; i<m_peri.size(); i++){
			//if(m_nome=="Mercurio") std::cout<<m_peri[i].angolo(m_peri[i-1])*Tterra/m_TT<<std::endl;
			m_histos[7]->Fill(m_peri[i].angolo(m_peri[i-1])*Tterra/m_TT);
		}
		//if(m_nome=="Mercurio") for(auto p: m_peri) std::cout<<p<<std::endl;
	}
}

TH1I* corpo::getisto(uint32_t index)
{ 
  //if( index >= m_histos.size()+m_graps.size() )
    //return NULL;            // Ritorna NULL se l'indice è fuori range
  if(index >= m_histos.size())
  	return NULL;
  	//return m_graps[index-m_histos.size()];
  return m_histos[index];   // Altrimenti l'elemento a indice = index
}

TGraph* corpo::getgraf(uint32_t index)
{ 
  if( index >= m_graps.size() )
    return NULL;            // Ritorna NULL se l'indice è fuori range
  return m_graps[index];   // Altrimenti l'elemento a indice = index
}

//crea gli istogrammmi
void corpo::inizia(){
  std::cout << "Inizializzo gli istogrammi di " << m_nome << std:: endl;
  int numBins = 400;
  std::string s;
 
  // Variabili utili per riempire gli istogrammi
  double d    = (double)m_pos.modulo();
  float v    = m_vel.modulo();

  // NB: le frazioni usate nei costruttori qui sotto servono per centrare
  //     l'istogramma rispetto al valore atteso - commento gli istogrammi inutili per simulazioni lunghe
  
  // Histo 0
  s = m_nome + ": Distanza dal Sole";
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";Distanza [m];Conteggi").c_str(), numBins, d*3/5, d*115/100) ) );

  // Histo 1
  s = m_nome + ": X vs Y";
  if(m_nome!="Sole") m_histos.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";x [m];y [m]").c_str(), numBins, -d*3/2, d*3/2,
                                                         numBins, -d*3/2, d*3/2) ) );
  else m_histos.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";x [m];y [m]").c_str(), numBins, -1e10, 1e10,
                                                         numBins, -5e11, 1e3) ) );  
   
  // Histo 2
  s = m_nome + ": |vel| vs |dist dal Sole|";
  if(m_nome!="Sole") m_histos.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";Distanza [m];|v| [m/s]").c_str(), numBins, d*60/100, d*115/100,  //mod lim in base a dati
                                                         numBins, v*85/100, v*155/100) ) );
  else m_histos.push_back( reinterpret_cast<TH1I*>
                        ( new TH2I(s.c_str(), (s+";Distanza [m];|v| [m/s]").c_str(), numBins, 0, 5e9,
                                                         numBins, 0, 5) ) );                                                
                                                         
  // Histo 3
  s = m_nome + ": |V|";
  if(m_nome=="Sole"){  
	  m_histos.push_back(
	    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";|v| [m/s];Conteggi").c_str(), numBins, 0, 5 ) ) );
  }
  else m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";|v| [m/s];Conteggi").c_str(), numBins, v*4/5, v*8/5 ) ) );
    
  // Histo 4
  s = m_nome + ": momento angolare (r x mv)";
  double L_att=m_L.modulo();
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";|L| [kg*m^2/s];Conteggi").c_str(), numBins, L_att*9/10, L_att*11/10) )
  );
  m_histos[4]->GetXaxis()->SetNdivisions(4, 2, 0, kFALSE);
  
  // Histo 5
  s = m_nome + ": eccentricita' secondo modo"; //NB e' negativa!!!
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";Eccentricita';Conteggi").c_str(), 2*numBins, 0, 0.25 ) ) );

  // Histo 6
  s = m_nome + ": inclinazione orbita"; //NB e' negativa!!!
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";Inclinazione [gradi];Conteggi").c_str(), 2*numBins, -m_teta*3/2, m_teta*3/2 ) ) );
  
  //Histo 7
  s = m_nome + ": precessione"; //NB e' negativa!!!
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";Precessione [arcsec/anno];Conteggi").c_str(), numBins, 0, 1 ) ) );  
    
  // Histo 8    prima inizializzo soo l'istograma rispetto al sole, poi metto i lsto dopo aver aggiunto tutto 
  s = m_nome + ": energia meccanica solo col sole"; 
  if(m_nome!="Sole"){  
	  double Epot=-G*1.98e30*m_massa/(double)m_pos.modulo();
	  double E=m_Ek+Epot; 
	  m_histos.push_back(
	    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, E*51/50, E*49/50) ) );
  }    
  else{
	  m_histos.push_back(
	    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, 0, 1e33 ) ) );
  } //per il sole rispetto a sé considero solo l'energia cineitca
  m_histos[8]->GetXaxis()->SetNdivisions(4, 2, 0, kFALSE);
  
  //l'istogramma dell'enegia mecanica completa lo aggiubngo dopo che ho messo tutti i pianeti
  
}

//questo va richiamato nel main, dopo aver aggiuntob tutti i corpi - Histo 12
void corpo::istEmec(std::vector<corpo*> cc){
	double E=m_Ek+m_Ep;
	std::cout<<E<<std::endl;
	int numBins = 400;
	std::string s=m_nome + ": energia meccanica";
	m_histos.push_back(
	    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, E*51/50, E*49/50) ) );
	m_histos[9]->GetXaxis()->SetNdivisions(4, 2, 0, kFALSE);

  	//aggiungo quì inizializzazione di m_sap per la precessione
  	m_sap=cc[0]->P0();
    m_s0=cc[0]->P0();	
  	//std::cout<<m_sap<<std::endl;
}

void corpo::fillgraf(std::string aa, std::string bb){
    // grafico d/t
	std::string s = m_nome + ": distanza dal sole in fz del tempo(anni)"; 
	std::string file="dist_sole_"+aa;
  	TGraph * lol = new TGraph(file.c_str(), bb.c_str());
    lol->SetName(s.c_str());
    s+=";Anni;Distanza (m)";
    lol->SetTitle(s.c_str());
    m_graps.push_back(
    	reinterpret_cast<TGraph*> ( lol ) );

    //grafico inmc/t
    s = m_nome + ": inclinazione orbita in fz del tempo(anni)"; 
	file="incl_"+aa;
  	TGraph * gtg = new TGraph(file.c_str(), bb.c_str());
    gtg->SetName(s.c_str());
    s+=";Anni;Inclinazione [gradi]";
    gtg->SetTitle(s.c_str());
    m_graps.push_back(
    	reinterpret_cast<TGraph*> ( gtg ) );
}
