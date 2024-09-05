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

void corpo::muovi(std::vector<corpo*> cc, unsigned int dt, uint32_t mode){
	switch(mode){
		case 0: //eulero 0
		{
			vettore a=acc(cc);
			vettore dv=a*dt;
			m_vel=m_vel+dv;
			vettore dp=m_vel*dt;
			m_pos=m_pos+dp;
			break;
		}
		case 1:
		{  //eulero modificato
			vettore v0=m_vel;
			vettore a0=acc(cc);
			vettore p0=m_pos;
			vettore dp=v0*dt;
			m_pos=m_pos+dp;
			vettore a2=acc(cc);
			vettore vv=((a0 + a2)*0.5)*dt; //quì medi accelerazioe e velocità
			m_vel = m_vel + vv;
			vettore pp=((m_vel + v0)*0.5)*dt;
			m_pos=p0+pp;
			break;
		}
		case 2:
		{  //velocity verlet
			vettore v0=m_vel;
			vettore a0=acc(cc);
			vettore dp= v0*dt;
			vettore dp2 = a0*dt*dt/2;
			m_pos=m_pos+dp+dp2;
			vettore a2=acc(cc);
			vettore vv=((a0 + a2)*0.5)*dt;  //quì fai media tra acc iniziale e finale
			m_vel = m_vel + vv;
			break;
		}
		case 4:
		{  //runge-kutta ordine 2 eulero-richardson sballato
			vettore a0=acc(cc);
			vettore dp= m_vel*dt/2;
			m_pos=m_pos+dp;
			vettore dv=a0*dt/2;
			vettore v1=m_vel+dv;
			vettore a2=acc(cc);
			dp = v1*dt;
			m_pos=m_pos + dp;
			vettore vv= a2*dt; //quì usi accelerazione calcolata in punto intermedio
			m_vel = m_vel + vv;
			break;

		  //runge-kutta ordine 2 eulero-richardson - è ancora sballato, la luna va a puttane
			vettore a=acc(cc); //accelerazione iniiale
			vettore dv=a*dt/2;
			vettore v12 = m_vel+dv; //velocità intermedia
			vettore dp= m_vel*dt;
			vettore dp2 = a*dt*dt/2;
			vettore p1=m_pos+dp+dp2; //posizione finale approssimata
			dp= m_vel*dt/2;
			dp2 = a*dt*dt/8;
			m_pos=m_pos+dp+dp2; //posizione intermedia	- serveper calcolare a
			a=acc(cc); //accelerazione intermedia		
			dp=v12*dt/2;
			dp2 = a*dt*dt/8;
			vettore p2=m_pos+dp+dp2; //posizione finale approssimata 2
			m_pos=p2*2-p1; //posizioone finale
			dv=a*dt;
			m_vel=m_vel+dv; //velocità finale non son molto convinto ma ok
			break;
			
			//sempre lo stesso ma non va ancora
			vettore k1v=acc(cc)*dt/2;
			vettore k1x = m_vel*dt/2;
			vettore p0=m_pos;
			m_pos=m_pos+k1x;
			vettore k2v=acc(cc)*dt;
			vettore k2x=(m_vel+k1v)*dt;
			m_pos=p0+k2x;
			m_vel=m_vel+k2v;
		}
		case 5:
		{  //Yoshida 4th order
			//step1
			vettore dp=m_vel*dt*W1/2;
			m_pos = m_pos+dp;
			vettore a=acc(cc);
			vettore dv= a*dt*W1;
			m_vel=m_vel+dv;
			//step2
			dp = m_vel*dt*(W0+W1)/2;
			m_pos= m_pos+dp;
			a = acc(cc);
			dv = a*dt*W0;
			m_vel= m_vel+dv;
			//step 3
			dp = m_vel*dt*(W0+W1)/2;
			m_pos= m_pos+dp;
			a = acc(cc);
			dv = a*dt*W1;
			m_vel= m_vel+dv;
			//step 4
			dp = m_vel*dt*W1/2;
			m_pos= m_pos+dp;	
			break;
		}
		case 3:
		{  //runge-kutta ordine 4
			//step 1
			vettore k1v=acc(cc)*dt/2;
			vettore k1x=m_vel*dt/2;
			vettore p0=m_pos;
			vettore v0=m_vel;
			//step 2
			m_pos=p0+k1x;
			vettore k2v=acc(cc)*dt/2;
			m_vel=v0+k1v;
			vettore k2x=m_vel*dt/2;
			//step 3
			m_pos=p0+k2x;
			vettore k3v=acc(cc)*dt;
			m_vel=v0+k2v;
			vettore k3x=m_vel*dt;
			//step 4
			m_pos=p0+k3x;
			vettore k4v=acc(cc)*dt;
			m_vel=v0+k3v;
			vettore k4x=m_vel*dt;			
			//finale
			k1v=k1v/3;
			k2v=k2v*4/6;
			k3v=k3v/3;
			k4v=k4v/6;
			m_vel=v0+k1v+k2v+k3v+k4v;
			k1x=k1x/3;
			k2x=k2x*4/6;
			k3x=k3x/3;
			k4x=k4x/6;
			m_pos=p0+k1x+k2x+k3x+k4x;
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
			k1v=k1v/2;
			k2v=k2v/2;
			m_vel=v0+k1v+k2v;
			k1x=k1x/2;
			k2x=k2x/2	;
			m_pos=p0+k1x+k2x;						
			break;
		}
		default:
			std::cout<<"invaild mode";
			break;
	}
}

void corpo::evolvidT(std::vector<corpo*> cc, unsigned int dt, uint32_t mode, uint32_t j){
	corpo *sole=cc[0];
	corpo *terra=cc[3];
	if(m_nome=="Sole" && j<2) m_app=vettore(0,0,0);
	else m_app=m_pos0;	
	m_pos0=m_pos; //devo salvare la posizione che il corpo aveva prima dell'evoluzione, soprattutto per il sole, perché serve per il calcolo della precessione
	
	muovi(cc,dt,mode);
	
	//seleziono sole per raccoliere dati rispetto a lui
	vettore sp=sole->P();
	//float sy=sp.y();
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
  	
  	vettore Ls=ds*m_vel*m_massa; //per calcolare l'eccentricità devo usare il momento angolare rispetto al sole, non rispetto all'origine
	double h2  = Ls.modulo()*Ls.modulo();
	//if(j%131400==0) std::cout<<h2<<std::endl;
	
	/*
	double h = m_L.modulo()*m_L.modulo();
	double num = 2 * h * Emec;
	double e = sqrt(1+num/den);  //eccentricità1
	//if(j%131400==0) std::cout<<e<<std::endl;
	*/
	
	double num2 = 2 * h2 * E;
	double e2 = sqrt(1+num2/den);  //eccentricità2
	//if(j%131400==0) std::cout<<e2<<std::endl;
	
	/*
	//calcolo eccentricità valutando L non rispetto al sole, ne all'origine, ma rispetto all'origine traslata nelle y come il sole
	//così va solo pk nella sim tutto trasla verso l'alto, ma dipende dai dati iniziali
	vettore yy(0,sy,0);
	vettore dc=m_pos-yy;
	vettore LL=dc*m_vel*m_massa;
	double h3=LL.modulo() * LL.modulo();
	double num3 = 2 * E * h3;
	double e3 = sqrt(1+num3/den);  //eccentricità3
	std::cout<<yy<<std::endl;	*/

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
//	m_histos[2]->Fill( m_pos.modulo(), m_vel.modulo() ); // Vel vs dist
	m_histos[2]->Fill( dSole , m_vel.modulo() ); 		// Vel vs dist dal sole
	m_histos[3]->Fill(m_vel.modulo());                  // Modulo velocità
	m_histos[4]->Fill( m_L.modulo() );                           // Momento angolare
	m_histos[5]->Fill( e2 );                              // Eccentricità modo 2
	m_histos[6]->Fill( m_teta );                              // inclinazione orbita
	m_histos[8]->Fill( E );                           // Enercia solo sole 
	m_histos[9]->Fill( Emec );                           // Enercia meccanica
	
//	m_histos[8]->Fill( m_pos.modulo() );                 // Dist da orig
//	m_histos[9]->Fill( m_pos.x() );                      // X
//	m_histos[10]->Fill( m_pos.y() );                      // Y
//	m_histos[11]->Fill( d0 );                    // Dist da pos iniz
//	m_histos[12]->Fill( e );                              // Eccentricità
//	m_histos[7]->Fill( e3 );                              // Eccentricità modo 3

//	m_graps[0]->AddPoint( j, dSole );  
	//std::cout<<m_pos<<std::endl;
	
	//raccoglo dati perielio
	if(m_nome!="Sole"){
		/*
		float d_cfr=(m_app-sole->P0()).modulo();
		if(dSole<d_cfr) m_app=m_pos;
		uint64_t Tstep = m_TT *24*3600 / dt ;
		if((j+1)%Tstep == 0){
			if(m_nome=="Mercurio") std::cout<<m_app-sp<<std::endl;
			m_peri.push_back(m_app-sp);
			m_app=m_pos0; //riinizializzo il vettore di confronto
		}
		*/
		
		vettore sp0=sole->P0();
		vettore sap=sole->AP();
		float d_media=(m_pos0-sp0).modulo();
		float d_pre=(m_app-sap).modulo();
		if(d_media<d_pre && d_media<dSole){
			m_peri.push_back(m_pos0-sp0);
			if(m_nome=="Mercurio"){
				std::cout<<d_pre<<" "<<d_media<<" "<<dSole<<std::endl;
				std::cout<<m_pos0-sp0<<std::endl;
			}
		}
		
	}
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
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";Precessione [arcsec/anno];Conteggi").c_str(), numBins, 0, 60 ) ) ); 
   
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
 
  //l'istogramma dell'enegia mecanica completa lo aggiubngo dopo che ho messo tutti i pianeti
 
  /* grafico d/t
  s = m_nome + ": distanza dal sole in fz del tempo(step)";
  m_graps.push_back(
    reinterpret_cast<TGraph*> ( new TGraph() ) );
    m_graps[0]->SetName(s.c_str());
    m_graps[0]->SetTitle(s.c_str());
 
   /*
     // Histo 0
  s = m_nome + ": Distanza da (0,0)";
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), s.c_str(), numBins, d*4/5, d*6/5) ) );
   
   
  // Histo 2
  s = m_nome + ": X";
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), s.c_str(), numBins, -d*3/2, d*3/2) ) );

  // Histo 3
  s = m_nome + ": Y";
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), s.c_str(), numBins, -d*3/2, d*3/2) ) );
   
  // Histo 4
  s = m_nome + ": dist. pos. iniziale";
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), s.c_str(), numBins, -d*0.01, 2.1*d) ) );
  
    // Histo 11
  s = m_nome + ": eccentricita'"; //NB e' negativa!!!
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), s.c_str(), numBins, 0, 1.1 ) ) );
   
  // Histo 12
  s = m_nome + ": eccentricita' terzo modo"; //NB e' negativa!!!
  m_histos.push_back(
    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), s.c_str(), 2*numBins, 0, 0.25 ) ) );
  */

}

//questo va richiamato nel main, dopo aver aggiuntob tutti i corpi - Histo 12
void corpo::istEmec(std::vector<corpo*> cc){
	double E=m_Ek+m_Ep;
	std::cout<<E<<std::endl;
	int numBins = 400;
	std::string s=m_nome + ": energia meccanica";
	m_histos.push_back(
	    reinterpret_cast<TH1I*> ( new TH1I(s.c_str(), (s+";E [J];Conteggi").c_str(), numBins, E*51/50, E*49/50) ) );
}

void corpo::fillgraf(std::string aa, std::string bb){
  /*/ grafico d/t
	std::string s = m_nome + ": distanza dal sole in fz del tempo(step)";
	std::string file=m_nome+"_dist_sole.txt";
  	TGraph * lol = new TGraph((aa+file).c_str());
    lol->SetName(s.c_str());
    lol->SetTitle(s.c_str());
    m_graps.push_back(
    	reinterpret_cast<TGraph*> ( lol ) );
    */
    	
    // grafico d/t
	std::string s = m_nome + ": distanza dal sole in fz del tempo(anni)";
	std::string file="dist_sole_"+aa;
  	TGraph * lol = new TGraph(file.c_str(), bb.c_str());
    lol->SetName(s.c_str());
    s+=";Step;Distanza (m)";
    lol->SetTitle(s.c_str());
    m_graps.push_back(
    	reinterpret_cast<TGraph*> ( lol ) );

    //grafico inmc/t
    s = m_nome + ": incli8nazione orbita in fz del tempo(anni)";
	file="incl_"+aa;
  	TGraph * gtg = new TGraph(file.c_str(), bb.c_str());
    gtg->SetName(s.c_str());
    s+=";Step;Inclinazione [gradi]";
    gtg->SetTitle(s.c_str());
    m_graps.push_back(
    	reinterpret_cast<TGraph*> ( gtg ) );
}

void corpo::precessione(float Tterra){
	if(m_nome=="Sole")m_histos[7]->Fill(0);
	else{
		for(int i=1; i<m_peri.size(); i++){
			if(m_nome=="Mercurio") std::cout<<m_peri[i].angolo(m_peri[i-1])*Tterra/m_TT<<std::endl;
			m_histos[7]->Fill(m_peri[i].angolo(m_peri[i-1])*3600*Tterra/m_TT);
		}
	}
}
