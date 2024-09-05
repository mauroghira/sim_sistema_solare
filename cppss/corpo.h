#include "vec.h"
#include<cstring>
#include<vector>
#include<iostream>
#include<TH1I.h>
#include <TH2I.h>
#include<TGraph.h>
#include<stdlib.h>
#ifndef CORPO
#define CORPO

const double G=6.67e-11;
const float C=3e8;
const float RS=2.95e3;
const float BETA = 0;
const float W0 = - std::pow(2, 1/3)/(2-std::pow(2,1/3));
const float W1 = 1/(2-std::pow(2,1/3));

class corpo{
    protected:
	double m_massa;
	std::string m_nome;
	vettore m_pos;
	vettore m_vel;
	//vettore m_acc;
	std::vector<TH1I*> m_histos;
	std::vector<TGraph*> m_graps;
	vettore m_pos0;
	double m_Ek;
	double m_Ep;
	vettore m_L;
	float m_teta;
	float m_TT;
	std::vector<vettore> m_peri;
	vettore m_app;
	vettore m_sap;
	vettore m_s0; //serve per la v2 del calolo della precessione, ma va anche l'altro
	
    public:
	corpo();
	corpo(std::string n, double m, vettore r, vettore v, float Torb, float t=0);
	void leggi();
	void ass(std::string n, double m, vettore r, vettore v, float Torb, float t=0);
	virtual void evolvidT(std::vector<corpo*> cc, unsigned int dt, uint32_t mode, uint64_t j);
	void muovi(std::vector<corpo*> cc, unsigned int dt, uint32_t mode);
	double MASS(){return m_massa;}
	vettore V(){return m_vel;}
	vettore P(){return m_pos;}
	std::string NAME(){return m_nome;}
	TH1I* getisto(uint32_t i);
	uint32_t numHistos(){return m_histos.size();}
	void inizia(); //crea istogrammi
	vettore acc(std::vector<corpo*> &cc);
	double ECIN(){return m_Ek;}
	double EPOT(){return m_Ep;}
	double EMEC(){return m_Ek+m_Ep;}
	vettore LA(){return m_L;}
	void modE(std::vector<corpo*> cc);
	void istEmec(std::vector<corpo*> cc);
	
	vettore P0(){return m_pos0;}
	vettore AP(){return m_app;}
	vettore SAP(){return m_sap;}
	vettore S0(){return m_s0;}
	void modAPP(const vettore a){m_app=a;}
	void modSAP(const vettore a){m_sap=a;}
	void precessione(float TTerra);
	float period(){return m_TT;}
	float incl(){return m_teta;}
	TGraph* getgraf(uint32_t i);
	uint32_t numgraf(){return m_graps.size();}
	void fillgraf(std::string aa, std::string bb="%lg %lg");
};
std::ostream& operator<<(std::ostream &stream, corpo& c);
#endif
