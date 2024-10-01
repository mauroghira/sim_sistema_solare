#include "vec.h"
#include "corpo.h"
#include<vector>
#include<cstring>
#include <TCanvas.h>
#include <TView.h>
#include <TMarker.h>
#include <time.h>
#ifndef SISTEMA
#define SISTEMA
class sistema{
	std::vector<corpo*> m_corpi;
	unsigned int m_dT;
	float m_T;
	std::vector<TH1I*> m_ist;	
	std::string m_inc; //pezzo di stringa per file dei dati per i grafici
	
    public:
	sistema(float T, unsigned int dt, std::string config, std::string of);
	void input(std::string config, std::string of);
	sistema();
	void add(corpo* c);
	void leggi(std::string config);  //leggere dati da file
	void evodt(uint32_t mode, uint64_t j);
	void evo(uint32_t mode, int st=0);
	void print();
	TH1I* getThisHisto(std::string nomeCorpo, uint32_t indice);
	void PrintHistos();
	void modT(float s){m_T=s;}
	TH1I* getist(uint32_t i);
	uint32_t numist(){return m_ist.size();}
	void ins();  //istogrammi
	
	TGraph* getThisGraph(std::string nomeCorpo, uint32_t indice);
	int select(std::string p, int v);
	void output(std::string file);
	void mkGraf(std::string pla="");
	void modSTR(std::string out){m_inc=out;}
	void clean();
	void savehist(std::string out);
	std::vector<double> anal();
	float tempo(){return m_T;}
	std::vector<float> times(std::vector<float>* teta);
};
#endif
