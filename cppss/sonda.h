#include "vec.h"
#include "corpo.h"
#ifndef SONDA
#define SONDA
class sonda: public corpo{
	corpo* m_origine;
	int m_lancio;
	vettore m_Fmot;
	int m_spinta;
    public:
	virtual void evolvidT(std::vector<corpo*> cc, float dt);
	sonda();
	sonda(std::string n, float m, vettore r, vettore v, vettore F, int Tsp, int data, corpo* ori);
}
#endif
