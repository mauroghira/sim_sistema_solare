#include "vec.h"
#include<cstring>
#include<vector>
#include<cmath>
#include "corpo.h"
#include<iostream>

int main(){
	corpo x;
	x.leggi();
	std::cout<<x.MASS()<<x.NAME()<<x.V()<<x.P()<<std::endl;
	vettore v(1,0,1), p(0,0,0);
	std::cout<<p<<std::endl;
	corpo y("lol",100000000000,p,v);
	std::cout<<y.MASS()<<y.NAME()<<y.V()<<y.P()<<std::endl;	
	std::vector<corpo*> lista;
	lista.push_back(&x);
	lista.push_back(&y);
	x.modE(lista);
	x.istEmec(lista);
	x.evolvidT(lista, 60, 1);
	std::cout<<x.MASS()<<x.NAME()<<x.V()<<x.P()<<std::endl;
	return 0;
}
