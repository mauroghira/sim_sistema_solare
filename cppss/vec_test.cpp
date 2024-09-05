#include<iostream>
#include<cstdlib>
#include<cmath>
#include "vec.h"

int main(){
	float k;
	vettore c;
	c.leggi();
	std::cout<<"c="<<c<<std::endl;
/*	c.sferiche();
	c.cilindriche();
	c.sfToCar();
	std::cout<<c<<std::endl;
	vettore t=c;
	std::cout<<t<<std::endl;
	std::cout<<"fattore di riscalamento: ";
	std::cin>>k;
	vettore s=c/k;
	std::cout<<"kc= "<<s<<std::endl;*/
	vettore b(1,0,1);
/*	std::cout<<"b="<<b<<std::endl;
	vettore a=b+c;
	std::cout<<"b+c="<<a<<std::endl;
	vettore e=b-c;
	std::cout<<"b-c="<<e<<std::endl;
	vettore p=b*c;
	std::cout<<"b vettor c="<<p<<std::endl;
	float q=b/c;
	std::cout<<"b scalar c="<<q<<std::endl; */
	std::cout<<"angolo: "<<c.angolo(b)<<std::endl;
	return 0;
}
