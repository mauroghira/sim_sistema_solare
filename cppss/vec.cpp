#include<iostream>
#include<cstdlib>
#include<cmath>
#include "vec.h"

vettore::vettore(){
	m_x=0;
	m_y=0;
	m_z=0;
}
vettore::vettore(double r, double i, double z){
	m_x=r;
	m_y=i;
	m_z=z;
}
void vettore::leggi(){
	std::cout<<"componente x: ";
	std::cin>>m_x;
	std::cout<<"componente y: ";
	std::cin>>m_y;
	std::cout<<"componente z: ";
	std::cin>>m_z;
}

void vettore::sfToCar(){
	double theta=m_z*M_PI/180;
	double phi=m_y*M_PI/180;
	m_z=m_x*cos(theta);
	m_y=m_x*sin(theta)*sin(phi);
	m_x=m_x*sin(theta)*cos(phi);
}

vettore vettore::operator+(const vettore &c){
	vettore s(m_x+c.m_x, m_y+c.m_y, m_z+c.m_z);
	return s;
}
vettore vettore::operator-(const vettore &c){
	vettore s(m_x-c.m_x, m_y-c.m_y, m_z-c.m_z);
	return s;
}

vettore vettore::operator*(const float k){
	vettore s(k*m_x, k*m_y, k*m_z);
	return s;
}
vettore vettore::operator*(const vettore &c){
	float x, y, z;
	x=m_y*c.m_z-m_z*c.m_y;
	y=-m_x*c.m_z+m_z*c.m_x;
	z=m_x*c.m_y-m_y*c.m_x;
	vettore s(x, y, z);
	return s;
}

double  vettore::operator/(const vettore &c){
	return m_x*c.m_x+m_y*c.m_y+m_z*c.m_z;
}
vettore vettore::operator/(const float k){
	vettore s(m_x/k, m_y/k, m_z/k);
	return s;
}

double vettore::x() const{
	return m_x;
}
double vettore::y() const{
	return m_y;
}
double vettore::z() const{
	return m_z;
}

float vettore::phi(){
	float t;
	if(m_x>0)t=atan(m_y/m_x);
	else if(m_x<0) t=3.141592658979-atan(m_y/(-m_x));
	else if (m_y!=0) t=abs(m_y)*3.141592658979/(2*m_y);
	else t=0;
	return t;
}
float vettore::teta(){
	float t;
	if (modulo()!=0) t=acos(m_z/modulo());
	else t=0;
	return t;
}
void vettore::cilindriche(){
	double m=sqrt(m_x*m_x+m_y*m_y);
	float t=phi();
	std::cout<<"coordinate cilindriche di c: rho="<<m<<", phi="<<t<<", z="<<m_z<<std::endl;
}
void vettore::sferiche(){
	double m=modulo();
	float t=teta();
	float f=phi();
	std::cout<<"coordinate sferiche di c: R="<<m<<", phi="<<f<<", teta="<<t<<std::endl;
}

double vettore::angolo(const vettore &c){
	double d=operator/(c);
	double a=modulo()*c.modulo();
	return acos(d/a)*180/M_PI;
}

std::ostream& operator<<(std::ostream &stream, const vettore c){
	/*float a=c.x(), b=c.y(), d=c.z();
	stream <<"("<<a<<", "<<b<<", "<<d<<")";*/
	stream <<"("<<c.m_x<<", "<<c.m_y<<", "<<c.m_z<<")";
	return stream;
}
