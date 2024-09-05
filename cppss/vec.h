#include<iostream>
#include<cmath>
#ifndef VETTORE
#define VETTORE
class vettore{
	private:
		double m_x, m_y, m_z;
		
	public:
		vettore(double x, double y, double z);
		vettore();
		void leggi();
		//~complessso();serve per ogggetti dinamici
		vettore operator+(vettore &c);
		double operator/(vettore &c); //prodotto scalare
		vettore operator/(float k);
		vettore operator*(float k);  //prodotto vettoriale
		vettore operator*(vettore &c);
		vettore operator-(vettore &c);
		float teta();
		float phi();
		double x() const;
		double y() const;
		double z() const;
		void sferiche();
		void cilindriche();
		double modulo(){return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);}
		friend std::ostream& operator<<(std::ostream &stream, const vettore c);
		void sfToCar();
		float angolo(vettore &c);
};
std::ostream& operator<<(std::ostream &stream, const vettore c);
//float angolo(vettore &b, vettore &c);
#endif
