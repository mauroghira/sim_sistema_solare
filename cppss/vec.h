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
		vettore operator+(const vettore &c);
		double operator/(const vettore &c); //prodotto scalare
		vettore operator/(const float k);
		vettore operator*(const float k);  //prodotto vettoriale
		vettore operator*(const vettore &c);
		vettore operator-(const vettore &c);
		float teta();
		float phi();
		double x() const;
		double y() const;
		double z() const;
		void sferiche();
		void cilindriche();
		double modulo() const {return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);}
		friend std::ostream& operator<<(std::ostream &stream, const vettore c);
		void sfToCar();
		float angolo(const vettore &c);
};
std::ostream& operator<<(std::ostream &stream, const vettore c);
//float angolo(const vettore &b, vettore &c);
#endif
