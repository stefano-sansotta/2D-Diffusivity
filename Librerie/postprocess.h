//Classe con alcune funzione per post process
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#define DST(X,Y)  abs(X-Y)
#define DST2(X,Y) pow((X-Y),2)

using namespace std;


class postpro{

public:

void Diffusion2(float**,float**,float*,int);
void ModPos(float**,float**,float**,float,float,float,float,float,float,int);
float distanza(float,float,float,float,float,float);
float distanza(float* const,float* const);
float distanza2(float,float,float,float,float,float);
void Normalize(float* const,const int&);
void Normalize(float* const,const int&,const float&);
void Intercept(const float* const,const  float* const,float*,const int&);
void Intercept(const float* const,const  float* const,float*,const int&,const int&);
float Distance_Min(float** const,const float&,const float&,const float&,const int&);

};

class postpro2{

public:

void BordCheck(float**,float**,float**,float,float,float,float,int,int);

};

//Funzione per controllare bordi
void BordCheck(float** Atm, float** Stp, float** AtmO, float Bord, float limit, float max, float min, int i, int at)
{
float Ram;

if (DST2(Atm[i][at],Stp[i][at]) > limit)
{
if (DST(Stp[i][at],max) < DST(Stp[i][at],min))
{
Ram = Atm[i][at];
Atm[i][at] = Stp[i][at] + (Atm[i][at]-AtmO[i][at]);
if (DST2(Atm[i][at],Stp[i][at]) > limit)
{
Atm[i][at] = Ram + Bord;
}
}

else
{
Ram = Atm[i][at];
Atm[i][at] = Stp[i][at] - (Atm[i][at]-AtmO[i][at]);
}
if (DST2(Atm[i][at],Stp[i][at]) > limit)
{
Atm[i][at] = Ram - Bord;
}
}
}


//Funzione per calcolare diffusion coeff 
void Diffusion2(float** Zero, float** Point, float* Result, int poin)
{
int i = 0, j = 0;
int Xpoin = 0;
int Ypoin = 0;
int Zpoin = 0;
float diffx = 0,diffy = 0,diffz = 0;

while ( i < poin )
{
if (Zero[j][0] == Point[i][0])
{
diffx += DST(Point[i][1],Zero[j][1]);
Xpoin++;

diffy += DST(Point[i][2],Zero[j][2]);
Ypoin++;

diffz += DST(Point[i][3],Zero[j][3]);
Zpoin++;

j++;
i++;
}

else
{j++;}
}

if (Xpoin > 0)
{
diffx /= Xpoin;
Result[0] = pow(diffx,2);
}
if (Ypoin > 0)
{
diffy /= Ypoin;
Result[1] = pow(diffy,2);
}
if (Zpoin > 0)
{
diffz /= Zpoin;
Result[2] = pow(diffz,2);
}
}  //End function



//Funzione per correggere la posizione degli atomi in un file xyz
//Considerando le condizioni al contorno periodiche
void ModPos(float** Atm, float** Stp, float** AtmO, float Xmin, float Xmax, float Ymin, float Ymax, float Zmin, float Zmax, int atom) {

	float Bord[3],limit[3];
	Bord[0]  =  DST(Xmax,Xmin);
	Bord[1]  =  DST(Ymax,Ymin);
	Bord[2]  =  DST(Zmax,Zmin);
	limit[0] =  pow(Bord[0]/2,2);
	limit[1] =  pow(Bord[1]/2,2);
	limit[2] =  pow(Bord[2]/2,2);
	int i = 0;

	while (i < atom) {
		BordCheck(Atm,Stp,AtmO,Bord[0],limit[0],Xmax,Xmin,i,1);
		BordCheck(Atm,Stp,AtmO,Bord[1],limit[1],Ymax,Ymin,i,2);
		BordCheck(Atm,Stp,AtmO,Bord[2],limit[2],Zmax,Zmin,i,3);
		i++;
	}
}

//funzione per il calcolo distanza tra 2 atomi
float distanza(float A, float A1, float B, float B1, float C,float C1) { 
	
	
	float result = sqrt(pow(A-A1,2)+pow(B-B1,2)+pow(C-C1,2));
	return result;

}

float distanza(float* const A1, float* const A2) {
	
	float result;
	result = sqrt(pow((A1[0]-A2[0]),2) + pow((A1[1]-A2[1]),2) + pow((A1[2]-A2[2]),2));
	return result;
}

//funzione per il calcolo distanza tra 2 atomi senza eseguire la radice quadrata
//in modo tale da risparmiare sul costo computazionale
float distanza2(float A, float A1, float B, float B1, float C,float C1) { 
	
	const float result = pow(A-A1,2)+pow(B-B1,2)+pow(C-C1,2);
	return result;
}






/*********************************************
Function to normalize line to the last 20 
points to have 1
**********************************************/
void Normalize(float* const Y, const int& righe) {

	const int venti = 20;

	enum {uno=1};

	float media = 0;

	for (int i = righe - venti; i < righe; i++)
		media += Y[i];

	media /= venti;
	media = uno/media;

	for (int i = 0; i < righe; i++)
		Y[i] = Y[i]*media;

}




/*********************************************
Function to normalize line to the last 20 
points to have given number
**********************************************/
void Normalize(float* const Y, const int& righe, const float& uno)
{

const int venti = 20;

float media = 0;

for (int i = righe - venti; i < righe; i++)
media += Y[i];
media /= venti;
media = uno/media;

for (int i = 0; i < righe; i++)
Y[i] = Y[i]*media;

}



/*********************************************
Function to calculate the intercept of a line
**********************************************/
inline void Intercept(const float* const X, const  float* const Y, float* Par, const int& time) {

	float Xm=0,Ym=0;
	for (int i = 0; i < time; i++) {
		Xm += X[i];
		Ym += Y[i];
	}
	Xm /= time;                    //Punto X medio
	Ym /= time;                    //Punto Y medio

	float num=0,den=0;
	for (int i = 0; i < time; i++) {
		num += (X[i]-Xm)*(Y[i]-Ym);
		den += pow((X[i]-Xm),2);
	}
	Par[1] = num/den;
	Par[0] = Ym-Par[1]*Xm;
}



/*****************************************************************
Function to calculate the intercept of a line with a minimum value 
*****************************************************************/
inline void Intercept(const float* const X, const  float* const Y, float* Par, const int& time, const int& MIN) {

	float Xm=0,Ym=0;
	for (int i = MIN; i < time; i++) {
		Xm += X[i];
		Ym += Y[i];
	}

	Xm /= time;                    //Punto X medio
	Ym /= time;                    //Punto Y medio

	float num=0,den=0;
	for (int i = 0; i < time; i++) {
		num += (X[i]-Xm)*(Y[i]-Ym);
		den += pow((X[i]-Xm),2);
	}

	Par[1] = num/den;
	Par[0] = Ym-Par[1]*Xm;
}


/*****************************************************************
  Function to calculate the intercept of a line with a minimum 
  and a maximum value
*****************************************************************/
inline void Intercept(const float* const X, const  float* const Y, float* Par, const int& time, const int& MIN, const int& MAX) {

	int End = 0;
	
	if (MAX < time)
		End = MAX;
	else 
		End = time;

	float Xm=0,Ym=0;
	for (int i = MIN; i < End; i++) {
		Xm += X[i];
		Ym += Y[i];
	}

	Xm /= End-MIN;                    //Punto X medio
	Ym /= End-MIN;                    //Punto Y medio

	float num=0,den=0;
	for (int i = MIN; i < End; i++) {
		num += (X[i]-Xm)*(Y[i]-Ym);
		den += pow((X[i]-Xm),2);
	}

	Par[1] = num/den;
	Par[0] = Ym-Par[1]*Xm;
}




/*********************************************
  Function to calculate the minimum distance 
  between one atom and a surface
*********************************************/
inline float Distance_Min(float** const Surface, const float& X, const float& Y, const float& Z, const int& Si ) {

        float DistMin = 5000;

        for (int i = 0; i < Si; i++) {
                const float Dist = distanza2(Surface[i][0],X,Surface[i][1],Y,Surface[i][2],Z);
                if (Dist < DistMin) DistMin = Dist;
	}

	return sqrt(DistMin);

}








