//Lista di funzioni disponibili per calcolare l'intercetta di alcuni punti
#include <iostream>
#include <cmath>
#include <cstdio>


using namespace std;

class function{


public:


float Function1(float,float,float,float);
float DerFunc1(string,float);
float Function2(float,float,float,float);
float DerFunc2(string,float);
float Function3(float,float,float,float);
float DerFunc3(string,float);

};



//Funzione di prova da interpolare
float Function1(float A, float B , float C, float r)
{
return A+B/pow(r,2)+C/r;
}



//Derivate della prima funzione
float DerFunc1(string tipo, float r)
{
if (tipo == "A" || tipo == "a")
{
return 1;
}
else if (tipo == "B" || tipo == "b")
{
return (1/(pow(r,2)));
}
else if (tipo == "C" || tipo == "c")
{
return 1/(r);
}
else
{
cout << "Impossibile fare la derivata rispetto a " << tipo <<endl;
exit(EXIT_FAILURE);
}
}


//Funzione
float Function2(float A, float B, float r)
{
return A+B/pow(r,2);
}


//Derivata
float DerFunc2(string tipo, float r)
{
if (tipo == "A" || tipo == "a")
{
return 1;
}
else if (tipo == "B" || tipo == "b")
{
return (1/(pow(r,2)));
}
else 
{
cout << "Impossibile fare la derivata rispetto a " << tipo <<endl;
exit(EXIT_FAILURE);
}
}




//Funzione3
float Function3(float A, float B, float C, float r)
{
return A+B*exp(-1/r)+C/pow(r,2);
}


//Derivata3 
float DerFunc3(string tipo, float r)
{
if (tipo == "A" || tipo == "a")
{
return 1;
}
else if (tipo == "B" || tipo == "b")
{
return exp(-1/r);;
}
else if (tipo == "C" || tipo == "c")
{
return 1/pow(r,2);
}
else
{
cout << "Impossibile fare la derivata rispetto a " << tipo <<endl;
exit(EXIT_FAILURE);
}
}


























