//Test funzione Determinante
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace std;


class Calc{

public:

float  DetMat(float**,int);
void   MatrixVecMol(float**,float*,float*,int);
void   TraspostaQ(float**,int);
void   MatrixScalMol(float**,int,float);


};


class Calc2{

public:

int    Fat(int);
void   Conj(float**,float**,int,int,int);
float  DetMatB(float**);
float  DetMatC(float**);


};


class Calc3{

public:

void   CoMatrix(float**,float**,int);

};



//Funzione per calcolare determinante matrice quadrata bidimensionale
float DetMatB(float** Mat)
{
float det = 0;

det = Mat[0][0]*Mat[1][1]-Mat[0][1]*Mat[1][0];

return det;
}


//Funzione per calcolare determinante matrice 3X3
float DetMatC(float** Mat)
{
float det = 0;

det = Mat[0][0]*Mat[1][1]*Mat[2][2]+Mat[0][1]*Mat[1][2]*Mat[2][0]+Mat[0][2]*Mat[1][0]*Mat[2][1]-(Mat[0][2]*Mat[1][1]*Mat[2][0]+Mat[0][1]*Mat[1][0]*Mat[2][2]+Mat[0][0]*Mat[1][2]*Mat[2][1]);

return det;
}



//Funzione per calcolare la funzione congiunta
void Conj(float** Mat ,float** MatEq ,int A ,int i, int k)
{

int l = 0, m = 0;
bool bol = true;
for (int j = 0; j < (A-1); j++)
{
for (int h = 0; h < (A-1); h++)
{


if (j == i && l < (A-1) && bol == true)
{
l++;
bol = false;
}
if (h == k && m < (A-1))
{m++;}


MatEq[j][h] = Mat[l][m];
if (m < A-1)
{m++;}
}
bol = true;
l++;
m = 0;
}

}



//Funzione per trovare numero fattoriale
int Fat(int A)
{

unsigned int result = 1;
if (A == 0)
{
}

else
{
for (int i = A; A > 0; A--)
{
result *= A;
}
}

return result;

}



//Funzione di supporto per calcolare il determinante di una matrice quadrata
float DetMat(float** Mat,int A)
{
float det = 0;

if ( A < 2 )
{det = 0;}

else if ( A == 2 )
{det = DetMatB(Mat);}

else if ( A == 3 )
{det = DetMatC(Mat);}

else
{

float** NewMat= (float**) malloc((A-1)*sizeof(float*));
for (int i = 0; i < (A-1); i++)
{NewMat[i] = (float*) malloc((A-1)*sizeof(float));}

int index = 0;

while (index < A)
{

//Conj(Mat,NewMat,A,index,0);
//pos[index] = (pow(-1,(index+2)))*Mat[index][0]*DetMat(Conj,(A-1));


index++;
}

}

return det;
}



//Funzione per calcolare il prodotto tra una matrice quadrata ed una colonna 
//ed una colonna di dimensione equivalente
void MatrixVecMol(float**Mat, float* Vec, float* Res ,int A)
{
for (int i = 0; i < A; i++)
{
for (int k = 0; k < A; k++)
{
Res[i] += Mat[i][k]*Vec[k];
}
}
}




//Funzione per fare la trasposta di una matrice quadrata
void   TraspostaQ(float** Mat, int A)
{
float RAM[A][A];


for (int i = 0; i < A; i++)
{
for (int k = 0; k < A; k++)
{
RAM[i][k] = Mat[k][i];
}
}
for (int i = 0; i < A; i++)
{
for (int k = 0; k < A; k++)
{
Mat[i][k] = RAM[i][k];
}
}
}


//Funzione per fare moltiplicazione tra scalare e matrice
void   MatrixScalMol(float** Mat, int A, float scalare)
{
for (int i = 0; i < A; i++)
{
for (int k = 0; k < A; k++)
{
Mat[i][k] *= scalare;
}
}
}



//Funzione per calcolare la matrice dei cofattori
void   CoMatrix(float** Mat,float** Cof, int A)
{

TraspostaQ(Mat,A);

if (A > 2)
{

float** MatEq = (float**) malloc((A-1)*sizeof(float*));
for (int i = 0; i < (A-1); i++)
{MatEq[i] = (float*) malloc((A-1)*sizeof(float));}

for (int i = 0; i < A; i++)
{
for (int k = 0; k < A; k++)
{

Conj(Mat,MatEq,A,i,k);

if (((i+1)+(k+1))%2==0)
{Cof[i][k] = DetMat(MatEq,A-1);}
else
{Cof[i][k] = DetMat(MatEq,A-1)*(-1);}

}
}
}

else
{

Cof[0][0]  =  Mat[1][1];
Cof[0][1]  =  (-1)*Mat[1][0];
Cof[1][0]  =  (-1)*Mat[0][1];
Cof[1][1]  =  Mat[0][0];

}

}












