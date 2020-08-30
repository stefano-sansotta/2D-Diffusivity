#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#define ATM 10
#define BND 4
#define ANG 5

using namespace std;

class KawZahn {

public:
  
void RemoveAtom(string,string,int,int,int);

};





//Funzione per rimuovere gli atomi
void RemoveAtom(string nome, string nome2, int Atype, int Btype, int Antype)
{


const int atoms = ATOMI(nome);
const int bonds = BOND(nome);
const int angle = ANGLE(nome);

string coord  = "Coordinate";
FILECOORD(nome,coord);

float** Atoms = (float**) malloc(atoms*sizeof(float*));
for (int i = 0; i < atoms; i++)
{Atoms[i] = (float*) malloc(ATM*sizeof(float));}
Fill(Atoms,atoms,ATM,coord);
remove(coord.c_str());

float** Bonds;
float** Angle;

if (bonds > 0)
{
string legami = "Bonds";
FILEBOND(nome,legami);

Bonds = (float**) malloc(bonds*sizeof(float*));
for (int i = 0; i < bonds; i++)
{Bonds[i] = (float*) malloc(BND*sizeof(float));}
Fill(Bonds,bonds,BND,legami);
remove(legami.c_str());
}

if (angle > 0)
{
string angoli = "Angles";
FILEANGLE(nome,angoli);

Angle = (float**) malloc(angle*sizeof(float*));
for (int i = 0; i < angle; i++)
{Angle[i] = (float*) malloc(ANG*sizeof(float));}
Fill(Angle,angle,ANG,angoli);
remove(angoli.c_str());
}

float** Atoms2 = (float**) malloc((atoms-1)*sizeof(float*));
for (int i = 0; i < (atoms-1); i++)
{Atoms2[i] = (float*) malloc(ATM*sizeof(float));}

int j = 0;
for (int i = 0; i < (atoms-1); i++)
{
if (Atoms[j][1] != 10000)
{
for (int k = 0; k < ATM; k++)
{

Atoms2[i][k] = Atoms[j][k];

}
}
j++;
}

float* Box  = (float*) malloc(6*sizeof(float));
float* Mass = (float*) malloc(Atype*sizeof(float));

BOXB(nome,Box);
GetMass(nome,Mass,Atype);

FILEDATA((atoms-1),Atype,bonds,Btype,angle,Antype,Atoms2,Bonds,Angle,Box,Mass,nome2);

}



















