#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#define ATM 10
#define BND 4
#define ANG 5
#define DHL 6

using namespace std;


class Lammps{

public:
  
int ATOMI(string);
int BOND(string);
int ANGLE(string);
int DIHEDRAL(string);

int ATOMT(string);
int BONDT(string);
int ANGLET(string);
int DIHEDRALT(string);

int StepNumbXYZ(const string&);
int AtomNumbXYZ(const string&);
int AtomTypeNumbXYZ(const string&,const string&);
int AtomNumb(const string&,const string&,const int&);
//int AtomTypeCounterXYZ(const string&);

float DistanzaTot(float**,int,int,int);
float DistanzaTot2(float**,int,float,float,float);

template <typename MN>
float MinDistAtmSup(MN** const,MN* const,const float&,const int&);

template <typename MN>
float MinDistAtmSup(MN** const, MN* const, const int&);

template <typename MN>
void CoordMat(const string&, MN** const);

template <typename MN>
void BondsMat(const string&, MN** const);

template <typename MN>
void AngleMat(const string&, MN** const);

template <typename MN>
void DihedralMat(const string&, MN** const);

void CenterMol(float** const,float* const,const int);

void CenterMolLammps(float** const,float* const,const int);

void FILECOORD(const string&,const string&);
void FILEBOND(string,string);
void FILEANGLE(string,string);
void BOXB(string,float*);
void FILEDATA(int,int,int,int,int,int,float**,float**,float**,float*,float*,string);
void GetMass(string,float*,int);
void StartBatch24(string);
void StartBatch2(string);

template <typename B>
void BorderCalculate(const string&,B*);

template <typename TM>
void TraslateM(TM** const, TM* const, TM* const,const double&,const int&);
//void Traslate(const string&,const float* const ,const float* const,const float&);

template <typename DC>
int DistCheckM(DC** const,DC* const,const float&,const int&);


};




//Funzione per conoscere quanti atomi ci sono nel datafile di Lammps
int ATOMI(string nome)
{
const string str = "atoms";
ifstream file;
file.open(nome.c_str());
const int posizione = FIND(nome,str);
file.seekg(posizione-1);
int A;
file >> A;
return A;
}


//Funzione per trovare il numero dei legami in un file data di Lammps
int BOND(string nome)
{
const string str = "bonds";
int posizione = FIND(nome,str);
ifstream file;
file.open(nome.c_str());
file.seekg(posizione-1);
int A;
file >> A;
file.close();
return A;
}



//Funzione per trovare il numero di angoli in un file data di Lammps
int ANGLE(string nome) {
	
	string str = "angles";
	int posizione = FIND(nome,str);
	ifstream file;
	file.open(nome.c_str());
	file.seekg(posizione-1);
	int A;
	file >> A;
	file.close();
	return A;
}


//Funzione per calcolare il numero di dihedral nel data file di lammps
int DIHEDRAL(string nome)  {

	string str = "dihedrals";
	int pos = FIND(nome,str);
	ifstream file;
	file.open(nome.c_str());
	file.seekg(pos-1);
	int H;
	file >> H;
	file.close();
	return H;

}
  



//Funzione per trasferire Atomi da datafile ad un file chiamato Coordinate
void FILECOORD(const string& data,const  string& nome)
{
const int Atomi = ATOMI(data);
string str = "Atoms";
int posizione = FIND(data,str);
ifstream file;
file.open(data.c_str());
file.seekg(posizione);
string line;
for (int i = 0; i < 2; i++)
{getline(file,line);}
ofstream fine;
fine.open(nome.c_str());
float A,B,C,D,E,F,G,H,I,J;
for (int i = 0; i < Atomi; i++)
{
file >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
fine << A << " " << B << " " << C << " " << D << " " << E << " " << F << " " << G << " " << H << " " << I << " " << J <<endl;
}
file.close();
fine.close();
}




//Funzione per trasferire Bonds da datafile ad un file chiamato Bonds
void FILEBOND(string data, string nome)
{
int Bond = BOND(data);
string str = "Bonds";
int posizione = FIND(data,str);
ifstream file;
file.open(data.c_str());
file.seekg(posizione);
string line;
for (int i = 0; i < 2; i++)
{getline(file,line);}
ofstream fine;
fine.open(nome.c_str());
float A,B,C,D;
for (int i = 0; i < Bond; i++)
{
file >> A >> B >> C >> D;
fine << A << " " << B << " " << C << " " << D <<endl;
}
file.close();
fine.close();
}




//Funzione per trasferire angoli da datafile ad un file chiamato Angle
void FILEANGLE(string data, string nome)
{
int Angle = ANGLE(data);
string str = "Angles";
int posizione = FIND(data,str);
ifstream file;
file.open(data.c_str());
file.seekg(posizione);
string line;
for (int i = 0; i < 2; i++)
{getline(file,line);}
ofstream fine;
fine.open(nome.c_str());
float A,B,C,D,E;
for (int i = 0; i < Angle; i++)
{
file >> A >> B >> C >> D >> E;
fine << A << " " << B << " " << C << " " << D << " " << E <<endl;
}
file.close();
fine.close();
}





//Funzione per estrarre box dimensions dal file e salvarle in un vettore
void BOXB(string data, float* V)
{
float A,B;
string str = "xlo";
int posizione = FIND(data,str);
ifstream file;
string line;
file.open(data.c_str());
file.seekg(posizione-1);
file >> A >> B;
V[0] = A;
V[1] = B;
getline(file,line);
file >> A >> B;
V[2] = A;
V[3] = B;
getline(file,line);
file >> A >> B;
V[4] = A;
V[5] = B;
}




//Funzione per creare file data e scriverci sopra coordinate bonds e angoli
void FILEDATA(int Atoms, int Atype, int Bond, int Btype, int Angle, int Antype, int Dihedral, int Dtype, float** At, float** Bn, float** An, float** Dh, float* V, float* M, string nome) { 



	FILE* file;
	file = fopen(nome.c_str(), "w" );
	fprintf(file,"Lammps datafile created with library created by Stefano Sansotta PhD student at FAU Erlangen Nuremberg\n\n");
	fprintf(file," %i atoms\n",Atoms);
	fprintf(file," %i atom types\n",Atype);
	fprintf(file," %i bonds\n",Bond);
	fprintf(file," %i bond types\n",Btype);
	fprintf(file," %i angles\n",Angle);
	fprintf(file," %i angle types\n",Antype);
        fprintf(file," %i dihedrals\n",Dihedral);
        fprintf(file," %i dihedral types\n",Dtype);
	fprintf(file,"\n %.16e %.16e xlo xhi\n",V[0],V[1]);
	fprintf(file," %.16e %.16e ylo yhi\n",V[2],V[3]);
	fprintf(file," %.16e %.16e zlo zhi\n",V[4],V[5]);
	fprintf(file,"\n Masses\n\n");

	for (int i = 0; i < Atype; i++) {
		fprintf(file,"%i %.4f\n",i+1,M[i]);
	}
	
	fprintf(file,"\n Atoms\n\n");
	
	for (int i = 0; i < Atoms; i++) {
		fprintf(file,"%.0f %.0f %.0f %.16e %.16e %.16e %.16e %.0f %.0f %.0f\n",At[i][0],At[i][1],At[i][2],At[i][3],At[i][4],At[i][5],At[i][6],At[i][7],At[i][8],At[i][9]);
	}

	if (Bond > 0) {
		fprintf(file,"\n Bonds\n\n");
		for (int i = 0; i < Bond; i++)	{
			for (int k = 0; k < BND; k++) {
				(k == (BND-1) ? fprintf(file,"%.0f\n",Bn[i][k]) : fprintf(file,"%.0f ",Bn[i][k]));
			}
		}
	}

	if (Angle > 0) {
		fprintf(file,"\n Angles\n\n");
		for (int i = 0; i < Angle; i++) {
			for (int k = 0; k < ANG; k++) {
				(k == (ANG-1) ? fprintf(file,"%.0f\n",An[i][k]) : fprintf(file,"%.0f ",An[i][k]));
			}
		}
	}

	if (Dihedral > 0) { 
		fprintf(file,"\n Dihedrals\n\n");
		for (int i = 0; i < Dihedral; i++) { 
			for (int k = 0; k < DHL; k++) {
				(k == (DHL-1) ? fprintf(file,"%.0f\n",Dh[i][k]) : fprintf(file,"%.0f ",Dh[i][k]));
			}
		}
	}

	fclose(file);
}


//Funzione per estrarre Masses da un datafile Lammps 
void GetMass(string data, float* V, int Atype)
{
float A,B;
string str = "Masses";
int posizione = FIND(data,str);
ifstream file;
string line;
file.open(data.c_str());
file.seekg(posizione);
for (int i = 0; i < 2; i++)
{getline(file,line);}
for (int i = 0; i < Atype; i++)
{
file >> A >> B;
V[i] = B;
}
}


//Funzione per calcolare la distanza tra 2 tipi di atomi in una simulazione
float DistanzaTot(float** P, int r, int A, int B) 
{

//Qui inizia il loop per il calcolo della distanza minima tra ogni singolo atomo di Si e di O
float D=5000;
int i = 0;

while (i < r)
{

// Inizia il loop quando vede un ossigeno 
if (P[i][0] == A)
{
float X,Y,Z;
X = P[i][1];
Y = P[i][2];
Z = P[i][3];

//Qui inizio il calcolo della distanza 
for (int j = 0; j < r; j++)
{

//controllo che la distanza sia calcolata con atomi di Silicio
if (P[j][0] == B)
{
float ram;
float X1,Y1,Z1;

X1 = P[j][1];
Y1 = P[j][2];
Z1 = P[j][3];

//Qui richiamo la funzione per il calcolo della distanza e quella del minimo.
ram = distanza2(X,X1,Y,Y1,Z,Z1);
D = fmin(D,ram);
}
}
}
i++;
}
return D;
}




//Funzione per calcolare la distanza tra un atomo ed un altro gruppo di atomi 
float DistanzaTot2(float** P, int r, float X, float Y, float Z) 
{

//Qui inizia il loop per il calcolo della distanza minima tra ogni singolo atomo di Si e di O
float D=5000;
int i = 0;

while (i < r)
{

//Qui inizio il calcolo della distanza 
for (int j = 0; j < r; j++)
{

float ram;
float X1,Y1,Z1;

X1 = P[j][1];
Y1 = P[j][2];
Z1 = P[j][3];

if (X1 == X && Y1 == Y && Z1 == Z)
{ram = D;}

else 
{ram = distanza(X,X1,Y,Y1,Z,Z1);}

D = fmin(D,ram);

}

i++;

}

return D;
}


int ATOMT(string nome)
{
string str = "atom types";
int posizione = FIND(nome,str);
ifstream file;
file.open(nome.c_str());
file.seekg(posizione-1);
int A;
file >> A;
file.close();
return A;
}




int BONDT(string nome)
{
string str = "bond types";
int posizione = FIND(nome,str);
ifstream file;
file.open(nome.c_str());
file.seekg(posizione-1);
int A;
file >> A;
file.close();
return A;
}




//Funzione per calcolare quanti tipi di angoli ci sono nell'input file
int ANGLET(string nome)
{
string str = "angle types";
int posizione = FIND(nome,str);
ifstream file;
file.open(nome.c_str());
file.seekg(posizione-1);
int A;
file >> A;
file.close();
return A;
}


//Funzione per contare il numero di dihedrals type nel data file di lammps
int DIHEDRALT(string nome)  {
	string str = "dihedral types";
	int posizione = FIND(nome,str);
	ifstream file;
	file.open(nome.c_str());
	file.seekg(posizione-1);
	int A;
	file >> A;
	file.close();
	return A;
}




//Funzione per iniziare un job in zahn server
void StartBatch24(string start)
{

string batch  =  "salloc mpirun /raid/bin/lmp_ubuntu_update -N1 --exclusive -p 24h -J ";
string in     =  " -in ";
string log    =  " -log log.";
string rest   = " -screen none &> out.";
string command;
command.append(batch);
command.append(start);
command.append(in);
command.append(start);
command.append(log);
command.append(start);
command.append(rest);
command.append(start);

system(command.c_str());

}




//Funzione per far partire simulazioni nel 2h que
void StartBatch2(string start)
{

string batch  =  "salloc mpirun /raid/bin/lmp_ubuntu_update -N1 --exclusive -p 2h -J ";
string in     =  " -in ";
string log    =  " -log log.";
string rest   = " -screen none &> out.";
string command;
command.append(batch);
command.append(start);
command.append(in);
command.append(start);
command.append(log);
command.append(start);
command.append(rest);
command.append(start);

system(command.c_str());

}





//Funzione per calcolare i bordi da file xyz
template <typename B>
inline void BorderCalculate(const string& nome, B* Brd)
{

ifstream file;
file.open(nome.c_str());

int atoms;
file >> atoms;

string line;
enum {dos=2};
for (int i = 0; i < dos; i++)  getline(file,line);

enum {vier=4};
float** Crd = (float**) malloc(atoms*sizeof(float*));
for (int i = 0; i < atoms; i++)
Crd[i] = (float*) malloc(vier*sizeof(float));
MatrixIn(Crd,atoms,vier);

string A;
for (int i = 0; i < atoms; i++)
{

file >> A;
file >> Crd[i][1];
file >> Crd[i][2];
file >> Crd[i][3];
}
file.close();

float* const X = (float*) malloc(atoms*sizeof(float));
float* const Y = (float*) malloc(atoms*sizeof(float));
float* const Z = (float*) malloc(atoms*sizeof(float));
Insert(Crd,atoms,X,1);
Insert(Crd,atoms,Y,2);
Insert(Crd,atoms,Z,3);

for (int i = 0; i < atoms; i++)
free(Crd[i]);
free(Crd);

Brd[0] = *min_element(X, X+atoms);
Brd[1] = *max_element(X, X+atoms);
Brd[2] = *min_element(Y, Y+atoms);
Brd[3] = *max_element(Y, Y+atoms);
Brd[4] = *min_element(Z, Z+atoms);
Brd[5] = *max_element(Z, Z+atoms);

free(X);
free(Y);
free(Z);

}



/****************************************************
Function to count the number of steps in a xyz file
****************************************************/
inline int StepNumbXYZ(const string& nome)
{

int Steps = 0;
ifstream file;
string line;
const string str = "Timestep:";
file.open(nome.c_str());
int atoms;
file >> atoms;

while (getline(file,line))
{
if (line.find(str, 0) != string::npos)
{
Steps++;
//const int POS = file.tellg();
//file.seekg((POS+atoms));
}
}
file.close();

return Steps;

}



/****************************************************
Function to count the number of atoms in a xyz file
****************************************************/
int AtomNumbXYZ(const string& nome)
{

int atoms;
ifstream file;
file.open(nome.c_str());
file >> atoms;
file.close();

return atoms;

}





/******************************************************
Function to calculate the atoms of one type in xyz file
******************************************************/
int AtomTypeNumbXYZ(const string& nome, const string& atomo)
{

int counter = 0;
int atoms;
string line,A;
float B,C,D;

ifstream file;
file.open(nome.c_str());
file >> atoms; 

getline(file,line);
getline(file,line);

for (int i = 0; i < atoms; i++)
{
file >> A >> B >> C >> D;
if (A == atomo)
counter++;
}
file.close();

return counter;

}

/********************************************************************
 Funzione per calcolare il numero di atomi di un tipo presenti in un 
          file xyz. Si considera solo il primo Step
********************************************************************/
int AtomNumb(const string& atom, const string& nome, const int& atoms) {

	int Atomo = 0;
	string A;
	float B,C,D;
	ifstream file;
	file.open(nome.c_str());
	string line;

	for (int i = 0; i < 2; i++)
		getline(file,line);

	for (int i = 0; i < atoms; i++) {
		file >> A >> B >> C >> D;
		if ( A == atom ) Atomo++;
	}
	file.close();

	return Atomo;
}




/********************************************************
Function to change positions atoms directly in the matrix
********************************************************/
template <typename TM>
void TraslateM(TM** const NH3, TM* const BoxB, TM* const Arg, const double& Toll, const int& righe)
{

//Bordo X inferiore
if ( abs(BoxB[0] - Arg[0]) < Toll )
{
for (int i = 0; i < righe; i++)
{
if (NH3[i][1] > (BoxB[1] - Toll)) NH3[i][1] -= BoxB[1];
}
}

//Bordo X Superiore
if ( abs(BoxB[1] - Arg[0]) < Toll )
{
for ( int i = 0; i < righe; i++)
{
if (NH3[i][1] < (BoxB[0] + Toll)) NH3[i][1] += BoxB[1];
}
}

//Bordo Y Inferiore
if ( abs(BoxB[2] - Arg[1]) < Toll )
{
for (int i = 0; i < righe; i++)
{
if (NH3[i][2] > (BoxB[3] - Toll)) NH3[i][2] -= BoxB[3];
}
}

//Bordo Y Superiore
if ( abs(BoxB[3] - Arg[1]) < Toll )
{
for (int i = 0; i < righe; i++)
{
if ( NH3[i][2] < (BoxB[2] + Toll) ) NH3[i][2] += BoxB[3];
}
}

//Limite Z inferiore
if ( abs(BoxB[4] - Arg[2]) < Toll)
{
for (int i = 0; i < righe; i++)
{
if (NH3[i][3] > (BoxB[5] - Toll)) NH3[i][3] -= BoxB[5];
}
}

//Limite Z Superiore
if ( abs(BoxB[5] - Arg[2] ) < Toll )
{
for (int i = 0; i < righe; i++)
{
if ( NH3[i][3] < (BoxB[4] + Toll) ) NH3[i][3] += BoxB[5];
}
}

}





/******************************************************
Function that gives back the minimum distance between 
the atoms in one vector and the other atoms in a matrix
the vector has only 3 slot and the vector has 4 slot
for each culumn the first one is a number
******************************************************/
template <typename MN>
float MinDistAtmSup(MN** const Sup, MN* const Tmp, const float& distmax, const int& righe)
{

float Dist = 5000;
float ram;

for (int i = 0; i < righe; i++)
{

if (Sup[i][3] > distmax)
{
ram = distanza(Sup[i][1],Tmp[0],Sup[i][2],Tmp[1],Sup[i][3],Tmp[2]);
Dist = fmin(Dist,ram);
}

}

return Dist;

}



/******************************************************
Function that gives back the minimum distance between 
the atoms in one vector and the other atoms in a matrix
the vector has only 3 slot and the vector has 4 slot
for each culumn the first one is a number
******************************************************/
template <typename MN>
float MinDistAtmSup(MN** const Sup, MN* const Tmp, const int& righe)
{

float Dist = 5000;
float ram;

for (int i = 0; i < righe; i++)
{

ram = distanza(Sup[i][1],Tmp[0],Sup[i][2],Tmp[1],Sup[i][3],Tmp[2]);
Dist = fmin(Dist,ram);

}

return Dist;

}


/******************************************************
Function to count the number of different atoms present 
in a XYZ file 
******************************************************/
int AtomTypesCounter(const string& nome)
{

string Atomi[20];
int atomi=0;

ifstream file;
file.open(nome.c_str());
int atoms;
file >> atoms;

string A;
float B,C,D;
bool bol = false;

file >> A;

for (int i = 0; i < atoms; i++)
{

file >> A >> B >> C >> D;

if (i == 0)
{
Atomi[atomi].append(A);
if (atomi < 1000)
atomi++;
}
if (i > 0)
{
for (int k = 0; k < atomi; k++)
{
if (Atomi[k] == A) bol = true;
}

if (bol == true) bol = false;
else
{
Atomi[atomi] = A;
bol = false;
if (atomi < 1000) atomi++;
}
}
}
file.close();

return atomi;

}




/**************************************************
Function to check how many molecule are around an
atom at a specific distance
**************************************************/
template <typename DC>
int DistCheckM(DC** const Atm, DC* const Arg, const float& dist, const int& atoms)
{

int counter = 0;

for (int i = 0; i < atoms; i++)
if ( distanza(Atm[i][1],Arg[0],Atm[i][2],Arg[1],Atm[i][3],Arg[2]) <= dist ) counter++;

return counter;

}



/***************************************************/
/*    Function that put all the data about the     */
/*coordinates in a lammps data file inside a matrix*/
/***************************************************/
template <typename MN>
void CoordMat(const string& nome, MN** const Mat)
{

const string str = "Atoms";

const int Atomi = ATOMI(nome);
const int posizione = FIND(nome,str);

string line;

ifstream file;
file.open(nome.c_str());
file.seekg(posizione);

for (int i = 0; i < 2; i++) getline(file,line);

for (int i = 0; i < Atomi; i++)
file >> Mat[i][0] >> Mat[i][1]  >> Mat[i][2]  >> Mat[i][3] >> Mat[i][4]  >> Mat[i][5]  >> Mat[i][6]  >> Mat[i][7] >> Mat[i][8] >> Mat[i][9];
file.close();

}




/***************************************************/
/*    Function that put all the data about the     */
/*   bonds in a lammps data file inside a matrix   */
/***************************************************/
template <typename MN>
void BondsMat(const string& nome, MN** const Mat)
{

const string str = "Bonds";

const int Bond = BOND(nome);
const int pos  = FIND(nome,str);

string line;

ifstream file;
file.open(nome.c_str());
file.seekg(pos);

for (int i = 0; i < 2; i++) getline(file,line);

for (int i = 0; i < Bond; i++) 
file >> Mat[i][0] >> Mat[i][1] >> Mat[i][2] >> Mat[i][3];
file.close();

}



/***************************************************/
/*    Function that put all the data about the     */
/*  angles in a lammps data file inside a matrix   */
/***************************************************/
template <typename MN>
void AngleMat(const string& nome, MN** const Mat)
{

const string str = "Angles";

const int Angle = ANGLE(nome);
const int   Pos = FIND(nome,str);

string line;

ifstream file;
file.open(nome.c_str());
file.seekg(Pos);

for (int i = 0; i < 2; i++) getline(file,line);

for (int i = 0; i < Angle; i++)
file >> Mat[i][0] >> Mat[i][1] >> Mat[i][2] >> Mat[i][3] >> Mat[i][4];
file.close();

}


/***************************************************/
/*    Function that put all the data about the     */
/* dihedral in a lammps data file inside a matrix  */
/***************************************************/
template <typename MN>
void DihedralMat(const string& nome, MN** const Mat)  {

	const string str = "Dihedrals";

	const int Dihedrals = DIHEDRAL(nome);
	const int   Pos = FIND(nome,str);

	string line;

	ifstream file;
	file.open(nome.c_str());
	file.seekg(Pos);

	for (int i = 0; i < 2; i++) getline(file,line);

	for (int i = 0; i < Dihedrals; i++)
		file >> Mat[i][0] >> Mat[i][1] >> Mat[i][2] >> Mat[i][3] >> Mat[i][4] >> Mat[i][5];
	file.close();

}


/***********************************/
/* Function to calculate the center*/
/*        of a molecule            */
/***********************************/
void CenterMol(float** const Mol, float* const Vec, const int NAt)
{

float Xc = 0, Yc = 0, Zc = 0;

for (int i = 0; i < NAt; i++)
{

Xc += Mol[i][0];
Yc += Mol[i][1];
Zc += Mol[i][2];

}

Vec[0] = Xc/NAt;
Vec[1] = Yc/NAt;
Vec[2] = Zc/NAt;

}


/***********************************/
/* Function to calculate the center*/
/* of a molecule from lammps data  */
/***********************************/
void CenterMolLammps(float** const Mol, float* const Vec, const int NAt)
{

float Xc = 0, Yc = 0, Zc = 0;

for (int i = 0; i < NAt; i++)
{

Xc += Mol[i][4];
Yc += Mol[i][5];
Zc += Mol[i][6];

}

Vec[0] = Xc/NAt;
Vec[1] = Yc/NAt;
Vec[2] = Zc/NAt;

}



