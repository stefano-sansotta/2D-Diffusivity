#ifndef EDIT_H
#define EDIT_H


#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <stdio.h>


using namespace std;


class Edit{

public:


template <typename MM>
inline void InsertLine(MM** const,MM* const,const int&,const int&);

template <typename R1>
void Insert(R1** const,const int&,R1* const,const int&);

void AppendM(const string&,const float** const,const int&,const int&);

void VERIF(const string&);

template <typename M1>
void Fill(M1** const,const int&,const int&,const string&);

template <typename V>
void Fill(V* const,const int&,const string&);

template <typename C>
void CopyM(C** const,C** const,const int&,const int&);

template <typename T>
void PrintM(T**,const int&,const int&);

template <typename M>
void MatrixIn(M** const,const int&,const int&);

template <typename HH>
void VectorIn(HH* const,const int&);

void VectorInG(float*,int,float);

template <typename D>
void File(D**,int,int,string);

template <typename D1>
void FileE(D1**,const int&,const int&,const string&);


void MatrixOrderMin(float** const&,const int&,const int&,const int&);
void MatrixOrderMax(float** const&,const int&,const int&,const int&);

void mkdirf(float);
void cpf(string,float);
int  FIND(const string&,const string&);
int  Righe(const string&);
int  Foldercounter(const string&);
int  Foldercounter();


};

#endif


//Funzione di ferifica se il esiste
void VERIF(const string& nome) {

	FILE* file;
	file = fopen(nome.c_str(), "r" );
	if (file == NULL) {
		cout << "Errore: impossibile aprire il file " << nome <<endl;
		exit(EXIT_FAILURE);
	}
	fclose(file);

} //End function 



//Funzione per riempire matrice
template <typename M1>
void Fill(M1** const Mat, const int& R, const int& C, const string& name)
{
ifstream file;
file.open(name.c_str());
for (int i = 0; i < R; i++)
{
for (int k = 0; k < C; k++)
{
file >> Mat[i][k];
}
}
file.close();
}


//Funzione per riempire vettore
template <typename V>
void Fill(V* const Mat, const int& R, const string& name) {

	ifstream file;
	file.open(name.c_str());
	for (int i = 0; i < R; i++) {
		file >> Mat[i];
	}
	file.close();

}



//Funzione per copiare 2 matrici
template <typename C>
void CopyM(C** const Mat, C** const Mat2, const int& A, const int& B)
{
for (int i = 0; i < A; i++)
{
for (int k = 0; k < B; k++)
{
Mat2[i][k] = Mat[i][k];
}
}
}


//Funzione per stampare la matrice
template <typename T>
void PrintM(T** M, const int& A, const int& B)
{

for (int i = 0; i < A; i++)
{
for (int k = 0; k < B; k++)
{
cout << M[i][k] << " ";
}
cout <<endl;
}
}



//Funzione per appendere data in una matrice in un altro file
void AppendM(const string& data, const float** const Mat, const int& size, const int& C)
{
ofstream file;
file.open(data.c_str(), ios::app);
for (int i = 0; i < size; i++)
{
for (int k = 0; k < C; k++)
{
(k == (C-1) ? file << Mat[i][k] <<endl : file << Mat[i][k] << " " );
}
}
file.close();
}



//Funzione per inserire un colonna di un array in un array monodimensionale
template <typename R1>
void Insert(R1** const Mat,const int& R, R1* const vet, const int& val)
{
for (int i = 0; i < R; i++)
{
vet[i] = Mat[i][val];
}
}




//Funzione per trovare una parola all'interno di un file e riporta la posizione iniziale della riga dove si trova il file
int FIND(const string& nome, const string& parola)
{
ifstream file;
file.open(nome.c_str());
string line;
int posizione = 0;
while (getline(file,line))
{
posizione = file.tellg();
int correttore = 0;
if (line.find(parola,0) != string::npos)
{
correttore = line.size();
posizione -= correttore;
break;
}
}
return posizione;
}



//Funzione per calcolare il numero delle righe
float Righe(const string& nome)
{
string line;
ifstream file;
file.open(nome.c_str());
int l = 0;

while (getline(file,line))
{
if (line.empty() == false)
{
l++;
}
}
file.close();
return l;
}



//Funzione per inizializzare a 0 una matrice
template <typename M>
void MatrixIn(M** const MAT, const int& A , const int& B) {

	for (int i = 0; i < A; i++)
		for (int k = 0; k < B; k++)
			MAT[i][k] = 0;
}




//Funzione per inizializzare vettore
template <typename HH>
void VectorIn(HH* const Vec,const int& A) {

	for (int i = 0; i < A; i++) 
		Vec[i] = 0;
}






//Funzione per inizializzare vettore con valore a piacere
void VectorInG(float* Vec,int A, float B)
{
for (int i = 0; i < A; i++)
{
Vec[i] = B;
}
}




//Funzione per scrivere nel file 
template <typename D>
void File(D** Mat, int size, int C, string nome)
{
FILE* file;
file = fopen(nome.c_str(), "w");
for (int i = 0; i < size; i++)
{
for (int k = 0; k < (C); k++)
{
k == (C-1) ? fprintf(file,"%.5f\n",Mat[i][k]) : fprintf(file,"%.5f ",Mat[i][k]);
}
}
fclose(file);
}


//Funzione per scrivere nel file 
template <typename D1>
void FileE(D1** Mat, const int& size, const int& C, const string& nome)
{
FILE* file;
file = fopen(nome.c_str(), "w");
for (int i = 0; i < size; i++)
{
for (int k = 0; k < (C); k++)
{
k == (C-1) ? fprintf(file,"%.5e\n",Mat[i][k]) : fprintf(file,"%.5e ",Mat[i][k]);
}
}
fclose(file);
}








//Funzione per creare un folder che ha per nome un numero
void mkdirf(float folder)
{
string command;
ostringstream ram;
ram << folder;
const string mkdir = "mkdir ";
command.append(mkdir);
command.append(ram.str());
system(command.c_str());
}


//Funzione per copiare un file(string in una cartella con un numero(folder)
void cpf(string file, float folder)
{
string command;
ostringstream ram;
ram << folder;
const string cp     = "cp ";
const string space  = " ";
command.append(cp);
command.append(file);
command.append(space);
command.append(ram.str());
system(command.c_str());
}



//Funzione per calcolare il numero di cartelle in un path
int Foldercounter(const string& path)
{
const string ls    = "ls -l ";
const string grep  = " | grep -c ^d";
const string slash = "/";
string command;

command.append(ls);
command.append(path);
command.append(slash);
command.append(grep);

FILE* value = popen(command.c_str(), "r");

enum {cinco=5};
char valore[cinco];
fgets(valore,cinco,value);

pclose(value);

const int value2 = atoi(valore);

return value2;
}


/***********************************************
 Function that calculates the number of folder
 in he current folder
***********************************************/
int Foldercounter()
{

//const string path  = system("pwd");
const string ls    = "ls -l ";
const string grep  = " | grep -c ^d";
const string slash = "/";
string command;

command.append(ls);
//command.append(path);
//command.append(slash);
command.append(grep);

FILE* value = popen(command.c_str(), "r");

enum {cinco=5};
char valore[cinco];
fgets(valore,cinco,value);

pclose(value);

const int value2 = atoi(valore);

return value2;
}




/***********************************************
 Function to insert a line of an array into a 
 vector
***********************************************/
template <typename MM>
inline void InsertLine(MM** const Mat, MM* const Vec, const int& col, const int& rig)
{

for (int i = 0; i < col; i++)
Vec[i] = Mat[rig][i];

}


/******************************************/
/* Order a matrix from the smaller value  */
/* the bigger one in respect of one culumn*/
/******************************************/
void MatrixOrderMin(float** const& Atm, const int& righe, const int& col, const int& ord)
{

const int uno = 1;

float* Vec = (float*) malloc(col*sizeof(float));

for (int i = 0; i < righe-uno; i++)
{
for (int k = i+1; k < righe; k++)
{
if (Atm[i][ord] > Atm[k][ord])
{
for (int j = 0; j < col; j++) Vec[j] = Atm[i][j];
for (int j = 0; j < col; j++) Atm[i][j] = Atm[k][j];
for (int j = 0; j < col; j++) Atm[k][j] = Vec[j];
}
}
}

free(Vec);
}



/*******************************************/
/* Order a matrix from the bigger  value   */
/* the smaller one in respect of one culumn*/
/*******************************************/
void MatrixOrderMax(float** const& Atm, const int& righe, const int& col, const int& ord)
{

const int uno = 1;

float* Vec = (float*) malloc(col*sizeof(float));

for (int i = 0; i < righe-uno; i++)
{
for (int k = i+1; k < righe; k++)
{
if (Atm[i][ord] < Atm[k][ord])
{
for (int j = 0; j < col; j++) Vec[j] = Atm[i][j];
for (int j = 0; j < col; j++) Atm[i][j] = Atm[k][j];
for (int j = 0; j < col; j++) Atm[k][j] = Vec[j];
}
}
}

free(Vec);
}





