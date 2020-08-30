#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>


using namespace std;


class Gauss{

	int PointGauss(const string);
	void GausstoVec(float*,float*,const int,const int,float,float,float,const string);



};


/***************************************
  Function to extract number of steps
     from gaussian output file
***************************************/
int PointGauss(const string nome) {

	int point;
	const string str = "A total of";
	int posizione = FIND(nome,str);
	ifstream file (nome.c_str());
	file.seekg(posizione-1);
	string A,B,C;
	file >> A >> B >> C >> point;
	file.close();
	
	return point;
}




/***********************************************
   Function to get points from gaussian output
   file
***********************************************/
void GausstoVec(float* X, float* Y, const int min, const int max, float diff, float Ener0, float conv, const string nome)
{

	const string str  = "Summary of the potential surface scan";
	ifstream file (nome.c_str());
	string line;
	while (getline(file,line)) {
		float A,B,C,D;
		if (line.find(str,0) != string::npos) {
			for (int i = 0; i < 2; i++) getline(file,line);
			for (int i = min; i < max; i++)	{
				file >> A >> B >> C >> D;
				X[i] = B - diff;
				Y[i] = (D - Ener0)*conv;
			}
		}
	}
	file.close();

}




















