/****************************************************

  Script to calculate the mean square displacement 
  from an xyz file. The MSD is calculated for each 
  direction ( x, y, z) and in a specific range of 
  distance from a particle of a surface.

  14/September/2018
  Developed by S. Sansotta
****************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include "Librerie/postprocess.h"
#include "Librerie/fileedit.h"
#include "Librerie/lammpslib.h"


using namespace std;


/****************************************************
#                                                   #
#                LIST OF FUNCTIONS                  #
#                                                   #
****************************************************/

inline void Atoms_Mover(float** const,const float* const,const float&,const int&);
inline void Distance_Check(float** const,float** const,float** const,const int&,const int&,const int&,const int&);
inline void ZeroAdjuster(float** const,float** const,float** const,const int&);
inline void BordersCheck(float** const,float** const,float* const,const int&);
inline void MSD_average(float** const,float* const,const int&,const int&);

/****************************************************
#                                                   #
#                       MAIN                        #
#                                                   #
****************************************************/



int main(int argc, char** argv) {


	const string nome = argv[1];            //Name of the xyz file 
	const float DMin  = atof(argv[2]);      //Minimum distance at which i consider the atom
	const float DMax  = atof(argv[3]);      //Maximum "          "        "       "
	const string Atom = argv[4];            //Atom to take under consideration
	const string Surface_Atom = argv[5];    //Surface atom
	VERIF(nome);

	enum {uno=1};
	enum {drei=3};
	enum {vier=4};
	enum {sex=6};
	enum {sieben=7};
	enum {otto=8};

	const string str          = "Timestep:";

	const string MSDX = "DiffusioneX.txt";
	const string MSDY = "DiffusioneY.txt";
	const string MSDZ = "DiffusioneZ.txt";

	const int Atoms = AtomNumbXYZ(nome);
	const int Steps = StepNumbXYZ(nome); 

	const int O2 = AtomNumb(Atom,nome,Atoms);
	const int Si = AtomNumb(Surface_Atom,nome,Atoms);

	if ( O2 <= 0 ) {
		cout << "Attenzione il numero di atomi selezionati e uguale a zero" <<endl;
		exit(EXIT_FAILURE);
	}
	
	const float ZMin = -5;

	float* BoxB = (float*) malloc(sex*sizeof(float));
	BorderCalculate(nome,BoxB);

	float* MSD_Final_X = (float*) malloc(Steps*sizeof(float));
	float* MSD_Final_Y = (float*) malloc(Steps*sizeof(float));
	float* MSD_Final_Z = (float*) malloc(Steps*sizeof(float));

	VectorIn(MSD_Final_X,Steps);
	VectorIn(MSD_Final_Y,Steps);
	VectorIn(MSD_Final_Z,Steps);

	float** Point_Zero = (float**) malloc(O2*sizeof(float*));
	for (int i = 0; i < O2; i++)
		Point_Zero[i] = (float*) malloc(otto*sizeof(float));

	float** Point_N = (float**) malloc(O2*sizeof(float*));
        for (int i = 0; i < O2; i++)
                Point_N[i] = (float*) malloc(otto*sizeof(float));

	float** Point_NM1 = (float**) malloc(O2*sizeof(float*));
        for (int i = 0; i < O2; i++)
                Point_NM1[i] = (float*) malloc(otto*sizeof(float));

	float** Surface = (float**) malloc(Si*sizeof(float*));
	for (int i = 0; i < Si; i++)
		Surface[i] = (float*) malloc(drei*sizeof(float));

	float** MSD_X = (float**) malloc(O2*sizeof(float*));
	for (int i = 0; i < O2; i++)
		MSD_X[i] = (float*) malloc(Steps*sizeof(float));

	float** MSD_Y = (float**) malloc(O2*sizeof(float*));
        for (int i = 0; i < O2; i++)
                MSD_Y[i] = (float*) malloc(Steps*sizeof(float));

	float** MSD_Z = (float**) malloc(O2*sizeof(float*));
        for (int i = 0; i < O2; i++)
                MSD_Z[i] = (float*) malloc(Steps*sizeof(float));

	MatrixIn(MSD_X,O2,Steps);
	MatrixIn(MSD_Y,O2,Steps);
	MatrixIn(MSD_Z,O2,Steps);

	const float Size_X = BoxB[1] - BoxB[0];
        const float Size_Y = BoxB[3] - BoxB[2];
	
	FILE* file_X;
	FILE* file_Y;
	FILE* file_Z;
	ifstream file;
	file.open(nome.c_str());

	string line;
	string A;
	float B,C,D;

	int counter_atom = 0;
	int counter_surface = 0;
	int counter_steps = 0;
	int counter_atom_range = 0;
	
	while (getline(file,line)) {
		if (line.find(str, 0) != string::npos) {
			/*For cicle to fill the arrays*/
			for ( int i = 0; i < Atoms; i++) {
				file >> A >> B >> C >> D;	
	
				if ( A == Atom && counter_steps == 0) {
					Point_Zero[counter_atom][0] = B;
					Point_Zero[counter_atom][1] = C;
					Point_Zero[counter_atom][2] = D;				
					Point_Zero[counter_atom][3] = 0;
					Point_Zero[counter_atom][4] = 0;
					Point_Zero[counter_atom][5] = 0;
					Point_Zero[counter_atom][6] = 0;
					Point_Zero[counter_atom][7] = 0;
					counter_atom++;

				} else if ( A == Atom && counter_steps > 0 ) {
					Point_N[counter_atom][0] = B;
					Point_N[counter_atom][1] = C;
					Point_N[counter_atom][2] = D;
					counter_atom++;
	
				} else if ( A == Surface_Atom ) { 
					Surface[counter_surface][0] = B;
					Surface[counter_surface][1] = C;
					Surface[counter_surface][2] = D;
					counter_surface++;
	
				}
			
			}

			/*Setting counters to zero*/
			counter_atom = 0;
			counter_surface = 0;


			/*Initialization*/
			if ( counter_steps == 0 ) 
				Atoms_Mover(Point_Zero,BoxB,ZMin,O2);
				cout << "Step 0" <<endl;
			else {
				Atoms_Mover(Point_N,BoxB,ZMin,O2);
				Distance_Check(Point_N,Point_NM1,Surface,O2,Si,DMin,DMax);
			}
		

			/*Calculation of the MSD*/
			if ( counter_steps == 1 ) {
  			        ZeroAdjuster(Point_Zero,Point_N,Point_NM1,O2);
				BordersCheck(Point_N,Point_Zero,BoxB,O2);
				for (int i = 0; i < O2; i++) {
					if (Point_N[i][3] == uno) {
						MSD_X[i][uno] = 0;
						MSD_Y[i][uno] = 0;
						MSD_Z[i][uno] = 0;
					}
				}


			} else if ( counter_steps > 1 ) {
				ZeroAdjuster(Point_Zero,Point_N,Point_NM1,O2);
				BordersCheck(Point_N,Point_NM1,BoxB,O2);
				for (int i = 0; i < O2; i++) {
					if ( Point_N[i][3] == uno || Point_NM1[i][3] == uno ) {
						if ( Point_Zero[i][3] == uno && Point_N[i][4] > uno ) {
						const int real_Step = Point_N[i][4];
						const float Delta_X = pow((Point_N[i][0] + Point_N[i][5]*Size_X) - Point_Zero[i][0],2);
						const float Delta_Y = pow((Point_N[i][1] + Point_N[i][6]*Size_Y) - Point_Zero[i][1],2);
						const float Delta_Z = pow(Point_N[i][2] - Point_Zero[i][2],2);
						MSD_X[i][real_Step] = Delta_X;
                                                MSD_Y[i][real_Step] = Delta_Y;
                                                MSD_Z[i][real_Step] = Delta_Z;
						}
					} 
				}
			}	
		
			/*Change N = N-1*/
			for (int i = 0; i < O2; i++)
				for (int j = 0; j < otto; j++)
					Point_NM1[i][j] = Point_N[i][j];
		
			counter_steps++;
		}		


	}
	file.close();
	
	MSD_average(MSD_X,MSD_Final_X,O2,Steps);
	MSD_average(MSD_Y,MSD_Final_Y,O2,Steps);
	MSD_average(MSD_Z,MSD_Final_Z,O2,Steps);

	
	ofstream filex,filey,filez;
	filex.open(MSDX.c_str());
	filey.open(MSDY.c_str());
	filez.open(MSDZ.c_str());

	for ( int i = 0; i < Steps; i++) {
		filex << MSD_Final_X[i] <<endl;
		filey << MSD_Final_Y[i] <<endl;
		filez << MSD_Final_Z[i] <<endl;
	}
	filex.close();
	filey.close();
	filez.close();
	



	return 0;

}


/****************************************************
#                                                   #
#                    END MAIN                       #
#                                                   #
****************************************************/



/*******************************************
 Function to check if the atoms of a matrix 
 (8 culumn) are inside a specific renge of 
 distance from a surface/particle
*******************************************/
inline void Distance_Check(float** const Atoms, float** const Atoms2, float** const Surface, const int& O2, const int& Si, const int& Min, const int& Max) {

	const int Uno = 1;

	for (int i = 0; i < O2; i++) {
		const float Dist = Distance_Min(Surface,Atoms[i][0],Atoms[i][1],Atoms[i][2],Si);
		if ( Dist >= Min && Dist <= Max ) {
			Atoms[i][3] = Uno;

		} else Atoms[i][3] = 0;
	}

	for (int i = 0; i < O2; i++) 
		if ( Atoms[i][3] == Uno || Atoms2[i][3] == Uno )
			Atoms[i][4]++;
	

}

/***********************************************
  Function to move the lower atoms to the upper 
  side of the simulation box
***********************************************/
inline void Atoms_Mover(float** const Atoms, const float* const BoxB, const float& ZMin, const int& O2) {

	const float Z_Size = BoxB[5] - BoxB[4];

	for (int i = 0; i < O2; i++) {
		if (Atoms[i][2] < ZMin) {
			const float Difference = abs(ZMin - Atoms[i][2]);
			Atoms[i][2] += Z_Size + Difference;
		}
	}
}





/*************************************************
  Function to change the Point_Zero once the 
  molecule is back inside again or out...
*************************************************/
inline void ZeroAdjuster(float** const Point_Zero, float** const Point_N, float** const Point_NM1, const int& O2) {

	for (int i = 0; i < O2; i++) {
		if (Point_N[i][3] == 1 && Point_NM1[i][3] == 0 && Point_Zero[i][3] == 0) {
			Point_Zero[i][0] = Point_N[i][0];
			Point_Zero[i][1] = Point_N[i][1];
			Point_Zero[i][2] = Point_N[i][2];
			Point_Zero[i][3] = Point_N[i][3];

		} else if (Point_N[i][3] == 0 && Point_NM1[i][3] == 0) {
			Point_N[i][4] = Point_NM1[i][4] = 0;
			Point_N[i][5] = Point_NM1[i][5] = 0;
			Point_N[i][6] = Point_NM1[i][6] = 0;
			Point_N[i][7] = Point_NM1[i][7] = 0;
			Point_Zero[i][3] = 0;
		}
	}
}





/*****************************************************
  THe function check the if the molecules pass the 
  borders and change the flags in the Point_N array
*****************************************************/
inline void BordersCheck(float** const Point_N, float** const Point_NM1, /* float** const Point_Zero, */ float* const BoxB, const int& O2) {

	const float Size_X = BoxB[1] - BoxB[0];
	const float Size_Y = BoxB[3] - BoxB[2];

	const int Tollerance = 5;
	const float Error = 0.5;
	
	const float Border_X_L = BoxB[0] - Error;
	const float Border_X_H = BoxB[1] + Error;
        const float Border_Y_L = BoxB[2] - Error;
        const float Border_Y_H = BoxB[3] + Error;

	for (int i = 0; i < O2; i++) {
		if ( (Point_N[i][3] == 1 || Point_NM1[i][3] == 1) && Point_N[i][4] > 1) {
			/* Part for the bordes in X direction */
			if ( (Point_N[i][0] + (Point_N[i][5]*Size_X)) - (Point_NM1[i][0] + (Point_NM1[i][5]*Size_X)) > Size_X/2 ) {
				if ( Point_N[i][0] < (Border_X_L + Tollerance) && Point_N[i][0] > Border_X_L ) { 
					Point_N[i][5]--;
				}
				else if ( Point_N[i][0] > (Border_X_H - Tollerance) && Point_N[i][0] < Border_X_H ) {
					Point_N[i][5]--;
				}
			} else if ( (Point_N[i][0] + (Point_N[i][5]*Size_X)) - (Point_NM1[i][0] + (Point_NM1[i][5]*Size_X)) < -Size_X/2 ) {
				if ( Point_N[i][0] <= (Border_X_L + Tollerance) && Point_N[i][0] >= Border_X_L ) {
					Point_N[i][5]++;
				}
				else if ( Point_N[i][0] > (Border_X_H - Tollerance) && Point_N[i][0] < Border_X_H ) {
					Point_N[i][5]++;
				}
			}
			
			/* Part for the borders in Y direction */ 
			if ( (Point_N[i][1] + (Point_N[i][6]*Size_Y)) - (Point_NM1[i][1] + (Point_NM1[i][6]*Size_Y)) > Size_Y/2 ) {
				if ( Point_N[i][1] < (Border_Y_L + Tollerance) && Point_N[i][1] > Border_Y_L ) 
					Point_N[i][6]--;
				else if ( Point_N[i][1] > (Border_Y_H - Tollerance) && Point_N[i][1] < Border_Y_H )
					Point_N[i][6]--;
			} else if ( (Point_N[i][1] + (Point_N[i][6]*Size_Y)) - (Point_NM1[i][1] + (Point_NM1[i][6]*Size_Y)) < -Size_Y/2 ) {
				if ( Point_N[i][1] < (Border_Y_L + Tollerance) && Point_N[i][1] > Border_Y_L )
					Point_N[i][6]++;
				else if ( Point_N[i][1] > (Border_Y_H - Tollerance) && Point_N[i][1] < Border_Y_H )
					Point_N[i][6]++;
			}
		}
	}
}





/********************************************************
  Function to make an average through a matrix with 
  several msd point and save it into an single dimension
  array
********************************************************/
inline void MSD_average(float** const MSD, float* const MSD_Finale, const int& O2, const int& Steps) {

	int counter = 0;

	for (int i = 0; i < Steps; i++) {
		for (int j = 0; j < O2; j++) {
			if ( MSD[j][i] > 0 ) {
				MSD_Finale[i] += MSD[j][i];
//				if (MSD[j][i] > 3) cout << MSD[j][i] << "  " << j << "  " << i <<endl;
				counter++;
			}
		}
		if (counter > 0) MSD_Finale[i] /= counter;
		counter = 0;
	}


}








