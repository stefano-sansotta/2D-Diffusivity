/************************************************************

   Library with some helpfull function to make graphs
   using root libraries

   25/October/2016

   released by S. Sansotta
***************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFrame.h"
#include "TAttLine.h"


using namespace std;


class Graph {


public:

//TApplication myapp("app",NULL,NULL);
void SAVEGraph(TCanvas*,const float&,const string&);


void IntMaker(TH1F*,TH1F*,const int&);

void NormalizeHist(TH1F*,const double&);



};



/******************************************
  Function to save graphs
******************************************/
void SAVEGraph(TCanvas* tela, const float& conversion, const string& nome)
{

const int width  = tela->GetWindowWidth();
const int height = tela->GetWindowHeight();
TCanvas* tela2 = new TCanvas("temp", "",width*conversion,height*conversion);
tela->DrawClonePad();
tela2->Draw();
tela2->SaveAs(nome.c_str());
tela2->Close();

}





/*********************************************
 Function to fill and integral
**********************************************/
void IntMaker(TH1F* Int, TH1F* hist,const int& MAX) {

	const int uno   = 1;
	const int dos   = 15;
	const int PUNTI = 200;

	int sum = 0;

	for (int i = 1; i < PUNTI + uno; i++)	{
		sum += hist -> GetBinContent(i);
		Int -> SetBinContent(i,dos*sum/MAX);
	}

}







/*********************************************
 Function to normalize a Histogram
**********************************************/
void NormalizeHist(TH1F *hist, const double& norm) {

	double A = 0;

	for (int i = 0; i < 9; i++) 	{
		A = hist -> GetBinContent(i);
		hist -> SetBinContent(i,A/norm);
	}

}


