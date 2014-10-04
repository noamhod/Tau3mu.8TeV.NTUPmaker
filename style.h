#ifndef TauLFVCommonTools_style_H
#define TauLFVCommonTools_style_H

namespace
{

#include "rawROOT.h"

void rgbPalette(Double_t r, Double_t g, Double_t b, Int_t nb=50)
{
	const UInt_t Number = 3;
	Double_t Red[Number]    = { r, 0.0, 0.0 };
	Double_t Green[Number]  = { g, 0.0, 0.0 };
	Double_t Blue[Number]   = { b, 0.0, 0.0 };
	Double_t Length[Number] = { 0.1, 0.5, 1.0 };
	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
}
void gryPalette(Int_t nb=50)
{
	const UInt_t Number = 2;
	Double_t Red[Number]    = { 1.0, 0.0 };
	Double_t Green[Number]  = { 1.0, 0.0 };
	Double_t Blue[Number]   = { 1.0, 0.0 };
	Double_t Length[Number] = { 0.0, 1.0,};
	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
}

// TExec* exeRed   = new TExec("exeRed",   "rgbPalette(0.7,0,0,50);");
// TExec* exeGreen = new TExec("exeGreen", "rgbPalette(0,0.5,0,50);");
// TExec* exeBlue  = new TExec("exeBlue",  "rgbPalette(0,0.2,0.9,50);");
// TExec* exeGray  = new TExec("exeGray",  "rgbPalette(1,1,1,50);");
// // TExec* exeGray  = new TExec("exeGray",  "gryPalette(50);");

void style(bool dostat=true, bool setpalette=true)
{
	// gStyle->Reset();

	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.13);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadLeftMargin(0.12);
	Int_t font=42;
	Double_t tsize=0.04;
	gStyle->SetTextFont(font);
	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font,"x");
	gStyle->SetTitleFont(font,"x");
	gStyle->SetLabelFont(font,"y");
	gStyle->SetTitleFont(font,"y");
	gStyle->SetLabelFont(font,"z");
	gStyle->SetTitleFont(font,"z");
	gStyle->SetLabelSize(tsize,"x");
	gStyle->SetTitleSize(tsize,"x");
	gStyle->SetLabelSize(tsize,"y");
	gStyle->SetTitleSize(tsize,"y");
	gStyle->SetLabelSize(tsize,"z");
	gStyle->SetTitleSize(tsize,"z");
	if(dostat)
	{
		cout << "STYLE WITH STATBOX" << endl;
		gStyle->SetStatX(0.82);
		gStyle->SetStatY(0.85);
		gStyle->SetStatFont(42);
		gStyle->SetStatFontSize(0.05);
		gStyle->SetOptStat(1110); // The default is (000001111) which means that rms, mean, total of entries plus the histogram are printed. If less than 9 parameters are given as the mode, they are read starting from the end (i.e. histogram name).
		gStyle->SetStatW(0.27);
		gStyle->SetStatH(0.27);
	}
	else
	{
		cout << "STYLE WITHOUT STATBOX" << endl;
		gStyle->SetStatColor(0);
		gStyle->SetStatBorderSize(0);
		gStyle->SetStatColor(0);
		gStyle->SetStatX(0);
		gStyle->SetStatY(0);
		gStyle->SetStatFont(42);
		gStyle->SetStatFontSize(0);
		gStyle->SetOptStat(0);
		gStyle->SetStatW(0);
		gStyle->SetStatH(0);
	}
	gStyle->SetTitleX(0.12);
	gStyle->SetTitleY(1);

	if(setpalette) gStyle->SetPalette(1);

	gStyle->SetTitleX(0.25); //title X location 
	gStyle->SetTitleY(0.94); //title Y location 
	gStyle->SetTitleW(0.5); //title width 
	gStyle->SetTitleH(0.05); //title height
	gStyle->SetTitleBorderSize(0);

	gStyle->cd();
}

void minimalstyle()
{
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.13);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadLeftMargin(0.12);
	Int_t font=42;
	Double_t tsize=0.04;
	gStyle->SetTextFont(font);
	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font,"x");
	gStyle->SetTitleFont(font,"x");
	gStyle->SetLabelFont(font,"y");
	gStyle->SetTitleFont(font,"y");
	gStyle->SetLabelFont(font,"z");
	gStyle->SetTitleFont(font,"z");
	gStyle->SetLabelSize(tsize,"x");
	gStyle->SetTitleSize(tsize,"x");
	gStyle->SetLabelSize(tsize,"y");
	gStyle->SetTitleSize(tsize,"y");
	gStyle->SetLabelSize(tsize,"z");
	gStyle->SetTitleSize(tsize,"z");
	gStyle->SetStatColor(0);
	gStyle->SetStatBorderSize(0);
	gStyle->SetStatColor(0);
	gStyle->SetStatX(0);
	gStyle->SetStatY(0);
	gStyle->SetStatFont(42);
	gStyle->SetStatFontSize(0);
	gStyle->SetOptStat(0);
	gStyle->SetStatW(0);
	gStyle->SetStatH(0);
	gStyle->SetTitleX(0.25); //title X location 
	gStyle->SetTitleY(0.94); //title Y location 
	gStyle->SetTitleW(0.5); //title width 
	gStyle->SetTitleH(0.05); //title height
	gStyle->SetTitleBorderSize(0);
}

}

#endif
