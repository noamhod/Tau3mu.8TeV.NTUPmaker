#ifndef types_H
#define types_H

namespace
{
#include "rawStd.h"
#include "rawROOT.h"

typedef map<int, const char*>   TMapiP2cc;
typedef multimap<int, int>      TMapii;
typedef multimap<double, int>   TMMapdi;
typedef multimap<float, unsigned int> TMMapfui;
typedef multimap<int, string>   TMapis;
typedef map<int, double>        TMapid;
typedef map<double, double>     TMapdd;
typedef map<double, vector<int> > TMapdvi;
typedef map<float,  vector<int> > TMapfvi;
typedef vector<TLorentzVector*> TVectorP2VL;
typedef map<TString, TLorentzVector*> TMapP2VL;
typedef map<TString, vector<TLorentzVector>* > TMapP2vVL;
typedef map<TString, TLorentzVector > TMapVL;
typedef map<TString, vector<vector<TLorentzVector> >* > TMapP2vvVL;
typedef map<TString, vector<int> >             TMapTSvi;
typedef map<TString, vector<int>* >            TMapP2vi;
typedef map<TString, TMapP2vi>                 TMap2P2vi;
typedef map<TString, vector<float>* >          TMapP2vf;
typedef map<TString, vector<string>* >         TMapP2vs;
typedef map<TString, TMapP2vf>                 TMap2P2vf;
typedef map<TString, vector<unsigned int>* >   TMapP2vui;
typedef map<TString, vector<unsigned int> >    TMapvui;
typedef map<TString, vector<vector<unsigned int> >* > TMapP2vvui;
typedef map<TString, vector<vector<int> >* >   TMapP2vvi;
typedef map<TString, vector<vector<float> >* > TMapP2vvf;
typedef map<string, TObject*>   TMapSP2TOBJ;
typedef map<TString, TObject*>  TMapTSP2TOBJ;
typedef map<string, TH1*>       TMapSP2TH1;
typedef map<string, vector<TH1*> > TMapSvTH1;
typedef map<string, TH2*>       TMapSP2TH2;
typedef map<string, TGraph*>    TMapSP2TGraph;
typedef map<int, TString>       TMapiTS;
typedef map<TString, TMapiTS>   TMap2iTS;
typedef map<TString, float>     TMapTSf;
typedef map<TString, TMapTSf>   TMap2TSf;
typedef map<TString, string>    TMapTSs;
typedef map<TString, int>       TMapTSi;
typedef map<TString, TMapTSi>   TMap2TSi;
typedef map<TString, bool>      TMapTSb;
typedef map<TString, TMapTSb>   TMap2TSb;
typedef map<TString, unsigned int> TMapTSui;
typedef map<TString, TMapTSui>     TMap2TSui;
typedef map<TString, TMapTSui>     TMap2TSui;
typedef map<unsigned int, TString> TMapuiTS;
typedef map<unsigned int, float>   TMapuif;
typedef map<TString, TString>   TMapTSTS;
typedef map<TString, TFile*>    TMapTSP2TF;
typedef map<TString, TH1*>      TMapTSP2TH1;
typedef map<TString, TH2*>      TMapTSP2TH2;
typedef map<TString, THStack*>  TMapTSP2THS;
typedef map<TString, TGraph*>   TMapTSP2TGraph;
typedef map<TString, TLine*>    TMapTSP2TLINE;
typedef map<TString, TTree*>    TMapTSP2TTREE;
typedef map<TString, TChain*>   TMapTSP2TCHAIN;
typedef map<TString, double>    TMapTSd;
typedef map<TString, TMapTSd>   TMap2TSd;
typedef map<TString, Color_t>   TMapTSTC;
typedef map<TString, Int_t>     TMapTSI;
typedef map<string, TH1D*>      TMapSP2TH1D;
typedef map<TString, TH1D*>     TMapTSP2TH1D;
typedef map<string, TH2D*>      TMapSP2TH2D;
typedef map<string, TCanvas*>   TMapSP2TCNV;
typedef map<TString, TCanvas*>  TMapTSP2TCNV;
typedef map<TString, TGraphAsymmErrors*> TMapTSP2TGAE;
typedef map<TString, TLegend*>  TMapTSP2TLeg;
typedef map<string, double>     TMapsd;
typedef map<string, int>        TMapsi;
typedef map<string, bool>       TMapsb;
typedef map<double, string>     TMapds;
typedef map<double, int>        TMapdi;
typedef map<float, int>         TMapfi;
typedef map<float, unsigned int> TMapfui;
typedef map<unsigned int, vector<float> >  TMapuivf;
typedef map<unsigned int, vector<int> >  TMapuivi;
typedef map<TString, vector<float> >  TMapTSvf;
typedef map<string, vector<double> >  TMapsvd;
typedef map<string, vector<double>* > TMapsP2vd;
typedef map<string, vector<float>* >   TMapsP2vf;
typedef map<TString, vector<float>* >  TMapTSP2vf;
typedef map<TString, vector<string>* > TMapTSP2vs;
typedef map<string, string> TMapss;
typedef map<string, vector<string> > TMapsvs;
typedef map<string, vector<string>* > TMapsP2vs;
typedef map<int, vector<string>* >    TMapiP2vs;


// typedef map<TString, int>                    TMapTSi;
typedef map<TString, vector<int>* >             TMapTSP2vi;
typedef map<TString, vector<vector<int> >* >    TMapTSP2vvi;
// typedef map<TString, double>                 TMapTSd;
typedef map<TString, vector<double>* >          TMapTSP2vd;
typedef map<TString, vector<vector<double> >* > TMapTSP2vvd;
typedef map<TString, vector<string>* >          TMapTSP2vs;
typedef map<TString, vector<vector<string> >* > TMapTSP2vvs;


typedef multimap<float, float > TMultimapff;
typedef vector<TCanvas*> TVecCanvas;
typedef complex<double> dcomplex;
typedef complex<int>    icomplex;
static dcomplex Im = dcomplex(0.,1.);



class fermion
{
	public:
		fermion(string Name, unsigned int PDGid, double Mass, double Charge, double WeakI3)
		{
			this->name = Name;
			this->id = PDGid;
			this->mass = Mass;
			this->charge = Charge;
			this->I3 = WeakI3;
			if(fabs(this->charge)<1.) this->Nc = 3;
			else                      this->Nc = 1;
		}
		~fermion();

	public:
		string name;
		unsigned int id;
		unsigned int Nc;
		double mass; // MeV
		double charge; // in units of the proton's charge
		double I3;
};
typedef map<unsigned int, fermion*> ui2fermion;
typedef map<string, fermion*>       s2fermion;


class GRPX
{
	public:
		GRPX(int ord, TString lbl,
				 Color_t fc, Int_t fs,
				 Color_t lc, Int_t ls,  Width_t lw,
				 Color_t mc, Int_t mst, Float_t mss)
		{
			order = ord;
			label = lbl;
			
			filcolor = fc;
			filstyle = fs;
			
			lincolor = lc;
			linstyle = ls;
			linwidth = lw;
			
			mrkcolor = mc;
			mrkstyle = mst;
			mrksize  = mss;
		}
		~GRPX();
	public:
		int     order; // order of appearance
		TString label;
		
		Color_t filcolor;
		Int_t   filstyle;
		
		Color_t lincolor;
		Int_t   linstyle;
		Width_t linwidth;
		
		Color_t mrkcolor;
		Int_t   mrkstyle;
		Float_t mrksize;
};
typedef map<TString, GRPX*> TMapTS2GRPX;

}

#endif

