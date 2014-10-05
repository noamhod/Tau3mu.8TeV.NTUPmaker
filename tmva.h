/*
	usage:
	1. call initTMVA (on init only)
	2. call bookMVAvars (on init only)
	3. call setMVAvars per event
	4. call setMVAspect per event
	4. call getMVAscore per event
*/

#include "rawStd.h"
#include "rawROOT.h"
#include "types.h"
#include "enums.h"
#include "logs.h"

// #if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
// #endif

using namespace TMVA;

// #ifdef __CINT__
// 	gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
// #endif

namespace
{

TMVA::Reader* reader;    // The Reader object
TMapTSf mva_variables;
TMapTSf mva_spectators;
TMapiTS mva_var_order;
TMapiTS mva_spe_order;
TString methodName;


void addVariable(TString name, int& order)
{
	mva_variables.insert(make_pair(name,-999.));
	mva_var_order.insert(make_pair(order,name));
	order++;
}
void addSpectator(TString name, int& order)
{
	mva_spectators.insert(make_pair(name,-999.));
	mva_spe_order.insert(make_pair(order,name));
	order++;
}

void bookMVAvars()
{
	mva_variables.clear();
	mva_spectators.clear();
	int order = -1;
	
	//// MVA variables
	order = -1;
	addVariable("vtx_pt",order);
	addVariable("vtx_mOS1",order);
	addVariable("vtx_mOS2",order);
	addVariable("jet_pt1",order);
	addVariable("jet_pt2",order);
	addVariable("jet_dphi3muJ1",order);
	addVariable("jet_dphiJ1J2",order);
	addVariable("vtx_isolation",order);
	addVariable("vtx_pval",order);
	addVariable("vtx_lxy",order);
	addVariable("vtx_a0xy",order);
	addVariable("vtx_cosTxy",order);
	addVariable("met_reffinal_et",order);
	addVariable("met_reffinal_mT",order);
	addVariable("met_reffinal_dPhi3mu",order);
	
	
	//// MVA spectators
	order = -1;
	addSpectator("EF_3mu4T",order);
	addSpectator("EF_3mu6",order);
	addSpectator("EF_2mu13",order);
	addSpectator("EF_mu18_tight_mu8_EFFS",order);
	addSpectator("EF_mu18_tight_2mu4_EFFS",order);
	addSpectator("EF_2mu8_EFxe30_tclcw",order);
	addSpectator("EF_mu24_tight_EFxe40",order);
	addSpectator("EF_mu24i_tight",order);
	addSpectator("EF_mu36_tight",order);
	
	addSpectator("vtx_code",order);
	addSpectator("vtx_charge",order);
	addSpectator("vtx_mass",order);
	addSpectator("vtx_mSS",order);
	addSpectator("vtx_rapidity",order);
	addSpectator("vtx_chi2",order);
	addSpectator("vtx_ndf",order);
	addSpectator("vtx_chi2ndf",order);
	addSpectator("vtx_lxySig",order);
	addSpectator("vtx_a0",order);
	addSpectator("vtx_cosT",order);
	addSpectator("vtx_dpt12",order);
	addSpectator("vtx_dpt23",order);
	addSpectator("vtx_dpt13",order);
	addSpectator("vtx_ptfrac12",order);
	addSpectator("vtx_ptfrac23",order);
	addSpectator("vtx_ptfrac13",order);
	addSpectator("vtx_dRmax",order);
	addSpectator("vtx_dRmin",order);
	
	addSpectator("mu_order1",order);
	addSpectator("mu_order2",order);
	addSpectator("mu_order3",order);
	addSpectator("mu_type1",order);
	addSpectator("mu_type2",order);
	addSpectator("mu_type3",order);
	addSpectator("mu_pt1",order);
	addSpectator("mu_pt2",order);
	addSpectator("mu_pt3",order);
	addSpectator("mu_eta1",order);
	addSpectator("mu_eta2",order);
	addSpectator("mu_eta3",order);
	addSpectator("mu_phi1",order);
	addSpectator("mu_phi2",order);
	addSpectator("mu_phi3",order);
	addSpectator("mu_sctangsig1",order);
	addSpectator("mu_sctangsig2",order);
	addSpectator("mu_sctangsig3",order);
	addSpectator("mu_sctngbsig1",order);
	addSpectator("mu_sctngbsig2",order);
	addSpectator("mu_sctngbsig3",order);
	addSpectator("mu_pbalsig1",order);
	addSpectator("mu_pbalsig2",order);
	addSpectator("mu_pbalsig3",order);
	addSpectator("mu_chi2trkfit1",order);
	addSpectator("mu_chi2trkfit2",order);
	addSpectator("mu_chi2trkfit3",order);
	addSpectator("mu_ndftrkfit1",order);
	addSpectator("mu_ndftrkfit2",order);
	addSpectator("mu_ndftrkfit3",order);
	addSpectator("mu_chi2ndftrkfit1",order);
	addSpectator("mu_chi2ndftrkfit2",order);
	addSpectator("mu_chi2ndftrkfit3",order);
	addSpectator("mu_pvaltrkfit1",order);
	addSpectator("mu_pvaltrkfit2",order);
	addSpectator("mu_pvaltrkfit3",order);
	addSpectator("mu_srcqoverp1",order);
	addSpectator("mu_srcqoverp2",order);
	addSpectator("mu_srcqoverp3",order);
	addSpectator("mu_trkqoverp1",order);
	addSpectator("mu_trkqoverp2",order);
	addSpectator("mu_trkqoverp3",order);
	addSpectator("mu_ptfrac1",order);
	addSpectator("mu_ptfrac2",order);
	addSpectator("mu_ptfrac3",order);
	addSpectator("mu_htTRThits1",order);
	addSpectator("mu_htTRThits2",order);
	addSpectator("mu_htTRThits3",order);
	
	addSpectator("jet_pt1",order);
	addSpectator("jet_pt2",order);
	addSpectator("jet_eta1",order);
	addSpectator("jet_eta2",order);
	addSpectator("jet_phi1",order);
	addSpectator("jet_phi2",order);
	addSpectator("jet_m1",order);
	addSpectator("jet_m2",order);
	addSpectator("jet_E1",order);
	addSpectator("jet_E2",order);
	addSpectator("jet_MV1w1",order);
	addSpectator("jet_MV1w2",order);
	addSpectator("jet_sumpt12",order);
	addSpectator("jet_dR3muJ1",order);
	addSpectator("jet_dRJ1J2",order);
}
void initTMVA(TString weightsfile, TString method)
{
	TMVA::Tools::Instance(); // This loads the library
	
	// Create the Reader object
	reader = new TMVA::Reader( "!Color:!Silent" );
	
	// Create a set of variables and declare them to the reader
	// the variable names MUST corresponds in name and type to those given in the weight file(s) used
	
	for(TMapiTS::iterator it= mva_var_order.begin() ; it!=mva_var_order.end() ; ++it)
	{
		TString name = it->second;
		reader->AddVariable(name,&mva_variables[name]);
	}
	
	// Spectator variables declared in the training have to be added to the reader, too
	for(TMapiTS::iterator it= mva_spe_order.begin() ; it!=mva_spe_order.end() ; ++it)
	{
		TString name = it->second;
		reader->AddSpectator(name,&mva_spectators[name]);
	}
	
	// Book MVA method
	methodName = method+" method";
	reader->BookMVA(methodName,weightsfile);
}
void setMVAvars(unsigned int vtx,
				TMapTSi& input_ints, TMapTSf& input_floats, 
				TMapTSP2vi& input_vints, TMapTSP2vf& input_vfloats)
{
	// cout << "MVA vars: vtx=" << vtx << endl;
	// cout << "MVA vars: input_vfloats->size()=" << input_vfloats["vtx_pt"]->size() << endl;
	
	
	mva_variables["vtx_pt"]               = input_vfloats["vtx_pt"]->at(vtx);
	mva_variables["vtx_mOS1"]             = input_vfloats["vtx_mOS1"]->at(vtx);
	mva_variables["vtx_mOS2"]             = input_vfloats["vtx_mOS2"]->at(vtx);
	mva_variables["jet_pt1"]              = input_vfloats["jet_pt1"]->at(vtx);
	mva_variables["jet_pt2"]              = input_vfloats["jet_pt2"]->at(vtx);
	mva_variables["jet_dphi3muJ1"]        = input_vfloats["jet_dphi3muJ1"]->at(vtx);
	mva_variables["jet_dphiJ1J2"]         = input_vfloats["jet_dphiJ1J2"]->at(vtx);
	mva_variables["vtx_isolation"]        = input_vfloats["vtx_isolation"]->at(vtx);
	mva_variables["vtx_pval"]             = input_vfloats["vtx_pval"]->at(vtx);
	mva_variables["vtx_lxy"]              = input_vfloats["vtx_lxy"]->at(vtx);
	mva_variables["vtx_a0xy"]             = input_vfloats["vtx_a0xy"]->at(vtx);
	mva_variables["vtx_cosTxy"]           = input_vfloats["vtx_cosTxy"]->at(vtx);
	mva_variables["met_reffinal_et"]      = input_floats["met_reffinal_et"];
	mva_variables["met_reffinal_mT"]      = input_vfloats["met_reffinal_mT"]->at(vtx);
	mva_variables["met_reffinal_dPhi3mu"] = input_vfloats["met_reffinal_dPhi3mu"]->at(vtx);
}
void setMVAspect(unsigned int vtx,
				TMapTSi& input_ints, TMapTSf& input_floats,
				TMapTSP2vi& input_vints, TMapTSP2vf& input_vfloats)
{
	// cout << "MVA spectators: vtx=" << vtx << endl;
	// cout << "MVA spectators: input_vints->size()=" << input_vints["vtx_code"]->size() << endl;
	// cout << "MVA spectators: input_vfloats->size()=" << input_vfloats["vtx_charge"]->size() << endl;
	
	mva_spectators["EF_3mu4T"]                = input_ints["EF_3mu4T"];
	mva_spectators["EF_3mu6"]                 = input_ints["EF_3mu6"];
	mva_spectators["EF_2mu13"]                = input_ints["EF_2mu13"];
	mva_spectators["EF_mu18_tight_mu8_EFFS"]  = input_ints["EF_mu18_tight_mu8_EFFS"];
	mva_spectators["EF_mu18_tight_2mu4_EFFS"] = input_ints["EF_mu18_tight_2mu4_EFFS"];
	mva_spectators["EF_2mu8_EFxe30_tclcw"]    = input_ints["EF_2mu8_EFxe30_tclcw"];
	mva_spectators["EF_mu24_tight_EFxe40"]    = input_ints["EF_mu24_tight_EFxe40"];
	mva_spectators["EF_mu24i_tight"]          = input_ints["EF_mu24i_tight"];
	mva_spectators["EF_mu36_tight"]           = input_ints["EF_mu36_tight"];
	
	mva_spectators["vtx_code"]     = input_vints["vtx_code"]->at(vtx);
	mva_spectators["vtx_charge"]   = input_vfloats["vtx_charge"]->at(vtx);
	mva_spectators["vtx_mass"]     = input_vfloats["vtx_mass"]->at(vtx);
	mva_spectators["vtx_mSS"]      = input_vfloats["vtx_mSS"]->at(vtx);
	mva_spectators["vtx_rapidity"] = input_vfloats["vtx_rapidity"]->at(vtx);
	mva_spectators["vtx_chi2"]     = input_vfloats["vtx_chi2"]->at(vtx);
	mva_spectators["vtx_ndf"]      = input_vfloats["vtx_ndf"]->at(vtx);
	mva_spectators["vtx_chi2ndf"]  = input_vfloats["vtx_chi2ndf"]->at(vtx);
	mva_spectators["vtx_lxySig"]   = input_vfloats["vtx_lxySig"]->at(vtx);
	mva_spectators["vtx_a0"]       = input_vfloats["vtx_a0"]->at(vtx);
	mva_spectators["vtx_cosT"]     = input_vfloats["vtx_cosT"]->at(vtx);
	mva_spectators["vtx_dpt12"]    = input_vfloats["vtx_dpt12"]->at(vtx);
	mva_spectators["vtx_dpt23"]    = input_vfloats["vtx_dpt23"]->at(vtx);
	mva_spectators["vtx_dpt13"]    = input_vfloats["vtx_dpt13"]->at(vtx);
	mva_spectators["vtx_ptfrac12"] = input_vfloats["vtx_ptfrac12"]->at(vtx);
	mva_spectators["vtx_ptfrac23"] = input_vfloats["vtx_ptfrac23"]->at(vtx);
	mva_spectators["vtx_ptfrac13"] = input_vfloats["vtx_ptfrac13"]->at(vtx);
	mva_spectators["vtx_dRmax"]    = input_vfloats["vtx_dRmax"]->at(vtx);
	mva_spectators["vtx_dRmin"]    = input_vfloats["vtx_dRmin"]->at(vtx);
	
	mva_spectators["mu_order1"]         = input_vints["mu_order1"]->at(vtx);
	mva_spectators["mu_order2"]         = input_vints["mu_order2"]->at(vtx);
	mva_spectators["mu_order3"]         = input_vints["mu_order3"]->at(vtx);
	mva_spectators["mu_type1"]          = input_vints["mu_type1"]->at(vtx);
	mva_spectators["mu_type2"]          = input_vints["mu_type2"]->at(vtx);
	mva_spectators["mu_type3"]          = input_vints["mu_type3"]->at(vtx);
	mva_spectators["mu_pt1"]            = input_vfloats["mu_pt1"]->at(vtx);
	mva_spectators["mu_pt2"]            = input_vfloats["mu_pt2"]->at(vtx);
	mva_spectators["mu_pt3"]            = input_vfloats["mu_pt3"]->at(vtx);
	mva_spectators["mu_eta1"]           = input_vfloats["mu_eta1"]->at(vtx);
	mva_spectators["mu_eta2"]           = input_vfloats["mu_eta2"]->at(vtx);
	mva_spectators["mu_eta3"]           = input_vfloats["mu_eta3"]->at(vtx);
	mva_spectators["mu_phi1"]           = input_vfloats["mu_phi1"]->at(vtx);
	mva_spectators["mu_phi2"]           = input_vfloats["mu_phi2"]->at(vtx);
	mva_spectators["mu_phi3"]           = input_vfloats["mu_phi3"]->at(vtx);
	mva_spectators["mu_sctangsig1"]     = input_vfloats["mu_sctangsig1"]->at(vtx);
	mva_spectators["mu_sctangsig2"]     = input_vfloats["mu_sctangsig2"]->at(vtx);
	mva_spectators["mu_sctangsig3"]     = input_vfloats["mu_sctangsig3"]->at(vtx);
	mva_spectators["mu_sctngbsig1"]     = input_vfloats["mu_sctngbsig1"]->at(vtx);
	mva_spectators["mu_sctngbsig2"]     = input_vfloats["mu_sctngbsig2"]->at(vtx);
	mva_spectators["mu_sctngbsig3"]     = input_vfloats["mu_sctngbsig3"]->at(vtx);
	mva_spectators["mu_pbalsig1"]       = input_vfloats["mu_pbalsig1"]->at(vtx);
	mva_spectators["mu_pbalsig2"]       = input_vfloats["mu_pbalsig2"]->at(vtx);
	mva_spectators["mu_pbalsig3"]       = input_vfloats["mu_pbalsig3"]->at(vtx);
	mva_spectators["mu_chi2trkfit1"]    = input_vfloats["mu_chi2trkfit1"]->at(vtx);
	mva_spectators["mu_chi2trkfit2"]    = input_vfloats["mu_chi2trkfit2"]->at(vtx);
	mva_spectators["mu_chi2trkfit3"]    = input_vfloats["mu_chi2trkfit3"]->at(vtx);
	mva_spectators["mu_ndftrkfit1"]     = input_vfloats["mu_ndftrkfit1"]->at(vtx);
	mva_spectators["mu_ndftrkfit2"]     = input_vfloats["mu_ndftrkfit2"]->at(vtx);
	mva_spectators["mu_ndftrkfit3"]     = input_vfloats["mu_ndftrkfit3"]->at(vtx);
	mva_spectators["mu_chi2ndftrkfit1"] = input_vfloats["mu_chi2ndftrkfit1"]->at(vtx);
	mva_spectators["mu_chi2ndftrkfit2"] = input_vfloats["mu_chi2ndftrkfit2"]->at(vtx);
	mva_spectators["mu_chi2ndftrkfit3"] = input_vfloats["mu_chi2ndftrkfit3"]->at(vtx);
	mva_spectators["mu_pvaltrkfit1"]    = input_vfloats["mu_pvaltrkfit1"]->at(vtx);
	mva_spectators["mu_pvaltrkfit2"]    = input_vfloats["mu_pvaltrkfit2"]->at(vtx);
	mva_spectators["mu_pvaltrkfit3"]    = input_vfloats["mu_pvaltrkfit3"]->at(vtx);
	mva_spectators["mu_srcqoverp1"]     = input_vfloats["mu_srcqoverp1"]->at(vtx);
	mva_spectators["mu_srcqoverp2"]     = input_vfloats["mu_srcqoverp2"]->at(vtx);
	mva_spectators["mu_srcqoverp3"]     = input_vfloats["mu_srcqoverp3"]->at(vtx);
	mva_spectators["mu_trkqoverp1"]     = input_vfloats["mu_trkqoverp1"]->at(vtx);
	mva_spectators["mu_trkqoverp2"]     = input_vfloats["mu_trkqoverp2"]->at(vtx);
	mva_spectators["mu_trkqoverp3"]     = input_vfloats["mu_trkqoverp3"]->at(vtx);
	mva_spectators["mu_ptfrac1"]        = input_vfloats["mu_ptfrac1"]->at(vtx);
	mva_spectators["mu_ptfrac2"]        = input_vfloats["mu_ptfrac2"]->at(vtx);
	mva_spectators["mu_ptfrac3"]        = input_vfloats["mu_ptfrac3"]->at(vtx);
	mva_spectators["mu_htTRThits1"]     = input_vints["mu_htTRThits1"]->at(vtx);
	mva_spectators["mu_htTRThits2"]     = input_vints["mu_htTRThits2"]->at(vtx);
	mva_spectators["mu_htTRThits3"]     = input_vints["mu_htTRThits3"]->at(vtx);
	
	mva_spectators["jet_pt1"]     = input_vfloats["jet_pt1"]->at(vtx);
	mva_spectators["jet_pt2"]     = input_vfloats["jet_pt2"]->at(vtx);
	mva_spectators["jet_eta1"]    = input_vfloats["jet_eta1"]->at(vtx);
	mva_spectators["jet_eta2"]    = input_vfloats["jet_eta2"]->at(vtx);
	mva_spectators["jet_phi1"]    = input_vfloats["jet_phi1"]->at(vtx);
	mva_spectators["jet_phi2"]    = input_vfloats["jet_phi2"]->at(vtx);
	mva_spectators["jet_m1"]      = input_vfloats["jet_m1"]->at(vtx);
	mva_spectators["jet_m2"]      = input_vfloats["jet_m2"]->at(vtx);
	mva_spectators["jet_E1"]      = input_vfloats["jet_E1"]->at(vtx);
	mva_spectators["jet_E2"]      = input_vfloats["jet_E2"]->at(vtx);
	mva_spectators["jet_MV1w1"]   = input_vfloats["jet_MV1w1"]->at(vtx);
	mva_spectators["jet_MV1w2"]   = input_vfloats["jet_MV1w2"]->at(vtx);
	mva_spectators["jet_sumpt12"] = input_vfloats["jet_sumpt12"]->at(vtx);
	mva_spectators["jet_dR3muJ1"] = input_vfloats["jet_dR3muJ1"]->at(vtx);
	mva_spectators["jet_dRJ1J2"]  = input_vfloats["jet_dRJ1J2"]->at(vtx);
}
Double_t getMVAscore()
{	
	return reader->EvaluateMVA(methodName);
}

}
