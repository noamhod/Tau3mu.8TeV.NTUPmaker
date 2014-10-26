/*
	usage:
	1. call flatout_init, flatout_book on initialize
	2. call flatout_clear on the beginning of each event
	3. fill the branches in your loop, e.g. flatout_doubles["vtx_pval"] = ...;
	4. call flatout_finit on finalization
*/

#include "rawStd.h"
#include "rawROOT.h"
#include "types.h"
#include "enums.h"
#include "logs.h"

namespace
{

TMapTSi    flatout_ints;
TMapTSf    flatout_floats;
TMapTSd    flatout_doubles;
TMapTSs    flatout_strings;
TMapTSP2vi flatout_vints;
TMapTSP2vf flatout_vfloats;
TMapTSP2vd flatout_vdoubles;
TMapTSP2vs flatout_vstrings;
TMapiTS    flatout_order;
TMapTSi    flatout_type;
int    flatout_idefault = -9999;
double flatout_ddefault = -9999.;
float  flatout_fdefault = -9999.;

TFile* flatout_file = NULL;
TTree* flatout_tree = NULL;


void flatout_finit(TDirectory* olddir)
{
	_INFO("Writing the the flatout tree");
	
	flatout_file->cd();
	flatout_tree->AutoSave();
	// flatout_tree->Write("", TObject::kOverwrite);
	flatout_file->Write();
	flatout_file->Close();
	delete flatout_file;
	olddir->cd();
}
void flatout_init(TString fname, TString tname, TDirectory* olddir)
{
	_INFO("Opening file for writing the flatout tree: "+(string)fname);
	TString fmode = (flatout_file==NULL) ? "RECREATE" : "UPDATE";
	if(flatout_file==NULL) // this will happen on the first call to flatout_init(...)
	{
		_INFO("Creating non-existing flatout file "+(string)fname);
		_INFO("for writing the flatout tree: "+(string)tname);
	}
	else
	{
		_INFO("Updating existing flatout file "+(string)fname);
		_INFO("for writing the flatout tree: "+(string)tname);
		flatout_finit(olddir); // first write the current tree
	}
	
	// this is universal
	flatout_file = new TFile(fname,fmode);
	flatout_file->cd();
	flatout_tree = new TTree("flatout_"+tname,"flatout_"+tname);
	olddir->cd();
}

void addBranch(TString name, int& order, int type)
{
	flatout_order.insert(make_pair(order,name));
	flatout_type.insert(make_pair(name,type));
	
	if     (type==INT)  flatout_ints.insert(make_pair(name, -9999));
	else if(type==FLT)  flatout_floats.insert(make_pair(name, -9999.));
	else if(type==DBL)  flatout_doubles.insert(make_pair(name, -9999.));
	else if(type==STR)  flatout_strings.insert(make_pair(name, ""));
	else if(type==VINT) flatout_vints.insert(make_pair(name, new vector<int>));
	else if(type==VFLT) flatout_vfloats.insert(make_pair(name, new vector<float>));
	else if(type==VDBL) flatout_vdoubles.insert(make_pair(name, new vector<double>));
	else if(type==VSTR) flatout_vstrings.insert(make_pair(name, new vector<string>));
	// else if(type==VVINT) flatout_doubles.insert(make_pair(name, new vector<vector<int> >));
	// else if(type==VVFLT) flatout_vvfloats.insert(make_pair(name, new vector<vector<float> >));
	// else if(type==VVDBL) flatout_doubles.insert(make_pair(name, new vector<vector<double> >));
	// else if(type==VVSTR) flatout_doubles.insert(make_pair(name, new vector<vector<string> >));
	else _FATAL("unsupported type of branch");
	
	order++;
}

void flatout_book(TDirectory* olddir)
{
	_INFO("Booking the flatout tree branches");
	
	int order = -1;
	
	// event info
	addBranch("evt_RunNumber",         order, INT);
	addBranch("evt_lbn",               order, INT);
	addBranch("evt_EventNumber",       order, INT);
	addBranch("evt_actualIntPerXing",  order, INT);
	addBranch("evt_averageIntPerXing", order, INT);
	
	// weights info
	addBranch("wgt_mcevt",      order, FLT);
	addBranch("wgt_shapeFONLL", order, FLT);
	addBranch("wgt_normFONLL",  order, FLT);
	addBranch("wgt_luminosity", order, FLT);
	addBranch("wgt_kfactor",    order, FLT);
	addBranch("wgt_dijets",     order, FLT);
	addBranch("wgt_pileup",     order, FLT);
	addBranch("wgt_total",      order, FLT);
	
	// trigger info
	addBranch("EF_3mu4T",                order, INT);
	addBranch("EF_3mu6",                 order, INT);
	addBranch("EF_3mu6_MSonly",          order, INT);
	addBranch("EF_2mu13",                order, INT);
	addBranch("EF_mu18_tight_mu8_EFFS",  order, INT);
	addBranch("EF_mu18_tight_2mu4_EFFS", order, INT);
	addBranch("EF_2mu8_EFxe30_tclcw",    order, INT);
	addBranch("EF_mu24_tight_EFxe40",    order, INT);
	addBranch("EF_mu24i_tight",          order, INT);
	addBranch("EF_mu36_tight",           order, INT);
	
	// vertex info
	addBranch("vtx_n",         order, INT);
	addBranch("vtx_type",      order, VSTR);
	addBranch("vtx_code",      order, VINT);
	addBranch("vtx_charge",    order, VFLT);
	addBranch("vtx_mass",      order, VFLT);
	addBranch("vtx_mOS1",      order, VFLT);
	addBranch("vtx_mOS2",      order, VFLT);
	addBranch("vtx_mSS",       order, VFLT);
	addBranch("vtx_mQuad4",    order, VFLT);
	addBranch("vtx_mQuad5",    order, VFLT);
	addBranch("vtx_mQuad6",    order, VFLT);
	addBranch("vtx_pt",        order, VFLT);
	addBranch("vtx_rapidity",  order, VFLT);
	addBranch("vtx_chi2",      order, VFLT);
	addBranch("vtx_ndf",       order, VFLT);
	addBranch("vtx_chi2ndf",   order, VFLT);
	addBranch("vtx_pval",      order, VFLT);
	addBranch("vtx_lxy",       order, VFLT);
	addBranch("vtx_lxySig",    order, VFLT);
	addBranch("vtx_a0",        order, VFLT);
	addBranch("vtx_a0xy",      order, VFLT);
	addBranch("vtx_cosT",      order, VFLT);
	addBranch("vtx_cosTxy",    order, VFLT);
	addBranch("vtx_tau",       order, VFLT);
	addBranch("vtx_ptfrac12",  order, VFLT);
	addBranch("vtx_ptfrac23",  order, VFLT);
	addBranch("vtx_ptfrac13",  order, VFLT);
	addBranch("vtx_dpt12",     order, VFLT);
	addBranch("vtx_dpt23",     order, VFLT);
	addBranch("vtx_dpt13",     order, VFLT);
	addBranch("vtx_dRmax",     order, VFLT);
	addBranch("vtx_dRmin",     order, VFLT);
	addBranch("vtx_isolation000", order, VFLT);
	addBranch("vtx_isolation001", order, VFLT);
	addBranch("vtx_isolation002", order, VFLT);
	addBranch("vtx_isolation003", order, VFLT);
	addBranch("vtx_isolation004", order, VFLT);
	addBranch("vtx_isolation005", order, VFLT);
	addBranch("vtx_isolation006", order, VFLT);
	addBranch("vtx_isolation007", order, VFLT);
	addBranch("vtx_isolation008", order, VFLT);
	addBranch("vtx_isolation009", order, VFLT);
	addBranch("vtx_isolation010", order, VFLT);
	addBranch("vtx_isolation012", order, VFLT);
	addBranch("vtx_isolation014", order, VFLT);
	addBranch("vtx_isolation016", order, VFLT);
	addBranch("vtx_isolation018", order, VFLT);
	addBranch("vtx_isolation020", order, VFLT);
	addBranch("vtx_isolation022", order, VFLT);
	addBranch("vtx_isolation024", order, VFLT);
	addBranch("vtx_isolation026", order, VFLT);
	addBranch("vtx_isolation028", order, VFLT);
	addBranch("vtx_isolation030", order, VFLT);
	
	// MET_RefFinal-based kinematics
	addBranch("met_reffinal_et",      order, FLT);
	addBranch("met_reffinal_phi",     order, FLT);
	addBranch("met_reffinal_mT",      order, VFLT);
	addBranch("met_reffinal_dPhi3mu", order, VFLT);
	
	// Muons
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_order"+x,          order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_type"+x,           order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_pt"+x,             order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_eta"+x,            order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_phi"+x,            order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_sctangsig"+x,      order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_sctngbsig"+x,      order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_pbalsig"+x,        order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_chi2trkfit"+x,     order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_ndftrkfit"+x,      order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_chi2ndftrkfit"+x,  order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_pvaltrkfit"+x,     order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_srcqoverp"+x,      order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_trkqoverp"+x,      order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_ptfrac"+x,         order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_pixeldEdx"+x,      order, VFLT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_isMedium"+x,       order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nPIXhits"+x,       order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nDeadPIX"+x,       order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nPIXholes"+x,      order, VINT); }  
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nSCThits"+x,       order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nDeadSCT"+x,       order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nSCTholes"+x,      order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nTRThits"+x,       order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nTRToutliers"+x,   order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_htTRThits"+x,      order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nUsedHitsdEdx"+x,  order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nMDThits"+x,                   order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nTGCPhiHits"+x,                order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nTGCEtaHits"+x,                order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nCSCPhiHits"+x,                order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nCSCEtaHits"+x,                order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nRPCPhiHits"+x,                order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nRPCEtaHits"+x,                order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nCSCEtaHoles"+x,               order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nCSCPhiHoles"+x,               order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nRPCEtaHoles"+x,               order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nRPCPhiHoles"+x,               order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nMDTHoles"+x,                  order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nTGCEtaHoles"+x,               order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nTGCPhiHoles"+x,               order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nOutliersOnTrack"+x,           order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_standardDeviationOfChi2OS"+x,  order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nPrecisionHits"+x,         order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nPhiLayers"+x,             order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nEtaPhiLayers"+x,          order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nPrecisionHoles"+x,        order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nEtaTriggerHoleLayers"+x,  order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nPhiHoleLayers"+x,         order, VINT); }
	for(int i=1 ; i<=6 ; ++i) { TString x = _s((float)i,0); addBranch("mu_nPrecisionOutliers"+x,     order, VINT); }
	
	
	// Jets
	addBranch("jet_pt1",       order, VFLT);
	addBranch("jet_pt2",       order, VFLT);
	addBranch("jet_pt3",       order, VFLT);
	addBranch("jet_pt4",       order, VFLT);
	addBranch("jet_eta1",      order, VFLT);
	addBranch("jet_eta2",      order, VFLT);
	addBranch("jet_eta3",      order, VFLT);
	addBranch("jet_eta4",      order, VFLT);
	addBranch("jet_phi1",      order, VFLT);
	addBranch("jet_phi2",      order, VFLT);
	addBranch("jet_phi3",      order, VFLT);
	addBranch("jet_phi4",      order, VFLT);
	addBranch("jet_m1",        order, VFLT);
	addBranch("jet_m2",        order, VFLT);
	addBranch("jet_m3",        order, VFLT);
	addBranch("jet_m4",        order, VFLT);
	addBranch("jet_E1",        order, VFLT);
	addBranch("jet_E2",        order, VFLT);
	addBranch("jet_E3",        order, VFLT);
	addBranch("jet_E4",        order, VFLT);
	addBranch("jet_MV1w1",     order, VFLT);
	addBranch("jet_MV1w2",     order, VFLT);
	addBranch("jet_MV1w3",     order, VFLT);
	addBranch("jet_MV1w4",     order, VFLT);
	addBranch("jet_vtxf1",     order, VFLT);
	addBranch("jet_vtxf2",     order, VFLT);
	addBranch("jet_vtxf3",     order, VFLT);
	addBranch("jet_vtxf4",     order, VFLT);
	addBranch("jet_sumpt12",   order, VFLT);
	addBranch("jet_dphi3muJ1", order, VFLT);
	addBranch("jet_dR3muJ1",   order, VFLT);
	addBranch("jet_dphiJ1J2",  order, VFLT);
	addBranch("jet_dRJ1J2",    order, VFLT);
	
	
	// addBranch("mc_channel_number",  order, INT);
	// addBranch("mc_event_number",    order, INT);
	// addBranch("mc_event_weight",    order, FLT);
	// addBranch("mc_pdgId",           order, VINT);
	// addBranch("mc_pt",              order, VFLT);
	// addBranch("mc_eta",             order, VFLT);
	// addBranch("mc_phi",             order, VFLT);
	// addBranch("mc_m",               order, VFLT);
	// addBranch("mc_charge",          order, VFLT);
	// addBranch("mc_prodvtx_x",       order, VFLT);
	// addBranch("mc_prodvtx_y",       order, VFLT);
	// addBranch("mc_prodvtx_z",       order, VFLT);
	// addBranch("mc_decayvtx_z",      order, VFLT);
	// addBranch("mc_decayvtx_z",      order, VFLT);
	// addBranch("mc_decayvtx_z",      order, VFLT);
	// /*
	// vector<int>*          mc_productionvtx_barcode;
	// vector<int>*          mc_decayvtx_barcode;
	// vector<int>*          mc_has_productionvtx;
	// vector<int>*          mc_has_decayvtx;
	// */
	
	
	
	
	// Attach the branches
	flatout_file->cd();
	for(TMapiTS::iterator it=flatout_order.begin() ; it!=flatout_order.end() ; ++it)
	{
		TString name = it->second;
		int type = flatout_type[name];
		if     (type==INT)   flatout_tree->Branch(name, &flatout_ints[name]);
		else if(type==FLT)   flatout_tree->Branch(name, &flatout_floats[name]);
		else if(type==DBL)   flatout_tree->Branch(name, &flatout_doubles[name]);
		else if(type==STR)   flatout_tree->Branch(name, &flatout_strings[name]);
		else if(type==VINT)  flatout_tree->Branch(name, &flatout_vints[name]);
		else if(type==VFLT)  flatout_tree->Branch(name, &flatout_vfloats[name]);
		else if(type==VDBL)  flatout_tree->Branch(name, &flatout_vdoubles[name]);
		else if(type==VSTR)  flatout_tree->Branch(name, &flatout_vstrings[name]);
		//else if(type==VVFLT) flatout_tree->Branch(name, &flatout_vvfloats[name]);
		else _FATAL("unsupported type of branch");
	}
	olddir->cd(); // go back to the old directory
}

// called on event basis
void flatout_clear() 
{
	for(TMapTSi::iterator    it=flatout_ints.begin()     ; it!=flatout_ints.end()     ; ++it) it->second = flatout_idefault;
	for(TMapTSf::iterator    it=flatout_floats.begin()   ; it!=flatout_floats.end()   ; ++it) it->second = flatout_fdefault;
	for(TMapTSd::iterator    it=flatout_doubles.begin()  ; it!=flatout_doubles.end()  ; ++it) it->second = flatout_ddefault;
	for(TMapTSs::iterator    it=flatout_strings.begin()  ; it!=flatout_strings.end()  ; ++it) it->second = "";
	for(TMapTSP2vi::iterator it=flatout_vints.begin()    ; it!=flatout_vints.end()    ; ++it) it->second->clear();
	for(TMapTSP2vf::iterator it=flatout_vfloats.begin()  ; it!=flatout_vfloats.end()  ; ++it) it->second->clear();
	for(TMapTSP2vd::iterator it=flatout_vdoubles.begin() ; it!=flatout_vdoubles.end() ; ++it) it->second->clear();
	for(TMapTSP2vs::iterator it=flatout_vstrings.begin() ; it!=flatout_vstrings.end() ; ++it) it->second->clear();
}

void flatout_fill(int nAll, TDirectory* olddir)
{
	flatout_file->cd();
	flatout_tree->Fill();

	if(nAll%10000==0 && nAll!=0)
	{
		flatout_tree->FlushBaskets();
		// flatout_tree->Write("", TObject::kOverwrite);
	}
	olddir->cd();
}

}
