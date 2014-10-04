/*
	usage:
	1. call flatout_init, flatout_book on initialize
	2. call flatout_clear on the beginning of each event
	3. fill the branches in your loop, e.g. flatout_doubles["vtx_pval"] = ...;
	4. call flatout_finit on finalization
*/

#include "TauLFVCommonTools/rawStd.h"
#include "TauLFVCommonTools/rawROOT.h"
#include "TauLFVCommonTools/types.h"
#include "TauLFVCommonTools/enums.h"
#include "TauLFVCommonTools/logs.h"

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
	addBranch("wgt_shapeFONLL", order, FLT);
	addBranch("wgt_normFONLL",  order, FLT);
	addBranch("wgt_luminosity", order, FLT);
	addBranch("wgt_kfactor",    order, FLT);
	addBranch("wgt_dijets",     order, FLT);
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
	addBranch("mu_order1",          order, VINT);
	addBranch("mu_order2",          order, VINT);
	addBranch("mu_order3",          order, VINT);
	addBranch("mu_type1",           order, VINT);
	addBranch("mu_type2",           order, VINT);
	addBranch("mu_type3",           order, VINT);
	addBranch("mu_pt1",             order, VFLT);
	addBranch("mu_pt2",             order, VFLT);
	addBranch("mu_pt3",             order, VFLT);
	addBranch("mu_eta1",            order, VFLT);
	addBranch("mu_eta2",            order, VFLT);
	addBranch("mu_eta3",            order, VFLT);
	addBranch("mu_phi1",            order, VFLT);
	addBranch("mu_phi2",            order, VFLT);
	addBranch("mu_phi3",            order, VFLT);
	addBranch("mu_sctangsig1",      order, VFLT);
	addBranch("mu_sctangsig2",      order, VFLT);
	addBranch("mu_sctangsig3",      order, VFLT);
	addBranch("mu_sctngbsig1",      order, VFLT);
	addBranch("mu_sctngbsig2",      order, VFLT);
	addBranch("mu_sctngbsig3",      order, VFLT);
	addBranch("mu_pbalsig1",        order, VFLT);
	addBranch("mu_pbalsig2",        order, VFLT);
	addBranch("mu_pbalsig3",        order, VFLT);
	addBranch("mu_chi2trkfit1",     order, VFLT);
	addBranch("mu_chi2trkfit2",     order, VFLT);
	addBranch("mu_chi2trkfit3",     order, VFLT);
	addBranch("mu_ndftrkfit1",      order, VFLT);
	addBranch("mu_ndftrkfit2",      order, VFLT);
	addBranch("mu_ndftrkfit3",      order, VFLT);
	addBranch("mu_chi2ndftrkfit1",  order, VFLT);
	addBranch("mu_chi2ndftrkfit2",  order, VFLT);
	addBranch("mu_chi2ndftrkfit3",  order, VFLT);
	addBranch("mu_pvaltrkfit1",     order, VFLT);
	addBranch("mu_pvaltrkfit2",     order, VFLT);
	addBranch("mu_pvaltrkfit3",     order, VFLT);
	addBranch("mu_srcqoverp1",      order, VFLT);
	addBranch("mu_srcqoverp2",      order, VFLT);
	addBranch("mu_srcqoverp3",      order, VFLT);
	addBranch("mu_trkqoverp1",      order, VFLT);
	addBranch("mu_trkqoverp2",      order, VFLT);
	addBranch("mu_trkqoverp3",      order, VFLT);
	addBranch("mu_ptfrac1",         order, VFLT);
	addBranch("mu_ptfrac2",         order, VFLT);
	addBranch("mu_ptfrac3",         order, VFLT);
	addBranch("mu_pixeldEdx1",      order, VFLT);
	addBranch("mu_pixeldEdx2",      order, VFLT);
	addBranch("mu_pixeldEdx3",      order, VFLT);
	addBranch("mu_isMedium1",       order, VINT);
	addBranch("mu_isMedium2",       order, VINT);
	addBranch("mu_isMedium3",       order, VINT);
	addBranch("mu_nPIXhits1",       order, VINT);
	addBranch("mu_nPIXhits2",       order, VINT);
	addBranch("mu_nPIXhits3",       order, VINT);
	addBranch("mu_nDeadPIX1",       order, VINT);
	addBranch("mu_nDeadPIX2",       order, VINT);
	addBranch("mu_nDeadPIX3",       order, VINT);
	addBranch("mu_nPIXholes1",      order, VINT);
	addBranch("mu_nPIXholes2",      order, VINT);
	addBranch("mu_nPIXholes3",      order, VINT);  
	addBranch("mu_nSCThits1",       order, VINT);
	addBranch("mu_nSCThits2",       order, VINT);
	addBranch("mu_nSCThits3",       order, VINT);
	addBranch("mu_nDeadSCT1",       order, VINT);
	addBranch("mu_nDeadSCT2",       order, VINT);
	addBranch("mu_nDeadSCT3",       order, VINT);
	addBranch("mu_nSCTholes1",      order, VINT);
	addBranch("mu_nSCTholes2",      order, VINT);
	addBranch("mu_nSCTholes3",      order, VINT);
	addBranch("mu_nTRThits1",       order, VINT);
	addBranch("mu_nTRThits2",       order, VINT);
	addBranch("mu_nTRThits3",       order, VINT);
	addBranch("mu_nTRToutliers1",   order, VINT);
	addBranch("mu_nTRToutliers2",   order, VINT);
	addBranch("mu_nTRToutliers3",   order, VINT);
	addBranch("mu_htTRThits1",      order, VINT);
	addBranch("mu_htTRThits2",      order, VINT);
	addBranch("mu_htTRThits3",      order, VINT);
	addBranch("mu_nUsedHitsdEdx1",  order, VINT);	
	addBranch("mu_nUsedHitsdEdx2",  order, VINT);	
	addBranch("mu_nUsedHitsdEdx3",  order, VINT);
	
	addBranch("mu_nMDThits1",                   order, VINT);
	addBranch("mu_nMDThits2",                   order, VINT);
	addBranch("mu_nMDThits3",                   order, VINT);
	addBranch("mu_nTGCPhiHits1",                order, VINT);
	addBranch("mu_nTGCPhiHits2",                order, VINT);
	addBranch("mu_nTGCPhiHits3",                order, VINT);
	addBranch("mu_nTGCEtaHits1",                order, VINT);
	addBranch("mu_nTGCEtaHits2",                order, VINT);
	addBranch("mu_nTGCEtaHits3",                order, VINT);
	addBranch("mu_nCSCPhiHits1",                order, VINT);
	addBranch("mu_nCSCPhiHits2",                order, VINT);
	addBranch("mu_nCSCPhiHits3",                order, VINT);
	addBranch("mu_nCSCEtaHits1",                order, VINT);
	addBranch("mu_nCSCEtaHits2",                order, VINT);
	addBranch("mu_nCSCEtaHits3",                order, VINT);
	addBranch("mu_nRPCPhiHits1",                order, VINT);
	addBranch("mu_nRPCPhiHits2",                order, VINT);
	addBranch("mu_nRPCPhiHits3",                order, VINT);
	addBranch("mu_nRPCEtaHits1",                order, VINT);
	addBranch("mu_nRPCEtaHits2",                order, VINT);
	addBranch("mu_nRPCEtaHits3",                order, VINT);
	addBranch("mu_nCSCEtaHoles1",               order, VINT);
	addBranch("mu_nCSCEtaHoles2",               order, VINT);
	addBranch("mu_nCSCEtaHoles3",               order, VINT);
	addBranch("mu_nCSCPhiHoles1",               order, VINT);
	addBranch("mu_nCSCPhiHoles2",               order, VINT);
	addBranch("mu_nCSCPhiHoles3",               order, VINT);
	addBranch("mu_nRPCEtaHoles1",               order, VINT);
	addBranch("mu_nRPCEtaHoles2",               order, VINT);
	addBranch("mu_nRPCEtaHoles3",               order, VINT);
	addBranch("mu_nRPCPhiHoles1",               order, VINT);
	addBranch("mu_nRPCPhiHoles2",               order, VINT);
	addBranch("mu_nRPCPhiHoles3",               order, VINT);
	addBranch("mu_nMDTHoles1",                  order, VINT);
	addBranch("mu_nMDTHoles2",                  order, VINT);
	addBranch("mu_nMDTHoles3",                  order, VINT);
	addBranch("mu_nTGCEtaHoles1",               order, VINT);
	addBranch("mu_nTGCEtaHoles2",               order, VINT);
	addBranch("mu_nTGCEtaHoles3",               order, VINT);
	addBranch("mu_nTGCPhiHoles1",               order, VINT);
	addBranch("mu_nTGCPhiHoles2",               order, VINT);
	addBranch("mu_nTGCPhiHoles3",               order, VINT);
	addBranch("mu_nOutliersOnTrack1",           order, VINT);
	addBranch("mu_nOutliersOnTrack2",           order, VINT);
	addBranch("mu_nOutliersOnTrack3",           order, VINT);
	addBranch("mu_standardDeviationOfChi2OS1",  order, VINT);
	addBranch("mu_standardDeviationOfChi2OS2",  order, VINT);
	addBranch("mu_standardDeviationOfChi2OS3",  order, VINT);
	/*	
	// Reference code
	for( std::set<int>::const_iterator it=sectors.begin();it!=sectors.end();++it )
	{
	        unsigned int x = (m_sectors->size()==0) ? 0 : m_sectors->size()-1; // protection
	        m_sectors->at(x).push_back(*it);
	}
	for( std::map<Muon::IMuonHitSummaryTool::StIndex,Muon::IMuonHitSummaryTool::HitSummary>::const_iterator it=stationLayers.begin();it!=stationLayers.end();++it )
	{
	        unsigned int x = (m_stationName->size()==0) ? 0 : m_stationName->size()-1; // protection
	        m_stationName->at(x).push_back(m_muonHitSummaryTool->stName(it->first));
	        m_nprecisionHits->at(x).push_back(it->second.nprecisionHits);
	        m_netaTriggerLayers->at(x).push_back(it->second.netaTriggerLayers);
	        m_nphiLayers->at(x).push_back(it->second.nphiLayers);
	        m_netaPhiLayers->at(x).push_back(it->second.netaPhiLayers);
	        m_nprecisionHoles->at(x).push_back(it->second.nprecisionHoles);
	        m_netaTriggerHoleLayers->at(x).push_back(it->second.netaTriggerHoleLayers);
	        m_nphiHoleLayers->at(x).push_back(it->second.nphiHoleLayers);
	        m_nprecisionOutliers->at(x).push_back(it->second.nprecisionOutliers);
	        m_nprecisionCloseHits->at(x).push_back(it->second.nprecisionCloseHits);
	}
	for( std::set<Muon::IMuonHitSummaryTool::PhiIndex>::const_iterator it=phiLayers.begin();it!=phiLayers.end();++it )
	{
	        unsigned int x = (m_phiName->size()==0) ? 0 : m_phiName->size()-1; // protection
	        m_phiName->at(x).push_back(m_muonHitSummaryTool->phiName(*it));
	}
	*/
	addBranch("mu_nPrecisionHits1",         order, VINT);
	addBranch("mu_nPrecisionHits2",         order, VINT);
	addBranch("mu_nPrecisionHits3",         order, VINT);
	addBranch("mu_nPhiLayers1",             order, VINT);
	addBranch("mu_nPhiLayers2",             order, VINT);
	addBranch("mu_nPhiLayers3",             order, VINT);
	addBranch("mu_nEtaPhiLayers1",          order, VINT);
	addBranch("mu_nEtaPhiLayers2",          order, VINT);
	addBranch("mu_nEtaPhiLayers3",          order, VINT);
	addBranch("mu_nPrecisionHoles1",        order, VINT);
	addBranch("mu_nPrecisionHoles2",        order, VINT);
	addBranch("mu_nPrecisionHoles3",        order, VINT);
	addBranch("mu_nEtaTriggerHoleLayers1",  order, VINT);
	addBranch("mu_nEtaTriggerHoleLayers2",  order, VINT);
	addBranch("mu_nEtaTriggerHoleLayers3",  order, VINT);
	addBranch("mu_nPhiHoleLayers1",         order, VINT);
	addBranch("mu_nPhiHoleLayers2",         order, VINT);
	addBranch("mu_nPhiHoleLayers3",         order, VINT);
	addBranch("mu_nPrecisionOutliers1",     order, VINT);
	addBranch("mu_nPrecisionOutliers2",     order, VINT);
	addBranch("mu_nPrecisionOutliers3",     order, VINT);
	
	
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
	addBranch("jet_sumpt12",   order, VFLT);
	addBranch("jet_dphi3muJ1", order, VFLT);
	addBranch("jet_dR3muJ1",   order, VFLT);
	addBranch("jet_dphiJ1J2",  order, VFLT);
	addBranch("jet_dRJ1J2",    order, VFLT);
	
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
