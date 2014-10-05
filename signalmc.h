/*
	usage:
		1. call getMC
		2. call getMU
		3. call getCAL
		4. call clearTruth
		5. call matchTruth
		6. get the output via: isTruthMatchedTriplet
*/

#include "rawStd.h"
#include "rawROOT.h"
#include "types.h"
#include "enums.h"
#include "constants.h"
#include "logs.h"

namespace
{

TMapdi muMapTruAll;
TMapdi muMapTruAcc;
TMapdi muMapRecMat;
TMapdi muMapTruMat;
TMapdi calMapRecMat;
TMapdi calMapTruMat;
TMapdi trkMapRecMat;
TMapdi trkMapTruMat;
TMapdi mucalMapTruMat;
int itau   = -1;
int inutau = -1;
int nrecmatmuons = 0;
int nrecmatcalos = 0;
int n_trumu_all    = 0;
int n_trumu_acc    = 0;
int n_recmu_mat    = 0;
int n_trumu_mat    = 0;
int n_reccal_mat   = 0;
int n_trucal_mat   = 0;
int n_trumucal_mat = 0;
vector<int> sorted_tru_mu_all(3,-1);    // sort truth muons by pt
vector<int> sorted_tru_mu_acc(3,-1);    // sort tru muons in acceptance
vector<int> sorted_tru_mu_mat(3,-1);    // sort tru muons reco-muon-matched
vector<int> sorted_tru_cal_mat(3,-1);   // sort tru muons reco-calomuon-matched
vector<int> sorted_tru_mucal_mat(3,-1); // sort tru reco-matched-muons || reco-matched-calomu 
vector<int> sorted_rec_mu_mat(3,-1);    // sort rec matched muons by pt
vector<int> sorted_rec_cal_mat(3,-1);   // sort rec matched calo muons by pt



vector<int>*          truth_pdgId    = 0;
vector<int>*          truth_status   = 0;
vector<int>*          truth_barcode  = 0;
vector<vector<int> >* truth_children = 0;
vector<vector<int> >* truth_parents  = 0;
vector<double>*       truth_pt       = 0;
vector<double>*       truth_eta      = 0;
vector<double>*       truth_phi      = 0;

vector<double>* mu_pt  = 0;
vector<double>* mu_eta = 0;
vector<double>* mu_phi = 0;
                        
vector<double>* cal_pt  = 0;
vector<double>* cal_eta = 0;
vector<double>* cal_phi = 0;


// identifiers
bool isNuTau(int pdgId) { return (pdgId==16); }
bool isTau(int pdgId)   { return (pdgId==15 || pdgId==-15); }
bool isMuon(int pdgId)  { return (pdgId==13 || pdgId==-13); } 
bool isW(int pdgId)     { return (pdgId==24 || pdgId==-24); }

// kinematic matching
double dRtru(double eta1, double phi1, double eta2, double phi2) { return sqrt((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2)); }
double dpTreltru(double pTtru, double pTrec)                     { return (pTtru>0.) ? fabs(pTtru-pTrec)/pTtru : +1.e20; }

// help
TLorentzVector getTruthMatchedTriplet(vector<int>& index, vector<double>* pt, vector<double>* eta, vector<double>* phi)
{
	TLorentzVector p1,p2,p3,p;
	p1.SetPtEtaPhiM(pt->at(index[0]),eta->at(index[0]),phi->at(index[0]),muonMassMeV);
	p2.SetPtEtaPhiM(pt->at(index[1]),eta->at(index[1]),phi->at(index[1]),muonMassMeV);
	p3.SetPtEtaPhiM(pt->at(index[2]),eta->at(index[2]),phi->at(index[2]),muonMassMeV);
	p = p1+p2+p3;
	return p;
}


// getters
void getMC(vector<int>* pdgId, vector<int>* status, vector<int>* barcode,
		   vector<vector<int> >* children, vector<vector<int> >* parents,
		   vector<double>* pt, vector<double>* eta, vector<double>* phi)
{
	truth_pdgId    = pdgId;
	truth_status   = status;
	truth_barcode  = barcode;
	truth_children = children;
	truth_parents  = parents;
	truth_pt       = pt;
	truth_eta      = eta;
	truth_phi      = phi;
}
void getMU(vector<double>* pt, vector<double>* eta, vector<double>* phi)
{
	mu_pt  = pt;
	mu_eta = eta;
	mu_phi = phi;
}
void getCAL(vector<double>* pt, vector<double>* eta, vector<double>* phi)
{
	cal_pt  = pt;
	cal_eta = eta;
	cal_phi = phi;
}


// information
void printEvent()
{
	unsigned int nparticles = truth_pdgId->size();
	for(unsigned int imc=0 ; imc<nparticles ; imc++)
	{
		bool istau = (isTau(truth_pdgId->at(imc)) && (truth_status->at(imc)==2 || truth_status->at(imc)==10902));
		if(!istau) continue;

		cout << "["<<imc<<"]pdgId=" << truth_pdgId->at(imc) << ",status=" << truth_status->at(imc) << ",barcode=" << truth_barcode->at(imc);
		int nchildren = truth_children->at(imc).size();
		cout << "\n\t children: ";
		for(int child=0 ; child<nchildren ; child++)
		{
			int ichild = truth_children->at(imc)[child];
			int child_pdgId  = truth_pdgId->at(ichild);
			int child_status = truth_status->at(ichild);
			int child_barcode = truth_barcode->at(ichild);
			cout << " ["<<ichild<<"]id=" << child_pdgId << ",st=" << child_status << ",bc=" << child_barcode;
		}
		int nparents = truth_parents->at(imc).size();
		cout << "\n\t parents: ";
		for(int parent=0 ; parent<nparents ; parent++)
		{
			int iparent = truth_parents->at(imc)[parent];
			int parent_pdgId  = truth_pdgId->at(iparent);
			int parent_status = truth_status->at(iparent);
			int parent_barcode = truth_barcode->at(iparent);
			cout << " ["<<iparent<<"]id=" << parent_pdgId << ",st=" << parent_status << ",bc=" << parent_barcode;
		}
		cout << endl;
	}
}

bool findInMap(TMapdi& m, int val)
{
	for(TMapdi::iterator it=m.begin() ; it!=m.end() ; ++it) if(it->second==val) return true;
	return false;
}
int getOrder(TMapdi& m, int val)
{
	int order = 0;
	for(TMapdi::reverse_iterator rit=m.rbegin() ; rit!=m.rend() ; ++rit)
	{
		if(rit->second==val) break;
		order++;
	}
	if(order<(int)m.size()) return order;
	return -1;
}
int getTruthSignalNeutrinoindex(int signalWindex)
{
	unsigned int nWchildren = truth_children->at(signalWindex).size();
	
	for(unsigned int imc=0 ; imc<nWchildren ; imc++)
	{
		int ichild = truth_children->at(signalWindex)[imc];
		
		bool isnutau = (isNuTau(truth_pdgId->at(ichild)));
		if(!isnutau) continue;
		
		return ichild;
	}
	return -1;
}
bool isInTriggerAcceptance(vector<int>& index, vector<double>* pt, vector<double>* eta, double met)
{
	TMapdi pt2idxMap;
	for(unsigned int i=0 ; i<index.size() ; i++)
	{
		int idx = index[i];
		if(idx<0) continue;
		pt2idxMap.insert(make_pair(pt->at(idx),idx));
	}
	
	unsigned int n = pt2idxMap.size();
	if(n<1) return false;
	
	int i1 = -1;
	int i3 = -1;
	int i2 = -1;
	TMapdi::reverse_iterator rit=pt2idxMap.rbegin();
	if(n>0) { i1 = rit->second; rit++; }
	if(n>1) { i2 = rit->second; rit++; }
	if(n>2)   i3 = rit->second;

	bool EF3mu4        = (n>2) ? (pt->at(i1)>4*GeV2MeV  && fabs(eta->at(i1))<2.5 && pt->at(i2)>4*GeV2MeV  && fabs(eta->at(i2))<2.5 && pt->at(i3)>4*GeV2MeV && fabs(eta->at(i3))<2.5) : 0;
	bool EF3mu6        = (n>2) ? (pt->at(i1)>6*GeV2MeV  && fabs(eta->at(i1))<2.5 && pt->at(i2)>6*GeV2MeV  && fabs(eta->at(i2))<2.5 && pt->at(i3)>6*GeV2MeV && fabs(eta->at(i3))<2.5) : 0;
	bool EF2mu13       = (n>1) ? (pt->at(i1)>13*GeV2MeV && fabs(eta->at(i1))<2.5 && pt->at(i2)>13*GeV2MeV && fabs(eta->at(i2))<2.5) : 0;
	bool EFmu18_mu8FS  = (n>1) ? (pt->at(i1)>18*GeV2MeV && fabs(eta->at(i1))<2.5 && pt->at(i2)>8*GeV2MeV  && fabs(eta->at(i2))<2.5) : 0;
	bool EFmu18_2mu4FS = (n>2) ? (pt->at(i1)>18*GeV2MeV && fabs(eta->at(i1))<2.5 && pt->at(i2)>4*GeV2MeV  && fabs(eta->at(i2))<2.5 && pt->at(i3)>4*GeV2MeV && fabs(eta->at(i3))<2.5) : 0;
	bool EFmu24        = (n>0) ? (pt->at(i1)>24*GeV2MeV && fabs(eta->at(i1))<2.5) : 0;
	bool EFmu36        = (n>0) ? (pt->at(i1)>36*GeV2MeV && fabs(eta->at(i1))<2.5) : 0;
	bool EFmu24_xe40   = (n>0) ? (pt->at(i1)>24*GeV2MeV && fabs(eta->at(i1))<2.5 && met>40*GeV2MeV) : 0;
	bool EF2mu8_xe30   = (n>1) ? (pt->at(i1)>8*GeV2MeV  && fabs(eta->at(i1))<2.5 && pt->at(i2)>8*GeV2MeV && fabs(eta->at(i2))<2.5 && met>30*GeV2MeV) : 0;

	return (EF3mu4 || EF3mu6 || EF2mu13 || EFmu18_mu8FS || EFmu18_2mu4FS || EFmu24 || EFmu36 || EFmu24_xe40 || EF2mu8_xe30);
}
bool isInTriggerAcceptance(vector<int>& index, vector<double>* pt, vector<double>* eta, double acc_pt, double acc_eta)
{
	int npassing = 0;
	for(unsigned int i=0 ; i<index.size() ; i++)
	{
		int idx = index[i];
		if(idx<0) continue;
		if(pt->at(idx)>acc_pt && fabs(eta->at(idx))<acc_eta) npassing++;
	}
	return (npassing>0);
}
bool isInAcceptance(double pt, double eta, double acc_pt=1.*GeV2MeV, double acc_eta=2.5)
{
	return (pt>acc_pt && fabs(eta)<acc_eta);
}
bool isMatched(double dr, double dptrel, double drmatch=0.005, double dptrelmatch=0.5)
{
	return (dr<drmatch && dptrel<dptrelmatch);
}
int match_tru2recmu(int itru, double drmatch=0.005, double dptrelmatch=0.5)
{
	for(int mu=0 ; mu<(int)mu_pt->size() ; mu++)
	{
		double dr = dRtru(mu_eta->at(mu),mu_phi->at(mu),truth_eta->at(itru),truth_phi->at(itru));
		double dptrel = dpTreltru(mu_pt->at(mu),truth_pt->at(itru));
		if(isMatched(dr,dptrel,drmatch,dptrelmatch)) return mu; 
	}
	return -1;
}
int match_tru2reccal(int itru, double drmatch=0.005, double dptrelmatch=0.5)
{	
	for(int cal=0 ; cal<(int)cal_pt->size() ; cal++)
	{	
		double dr = dRtru(cal_eta->at(cal),cal_phi->at(cal),truth_eta->at(itru),truth_phi->at(itru));
		double dptrel = dpTreltru(cal_pt->at(cal),truth_pt->at(itru));
		if(isMatched(dr,dptrel,drmatch,dptrelmatch)) return cal; 
	}
	return -1;
}


void clearTruth()
{
	truth_pdgId    = 0;
	truth_status   = 0;
	truth_barcode  = 0;
	truth_children = 0;
	truth_parents  = 0;
	truth_pt       = 0;
	truth_eta      = 0;
	truth_phi      = 0;
    
	mu_pt  = 0;
	mu_eta = 0;
	mu_phi = 0;
    
	cal_pt  = 0;
	cal_eta = 0;
	cal_phi = 0;
	
	itau   = -1;
	inutau = -1;
	
	nrecmatmuons = 0;
	nrecmatcalos = 0;
	
	n_trumu_all    = 0;
	n_trumu_acc    = 0;
	n_recmu_mat    = 0;
	n_trumu_mat    = 0;
	n_reccal_mat   = 0;
	n_trucal_mat   = 0;
	n_trumucal_mat = 0;
	
	muMapTruAll.clear();
	muMapTruAcc.clear();
	muMapRecMat.clear();
	muMapTruMat.clear();
	calMapRecMat.clear();
	calMapTruMat.clear();
	mucalMapTruMat.clear();

	for(int i=0 ; i<3 ; i++)
	{
		sorted_tru_mu_all[i]    = -1;
		sorted_tru_mu_acc[i]    = -1;
		sorted_tru_mu_mat[i]    = -1;
		sorted_tru_cal_mat[i]   = -1;
		sorted_tru_mucal_mat[i] = -1;
		sorted_rec_mu_mat[i]    = -1;
		sorted_rec_cal_mat[i]   = -1;
	}
}
void matchTruth()
{
	// truth particles
	unsigned int nparticles = truth_pdgId->size();
	for(unsigned int imc=0 ; imc<nparticles ; imc++)
	{
		bool istau = (isTau(truth_pdgId->at(imc)) && (truth_status->at(imc)==2 || truth_status->at(imc)==10902));
		if(!istau) continue;
		
		int nchildren = truth_children->at(imc).size();
		if(nchildren!=3) continue;
		
		int imu1 = truth_children->at(imc)[0];
		int imu2 = truth_children->at(imc)[1];
		int imu3 = truth_children->at(imc)[2];
		bool ismu1 = (isMuon(truth_pdgId->at(imu1)) && truth_status->at(imu1)==1);
		bool ismu2 = (isMuon(truth_pdgId->at(imu2)) && truth_status->at(imu2)==1);
		bool ismu3 = (isMuon(truth_pdgId->at(imu3)) && truth_status->at(imu3)==1);
		if((ismu1+ismu2+ismu3)!=3) continue;

		int nparents = truth_parents->at(imc).size();
		if(nparents!=1) continue;
		
		int iparent = truth_parents->at(imc)[0];
		if(!isW(truth_pdgId->at(iparent))) continue;
		
		inutau = getTruthSignalNeutrinoindex(iparent);
		itau = imc;
		
		/////////////////////////////////////////////
		
		// all truth signal muons
		muMapTruAll.insert(make_pair(truth_pt->at(imu1),imu1));
		muMapTruAll.insert(make_pair(truth_pt->at(imu2),imu2));
		muMapTruAll.insert(make_pair(truth_pt->at(imu3),imu3));

		// in acceptance signal muons
		bool ismu1acc = isInAcceptance(truth_pt->at(imu1),truth_eta->at(imu1));
		bool ismu2acc = isInAcceptance(truth_pt->at(imu2),truth_eta->at(imu2));
		bool ismu3acc = isInAcceptance(truth_pt->at(imu3),truth_eta->at(imu3));
		if(ismu1acc) muMapTruAcc.insert(make_pair(truth_pt->at(imu1),imu1));
		if(ismu2acc) muMapTruAcc.insert(make_pair(truth_pt->at(imu2),imu2));
		if(ismu3acc) muMapTruAcc.insert(make_pair(truth_pt->at(imu3),imu3));
		
		// match to reco muons
		int imumatch1 = match_tru2recmu(imu1);
		int imumatch2 = match_tru2recmu(imu2);
		int imumatch3 = match_tru2recmu(imu3);
		if(imumatch1>=0 && ismu1acc) muMapTruMat.insert(make_pair(truth_pt->at(imu1),imu1));
		if(imumatch2>=0 && ismu2acc) muMapTruMat.insert(make_pair(truth_pt->at(imu2),imu2));
		if(imumatch3>=0 && ismu3acc) muMapTruMat.insert(make_pair(truth_pt->at(imu3),imu3));
		
		// match to reco calo muons
		int icalmatch1 = match_tru2reccal(imu1);
		int icalmatch2 = match_tru2reccal(imu2);
		int icalmatch3 = match_tru2reccal(imu3);
		if(icalmatch1>=0 && ismu1acc) calMapTruMat.insert(make_pair(truth_pt->at(imu1),imu1));
		if(icalmatch2>=0 && ismu2acc) calMapTruMat.insert(make_pair(truth_pt->at(imu2),imu2));
		if(icalmatch3>=0 && ismu3acc) calMapTruMat.insert(make_pair(truth_pt->at(imu3),imu3));
		
		nrecmatmuons = ((imumatch1>=0)  + (imumatch2>=0)  + (imumatch3>=0));
		nrecmatcalos = ((icalmatch1>=0) + (icalmatch2>=0) + (icalmatch3>=0));

		// match to reco muons || calo muon
		if(ismu1acc && (imumatch1>=0 || icalmatch1>=0)) mucalMapTruMat.insert(make_pair(truth_pt->at(imu1),imu1));
		if(ismu2acc && (imumatch2>=0 || icalmatch2>=0)) mucalMapTruMat.insert(make_pair(truth_pt->at(imu2),imu2));
		if(ismu3acc && (imumatch3>=0 || icalmatch3>=0)) mucalMapTruMat.insert(make_pair(truth_pt->at(imu3),imu3));
		
		/////////////////////////////////////////////////////
		// allow only one signal W->tau->3mu per event... ///
		break; //////////////////////////////////////////////
		/////////////////////////////////////////////////////
	}
	
	/////////////////////////////////////////////
	// reco muons
	for(int mu=0 ; mu<(int)mu_pt->size() ; mu++)
	{
		// matched reco-to-signal muons
		for(TMapdi::iterator it=muMapTruAll.begin() ; it!=muMapTruAll.end() ; ++it)
		{
			int itru = it->second;
			
			if(!isInAcceptance(truth_pt->at(itru),truth_eta->at(itru))) continue;
			
			double dr = dRtru(mu_eta->at(mu),mu_phi->at(mu),truth_eta->at(itru),truth_phi->at(itru));
			double dptrel = dpTreltru(mu_pt->at(mu),truth_pt->at(itru));
			if(isMatched(dr,dptrel))
			{
				muMapRecMat.insert(make_pair(mu_pt->at(mu),mu));
				break;
			}
		}
	}
	
	/////////////////////////////////////////////
	// reco calons
	for(int cal=0 ; cal<(int)cal_pt->size() ; cal++)
	{
		// matched reco-to-signal muons
		for(TMapdi::iterator it=muMapTruAll.begin() ; it!=muMapTruAll.end() ; ++it)
		{
			int itru = it->second;
			
			if(!isInAcceptance(truth_pt->at(itru),truth_eta->at(itru))) continue;
			
			double dr = dRtru(cal_eta->at(cal),cal_phi->at(cal),truth_eta->at(itru),truth_phi->at(itru));
			double dptrel = dpTreltru(cal_pt->at(cal),truth_pt->at(itru));
			if(isMatched(dr,dptrel))
			{
				calMapRecMat.insert(make_pair(cal_pt->at(cal),cal));
				break;
			}
		}
	}
	
	n_trumu_all    = muMapTruAll.size();
	n_trumu_acc    = muMapTruAcc.size();
	n_recmu_mat    = muMapRecMat.size();	
	n_trumu_mat    = muMapTruMat.size();
	n_reccal_mat   = calMapRecMat.size();	
	n_trucal_mat   = calMapTruMat.size();
	n_trumucal_mat = mucalMapTruMat.size();
	
	// sort truth muons by pt
	TMapdi::reverse_iterator rit_tru_all=muMapTruAll.rbegin();
	if(n_trumu_all!=3) _FATAL("Didn't find 3 truth signal muons.");
	sorted_tru_mu_all[0] = rit_tru_all->second; rit_tru_all++;
	sorted_tru_mu_all[1] = rit_tru_all->second; rit_tru_all++;
	sorted_tru_mu_all[2] = rit_tru_all->second;

	// sort tru muons in acceptance
	for(TMapdi::reverse_iterator rit=muMapTruAcc.rbegin() ; rit!=muMapTruAcc.rend() ; ++rit)
	{
		int itru = rit->second;
		int order = getOrder(muMapTruAll,itru);
		if(order>=0) sorted_tru_mu_acc[order] = itru;
	}

	// sort tru muons reco-muon-matched
	for(TMapdi::reverse_iterator rit=muMapTruMat.rbegin() ; rit!=muMapTruMat.rend() ; ++rit)
	{
		int itru = rit->second;
		int order = getOrder(muMapTruAll,itru);
		if(order>=0) sorted_tru_mu_mat[order] = itru;
	}

	// sort tru muons reco-calomuon-matched
	for(TMapdi::reverse_iterator rit=calMapTruMat.rbegin() ; rit!=calMapTruMat.rend() ; ++rit)
	{
		int itru = rit->second;
		int order = getOrder(muMapTruAll,itru);
		if(order>=0) sorted_tru_cal_mat[order] = itru;
	}

	// sort tru reco-matched-muons || reco-matched-calomu 
	for(TMapdi::reverse_iterator rit=mucalMapTruMat.rbegin() ; rit!=mucalMapTruMat.rend() ; ++rit)
	{
		int itru = rit->second;
		int order = getOrder(muMapTruAll,itru);
		if(order>=0) sorted_tru_mucal_mat[order] = itru;
	}

	// sort rec matched muons by pt
	TMapdi::reverse_iterator rit_recmu_mat=muMapRecMat.rbegin();
	if(n_recmu_mat>0) sorted_rec_mu_mat[0] = rit_recmu_mat->second; rit_recmu_mat++;
	if(n_recmu_mat>1) sorted_rec_mu_mat[1] = rit_recmu_mat->second; rit_recmu_mat++;
	if(n_recmu_mat>2) sorted_rec_mu_mat[2] = rit_recmu_mat->second;
	
	// sort rec matched calo muons by pt
	TMapdi::reverse_iterator rit_reccal_mat=calMapRecMat.rbegin();
	if(n_reccal_mat>0) sorted_rec_cal_mat[0] = rit_reccal_mat->second; rit_reccal_mat++;
	if(n_reccal_mat>1) sorted_rec_cal_mat[1] = rit_reccal_mat->second; rit_reccal_mat++;
	if(n_reccal_mat>2) sorted_rec_cal_mat[2] = rit_reccal_mat->second;
	
	
	// counters
	// if(n_trumu_all==3)      nSignals++;
	// if(n_trumu_acc==3)      nAcceptance++;
	// if(n_recmu_mat==3)      nMatchedMuons++;
	// if(n_reccal_mat==3)     nMatchedCaloMuons++;
	// if(n_trumucal_mat==3 && n_trumu_mat==2) nMatchedMuonsCalos++;
}

bool isTruthMatchedTriplet(TString triplet_type, TLorentzVector& p)
{
	vector<int> index;
	vector<double>* pt;
	vector<double>* eta;
	vector<double>* phi;
	
	if     (triplet_type=="tru_all")           { index = sorted_tru_mu_all;     pt = truth_pt; eta = truth_eta; phi = truth_phi; }
	else if(triplet_type=="tru_in-acc")        { index = sorted_tru_mu_acc;     pt = truth_pt; eta = truth_eta; phi = truth_phi; }
	else if(triplet_type=="tru_muon-mat")      { index = sorted_tru_mu_mat;     pt = truth_pt; eta = truth_eta; phi = truth_phi; }
	else if(triplet_type=="tru_calo-mat")      { index = sorted_tru_cal_mat;    pt = truth_pt; eta = truth_eta; phi = truth_phi; }
	else if(triplet_type=="tru_muon/calo-mat") { index = sorted_tru_mucal_mat;  pt = truth_pt; eta = truth_eta; phi = truth_phi; }
	else if(triplet_type=="rec_muon-mat")      { index = sorted_rec_mu_mat;     pt = mu_pt;    eta = mu_eta;    phi = mu_phi;    }
	else if(triplet_type=="rec_calo-mat")      { index = sorted_rec_cal_mat;    pt = cal_pt;   eta = cal_eta;   phi = cal_phi;   }
	else _FATAL("unknown triplet type: "+(string)triplet_type);
	
	if(triplet_type.Contains("tru_"))
	{
		bool mc1 = (index[0]>=0);
		bool mc2 = (index[1]>=0);
		bool mc3 = (index[2]>=0);
		if((mc1+mc2+mc3)==3)
		{
			p = getTruthMatchedTriplet(index,pt,eta,phi);
			return true;
		}
	}
	return false;
}



}