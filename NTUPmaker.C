//////////////////////////////////////////////
//// run ./execute here and follow ///////////
//// the instructions ////////////////////////
//////////////////////////////////////////////

#include "../PileupReweighting/PileupReweighting/TPileupReweighting.h"
#include "../ApplyJetCalibration/ApplyJetCalibration/ApplyJetCalibration.h"
#include "../JetUncertainties/JetUncertainties/MultijetJESUncertaintyProvider.h"
#include "../ApplyJetResolutionSmearing/ApplyJetResolutionSmearing/ApplyJetSmearing.h"
#include "../JetUncertainties/JetUncertainties/JESUncertaintyProvider.h"
#include "../MissingETUtility/MissingETUtility/METUtility.h"
#include "../TileTripReader/TileTripReader/TTileTripReader.h"
#include "../BCHCleaningTool/BCHCleaningTool/BCHCleaningToolRoot.h"

#include "rawStd.h"
#include "rawROOT.h"
#include "types.h"
#include "enums.h"
#include "logs.h"
// #include "stle.h"
#include "histos.h"
#include "constants.h"
// #include "AtlasStyle.C"
#include "flatout.h"
#include "vertex.h"
#include "signalmc.h"
#include "tmva.h"

bool makepufile = false;
bool glob_isMC  = true;
bool glob_isSig  = true;

struct sources
{
	unsigned int   vtx;
	TString        type;
	TString        shorttype;
	TString        srcName[3];
	int            srcType[3];
	int            srcCode[3]; // Muons=MUONS / Muid=MUID, TPmuA=, TPmuB=, Calo=
	bool           isMuon[3];
	bool           isCalo[3];
	bool           isTPmu[3];
	bool           isTPa[3];
	bool           isTPb[3];
	int            trkIndex[3];
	int            srcIndex[3];
	int            srcOrder[3];
	TLorentzVector srcTlv[3];
	TLorentzVector p3body,p2body12,p2body13,p2body23;
	double q3body,q2body12,q2body13,q2body23;
};


void rgbPalette(Double_t r, Double_t g, Double_t b, Int_t nb=50)
{
	const UInt_t Number = 3;
	Double_t Red[Number]    = { r,   0.0, 0.0 };
	Double_t Green[Number]  = { g,   0.0, 0.0 };
	Double_t Blue[Number]   = { b,   0.0, 0.0 };
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
TExec* exeRed   = new TExec("exeRed",   "rgbPalette(0.7,0,0,50);");
TExec* exeGreen = new TExec("exeGreen", "rgbPalette(0,0.5,0,50);");
TExec* exeBlue  = new TExec("exeBlue",  "rgbPalette(0,0.2,0.9,50);");
TExec* exeGray  = new TExec("exeGray",  "rgbPalette(1,1,1,50);");

ofstream* ofstr = new ofstream("events.txt");
ofstream* ofstr1 = new ofstream("eventrecord.txt");


TString pdffilename;
TCanvas* c1;
double BRbelle  = 2.1e-8;
double BRfactor = 5e2;
TString sBRf    = "5e2";
TMapTSui counters;
TMapTSd  counters_weighted;
TMapTSui counters_evtvisited;
TMapiTS  counters_ordered;
TMapTSui triggerbits;
TMapuiTS triggerorder;
TMapTSd  iso;
TMapTSI  nevents;
TMapTSI  drawchannels;
TMapTSd  periodlumi;
TMapTSb  periodenable;
TMapTSb  binnedmcenable;
TMapTSd  weights;
double totalLumi;
TMapuivi isFourthTrack;
TString JetQuality = "VeryLooseBad";
vector<TLorentzVector> calibJets;
vector<bool>           badBCHJets;
vector<TLorentzVector> smearedJets;
vector<double>         calibJetsUnc;
vector<unsigned int>   calibJetsIndex;

vector<float>* AntiKt4LCTopoLets_calibrated_pt;
vector<float>* AntiKt4LCTopoLets_calibrated_E;
vector<float>* AntiKt4LCTopoLets_calibrated_m;
vector<float>* AntiKt4LCTopoLets_calibrated_eta;
vector<float>* AntiKt4LCTopoLets_calibrated_phi;

double isolation;
vector<double> iso10;
vector<double> iso20;
vector<double> iso30;
vector<double> iso40;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int VtxType;
float Pt1; float Eta1; float Phi1; int Code1; float ptFraction1;
float Pt2; float Eta2; float Phi2; int Code2; float ptFraction2; 
float Pt3; float Eta3; float Phi3; int Code3; float ptFraction3; 
float dPt12; double ptFraction12;
float dPt23; double ptFraction23;
float dPt13; double ptFraction13;
float SctAngSig1; float SctNgbSig1; float PbalSig1; 
float SctAngSig2; float SctNgbSig2; float PbalSig2; 
float SctAngSig3; float SctNgbSig3; float PbalSig3;
float Chi2TrkFit1; float NdfTrkFit1; float PvalTrkFit1; float QoverP1;
float Chi2TrkFit2; float NdfTrkFit2; float PvalTrkFit2; float QoverP2;
float Chi2TrkFit3; float NdfTrkFit3; float PvalTrkFit3; float QoverP3;
int NhtTRThits1; int NhtTRThits2; int NhtTRThits3;
float JetPt1; float JetEta1; float JetPhi1; float JetM1; float JetE1; float JetMV1w1;
float JetPt2; float JetEta2; float JetPhi2; float JetM2; float JetE2; float JetMV1w2;
float JetSumPt12; float JetdPhiJ1J2;
float JetdPhi3bodyJ1;
int nDijets; int nDijetsLoose; int nJet3muOverlaps; int nJet3muOverlapsLoose;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


string timestr;
string ftxtname;
double minpvalOS;
double mindchi2;
//TString basepath = "root://eosatlas//eos/atlas/user/h/hod/data/bphys/"; // all periods
TString basepath = "root://eosatlas//eos/atlas/user/h/hod/data/bphys_d3pd/"; // all periods
// TString basepath = "root://eosatlas//eos/atlas/user/h/hod/data/tests/overlaps/floatingOverlaps/bout/data/bphys/";
// TString basepath = "root://eosatlas//eos/atlas/user/h/hod/data/tests/overlaps/fixedOverlaps/bout/data/bphys/";
// TString basepath = "root://eosatlas//eos/atlas/user/h/hod/data/tests/overlaps/noOverlaps/bout/data/bphys/";
// TString basepath = "/tmp/hod/data/bphys/";
// TString basepath = "/stage1/hod/data/bphys/";
TString splitstr;
Int_t nChunks;
Int_t iChunk; 
Int_t lChunk; 
Int_t nChunksMax;
Int_t iFirstFile;
Int_t iLastFile;
Bool_t chainWasPrinted = false;
bool doBlind = true;
double mBlindMinGlob;
double mBlindMaxGlob;
TDirectory* olddir = gDirectory;
int idefault = -9999;
int ddefault = -9999.;
bool skim = false;
TString mastername;
TString foutname;
TString selectionMethod; // "CUTS:vertex,MET"; // "MVA:BDTG";
TString smethod;
bool doMVA;
TMapiTS tpmus;

// canvas
TCanvas* cnv = NULL;
vector<TVirtualPad*> pads;
int divx=-1;
int divy=-1;

Int_t  RunNumber,EventNumber,lbn;
vector<double>*          vtx_chi2;
vector<int>*             vtx_ndf;
vector<int>*             vtx_index;
vector<int>*             vtx_charge;
vector<double>*          vtx_mass;
vector<double>*          vtx_pt;
vector<double>*          vtx_rapidity;
vector<double>*          vtx_x;
vector<double>*          vtx_y;
vector<double>*          vtx_z;
vector<double>*          vtx_xErr;
vector<double>*          vtx_yErr;
vector<double>*          vtx_zErr;
vector<vector<double> >* vtx_reftrks_px;
vector<vector<double> >* vtx_reftrks_py;
vector<vector<double> >* vtx_reftrks_pz;
vector<vector<double> >* vtx_lxy;
vector<vector<double> >* vtx_lxyErr;
vector<vector<double> >* vtx_tau;
vector<vector<double> >* vtx_tauErr;
vector<vector<double> >* vtx_a0;
vector<vector<double> >* vtx_a0Err;
vector<vector<double> >* vtx_a0XY;
vector<vector<double> >* vtx_a0XYErr;
vector<vector<double> >* vtx_cosTheta;
vector<vector<double> >* vtx_cosThetaXY;
vector<vector<int> >*    vtx_trks_index;
vector<vector<int> >*    vtx_muons_index;
vector<vector<double> >* vtx_refPVx;
vector<vector<double> >* vtx_refPVy;
vector<vector<double> >* vtx_refPVz;
vector<vector<double> >* vtx_refPVcovxx;
vector<vector<double> >* vtx_refPVcovyy;
vector<vector<double> >* vtx_refPVcovzz;
vector<vector<int> >*    vtx_srcIndex;
vector<vector<int> >*    vtx_trkIndex;
vector<vector<string> >* vtx_srcName;

vector<double>* pv_x;
vector<double>* pv_y;
vector<double>* pv_z;
vector<double>* pv_xErr;
vector<double>* pv_yErr;
vector<double>* pv_zErr;
vector<int>*    pv_mu;
vector<double>* pv_chi2;
vector<double>* pv_ndf;
vector<int>*    pv_ntrk;
vector<int>*    pv_type;
vector<int>*    pv_index;

vector<string>* L1_trigger_name;
vector<string>* EF_trigger_name;

int                   mc_channel_number;
int                   mc_event_number;
double                mc_event_weight;
vector<int>*          mc_pdgId;
vector<int>*          mc_status;
vector<int>*          mc_barcode;
vector<int>*          mc_productionvtx_barcode;
vector<int>*          mc_decayvtx_barcode;
vector<int>*          mc_has_productionvtx;
vector<int>*          mc_has_decayvtx;
vector<double>*       mc_pt;
vector<double>*       mc_eta;
vector<double>*       mc_phi;
vector<double>*       mc_m;
vector<double>*       mc_charge;
vector<double>*       mc_productionvtx_x;
vector<double>*       mc_productionvtx_y;
vector<double>*       mc_productionvtx_z;
vector<double>*       mc_decayvtx_x;
vector<double>*       mc_decayvtx_y;
vector<double>*       mc_decayvtx_z;
vector<vector<int> >* mc_parents;
vector<vector<int> >* mc_children;

vector<double>* trks_chi2;
vector<double>* trks_qoverpErr;
vector<double>* trks_qoverp;
vector<double>* trks_pt;
vector<double>* trks_eta;
vector<double>* trks_px;
vector<double>* trks_py;
vector<double>* trks_pz;
vector<int>*    trks_ndf;
vector<int>*    trks_nBLayer;
vector<int>*    trks_nPix;
vector<int>*    trks_nSCT;
vector<int>*    trks_nTRT;
vector<int>*    trks_nTRTOutliers;
vector<int>*    trks_nHighThresholdTRTHits;
vector<int>*    trks_nPixHoles;
vector<int>*    trks_nSCTHoles;
vector<int>*    trks_nDeadPixels;
vector<int>*    trks_nDeadSCT;
vector<int>*    trks_expectBLayer;
vector<double>* trks_pixeldEdx;
vector<int>*    trks_nUsedHitsdEdx;

TMapTSP2vd tpmu_vd;
TMapTSP2vi tpmu_vi;
TMapTSP2vvi tpmu_vvi;
TMapTSP2vvs tpmu_vvs;

vector<int>* muons_charge;
vector<double>* muons_ptcone10;
vector<double>* muons_ptcone20;
vector<double>* muons_ptcone30;
vector<double>* muons_ptcone40;
vector<double>* muons_etcone10;
vector<double>* muons_etcone20;
vector<double>* muons_etcone30;
vector<double>* muons_etcone40;
vector<double>* muons_eLoss;
vector<double>* muons_eLossErr;
vector<double>* muons_px;
vector<double>* muons_py;
vector<double>* muons_pz;
vector<double>* muons_e;
vector<double>* muons_pt;
vector<double>* muons_eta;
vector<double>* muons_phi;
vector<float>*  muons_sctangsig;
vector<float>*  muons_sctngbsig;
vector<float>*  muons_pbalsig;
vector<bool>*   muons_isCombined;
vector<bool>*   muons_hasInDetTrackParticle;
vector<int>*    muons_inDetTrackIndex;
vector<int>*    muons_index;
vector<double>* muons_chi2;
vector<int>*    muons_ndf;
vector<int>*    muons_author;
vector<double>* muons_matchchi2;
vector<double>* muons_matchchi2ndf;
vector<double>* muons_phi_me;
vector<double>* muons_phi_ie;
vector<double>* muons_eta_me;
vector<double>* muons_eta_ie;
vector<double>* muons_pt_me;
vector<double>* muons_pt_ie;
vector<double>* muons_px_me;
vector<double>* muons_px_ie;
vector<double>* muons_py_me;
vector<double>* muons_py_ie;
vector<double>* muons_pz_me;
vector<double>* muons_pz_ie;
vector<double>* muons_e_me;
vector<double>* muons_e_ie;
vector<bool>*   muons_isLoose;
vector<bool>*   muons_isMedium;
vector<bool>*   muons_isTight;
vector<int>* muons_numberOfMdtHits;
vector<int>* muons_numberOfTgcPhiHits;
vector<int>* muons_numberOfTgcEtaHits;
vector<int>* muons_numberOfCscPhiHits;
vector<int>* muons_numberOfCscEtaHits;
vector<int>* muons_numberOfRpcPhiHits;
vector<int>* muons_numberOfRpcEtaHits;
vector<int>* muons_numberOfCscEtaHoles;
vector<int>* muons_numberOfCscPhiHoles;
vector<int>* muons_numberOfRpcEtaHoles;
vector<int>* muons_numberOfRpcPhiHoles;
vector<int>* muons_numberOfMdtHoles;
vector<int>* muons_numberOfTgcEtaHoles;
vector<int>* muons_numberOfTgcPhiHoles;
vector<int>* muons_numberOfOutliersOnTrack;
vector<int>* muons_standardDeviationOfChi2OS;
vector<vector<int> >* muons_sectors;
vector<vector<int> >* muons_nprecisionHits;
vector<vector<int> >* muons_nphiLayers;
vector<vector<int> >* muons_netaPhiLayers;
vector<vector<int> >* muons_nprecisionHoles;
vector<vector<int> >* muons_netaTriggerHoleLayers;
vector<vector<int> >* muons_nphiHoleLayers;
vector<vector<int> >* muons_nprecisionOutliers;
vector<vector<string> >* muons_stationName;
vector<vector<string> >* muons_phiName;

vector<int>* muid_charge;
vector<double>* muid_ptcone10;
vector<double>* muid_ptcone20;
vector<double>* muid_ptcone30;
vector<double>* muid_ptcone40;
vector<double>* muid_etcone10;
vector<double>* muid_etcone20;
vector<double>* muid_etcone30;
vector<double>* muid_etcone40;
vector<double>* muid_eLoss;
vector<double>* muid_eLossErr;
vector<double>* muid_px;
vector<double>* muid_py;
vector<double>* muid_pz;
vector<double>* muid_e;
vector<double>* muid_pt;
vector<double>* muid_eta;
vector<double>* muid_phi;
vector<float>*  muid_sctangsig;
vector<float>*  muid_sctngbsig;
vector<float>*  muid_pbalsig;
vector<bool>*   muid_isCombined;
vector<bool>*   muid_hasInDetTrackParticle;
vector<int>*    muid_inDetTrackIndex;
vector<int>*    muid_index;
vector<double>* muid_chi2;
vector<int>*    muid_ndf;
vector<int>*    muid_author;
vector<double>* muid_matchchi2;
vector<double>* muid_matchchi2ndf;
vector<double>* muid_phi_me;
vector<double>* muid_phi_ie;
vector<double>* muid_eta_me;
vector<double>* muid_eta_ie;
vector<double>* muid_pt_me;
vector<double>* muid_pt_ie;
vector<double>* muid_px_me;
vector<double>* muid_px_ie;
vector<double>* muid_py_me;
vector<double>* muid_py_ie;
vector<double>* muid_pz_me;
vector<double>* muid_pz_ie;
vector<double>* muid_e_me;
vector<double>* muid_e_ie;
vector<bool>*   muid_isLoose;
vector<bool>*   muid_isMedium;
vector<bool>*   muid_isTight;
vector<int>* muid_numberOfMdtHits;
vector<int>* muid_numberOfTgcPhiHits;
vector<int>* muid_numberOfTgcEtaHits;
vector<int>* muid_numberOfCscPhiHits;
vector<int>* muid_numberOfCscEtaHits;
vector<int>* muid_numberOfRpcPhiHits;
vector<int>* muid_numberOfRpcEtaHits;
vector<int>* muid_numberOfCscEtaHoles;
vector<int>* muid_numberOfCscPhiHoles;
vector<int>* muid_numberOfRpcEtaHoles;
vector<int>* muid_numberOfRpcPhiHoles;
vector<int>* muid_numberOfMdtHoles;
vector<int>* muid_numberOfTgcEtaHoles;
vector<int>* muid_numberOfTgcPhiHoles;
vector<int>* muid_numberOfOutliersOnTrack;
vector<int>* muid_standardDeviationOfChi2OS;
vector<vector<int> >* muid_sectors;
vector<vector<int> >* muid_nprecisionHits;
vector<vector<int> >* muid_nphiLayers;
vector<vector<int> >* muid_netaPhiLayers;
vector<vector<int> >* muid_nprecisionHoles;
vector<vector<int> >* muid_netaTriggerHoleLayers;
vector<vector<int> >* muid_nphiHoleLayers;
vector<vector<int> >* muid_nprecisionOutliers;
vector<vector<string> >* muid_stationName;
vector<vector<string> >* muid_phiName;


vector<int>*    calo_muons_charge;
vector<double>* calo_muons_ptcone10;
vector<double>* calo_muons_ptcone20;
vector<double>* calo_muons_ptcone30;
vector<double>* calo_muons_ptcone40;
vector<double>* calo_muons_etcone10;
vector<double>* calo_muons_etcone20;
vector<double>* calo_muons_etcone30;
vector<double>* calo_muons_etcone40;
vector<double>* calo_muons_eLoss;
vector<double>* calo_muons_eLossErr;
vector<double>* calo_muons_px;
vector<double>* calo_muons_py;
vector<double>* calo_muons_pz;
vector<double>* calo_muons_e;
vector<double>* calo_muons_pt;
vector<double>* calo_muons_eta;
vector<double>* calo_muons_phi;
vector<float>*  calo_muons_sctangsig;
vector<float>*  calo_muons_sctngbsig;
vector<float>*  calo_muons_pbalsig;
vector<bool>*   calo_muons_isCombined;
vector<bool>*   calo_muons_hasInDetTrackParticle;
vector<int>*    calo_muons_inDetTrackIndex;
vector<int>*    calo_muons_index;
vector<double>* calo_muons_chi2;
vector<int>*    calo_muons_ndf;
vector<int>*    calo_muons_author;
vector<double>* calo_muons_matchchi2;
vector<double>* calo_muons_matchchi2ndf;
vector<double>* calo_muons_phi_me;
vector<double>* calo_muons_phi_ie;
vector<double>* calo_muons_eta_me;
vector<double>* calo_muons_eta_ie;
vector<double>* calo_muons_pt_me;
vector<double>* calo_muons_pt_ie;
vector<double>* calo_muons_px_me;
vector<double>* calo_muons_px_ie;
vector<double>* calo_muons_py_me;
vector<double>* calo_muons_py_ie;
vector<double>* calo_muons_pz_me;
vector<double>* calo_muons_pz_ie;
vector<double>* calo_muons_e_me;
vector<double>* calo_muons_e_ie;
vector<bool>*   calo_muons_isLoose;
vector<bool>*   calo_muons_isMedium;
vector<bool>*   calo_muons_isTight;


int aux_nTriplets;
int aux_nDoublets;
int aux_nQuadruplets;
double aux_averageIntPerXing;
double aux_actualIntPerXing;


vector<int>*    jets_nTrks;
vector<int>*    jets_LArQuality;
vector<int>*    jets_isGood;
vector<int>*    jets_isBad;
vector<int>*    jets_isUgly;
vector<int>*    jets_eMfrac;
vector<double>* jets_px;
vector<double>* jets_py;
vector<double>* jets_pz;
vector<double>* jets_pt;
vector<double>* jets_et;
vector<double>* jets_e;
vector<double>* jets_m;
vector<double>* jets_eta;
vector<double>* jets_phi;
vector<double>* jets_flwgt_MV1;

int met_RefFinal_source;
double met_RefFinal_etX;
double met_RefFinal_etY;
double met_RefFinal_sumet;
double met_RefFinal_et;
double met_RefFinal_phi;

vector<vector<int> >*    muons_cal_subCalo_subCaloId;
vector<vector<double> >* muons_cal_subCalo_energyDeposited;
vector<vector<double> >* muons_cal_subCalo_muonEnergyLoss;
vector<int>*             muons_cal_energyLossType;
vector<int>*             muons_cal_caloMuonIdTag;
vector<double>*          muons_cal_caloLRLikelihood;
vector<double>*          muons_cal_fsrCandidateEnergy;
vector<double>*          muons_cal_etCore;


int    phys_RunNumber;
int    phys_EventNumber;
int    phys_mc_channel_number;
int    phys_mc_event_number;
float  phys_mc_event_weight;
int    phys_lbn;
float  phys_actualIntPerXing;
float  phys_averageIntPerXing;

int            AntiKt4LCTopoLets_n;
vector<float>* AntiKt4LCTopoLets_E;
vector<float>* AntiKt4LCTopoLets_pt;
vector<float>* AntiKt4LCTopoLets_m;
vector<float>* AntiKt4LCTopoLets_eta;
vector<float>* AntiKt4LCTopoLets_phi;
vector<float>* AntiKt4LCTopoLets_EtaOrigin;
vector<float>* AntiKt4LCTopoLets_PhiOrigin;
vector<float>* AntiKt4LCTopoLets_MOrigin;
vector<float>* AntiKt4LCTopoLets_WIDTH;
vector<float>* AntiKt4LCTopoLets_n90;
vector<float>* AntiKt4LCTopoLets_Timing;
vector<float>* AntiKt4LCTopoLets_LArQuality;
vector<float>* AntiKt4LCTopoLets_nTrk;
vector<float>* AntiKt4LCTopoLets_sumPtTrk;
vector<float>* AntiKt4LCTopoLets_OriginIndex;
vector<float>* AntiKt4LCTopoLets_HECQuality;
vector<float>* AntiKt4LCTopoLets_NegativeE;
vector<float>* AntiKt4LCTopoLets_AverageLArQF;
vector<float>* AntiKt4LCTopoLets_BCH_CORR_CELL;
vector<float>* AntiKt4LCTopoLets_BCH_CORR_DOTX;
vector<float>* AntiKt4LCTopoLets_BCH_CORR_JET;
vector<float>* AntiKt4LCTopoLets_BCH_CORR_JET_FORCELL;
vector<float>* AntiKt4LCTopoLets_ENG_BAD_CELLS;
vector<float>* AntiKt4LCTopoLets_N_BAD_CELLS;
vector<float>* AntiKt4LCTopoLets_N_BAD_CELLS_CORR;
vector<float>* AntiKt4LCTopoLets_BAD_CELLS_CORR_E;
vector<float>* AntiKt4LCTopoLets_NumTowers;
vector<float>* AntiKt4LCTopoLets_ootFracCells5;
vector<float>* AntiKt4LCTopoLets_ootFracCells10;
vector<float>* AntiKt4LCTopoLets_ootFracClusters5;
vector<float>* AntiKt4LCTopoLets_ootFracClusters10;
vector<int>*   AntiKt4LCTopoLets_SamplingMax;
vector<float>* AntiKt4LCTopoLets_fracSamplingMax;
vector<float>* AntiKt4LCTopoLets_hecf;
vector<float>* AntiKt4LCTopoLets_tgap3f;
vector<int>*   AntiKt4LCTopoLets_isUgly;
vector<int>*   AntiKt4LCTopoLets_isBadLooseMinus;
vector<int>*   AntiKt4LCTopoLets_isBadLoose;
vector<int>*   AntiKt4LCTopoLets_isBadMedium;
vector<int>*   AntiKt4LCTopoLets_isBadTight;
vector<float>* AntiKt4LCTopoLets_emfrac;
vector<float>* AntiKt4LCTopoLets_Offset;
vector<float>* AntiKt4LCTopoLets_EMJES;
vector<float>* AntiKt4LCTopoLets_EMJES_EtaCorr;
vector<float>* AntiKt4LCTopoLets_EMJESnooffset;
vector<float>* AntiKt4LCTopoLets_LCJES;
vector<float>* AntiKt4LCTopoLets_LCJES_EtaCorr;
vector<float>* AntiKt4LCTopoLets_emscale_E;
vector<float>* AntiKt4LCTopoLets_emscale_pt;
vector<float>* AntiKt4LCTopoLets_emscale_m;
vector<float>* AntiKt4LCTopoLets_emscale_eta;
vector<float>* AntiKt4LCTopoLets_emscale_phi;
vector<float>* AntiKt4LCTopoLets_jvtx_x;
vector<float>* AntiKt4LCTopoLets_jvtx_y;
vector<float>* AntiKt4LCTopoLets_jvtx_z;
vector<float>* AntiKt4LCTopoLets_jvtxf;
vector<float>* AntiKt4LCTopoLets_GSCFactorF;
vector<float>* AntiKt4LCTopoLets_WidthFraction;
vector<float>* AntiKt4LCTopoLets_e_PreSamplerB;
vector<float>* AntiKt4LCTopoLets_e_EMB1;
vector<float>* AntiKt4LCTopoLets_e_EMB2;
vector<float>* AntiKt4LCTopoLets_e_EMB3;
vector<float>* AntiKt4LCTopoLets_e_PreSamplerE;
vector<float>* AntiKt4LCTopoLets_e_EME1;
vector<float>* AntiKt4LCTopoLets_e_EME2;
vector<float>* AntiKt4LCTopoLets_e_EME3;
vector<float>* AntiKt4LCTopoLets_e_HEC0;
vector<float>* AntiKt4LCTopoLets_e_HEC1;
vector<float>* AntiKt4LCTopoLets_e_HEC2;
vector<float>* AntiKt4LCTopoLets_e_HEC3;
vector<float>* AntiKt4LCTopoLets_e_TileBar0;
vector<float>* AntiKt4LCTopoLets_e_TileBar1;
vector<float>* AntiKt4LCTopoLets_e_TileBar2;
vector<float>* AntiKt4LCTopoLets_e_TileGap1;
vector<float>* AntiKt4LCTopoLets_e_TileGap2;
vector<float>* AntiKt4LCTopoLets_e_TileGap3;
vector<float>* AntiKt4LCTopoLets_e_TileExt0;
vector<float>* AntiKt4LCTopoLets_e_TileExt1;
vector<float>* AntiKt4LCTopoLets_e_TileExt2;
vector<float>* AntiKt4LCTopoLets_e_FCAL0;
vector<float>* AntiKt4LCTopoLets_e_FCAL1;
vector<float>* AntiKt4LCTopoLets_e_FCAL2;
vector<vector<float> >* AntiKt4LCTopoLets_shapeBins;
vector<int>*  AntiKt4LCTopoLets_Nconst;
vector<vector<float> >* AntiKt4LCTopoLets_ptconst_default;
vector<vector<float> >* AntiKt4LCTopoLets_econst_default;
vector<vector<float> >* AntiKt4LCTopoLets_etaconst_default;
vector<vector<float> >* AntiKt4LCTopoLets_phiconst_default;
vector<vector<float> >* AntiKt4LCTopoLets_weightconst_default;
vector<float>* AntiKt4LCTopoLets_constscale_E;
vector<float>* AntiKt4LCTopoLets_constscale_pt;
vector<float>* AntiKt4LCTopoLets_constscale_m;
vector<float>* AntiKt4LCTopoLets_constscale_eta;
vector<float>* AntiKt4LCTopoLets_constscale_phi;
vector<float>* AntiKt4LCTopoLets_LArBadHVEnergy;
vector<float>* AntiKt4LCTopoLets_LArBadHVRatio;
vector<float>* AntiKt4LCTopoLets_flavor_weight_Comb;
vector<float>* AntiKt4LCTopoLets_flavor_weight_IP2D;
vector<float>* AntiKt4LCTopoLets_flavor_weight_IP3D;
vector<float>* AntiKt4LCTopoLets_flavor_weight_SV0;
vector<float>* AntiKt4LCTopoLets_flavor_weight_SV1;
vector<float>* AntiKt4LCTopoLets_flavor_weight_SV2;
vector<float>* AntiKt4LCTopoLets_flavor_weight_SoftMuonTagChi2;
vector<float>* AntiKt4LCTopoLets_flavor_weight_SecondSoftMuonTagChi2;
vector<float>* AntiKt4LCTopoLets_flavor_weight_JetFitterTagNN;
vector<float>* AntiKt4LCTopoLets_flavor_weight_JetFitterCOMBNN;
vector<float>* AntiKt4LCTopoLets_flavor_weight_MV1;
vector<float>* AntiKt4LCTopoLets_flavor_weight_MV2;
vector<float>* AntiKt4LCTopoLets_flavor_weight_GbbNN;
vector<int>*   AntiKt4LCTopoLets_flavor_truth_label;
vector<float>* AntiKt4LCTopoLets_flavor_truth_dRminToB;
vector<float>* AntiKt4LCTopoLets_flavor_truth_dRminToC;
vector<float>* AntiKt4LCTopoLets_flavor_truth_dRminToT;
vector<int>*   AntiKt4LCTopoLets_flavor_truth_BHadronpdg;
vector<float>* AntiKt4LCTopoLets_flavor_truth_vx_x;
vector<float>* AntiKt4LCTopoLets_flavor_truth_vx_y;
vector<float>* AntiKt4LCTopoLets_flavor_truth_vx_z;
vector<int>*   AntiKt4LCTopoLets_flavor_putruth_label;
vector<float>* AntiKt4LCTopoLets_flavor_putruth_dRminToB;
vector<float>* AntiKt4LCTopoLets_flavor_putruth_dRminToC;
vector<float>* AntiKt4LCTopoLets_flavor_putruth_dRminToT;
vector<int>*   AntiKt4LCTopoLets_flavor_putruth_BHadronpdg;
vector<float>* AntiKt4LCTopoLets_flavor_putruth_vx_x;
vector<float>* AntiKt4LCTopoLets_flavor_putruth_vx_y;
vector<float>* AntiKt4LCTopoLets_flavor_putruth_vx_z;
vector<float>* AntiKt4LCTopoLets_flavor_component_ip2d_pu;
vector<float>* AntiKt4LCTopoLets_flavor_component_ip2d_pb;
vector<int>*   AntiKt4LCTopoLets_flavor_component_ip2d_isValid;
vector<int>*   AntiKt4LCTopoLets_flavor_component_ip2d_ntrk;
vector<float>* AntiKt4LCTopoLets_flavor_component_ip3d_pu;
vector<float>* AntiKt4LCTopoLets_flavor_component_ip3d_pb;
vector<int>  * AntiKt4LCTopoLets_flavor_component_ip3d_isValid;
vector<int>  * AntiKt4LCTopoLets_flavor_component_ip3d_ntrk;
vector<int>  * AntiKt4LCTopoLets_flavor_component_jetprob_ntrk;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv1_pu;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv1_pb;
vector<int>*   AntiKt4LCTopoLets_flavor_component_sv1_isValid;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv2_pu;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv2_pb;
vector<int>*   AntiKt4LCTopoLets_flavor_component_sv2_isValid;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_pu;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_pb;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_pc;
vector<int>*   AntiKt4LCTopoLets_flavor_component_jfit_isValid;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfitcomb_pu;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfitcomb_pb;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfitcomb_pc;
vector<int>*   AntiKt4LCTopoLets_flavor_component_jfitcomb_isValid;
vector<int>*   AntiKt4LCTopoLets_flavor_component_gbbnn_nMatchingTracks;
vector<double>* AntiKt4LCTopoLets_flavor_component_gbbnn_trkJetWidth;
vector<double>* AntiKt4LCTopoLets_flavor_component_gbbnn_trkJetMaxDeltaR;
vector<int>*   AntiKt4LCTopoLets_flavor_component_jfit_nvtx;
vector<int>*   AntiKt4LCTopoLets_flavor_component_jfit_nvtx1t;
vector<int>*   AntiKt4LCTopoLets_flavor_component_jfit_ntrkAtVx;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_efrc;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_mass;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_sig3d;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_deltaPhi;
vector<float>* AntiKt4LCTopoLets_flavor_component_jfit_deltaEta;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_d0val;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_d0sig;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_z0val;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_z0sig;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_w2D;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_w3D;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_pJP;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_pJPneg;
vector<vector<int> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_grade;
vector<vector<int> >* AntiKt4LCTopoLets_flavor_component_ipplus_trk_isFromV0;
vector<int>* AntiKt4LCTopoLets_flavor_component_svp_isValid;
vector<int>* AntiKt4LCTopoLets_flavor_component_svp_ntrkv;
vector<int>* AntiKt4LCTopoLets_flavor_component_svp_ntrkj;
vector<int>* AntiKt4LCTopoLets_flavor_component_svp_n2t;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_mass;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_efrc;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_x;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_y;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_z;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_err_x;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_err_y;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_err_z;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_cov_xy;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_cov_xz;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_cov_yz;
vector<float>* AntiKt4LCTopoLets_flavor_component_svp_chi2;
vector<int>*  AntiKt4LCTopoLets_flavor_component_svp_ndof;
vector<int>*  AntiKt4LCTopoLets_flavor_component_svp_ntrk;
vector<int>*  AntiKt4LCTopoLets_flavor_component_sv0p_isValid;
vector<int>*  AntiKt4LCTopoLets_flavor_component_sv0p_ntrkv;
vector<int>*  AntiKt4LCTopoLets_flavor_component_sv0p_ntrkj;
vector<int>*  AntiKt4LCTopoLets_flavor_component_sv0p_n2t;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_mass;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_efrc;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_x;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_y;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_z;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_err_x;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_err_y;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_err_z;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_cov_xy;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_cov_xz;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_cov_yz;
vector<float>* AntiKt4LCTopoLets_flavor_component_sv0p_chi2;
vector<int>*   AntiKt4LCTopoLets_flavor_component_sv0p_ndof;
vector<int>*   AntiKt4LCTopoLets_flavor_component_sv0p_ntrk;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_softmuoninfo_muon_w;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_softmuoninfo_muon_pTRel;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_softmuoninfo_muon_dRJet;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_softmuonchi2info_muon_w;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_softmuonchi2info_muon_pTRel;
vector<vector<float> >* AntiKt4LCTopoLets_flavor_component_softmuonchi2info_muon_dRJet;
vector<float>* AntiKt4LCTopoLets_el_dr;
vector<int>*   AntiKt4LCTopoLets_el_matched;
vector<float>* AntiKt4LCTopoLets_mu_dr;
vector<int>*   AntiKt4LCTopoLets_mu_matched;
vector<float>* AntiKt4LCTopoLets_L1_dr;
vector<int>*   AntiKt4LCTopoLets_L1_matched;
vector<float>* AntiKt4LCTopoLets_L2_dr;
vector<int>*   AntiKt4LCTopoLets_L2_matched;
vector<float>* AntiKt4LCTopoLets_EF_dr;
vector<int>*   AntiKt4LCTopoLets_EF_matched;
vector<float>* AntiKt4LCTopoLets_LikeLihood_0;
vector<float>* AntiKt4LCTopoLets_ActiveArea;
vector<float>* AntiKt4LCTopoLets_ActiveAreaPx;
vector<float>* AntiKt4LCTopoLets_ActiveAreaPy;
vector<float>* AntiKt4LCTopoLets_ActiveAreaPz;
vector<float>* AntiKt4LCTopoLets_ActiveAreaE;
vector<float>* AntiKt4LCTopoLets_VoronoiArea;
vector<float>* AntiKt4LCTopoLets_VoronoiAreaPx;
vector<float>* AntiKt4LCTopoLets_VoronoiAreaPy;
vector<float>* AntiKt4LCTopoLets_VoronoiAreaPz;
vector<float>* AntiKt4LCTopoLets_VoronoiAreaE;
vector<float>* AntiKt4LCTopoLets_LowEtConstituentsFrac;
vector<float>* AntiKt4LCTopoLets_pt_truth;
vector<float>* AntiKt4LCTopoLets_IsoKR20Perp;
vector<float>* AntiKt4LCTopoLets_IsoKR20Par;
vector<float>* AntiKt4LCTopoLets_IsoKR20SumPt;
vector<float>* AntiKt4LCTopoLets_IsoDelta2Perp;
vector<float>* AntiKt4LCTopoLets_IsoDelta2Par;
vector<float>* AntiKt4LCTopoLets_IsoDelta2SumPt;
vector<float>* AntiKt4LCTopoLets_IsoFixedCone8Perp;
vector<float>* AntiKt4LCTopoLets_IsoFixedCone8Par;
vector<float>* AntiKt4LCTopoLets_IsoFixedCone8SumPt;
vector<float>* AntiKt4LCTopoLets_IsoFixedArea13Perp;
vector<float>* AntiKt4LCTopoLets_IsoFixedArea13Par;
vector<float>* AntiKt4LCTopoLets_IsoFixedArea13SumPt;
vector<float>* AntiKt4LCTopoLets_Iso6To88Perp;
vector<float>* AntiKt4LCTopoLets_Iso6To88Par;
vector<float>* AntiKt4LCTopoLets_Iso6To88SumPt;
vector<float>* AntiKt4LCTopoLets_KtDr;
vector<float>* AntiKt4LCTopoLets_Centroid_r;
vector<float>* AntiKt4LCTopoLets_nTrk_pv0_1GeV;
vector<float>* AntiKt4LCTopoLets_sumPtTrk_pv0_1GeV;
vector<float>* AntiKt4LCTopoLets_nTrk_allpv_1GeV;
vector<float>* AntiKt4LCTopoLets_sumPtTrk_allpv_1GeV;
vector<float>* AntiKt4LCTopoLets_nTrk_pv0_500MeV;
vector<float>* AntiKt4LCTopoLets_sumPtTrk_pv0_500MeV;
vector<float>* AntiKt4LCTopoLets_trackWIDTH_pv0_1GeV;
vector<float>* AntiKt4LCTopoLets_trackWIDTH_allpv_1GeV;
vector<float>* AntiKt4LCTopoLets_TrackMFindex;
vector<float>* AntiKt4LCTopoLets_TrackMFPt;
vector<float>* AntiKt4LCTopoLets_TruthMFindex;
vector<float>* AntiKt4LCTopoLets_TruthMFPt;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
float MET_RefFinal_etx;
float MET_RefFinal_ety;
float MET_RefFinal_phi;
float MET_RefFinal_et;
float MET_RefFinal_sumet;

float MET_RefMuon_etx;
float MET_RefMuon_ety;
float MET_RefMuon_phi;
float MET_RefMuon_et;
float MET_RefMuon_sumet;

float MET_RefTau_etx;
float MET_RefTau_ety;
float MET_RefTau_phi;
float MET_RefTau_et;
float MET_RefTau_sumet;

float MET_CellOut_Eflow_etx;
float MET_CellOut_Eflow_ety;
float MET_CellOut_Eflow_phi;
float MET_CellOut_Eflow_et;
float MET_CellOut_Eflow_sumet;

float MET_Track_etx;
float MET_Track_ety;
float MET_Track_phi;
float MET_Track_et;
float MET_Track_sumet;

float MET_Staco_etx;
float MET_Staco_ety;
float MET_Staco_phi;
float MET_Staco_et;
float MET_Staco_sumet;

float MET_Muid_etx;
float MET_Muid_ety;
float MET_Muid_phi;
float MET_Muid_et;
float MET_Muid_sumet;

float MET_RefGamma_etx;
float MET_RefGamma_ety;
float MET_RefGamma_phi;
float MET_RefGamma_et;
float MET_RefGamma_sumet;

float MET_RefEle_etx;
float MET_RefEle_ety;
float MET_RefEle_phi;
float MET_RefEle_et;
float MET_RefEle_sumet;

float MET_RefJet_etx;
float MET_RefJet_ety;
float MET_RefJet_phi;
float MET_RefJet_et;
float MET_RefJet_sumet;

float MET_SoftJets_etx;
float MET_SoftJets_ety;
float MET_SoftJets_phi;
float MET_SoftJets_et;
float MET_SoftJets_sumet;

float MET_Truth_NonInt_etx;
float MET_Truth_NonInt_ety;
float MET_Truth_NonInt_phi;
float MET_Truth_NonInt_et;
float MET_Truth_NonInt_sumet;

Int_t          trk_n;
vector<float>* trk_pt;
vector<float>* trk_eta;
vector<float>* trk_phi_wrtPV;

Int_t          mu_muons_n;
vector<float>* mu_muons_E;
vector<float>* mu_muons_pt;
vector<float>* mu_muons_m;
vector<float>* mu_muons_eta;
vector<float>* mu_muons_phi;
vector<float>* mu_muons_charge;
vector<float>* mu_muons_ms_phi;
vector<float>* mu_muons_ms_theta;
vector<float>* mu_muons_ms_qoverp;

Int_t          mu_staco_n;
vector<float>* mu_staco_E;
vector<float>* mu_staco_pt;
vector<float>* mu_staco_m;
vector<float>* mu_staco_eta;
vector<float>* mu_staco_phi;
vector<float>* mu_staco_charge;
vector<float>* mu_staco_ms_phi;
vector<float>* mu_staco_ms_theta;
vector<float>* mu_staco_ms_qoverp;

Int_t          mu_muid_n;
vector<float>* mu_muid_E;
vector<float>* mu_muid_pt;
vector<float>* mu_muid_m;
vector<float>* mu_muid_eta;
vector<float>* mu_muid_phi;
vector<float>* mu_muid_charge;
vector<float>* mu_muid_ms_phi;
vector<float>* mu_muid_ms_theta;
vector<float>* mu_muid_ms_qoverp;

Int_t          el_n;
vector<float>* el_E;
vector<float>* el_Et;
vector<float>* el_pt;
vector<float>* el_m;
vector<float>* el_eta;
vector<float>* el_phi;

Int_t          ph_n;
vector<float>* ph_E;
vector<float>* ph_Et;
vector<float>* ph_pt;
vector<float>* ph_m;
vector<float>* ph_eta;
vector<float>* ph_phi;

Int_t          tau_n;
vector<float>* tau_Et;
vector<float>* tau_pt;
vector<float>* tau_m;
vector<float>* tau_eta;
vector<float>* tau_phi;

int  el_MET_n;
vector<vector<float> > *el_MET_wpx;
vector<vector<float> > *el_MET_wpy;
vector<vector<float> > *el_MET_wet;
vector<vector<unsigned int> > *el_MET_statusWord;
int           ph_MET_n;
vector<vector<float> > *ph_MET_wpx;
vector<vector<float> > *ph_MET_wpy;
vector<vector<float> > *ph_MET_wet;
vector<vector<unsigned int> > *ph_MET_statusWord;
int           mu_staco_MET_n;
vector<vector<float> > *mu_staco_MET_wpx;
vector<vector<float> > *mu_staco_MET_wpy;
vector<vector<float> > *mu_staco_MET_wet;
vector<vector<unsigned int> > *mu_staco_MET_statusWord;
int           mu_muid_MET_n;
vector<vector<float> > *mu_muid_MET_wpx;
vector<vector<float> > *mu_muid_MET_wpy;
vector<vector<float> > *mu_muid_MET_wet;
vector<vector<unsigned int> > *mu_muid_MET_statusWord;
int           mu_MET_n;
vector<vector<float> > *mu_MET_wpx;
vector<vector<float> > *mu_MET_wpy;
vector<vector<float> > *mu_MET_wet;
vector<vector<unsigned int> > *mu_MET_statusWord;
int           tau_MET_n;
vector<vector<float> > *tau_MET_wpx;
vector<vector<float> > *tau_MET_wpy;
vector<vector<float> > *tau_MET_wet;
vector<vector<unsigned int> > *tau_MET_statusWord;
int           jet_antikt4LCtopo_MET_n;
vector<vector<float> > *jet_antikt4LCtopo_MET_wpx;
vector<vector<float> > *jet_antikt4LCtopo_MET_wpy;
vector<vector<float> > *jet_antikt4LCtopo_MET_wet;
vector<vector<unsigned int> > *jet_antikt4LCtopo_MET_statusWord;
int           cl_MET_n;
vector<vector<float> > *cl_MET_wpx;
vector<vector<float> > *cl_MET_wpy;
vector<vector<float> > *cl_MET_wet;
vector<vector<unsigned int> > *cl_MET_statusWord;
int           trk_MET_n;
vector<vector<float> > *trk_MET_wpx;
vector<vector<float> > *trk_MET_wpy;
vector<vector<float> > *trk_MET_wet;
vector<vector<unsigned int> > *trk_MET_statusWord;


float rhorhoKt3EM;
float rhorhoKt4EM;
float rhorhoKt3LC;
float rhorhoKt4LC;

Int_t                   musp_n;
vector<float>*          musp_eta;
vector<float>*          musp_phi;
vector<unsigned short>* musp_trigHits;
vector<unsigned short>* musp_innerHits;
vector<unsigned short>* musp_middleHits;
vector<unsigned short>* musp_outerHits;
vector<unsigned short>* musp_innerSegments;
vector<unsigned short>* musp_middleSegments;
vector<unsigned short>* musp_outerSegments;

Int_t          vxp_n;
vector<float>* vxp_x;
vector<float>* vxp_y;
vector<float>* vxp_z;
vector<int>*   vxp_type;
vector<float>* vxp_chi2;
vector<int>*   vxp_ndof;
vector<float>* vxp_px;
vector<float>* vxp_py;
vector<float>* vxp_pz;
vector<float>* vxp_E;
vector<float>* vxp_m;
vector<int>*   vxp_nTracks;
vector<float>* vxp_sumPt;



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

Root::TPileupReweighting* pileupTool;
void initializePileup()
{
	pileupTool = new Root::TPileupReweighting("pileuptool");
	if(makepufile)
	{
		_INFO("making pileup file");
		pileupTool->UsePeriodConfig("MC12b");
	}
	else
	{
		_INFO("reading pileup file");
		pileupTool->AddConfigFile("Wtaunu_3mu.prw.root");
		pileupTool->AddLumiCalcFile("ilumicalc_histograms_None_200842-215643.root");
	}
	pileupTool->Initialize();
}
float getPileupWeight()
{
	/////////////////////////////
	float pileup_weight = 1.; ///
	/////////////////////////////

	// NOTE (23/01/2013): A bug has been found in the d3pd making code,
	// causing all MC12 samples to have a few of the averageIntPerXing
	// values incorrectly set (some should be 0 but are set to 1).
	// The bug does not affect data. To resolve this, when reading this branch,
	// for both prw file generating and for when retrieving pileup weights,
	// you should amend the value with the following line of code:

	float averageIntPerXing_fixed = (lbn==1 && int(phys_averageIntPerXing+0.5)==1) ? 0. : phys_averageIntPerXing;
	if(makepufile) pileupTool->Fill(phys_RunNumber,phys_mc_channel_number,phys_mc_event_weight,averageIntPerXing_fixed);
	else           pileup_weight = pileupTool->GetCombinedWeight(phys_RunNumber,phys_mc_channel_number,averageIntPerXing_fixed);
	return pileup_weight;
}
void finalizePileup()
{
	pileupTool->WriteToFile("Wtaunu_3mu.prw.root");
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

TLorentzVector getCalibratedJet(unsigned int jet)
{
	if(calibJets.size()==0) _FATAL("calibJets size is 0");
	return calibJets[jet];
}

BCHTool::BCHCleaningToolRoot* BCH;
Root::TTileTripReader* TTR;
void initBCH(bool isData)
{
	TTR = new Root::TTileTripReader("TileTripReader");	
	TTR->setTripFile("/afs/cern.ch/user/h/hod/TileTripReader/data/CompleteTripList_2011-2012.root" );
	
	// initialize the tool, giving the path to the file to emulate masked modules (used in MC only)
	BCH = new BCHTool::BCHCleaningToolRoot("BCH");
	BCH->InitializeTool(isData, TTR, "/afs/cern.ch/user/h/hod/BCHCleaningTool/share/FractionsRejectedJetsMC.root");
}

bool isBadCalibJetBCH(unsigned int jet)
{
	int randomRunNumber = pileupTool->GetRandomRunNumber(phys_RunNumber,phys_averageIntPerXing);
	int randomLumiBlockNumber = 0;
	if(randomRunNumber>0)
	{
		randomLumiBlockNumber = pileupTool->GetRandomLumiBlockNumber(randomRunNumber);
		double calibrated_calo_jet_pt  = AntiKt4LCTopoLets_calibrated_pt->at(jet);
		double calibrated_calo_jet_eta = AntiKt4LCTopoLets_calibrated_eta->at(jet);
		double calibrated_calo_jet_phi = AntiKt4LCTopoLets_calibrated_phi->at(jet);
		double BCH_CORR_CELL           = AntiKt4LCTopoLets_BCH_CORR_CELL->at(jet);
		double emfrac                  = AntiKt4LCTopoLets_emfrac->at(jet);
		int isUp                       = 0; // isUp (MC only, default=0): variation due to flavor dependence, can be either -1 (sample = 100% light quark jets), 0 (default QCD composition), 1 (sample = 100% gluon jets)
		
		 bool isBad = false;
		
		// bool isInMaskedRegion = BCH->IsInMaskedRegion(randomRunNumber, randomLumiBlockNumber, calibrated_calo_jet_eta, calibrated_calo_jet_phi);
		
		// OR: The medium cut, which comes with associated systematics in MC
		BCH->SetSeed(phys_RunNumber); // For setting the random seed used within IsBadMediumBCH(). Make sure to use this when comparing to others to make sure the same events are selected
		isBad = BCH->IsBadMediumBCH(randomRunNumber, randomLumiBlockNumber,
									calibrated_calo_jet_eta, calibrated_calo_jet_phi,
									BCH_CORR_CELL, emfrac, calibrated_calo_jet_pt, isUp);
		// OR: The tight cut
		// BCH->IsBadTightBCH(const int run, const int lbn, const double eta, const double phi, 
		//                    const double BCH_CORR_CELL, const double emfrac, const double pt);
		if(isBad) return true;
	}
	return false;
}

bool isGoodJet(unsigned int jet, TString configuration = "VeryLooseBad")
{
	bool isGood = true;
	
	// first do BCH cleaning
	if(badBCHJets.size()>0 && badBCHJets[jet]) return false;
	
	// then the rest of usual jet cleaning
	double emf          = AntiKt4LCTopoLets_emfrac->at(jet);
	double hecf         = AntiKt4LCTopoLets_hecf->at(jet);
	double larq         = AntiKt4LCTopoLets_LArQuality->at(jet);
	double hecq         = AntiKt4LCTopoLets_HECQuality->at(jet);
	double time         = AntiKt4LCTopoLets_Timing->at(jet); //in ns
	double sumpttrk     = AntiKt4LCTopoLets_sumPtTrk_pv0_500MeV->at(jet); //in MeV, same as sumpttrk
	double eta          = AntiKt4LCTopoLets_emscale_eta->at(jet); //emscale Eta // AntiKt4LCTopoLets_eta->at(jet);  ?
	double pt           = (skim) ? AntiKt4LCTopoLets_pt->at(jet) : AntiKt4LCTopoLets_calibrated_pt->at(jet)/*getCalibratedJet(jet).Pt()*/; // should be calibrated jet pT in the analysis mode (!skim)
	double fmax         = AntiKt4LCTopoLets_fracSamplingMax->at(jet);
	double negE         = AntiKt4LCTopoLets_NegativeE->at(jet); //in MeV
	double AverageLArQF = AntiKt4LCTopoLets_AverageLArQF->at(jet);
	
	// -----------------------------------------------------------
	// Do the actual selection
	double chf=sumpttrk/pt;
	
	//=============================================================
	//VeryLoose cuts
	//=============================================================
	//Non-collision background & cosmics
	if(emf<0.05 && chf<0.05 && fabs(eta)<2)                                   isGood = false;
	if(emf<0.05 && fabs(eta)>=2)                                              isGood = false;
	if(fmax>0.99 && fabs(eta)<2)                                              isGood = false;
	//HEC spike
	if(fabs(negE/1000.)>60)                                                   isGood = false;
	if(hecf>0.5 && fabs(hecq)>0.5 && AverageLArQF/65535>0.8)                  isGood = false;
	//EM calo noise
	if(emf>0.95 && fabs(larq)>0.8 && fabs(eta)<2.8 && AverageLArQF/65535>0.8) isGood = false;
	if ("VeryLooseBad"==configuration) return isGood;
	
	//=============================================================
	//Loose cuts
	//=============================================================
	//Non-collision background & cosmics
	if(fabs(time)>25)                                  isGood = false;
	//HEC spike
	if(hecf>0.5 && fabs(hecq)>0.5)                     isGood = false;
	//EM calo noise
	if(emf>0.95 && fabs(larq)>0.8 && fabs(eta)<2.8)    isGood = false;
	if ("LooseBad"==configuration) return isGood;

	//=============================================================
	//Additionnal medium cuts
	//=============================================================
	//Non-collision background & cosmics
	if(fabs(time)>10)                                isGood = false;
	if(emf<0.05 && chf<0.1  && fabs(eta)<2)          isGood = false;
	if(emf>0.95 && chf<0.05 && fabs(eta)<2)          isGood = false;
	//HEC spike
	if(hecf>1-fabs(hecq))                            isGood = false;
	//EM calo noise
	if(emf>0.9 && fabs(larq)>0.8 && fabs(eta)<2.8)   isGood = false;
	if ("MediumBad"==configuration) return isGood;
	
	//=============================================================
	//Additionnal tight cuts
	//=============================================================
	//Non-collision background & cosmics
	if(emf<0.1 && chf<0.2 && fabs(eta)<2.5)          isGood = false;
	if(emf<0.1 && fabs(eta)>2.5 )                    isGood = false;
	if(emf>0.9 && chf<0.1 && fabs(eta)<2.5)          isGood = false;
	//EM calo noise
	if(fabs(larq)>0.95)                              isGood = false;
	if(emf>0.98 && fabs(larq)>0.05)                  isGood = false;
	if(chf<0.01 && fabs(eta)<2.5 )                   isGood = false;
	if ("TightBad"==configuration) return isGood;
	
	// We should never arrive here!
	_FATAL("Unknown configuration: "+(string)configuration+" in jet cleaning");

	return isGood;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

JetAnalysisCalib::JetCalibrationTool* JES;
void initJES(bool isdata, TString tag="Moriond", TString jetAlgo="AntiKt4TopoEM") 
{
	if(JES) delete JES;
	
	if(AntiKt4LCTopoLets_calibrated_pt)  delete AntiKt4LCTopoLets_calibrated_pt;
	if(AntiKt4LCTopoLets_calibrated_E)   delete AntiKt4LCTopoLets_calibrated_E;
	if(AntiKt4LCTopoLets_calibrated_m)   delete AntiKt4LCTopoLets_calibrated_m;
	if(AntiKt4LCTopoLets_calibrated_eta) delete AntiKt4LCTopoLets_calibrated_eta;
	if(AntiKt4LCTopoLets_calibrated_phi) delete AntiKt4LCTopoLets_calibrated_phi;
	AntiKt4LCTopoLets_calibrated_pt  = new vector<float>;
    AntiKt4LCTopoLets_calibrated_E   = new vector<float>;
    AntiKt4LCTopoLets_calibrated_m   = new vector<float>;
    AntiKt4LCTopoLets_calibrated_eta = new vector<float>;
    AntiKt4LCTopoLets_calibrated_phi = new vector<float>;

	TString JES_config_file = "../ApplyJetCalibration/data/CalibrationConfigs/";
	
	if(tag=="Moriond") JES_config_file += "JES_Full2012dataset_Preliminary_Jan13.config"; // Data/MC12 Configuration file
	else
	{
		if(isdata) JES_config_file += "JES_Full2012dataset_May2014.config";       // Data Configuration file
		else       JES_config_file += "JES_Full2012dataset_MC12b_May2014.config"; // MC12b Configuration file (Fullsim)
	}
	JES = new JetAnalysisCalib::JetCalibrationTool(jetAlgo,JES_config_file,isdata);
}

int findMuspContainer(float jet_eta, float jet_phi)
{
	double delR = 0;
	double delta_phi;
	double delta_eta;
	// Match mspn container to jet delR<0.4
	delR = 100;
	int index_musp=-1;
	for(unsigned int i=0 ; i<musp_phi->size(); i++)
	{
		delta_phi = fabs(jet_phi-musp_phi->at(i)); //calculate the distance in phi
		if(delta_phi>TMath::Pi()) delta_phi = (2*TMath::Pi()) - delta_phi; // always take the smaller angle (below 180°)
		delta_eta = jet_eta - musp_eta->at(i); // distance in eta
		if(sqrt(delta_phi*delta_phi + delta_eta*delta_eta) < delR)
		{
			delR = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
			index_musp = i;
		}
	}
	if(delR<0.4) return index_musp;
	return -1;
}
double getNsegments(unsigned int jet)
{
	double eta_det   = AntiKt4LCTopoLets_emscale_eta->at(jet);
	double phi       = AntiKt4LCTopoLets_emscale_phi->at(jet);
	double Nsegments = 0;  //see below
	int musp_index = findMuspContainer(eta_det, phi); //Detector jet eta and phi, obtained as above.
	if(musp_index>=0) Nsegments = musp_innerSegments->at(musp_index) + musp_outerSegments->at(musp_index) + musp_middleSegments->at(musp_index);
	return Nsegments;
}
int getNPV()
{
	int NPV=0; // count the number of vertices with 2 or more tracks
	// for(unsigned tracki=0 ; tracki<vxp_nTracks->size(); tracki++)
	for(unsigned tracki=0 ; tracki<pv_ntrk->size() ; tracki++)
	{
		// if(vxp_nTracks->at(tracki)>=2) NPV++;
		if(pv_ntrk->at(tracki)>=2) NPV++;
	}
	return NPV;
}
void clearCalibratedJets()
{
	if(!AntiKt4LCTopoLets_calibrated_pt) return;
	
	AntiKt4LCTopoLets_calibrated_pt->clear();
	AntiKt4LCTopoLets_calibrated_E->clear();
	AntiKt4LCTopoLets_calibrated_m->clear();
	AntiKt4LCTopoLets_calibrated_eta->clear();
	AntiKt4LCTopoLets_calibrated_phi->clear();
}
void addCalibratedJet(TLorentzVector Jet)
{
	if(!AntiKt4LCTopoLets_calibrated_pt) return;
	
	AntiKt4LCTopoLets_calibrated_pt->push_back(Jet.Pt());
	AntiKt4LCTopoLets_calibrated_E->push_back(Jet.E());
	AntiKt4LCTopoLets_calibrated_m->push_back(Jet.M());
	AntiKt4LCTopoLets_calibrated_eta->push_back(Jet.Eta());
	AntiKt4LCTopoLets_calibrated_phi->push_back(Jet.Phi());
}

TLorentzVector getJES(unsigned int jet, TString tag="Moriond")
{
	TLorentzVector Jet;
	
	if(tag=="Moriond")
	{
		double Eraw = AntiKt4LCTopoLets_emscale_E->at(jet);
		double eta  = AntiKt4LCTopoLets_emscale_eta->at(jet);
		double phi  = AntiKt4LCTopoLets_emscale_phi->at(jet);
		double m    = AntiKt4LCTopoLets_emscale_m->at(jet);
		double Ax   = AntiKt4LCTopoLets_ActiveAreaPx->at(jet);
		double Ay   = AntiKt4LCTopoLets_ActiveAreaPy->at(jet);
		double Az   = AntiKt4LCTopoLets_ActiveAreaPz->at(jet);
		double Ae   = AntiKt4LCTopoLets_ActiveAreaE->at(jet);
		double rho  = rhorhoKt4EM;
		
		// For the pile-up correction, we need mu and NPV(2+ tracks)
		double mu = phys_averageIntPerXing;
		int NPV = getNPV(); // count the number of vertices with 2 or more tracks
		
		Jet = JES->ApplyJetAreaOffsetEtaJES(Eraw,eta,phi,m,Ax,Ay,Az,Ae,rho,mu,NPV);
	}
	else
	{
		double Eraw       = AntiKt4LCTopoLets_emscale_E->at(jet);
		double eta_det    = AntiKt4LCTopoLets_emscale_eta->at(jet);
		double phi        = AntiKt4LCTopoLets_emscale_phi->at(jet);
		double m          = AntiKt4LCTopoLets_emscale_m->at(jet);
		double eta_origin = AntiKt4LCTopoLets_EtaOrigin->at(jet);
		double phi_origin = AntiKt4LCTopoLets_PhiOrigin->at(jet);
		double m_origin   = AntiKt4LCTopoLets_MOrigin->at(jet);
		double Ax         = AntiKt4LCTopoLets_ActiveAreaPx->at(jet);
		double Ay         = AntiKt4LCTopoLets_ActiveAreaPy->at(jet);
		double Az         = AntiKt4LCTopoLets_ActiveAreaPz->at(jet);
		double Ae         = AntiKt4LCTopoLets_ActiveAreaE->at(jet);
		double rho        = rhorhoKt4EM;       
		double fEM3       = (AntiKt4LCTopoLets_e_EMB3->at(jet)+AntiKt4LCTopoLets_e_EME3->at(jet))/AntiKt4LCTopoLets_emscale_E->at(jet);
		double fTile0     = (AntiKt4LCTopoLets_e_TileBar0->at(jet)+AntiKt4LCTopoLets_e_TileExt0->at(jet))/AntiKt4LCTopoLets_emscale_E->at(jet);
		double nTrk       = AntiKt4LCTopoLets_nTrk_pv0_1GeV->at(jet); 
		double trackWIDTH = AntiKt4LCTopoLets_trackWIDTH_pv0_1GeV->at(jet); //This variable may also be called AntiKt4LCTopoLets_trackWIDTH_pv0_1GeV
		double Nsegments = getNsegments(jet);
				
		// For the pile-up correction, we need mu and NPV(2+ tracks)
		double mu = phys_averageIntPerXing;
		int NPV = getNPV(); // count the number of vertices with 2 or more tracks
		
		Jet = JES->ApplyJetAreaOffsetOriginEtaJESGSC(Eraw,eta_det,phi,m,eta_origin,phi_origin,m_origin,Ax,Ay,Az,Ae,rho,trackWIDTH,nTrk,fTile0,fEM3,Nsegments,mu,NPV);
	}
	return Jet;
}

MultijetJESUncertaintyProvider* JUNCP;
void initJUN(TString tag="Moriond")
{
	if(JUNCP) delete JUNCP;
	
	if(tag=="Moriond")
	{
		TString multijetConfig = "JES_2012/Moriond2013/MultijetJES_2012.config";
		TString jesprovConfig  = "JES_2012/Moriond2013/InsituJES2012_14NP.config"; // Two other nominal options
		TString jetAlgorithm   = "AntiKt4LCTopo"; // One of {AntiKt4TopoEM,AntiKt4LCTopo,AntiKt6TopoEM,AntiKt6LCTopo}
		TString mcType         = "MC12a"; // Also accepts {Pythia8,MC12c,AFII} (note: MC12a==Pythia8)
		TString path           = "../JetUncertainties/share/"; // Path to the share directory
		JUNCP = new MultijetJESUncertaintyProvider(multijetConfig,jesprovConfig,jetAlgorithm,mcType,path);
	}
	else
	{
		TString base_path_for_config_files = "../JetUncertainties/share/";
		TString flavour_and_topology_uncertainties_config_file = "JES_2012/Final/MultijetJES_2012.config";
		TString baseline_uncertainties_uncertainties_config_file = "JES_2012/Final/InsituJES2012_AllNuisanceParameters.config";
		TString collection_name = "AntiKt4LCTopo"; //  AntiKt4TopoEM, AntiKt6TopoEM, AntiKt4TopoLC and AntiKt6TopoLC.
		TString MC_youre_running_on = "MC12a"; // "Pythia8";
		JUNCP = new MultijetJESUncertaintyProvider(flavour_and_topology_uncertainties_config_file,
											  	   baseline_uncertainties_uncertainties_config_file,
											  	   collection_name,
											  	   MC_youre_running_on,
											  	   base_path_for_config_files);
	}
}


double getJUNC(unsigned int jet, TLorentzVector Jet, TMapVL& JetShiftsUp, TMapVL& JetShiftsDwn, TString tag="Moriond")
{
	double total_shift = 0;
	JetShiftsUp.clear();
	JetShiftsDwn.clear();
	TLorentzVector JetShiftUp, JetShiftDwn;
	
	/////////////////////////
	// Baseline uncertainties
	int Ncomp = JUNCP->getNUncertaintyComponents();
	for(int icomp=0 ; icomp<Ncomp ; ++icomp)
	{
		TString compName = JUNCP->getComponentNames().at(icomp); // returns the name of the nuisance parameter
		TString compDesc = JUNCP->getComponentDescriptions().at(icomp); // returns a one-sentence description of the nuisance parameter
		compCategory categoryEnum = JUNCP->getComponentCategories().at(icomp); // returns a enumeration of the category for the nuisance parameter
		TString categoryName = JUNCP->getCategoryStringFromEnum(categoryEnum); // returns a string for the category name from the enumeration
		double shift = JUNCP->getRelUncertComponent(icomp, Jet.Pt(), Jet.Eta()); // nuisance paramter amplitude (with sign) ("relative uncertainty" of component) 
		// can also use  jUNCP->getRelUncertComponent(compName, Jet.Pt(), Jet.Eta());
		
		JetShiftUp  = Jet;
		JetShiftDwn = Jet;
		
		// now scale your jet (the full 4-vector)
		// Upward variation: Jet *= 1+shift
		// Downward variation: Jet *= 1-shift
		JetShiftUp *= (1.0+shift);
		JetShiftDwn *= (1.0-shift);
		
		total_shift += shift*shift;
		
		JetShiftsUp.insert(make_pair(compName,JetShiftUp));
		JetShiftsDwn.insert(make_pair(compName,JetShiftDwn));
		
		// print the uncertainty and information
		// printf("JES uncertainty source %2d: %s\n",icomp,compName.Data());  
		// printf("Description : %s\n",compDesc.Data());
		// printf("Category: %s\n",categoryName.Data());
		// printf("The uncertainty for (pT,eta) = (%.1f,%5.1f) is %.1f%%\n",Jet.Pt(),Jet.Eta(),unc*100);
	}
	
	///////////////////////
	// Pileup uncertainties
	double mu = phys_averageIntPerXing;
	int NPV = getNPV(); // count the number of vertices with 2 or more tracks

	double shift = 0;
	
	shift = JUNCP->getRelNPVOffsetTerm(Jet.Pt(),Jet.Eta(),NPV); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("NPVOffsetTerm", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("NPVOffsetTerm",JetShiftDwn));
	
	shift = JUNCP->getRelMuOffsetTerm(Jet.Pt(),Jet.Eta(),mu); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("MuOffsetTerm", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("MuOffsetTerm",JetShiftDwn));
	
	shift = JUNCP->getRelPileupPtTerm(Jet.Pt(),Jet.Eta(),NPV,mu); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("PileupPtTerm", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("PileupPtTerm",JetShiftDwn));
	
	shift = JUNCP->getRelPileupRhoTopology(Jet.Pt(),Jet.Eta()); total_shift += shift*shift;
	JetShiftUp  = Jet; JetShiftUp  *= (1.0+shift); JetShiftsUp.insert(make_pair("PileupRhoTopology", JetShiftUp));
	JetShiftDwn = Jet; JetShiftDwn *= (1.0-shift); JetShiftsDwn.insert(make_pair("PileupRhoTopology",JetShiftDwn));
	
	
	///////////////////////////////////////
	//// b-tagging / flavor uncertainty
	double Nsegments = getNsegments(jet);
	shift = JUNCP->getRelFlavorCompUncert(Jet.Pt(),Jet.Eta(),true);  total_shift += shift*shift; // ONLY FOR LIGHT JETS
	shift = JUNCP->getRelFlavorResponseUncert(Jet.Pt(),Jet.Eta());   total_shift += shift*shift; // ONLY FOR LIGHT JETS
	shift = JUNCP->getRelBJESUncert(Jet.Pt(),Jet.Eta());             total_shift += shift*shift; // ONLY FOR b-JETS
	
	
	if(tag=="Moriond")
	{
		// ????
	}
	else
	{
		shift = JUNCP->getRelPunchThroughUncert(Jet.Pt(),Jet.Eta(),Nsegments);  total_shift += shift*shift; // Punch throgh
	}
	
	return sqrt(total_shift);
}


JetSmearingTool* JER;
void initJER()
{
	TString JERinputFile = "/afs/cern.ch/user/h/hod/JetResolution/share/JERProviderPlots_2012.root";
	TString jetAlgo = "AntiKt4LCTopo";
	JER = new JetSmearingTool(jetAlgo,JERinputFile);
	JER->init();
}

void seedJER()
{
	JER->SetSeed(phys_RunNumber);
}

TLorentzVector getJER(TLorentzVector Jet)
{
	TLorentzVector smearedJet = Jet;
	// The below is systematic evaluation, and ONLY for MC
	// Smear the jet to match the MC resolution+1 sigma!
	JER->SmearJet_Syst(smearedJet); // FOR 7 TeV 2011 and 8 TeV 2012 MC12 Fullsim
	return smearedJet;
}


METUtility *METU;
void initMET()
{
	METU = new METUtility;
	// METU->setVerbosity(isVerbose);
	// METU->setSoftJetCut(20); // soft jet cut
	
	// // set pre-defined regions as above
	// m_util->setCaloRegion(region);
	// // region = METUtil::Central, METUtil::EndCap or METUtil::Forward
	// // or set eta cuts explicitly
	// m_util->setObjectEtaCut(etaLow, etaHigh);
}

METUtil::METObject getMETU()
{
	METU->reset();

	if(el_pt==NULL || el_MET_wet==NULL) _FATAL("");

	METU->setMETTerm(METUtil::RefMuon, MET_RefMuon_etx, MET_RefMuon_ety, MET_RefMuon_sumet); // must be done in addition to the staco code below
	METU->setMETTerm(METUtil::RefTau, MET_RefTau_etx, MET_RefTau_ety, MET_RefTau_sumet);
	METU->setMETTerm(METUtil::SoftTerms, MET_CellOut_Eflow_etx, MET_CellOut_Eflow_ety, MET_CellOut_Eflow_sumet);
	
	//METU->setElectronParameters(el_pt, el_eta, el_phi,
	//			    el_MET_wet, el_MET_wpx, el_MET_wpy, el_MET_statusWord);
	METU->setMETTerm(METUtil::RefEle, MET_RefEle_etx, MET_RefEle_ety, MET_RefEle_sumet);

	//METU->setPhotonParameters(ph_pt, ph_eta, ph_phi,
	//		 	  ph_MET_wet, ph_MET_wpx, ph_MET_wpy, ph_MET_statusWord);
	METU->setMETTerm(METUtil::RefGamma, MET_RefGamma_etx, MET_RefGamma_ety, MET_RefGamma_sumet);
	
	METU->setJetParameters(AntiKt4LCTopoLets_calibrated_pt,AntiKt4LCTopoLets_calibrated_eta,AntiKt4LCTopoLets_calibrated_phi,AntiKt4LCTopoLets_calibrated_E,
						   jet_antikt4LCtopo_MET_wet, jet_antikt4LCtopo_MET_wpx, jet_antikt4LCtopo_MET_wpy, jet_antikt4LCtopo_MET_statusWord);
	// METU->setMETTerm(METUtil::RefJet, MET_RefJet_etx, MET_RefJet_ety, MET_RefJet_sumet);

	METU->setMuonParameters(mu_staco_pt, mu_staco_eta, mu_staco_phi,
							mu_staco_MET_wet, mu_staco_MET_wpx, mu_staco_MET_wpy, mu_staco_MET_statusWord);
	// The default pt for muons is from the combined ID/MS track.
	// Spectro-only muons need to have the MS momentum set separately.
	// This version of the method uses commonly available D3PD branches.
	METU->setExtraMuonParameters(mu_staco_ms_qoverp, mu_staco_ms_theta, mu_staco_ms_phi, mu_staco_charge);

	
	//// Now get the MET (RefFinal)
	METUtil::METObject refFinal = METU->getMissingET(METUtil::RefFinal);
	
	// METUtil::METObject refFinal_test(MET_RefFinal_etx, MET_RefFinal_ety, MET_RefFinal_sumet);
	// bool check_refFinal = METU->checkConsistency(refFinal_test,METUtil::RefFinal);
	// if(check_refFinal) cout << "RefFinal checks out!" << endl;
	// else               cout << "RefFinal doesn't check out!" << endl;
	
 	return refFinal;
}

// void setMETUNC()
// {
// 	// Fill vectors here
// 	vector<float> jesUp, jesDown, jerUp, jerDown;
// 	unsigned int nJets = calibJetsUnc.size();
// 	for(unsigned int jet=0 ; jet<nJets ;  ++jet)
// 	{
// 		float jesUncUp(0.), jesUncDown(0.); // , jerShiftUp(0.), jerShiftDown(0.);
// 		
// 		// get these values from the jet uncertainties
// 		jesUncUp   = calibJetsUnc[jet];
// 		jesUncDown = calibJetsUnc[jet];
// 	
// 		jesUp.push_back(jesUncUp);
// 		jesDown.push_back(jesUncUp);
// 		// jerUp.push_back(jerShiftUp);
// 		// jerDown.push_back(jerShiftDown); 
// 	}
// 	METU->setObjectEnergyUncertainties(METUtil::Jets, jesUp, jesDown);
// 	// METU->setObjectResolutionShift(METUtil::Jets, jerUp, jerDown);
// 	
// 	
// 	// METU->setAverageIntPerXing(phys_averageIntPerXing);
// 	// METU->setMETTerm(METUtil::Truth, MET_Truth_NonInt_etx, MET_Truth_NonInt_ety, MET_Truth_NonInt_sumet);
// }
// 
// void setMETUNC()
// {
// 	METUtility::METObject refFinal_nominal = METU->getMissingET(METUtil::RefFinal);
// 	METUtility::METObject refFinal_JESUP   = METU->getMissingET(METUtil::RefFinal, METUtil::JESUp);
// 	METUtility::METObject refFinal_JESDOWN = METU->getMissingET(METUtil::RefFinal, METUtil::JESDown);
// }


void fillJetCalibrationHistos(unsigned int jet, TLorentzVector Jet, TLorentzVector smearedJet, TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double wgt)
{
	double pTreldiff  = Jet.Pt()/AntiKt4LCTopoLets_pt->at(jet) - 1.;
	double Ereldiff   = Jet.E()/AntiKt4LCTopoLets_E->at(jet) - 1.;
	double mreldiff   = Jet.M()/AntiKt4LCTopoLets_m->at(jet) - 1.;
	double etareldiff = Jet.Eta()/AntiKt4LCTopoLets_eta->at(jet) - 1.;
	double phireldiff = Jet.Phi()/AntiKt4LCTopoLets_phi->at(jet) - 1.;
	
	double pTJetUncalib  = AntiKt4LCTopoLets_pt->at(jet);
	double EJetUncalib   = AntiKt4LCTopoLets_E->at(jet);
	double mJetUncalib   = AntiKt4LCTopoLets_m->at(jet);
	double etaJetUncalib = AntiKt4LCTopoLets_eta->at(jet);
	double phiJetUncalib = AntiKt4LCTopoLets_phi->at(jet);
	
	// if(pTJetUncalib>15.*GeV2MeV && fabs(etaJetUncalib)<4.5 && isGoodJet(jet,JetQuality))
	// {
		histos[name+"_jet_calibration_pt_reldiff"]->Fill(pTreldiff,wgt);
		histos[name+"_jet_calibration_E_reldiff"]->Fill(Ereldiff,wgt);
		histos[name+"_jet_calibration_m_reldiff"]->Fill(mreldiff,wgt);
		histos[name+"_jet_calibration_eta_reldiff"]->Fill(etareldiff,wgt);
		histos[name+"_jet_calibration_phi_reldiff"]->Fill(phireldiff,wgt);
	
		histos2[name+"_jet_calibration_reldiff_vs_pt"]->Fill(pTJetUncalib,pTreldiff,wgt);
		histos2[name+"_jet_calibration_reldiff_vs_E"]->Fill(EJetUncalib,Ereldiff,wgt);
		histos2[name+"_jet_calibration_reldiff_vs_m"]->Fill(mJetUncalib,mreldiff,wgt);
		histos2[name+"_jet_calibration_reldiff_vs_eta"]->Fill(etaJetUncalib,etareldiff,wgt);
		histos2[name+"_jet_calibration_reldiff_vs_phi"]->Fill(phiJetUncalib,phireldiff,wgt);
	// }
	
	pTreldiff  = smearedJet.Pt()/Jet.Pt() - 1.;
	Ereldiff   = smearedJet.E()/Jet.E() - 1.;
	mreldiff   = smearedJet.M()/Jet.M() - 1.;
	etareldiff = smearedJet.Eta()/Jet.Eta() - 1.;
	phireldiff = smearedJet.Phi()/Jet.Phi() - 1.;
	
	double pTJetCalib  = Jet.Pt();
	double EJetCalib   = Jet.E();
	double mJetCalib   = Jet.M();
	double etaJetCalib = Jet.Eta();
	double phiJetCalib = Jet.Phi();
	
	// if(pTJetUncalib>15.*GeV2MeV && fabs(etaJetUncalib)<4.5 && isGoodJet(jet,JetQuality))
	// {
		histos[name+"_jet_resolution_pt_reldiff"]->Fill(pTreldiff,wgt);
		histos[name+"_jet_resolution_E_reldiff"]->Fill(Ereldiff,wgt);
		histos[name+"_jet_resolution_m_reldiff"]->Fill(mreldiff,wgt);
		histos[name+"_jet_resolution_eta_reldiff"]->Fill(etareldiff,wgt);
		histos[name+"_jet_resolution_phi_reldiff"]->Fill(phireldiff,wgt);
	
		histos2[name+"_jet_resolution_reldiff_vs_pt"]->Fill(pTJetCalib,pTreldiff,wgt);
		histos2[name+"_jet_resolution_reldiff_vs_E"]->Fill(EJetCalib,Ereldiff,wgt);
		histos2[name+"_jet_resolution_reldiff_vs_m"]->Fill(mJetCalib,mreldiff,wgt);
		histos2[name+"_jet_resolution_reldiff_vs_eta"]->Fill(etaJetCalib,etareldiff,wgt);
		histos2[name+"_jet_resolution_reldiff_vs_phi"]->Fill(phiJetCalib,phireldiff,wgt);
	// }
}

void fillMETCalibrationHistos(METUtil::METObject calibMET, TString name, TMapTSP2TH1& histos/*, TMapTSP2TH2& histos2*/, double wgt)
{
	double refFinal_et = calibMET.et();
	// double refFinal_etx = calibMET.etx();
	// double refFinal_ety = calibMET.ety();
	double refFinal_phi = calibMET.phi();
	//cout << " MET_RefFinal_et : " << MET_RefFinal_et << " -> " << refFinal_et << endl;
	//cout << " MET_RefFinal_phi: " << MET_RefFinal_phi << " -> " << refFinal_phi << endl;
	
	double relDiff_et  = (refFinal_et-MET_RefFinal_et)/MET_RefFinal_et; 
	double relDiff_phi = (refFinal_phi-MET_RefFinal_phi)/MET_RefFinal_phi; 

	histos[name+"_met_calibration_et_reldiff"]->Fill(relDiff_et,wgt);
	histos[name+"_met_calibration_phi_reldiff"]->Fill(relDiff_phi,wgt);
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

// void atlstyle()
// {
// 	gROOT->Reset();
// 	gROOT->ForceStyle();
// 	gROOT->LoadMacro("../../../TauLFVCommonTools/TauLFVCommonTools/AtlasStyle.C");
// 	SetAtlasStyle();
// }
int increment(int& counter)
{
	counter += 1;
	return counter;
}
int findInVector(vector<int>& v, int x)
{
	vector<int>::iterator it = find(v.begin(),v.end(),x);
	int index = (it!=v.end()) ? distance(v.begin(),it) : -1;
	return index;
}
int findInVector(vector<int>* v, int x)
{
	vector<int>::iterator it = find(v->begin(),v->end(),x);
	int index = (it!=v->end()) ? distance(v->begin(),it) : -1;
	return index;
}
TVector3 vtrk(double px, double py, double pz)
{
	TVector3 v;
	v.SetXYZ(px,py,pz);
	return v;
}
double phitrk(double px, double py, double pz)
{
	return vtrk(px,py,pz).Phi();
}
double qtrk(double qoverp)
{
	if(qoverp==0.) { _ERROR("a neutral track ?"); return 0.; }
	return (qoverp>0) ? +1. : -1.;
}
double mT2(double MET, double phiMET, double pTX, double phiX)
{
	return 2*pTX*MET*(1-TMath::Cos(phiX-phiMET));
}
double mT(double MET, double phiMET, double pTX, double phiX)
{
	return TMath::Sqrt(mT2(MET,phiMET,pTX,phiX));
}
float dPhi(float phi1, float phi2)
{
	return TVector2::Phi_mpi_pi(phi1-phi2);
}
float deltaR(float eta1, float phi1, float eta2, float phi2)
{
	float deta = eta1-eta2;
	float dphi = dPhi(phi1,phi2);
	return TMath::Sqrt(deta*deta + dphi*dphi);
}
void minmaxDeltaR(vector<float>* eta, vector<float>* phi, float& dRmin, float& dRmax)
{
	dRmin = +1.e20;
	dRmax = -1.e20;
	int n = eta->size();
	for(int i=0 ; i<n ; i++)
	{
		for(int j=i+1 ; j<n ; j++)
		{
			float dr = deltaR(eta->at(i),phi->at(i), eta->at(j),phi->at(j));
			dRmin = (dr<dRmin) ? dr : dRmin;
			dRmax = (dr>dRmax) ? dr : dRmax;
		}
	}
}
void minmaxDeltaR(float etaRef, float phiRef, vector<float>* eta, vector<float>* phi, float& dRmin, float& dRmax)
{
	dRmin = +1.e20;
	dRmax = -1.e20;
	int n = eta->size();
	for(int i=0 ; i<n ; i++)
	{
		float dr = deltaR(etaRef,phiRef,eta->at(i),phi->at(i));
		dRmin = (dr<dRmin) ? dr : dRmin;
		dRmax = (dr>dRmax) ? dr : dRmax;
	}
}
//////////////////////////////////////////////////////////////////////




TString getMasterTreeName()
{
	if     (mastername.Contains("muons")) return "MUONS_TRIPLET";
	else if(mastername.Contains("muid"))  return "MUID_TRIPLET";
	else _FATAL("The name has to contain triplet type - cannot be:"+(string)mastername);
	return "XXX";
}
TString getPrimaryTPTreeName()
{
	TString chainContains = "";
	if     (mastername.Contains("muons")) return tpmus[MUONSTPA];
	else if(mastername.Contains("muid"))  return tpmus[MUIDTPA];
	else _FATAL("The name has to contain triplet type");
	return "XXX";
}
TString getSecondaryTPTreeName()
{
	TString chainContains = "";
	if     (mastername.Contains("muons")) return tpmus[MUONSTPB]; // or 3
	else if(mastername.Contains("muid"))  return tpmus[MUIDTPB]; // or 3 
	else _FATAL("The name has to contain triplet type");
	return "XXX";
}
TString classifyTriplet(unsigned int vtx)
{
	// 1, "CombinedFitMuonParticles"
	// 2, "SegmentTagTrackParticles"
	// 3, "MuGirlRefittedTrackParticles"
	// 4, "MuidCombTrackParticles"
	// 5, "MuTagIMOTrackParticles"
	
	// TString tpmuA = getPrimaryTPTreeName();
	// TString tpmuB = getSecondaryTPTreeName();
	
	TString src1 = vtx_srcName->at(vtx)[0];
	TString src2 = vtx_srcName->at(vtx)[1];
	TString src3 = vtx_srcName->at(vtx)[2];
	
	TString type = src1+":"+src2+":"+src3;
	return type;
}
TString classifyTripletShort(unsigned int vtx)
{
	TString type = classifyTriplet(vtx);
	TString shorttype = "";
	
	if     (type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[MUONS])          shorttype = "3muons";
	
	else if(type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPA])       shorttype = "2muons1tpmuA";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONS])       shorttype = "2muons1tpmuA";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONS]+":"+tpmus[MUONS])       shorttype = "2muons1tpmuA";
	
	else if(type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPB])       shorttype = "2muons1tpmuB";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONS])       shorttype = "2muons1tpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONS]+":"+tpmus[MUONS])       shorttype = "2muons1tpmuB";
	
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPA])    shorttype = "1muons2tpmuA";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPA])    shorttype = "1muons2tpmuA";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONS])    shorttype = "1muons2tpmuA";
	
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPB])    shorttype = "1muons2tpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPB])    shorttype = "1muons2tpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONS])    shorttype = "1muons2tpmuB";
	
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPB])    shorttype = "1muons1tpmuA1tpmuB";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPA])    shorttype = "1muons1tpmuA1tpmuB";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPB])    shorttype = "1muons1tpmuA1tpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPA])    shorttype = "1muons1tpmuA1tpmuB";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONS])    shorttype = "1muons1tpmuA1tpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONS])    shorttype = "1muons1tpmuA1tpmuB";
	
	else if(type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[CALO])           shorttype = "2muons1calomu";
	else if(type==tpmus[MUONS]+":"+tpmus[CALO]+":"+tpmus[MUONS])           shorttype = "2muons1calomu";
	else if(type==tpmus[CALO]+":"+tpmus[MUONS]+":"+tpmus[MUONS])           shorttype = "2muons1calomu";
	
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]) shorttype = "0muons3tpmuA";
	
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]) shorttype = "0muons3tpmuB";
	
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPB]) shorttype = "0muons2tpmuA1tmpmuB";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPA]) shorttype = "0muons2tpmuA1tmpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]) shorttype = "0muons2tpmuA1tmpmuB";
	
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]) shorttype = "0muons1tpmuA2tmpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPB]) shorttype = "0muons1tpmuA2tmpmuB";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPA]) shorttype = "0muons1tpmuA2tmpmuB";
	
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]+":"+tpmus[CALO])     shorttype = "0muons2tpmuA1calomu";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[CALO]+":"+tpmus[MUONSTPA])     shorttype = "0muons2tpmuA1calomu";
	else if(type==tpmus[CALO]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPA])     shorttype = "0muons2tpmuA1calomu";
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[MUID])          shorttype = "3muid";
	
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[MUIDTPA])       shorttype = "2muid1tpmuA";
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPA]+":"+tpmus[MUID])       shorttype = "2muid1tpmuA";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUID]+":"+tpmus[MUID])       shorttype = "2muid1tpmuA";
	
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[MUIDTPB])       shorttype = "2muid1tpmuB";
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPB]+":"+tpmus[MUID])       shorttype = "2muid1tpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUID]+":"+tpmus[MUID])       shorttype = "2muid1tpmuB";
	
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPA])    shorttype = "1muid2tpmuA";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUID]+":"+tpmus[MUIDTPA])    shorttype = "1muid2tpmuA";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]+":"+tpmus[MUID])    shorttype = "1muid2tpmuA";
	
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPB])    shorttype = "1muid2tpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUID]+":"+tpmus[MUIDTPB])    shorttype = "1muid2tpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUIDTPB]+":"+tpmus[MUID])    shorttype = "1muid2tpmuB";
	
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPB])    shorttype = "1muid1tpmuA1tpmuB";
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPA])    shorttype = "1muid1tpmuA1tpmuB";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUID]+":"+tpmus[MUIDTPB])    shorttype = "1muid1tpmuA1tpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUID]+":"+tpmus[MUIDTPA])    shorttype = "1muid1tpmuA1tpmuB";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPB]+":"+tpmus[MUID])    shorttype = "1muid1tpmuA1tpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUIDTPA]+":"+tpmus[MUID])    shorttype = "1muid1tpmuA1tpmuB";
	
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[CALO])           shorttype = "2muid1calomu";
	else if(type==tpmus[MUID]+":"+tpmus[CALO]+":"+tpmus[MUID])           shorttype = "2muid1calomu";
	else if(type==tpmus[CALO]+":"+tpmus[MUID]+":"+tpmus[MUID])           shorttype = "2muid1calomu";
	
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]) shorttype = "0muid3tpmuA";
	
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPB]) shorttype = "0muid3tpmuB";
	
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPB]) shorttype = "0muid2tpmuA1tmpmuB";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPA]) shorttype = "0muid2tpmuA1tmpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]) shorttype = "0muid2tpmuA1tmpmuB";
	
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPB]) shorttype = "0muid1tpmuA2tmpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPB]) shorttype = "0muid1tpmuA2tmpmuB";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPA]) shorttype = "0muid1tpmuA2tmpmuB";
	
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]+":"+tpmus[CALO])     shorttype = "0muid2tpmuA1calomu";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[CALO]+":"+tpmus[MUIDTPA])     shorttype = "0muid2tpmuA1calomu";
	else if(type==tpmus[CALO]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPA])     shorttype = "0muid2tpmuA1calomu";

	else _FATAL("Triplet cannot be classified: "+(string)type);

	return shorttype;
}
int classifyTripletCode(TString shorttype)
{
	if     (shorttype=="3muons")              return MUONS3;
	else if(shorttype=="2muons1tpmuA")        return MUONS2TPA1;
	else if(shorttype=="2muons1tpmuB")        return MUONS2TPB1;
	else if(shorttype=="2muons1calomu")       return MUONS2CALO1;
	else if(shorttype=="1muons2tpmuA")        return MUONS1TPA2;
	else if(shorttype=="1muons2tpmuB")        return MUONS1TPB2;
	else if(shorttype=="1muons1tpmuA1tpmuB")  return MUONS1TPA1TPB1;
	else if(shorttype=="0muons3tpmuA")        return MUONS0TPA3;
	else if(shorttype=="0muons3tpmuB")        return MUONS0TPB3;
	else if(shorttype=="0muons2tpmuA1tmpmuB") return MUONS0TPA2TPB1;
	else if(shorttype=="0muons2tpmuA1calomu") return MUONS0TPA2CALO1;
	else if(shorttype=="0muons1tpmuA2tmpmuB") return MUONS0TPA1TPB2;
	
	else if(shorttype=="3muid")              return MUID3;
	else if(shorttype=="2muid1tpmuA")        return MUID2TPA1;
	else if(shorttype=="2muid1tpmuB")        return MUID2TPB1;
	else if(shorttype=="2muid1calomu")       return MUID2CALO1;
	else if(shorttype=="1muid2tpmuA")        return MUID1TPA2;
	else if(shorttype=="1muid2tpmuB")        return MUID1TPB2;
	else if(shorttype=="1muid1tpmuA1tpmuB")  return MUID1TPA1TPB1;
	else if(shorttype=="0muid3tpmuA")        return MUID0TPA3;
	else if(shorttype=="0muid3tpmuB")        return MUID0TPB3;
	else if(shorttype=="0muid2tpmuA1tmpmuB") return MUID0TPA2TPB1;
	else if(shorttype=="0muid2tpmuA1calomu") return MUID0TPA2CALO1;
	else if(shorttype=="0muid1tpmuA2tmpmuB") return MUID0TPA1TPB2;
	
	else _FATAL("Triplet cannot be classified");
	
	return -999;
}
int getSrcType(unsigned int vtx, unsigned int index)
{
	TString srcName = vtx_srcName->at(vtx)[index]; // this is not necessarily ordered
	int srcType = -999;
	for(TMapiTS::iterator it=tpmus.begin() ; it!=tpmus.end() ; ++it)
	{
		if(it->second==srcName) { srcType = it->first; break; }
	}
	if(srcType==-999) _FATAL("Source type could not be found");
	return srcType;
}
bool isMuonSrc(int srcType)
{
	return (srcType>CALO && srcType<=MUID);
}
bool isCaloSrc(int srcType)
{
	return (srcType==CALO);
}
bool isTPmuSrc(int srcType)
{
	return (srcType>MUID);
}
bool isTPaSrc(int srcType)
{
	return (srcType==MUONSTPA || srcType==MUIDTPA);
}
bool isTPbSrc(int srcType)
{
	return (srcType==MUONSTPB || srcType==MUIDTPB);
}
int getSrcCode(int srcType)
{
	if(isMuonSrc(srcType)) return MUON;
	if(isTPaSrc(srcType))  return TPA;
	if(isTPbSrc(srcType))  return TPB;
	if(isCaloSrc(srcType)) return CALO;
	return -1;
}
void getSrc(unsigned int vtx, sources& src)
{	
	src.vtx       = vtx;
	src.type      = classifyTriplet(vtx);
	src.shorttype = classifyTripletShort(vtx);

	//cout << "vtx=" << vtx << ", type=" << src.type << ", src.shorttype=" << src.shorttype << endl;
	for(int i=0 ; i<3 ; ++i)
	{
		src.srcName[i]  = vtx_srcName->at(vtx)[i];
		src.trkIndex[i] = vtx_trkIndex->at(vtx)[i];
		src.srcIndex[i] = vtx_srcIndex->at(vtx)[i];
		src.srcType[i]  = getSrcType(vtx,i);
		src.isMuon[i]   = isMuonSrc(src.srcType[i]);
		src.isCalo[i]   = isCaloSrc(src.srcType[i]);
		src.isTPmu[i]   = isTPmuSrc(src.srcType[i]);
		src.isTPa[i]    = isTPaSrc(src.srcType[i]);
		src.isTPb[i]    = isTPbSrc(src.srcType[i]);
		src.srcCode[i]  = getSrcCode(src.srcType[i]);
	
		//cout << "   src.srcName["<<i<<"]=" << src.srcName[i] << endl;
		//cout << "   src.trkIndex["<<i<<"]=" << src.trkIndex[i] << endl;
		//cout << "   src.srcIndex["<<i<<"]=" << src.srcIndex[i] << endl;
		//cout << "   src.srcType["<<i<<"]=" << src.srcType[i] << endl;
		//cout << "   src.isMuon["<<i<<"]=" << src.isMuon[i] << endl;
		//cout << "   src.isCalo["<<i<<"]=" << src.isCalo[i] << endl;
		//cout << "   src.isTPmu["<<i<<"]=" << src.isTPmu[i] << endl;
		//cout << "   src.isTPa["<<i<<"]=" << src.isTPa[i] << endl;
		//cout << "   src.isTPb["<<i<<"]=" << src.isTPb[i] << endl;
		//cout << "   src.srcCode["<<i<<"]=" << src.srcCode[i] << endl;
	}
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3;
	TMapdi srcMap;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV); src.srcTlv[0] = p1; double pT1 = p1.Pt(); srcMap.insert(make_pair(pT1,0));
	p2.SetXYZM(px2,py2,pz2,muonMassMeV); src.srcTlv[1] = p2; double pT2 = p2.Pt(); srcMap.insert(make_pair(pT2,1));
	p3.SetXYZM(px3,py3,pz3,muonMassMeV); src.srcTlv[2] = p3; double pT3 = p3.Pt(); srcMap.insert(make_pair(pT3,2));
	TMapdi::reverse_iterator rit=srcMap.rbegin();
	src.srcOrder[rit->second] = 1; ++rit;
	src.srcOrder[rit->second] = 2; ++rit;
	src.srcOrder[rit->second] = 3;
	
	src.p3body  = p1+p2+p3;
	src.p2body12 = p1+p2;
	src.p2body13 = p1+p3;
	src.p2body23 = p2+p3;
	
	// double q1 = (src.isMuon[0] || src.isCalo[0]) ? qtrk(trks_qoverp->at(src.trkIndex[0])) : qtrk(tpmu_vd[src1+"_qOverP"]->at(src.srcIndex[0]));
	// double q2 = (src.isMuon[1] || src.isCalo[1]) ? qtrk(trks_qoverp->at(src.trkIndex[1])) : qtrk(tpmu_vd[src2+"_qOverP"]->at(src.srcIndex[1]));
	// double q3 = (src.isMuon[2] || src.isCalo[2]) ? qtrk(trks_qoverp->at(src.trkIndex[2])) : qtrk(tpmu_vd[src3+"_qOverP"]->at(src.srcIndex[2]));
	src.q3body  = vtx_charge->at(vtx);
	src.q2body12 = qtrk(trks_qoverp->at(src.trkIndex[0]))+qtrk(trks_qoverp->at(src.trkIndex[1])); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles)
	src.q2body13 = qtrk(trks_qoverp->at(src.trkIndex[0]))+qtrk(trks_qoverp->at(src.trkIndex[2])); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles) 
	src.q2body23 = qtrk(trks_qoverp->at(src.trkIndex[1]))+qtrk(trks_qoverp->at(src.trkIndex[2])); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles) 
}
bool validateSrcChain(unsigned int vtx)
{
	TString type = classifyTripletShort(vtx);
	if(type.Contains("muid")  && !mastername.Contains("muid"))  _FATAL("type="+(string)type+" mastername="+(string)mastername);
	if(type.Contains("muons") && !mastername.Contains("muons")) _FATAL("type="+(string)type+" mastername="+(string)mastername);
	return true;
}
bool validatedVertexChain(unsigned int vtx)
{
	TString shorttype = classifyTripletShort(vtx);
	bool okchain1 = ((shorttype.Contains("muid")  && mastername.Contains("muid")));
	bool okchain2 = ((shorttype.Contains("muons") && mastername.Contains("muons")));
	if(!okchain1 && !okchain2) return false;
	return true;
}
bool isCaloCrack(unsigned int vtx, double etaCrack=0.1) // returns true if all calos belonging to this vertex are in the crack region. otherwise, return false
{
	sources src;
	getSrc(vtx,src);
	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
	int  otrk1  = src.srcOrder[0]; int  otrk2   = src.srcOrder[1]; int  otrk3   = src.srcOrder[2];
	// int  itrk1  = src.trkIndex[0]; int  itrk2   = src.trkIndex[1]; int  itrk3   = src.trkIndex[2];
	int  isrc1  = src.srcIndex[0]; int  isrc2   = src.srcIndex[1]; int  isrc3   = src.srcIndex[2];
	bool isCalo1 = src.isCalo[0];  bool isCalo2 = src.isCalo[1];   bool isCalo3 = src.isCalo[2];
	bool okcalo1 = (isCalo1 && fabs(calo_muons_eta->at(isrc1))<etaCrack  && otrk1!=3);
	bool okcalo2 = (isCalo2 && fabs(calo_muons_eta->at(isrc2))<etaCrack  && otrk2!=3);
	bool okcalo3 = (isCalo3 && fabs(calo_muons_eta->at(isrc3))<etaCrack  && otrk3!=3);
	int nfoundcalos = (isCalo1+isCalo2+isCalo3);
	if(okcalo1+okcalo2+okcalo3<nfoundcalos) return false;
	return true;
}
bool validatedVertex(unsigned int vtx)
{
	//// validate the chain
	TString shorttype = classifyTripletShort(vtx);
	if(!validatedVertexChain(vtx)) return false;
	
	//// restrictions on the calo muons (be only in the crack)
	bool iscalo   = shorttype.Contains("calo");
	// if(skim  && iscalo && !isCaloCrack(vtx)) return false; // Calo only in the crack (only in skim mode)
	// if(!skim && iscalo)                      return false; // NO calo in the analysis mode (!skim)
	if(iscalo)                                  return false; // NO calo ever !!!
	
	//// Do not allow SegmentTagTrackParticles (muons chain only)
	bool istagged = (shorttype.Contains("muons") && shorttype.Contains("tpmuB"));
	if(!skim && istagged) return false; // NO tagged in the skim level
	// if(istagged)             return false; // NO tagged ever !!!
	
	//// if you wish to keep only 3CBmuons vertices
	// bool is3mu = (shorttype.Contains("3muons") || shorttype.Contains("3muid"));
	// if(!is3mu) return false;
	
	return true;
}


//////////////////////////////


TLorentzVector getTlv3mu(unsigned int vtx)
{
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];	
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3,psum;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
	psum = p1+p2+p3;
	return psum;
}
double momentumBalanceSig(int idtrk, int tptrk, TString collection)
{	
	double diff_qoverp     = tpmu_vd[collection+"_qOverP"]->at(tptrk)-trks_qoverp->at(idtrk);
	double err_diff_qoverp = sqrt(tpmu_vd[collection+"_qOverPErr"]->at(tptrk) + trks_qoverpErr->at(idtrk));
	// Signed momentum difference significance
	return  (err_diff_qoverp>0.) ? diff_qoverp / err_diff_qoverp : 999.;
}
TString getMuonCollection(int i)
{
	TString name = "";
	switch(i)
	{		
		case MUONS: name = "Muons";                                      break;
		case MUID:  name = "MuidMuonCollection";                         break;
		default: _FATAL("there is no muon container with index "+_s(i)); break;
	}
	return name;
}
TString getCollection(int i)
{
	TString name = "";
	switch(i)
	{		
		case MUONSTPA: name = "CombinedFitMuonParticles";           break;
		case MUONSTPB: name = "SegmentTagTrackParticles";           break;
		case MUONSTPC: name = "MuGirlRefittedTrackParticles";       break;
		case MUIDTPA:  name = "MuidCombTrackParticles";             break;
		case MUIDTPB:  name = "MuTagIMOTrackParticles";             break;
		case MUIDTPC:  name = "MuGirlRefittedTrackParticles";       break;
		default: _FATAL("there is no container with index "+_s(i)); break;
	}
	return name;
}
void setTPmusPhi()
{
	for(int i=MUONSTPA ; i<=MUIDTPC ; ++i)
	{
		if(mastername.Contains("muons") && i>MUONSTPC) continue; // muons chain (see sourcecode enum)
		if(mastername.Contains("muid")  && i<MUIDTPA)  continue; // muid chains (see sourcecode enum)
		TString collection = getCollection(i);
		tpmu_vd[collection+"_phi"]->clear();
		unsigned int ntp = tpmu_vd[collection+"_pt"]->size();
		for(unsigned int tp=0 ; tp<ntp ; ++tp)
		{
			double px = tpmu_vd[collection+"_px"]->at(tp);
			double py = tpmu_vd[collection+"_py"]->at(tp);
			double pz = tpmu_vd[collection+"_pz"]->at(tp);
			tpmu_vd[collection+"_phi"]->push_back(vtrk(px,py,pz).Phi());
		}
	}
}
bool isMatchedTP2Muon(float dr,float dpt, float dq)
{
	return (dr<0.001 && dq==0. && dpt<0.01);
}
int getorder(TMapdi& m, int val)
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
vector<int> getSignalMuons(bool isMC, TString type) // type="All" or "Acc"
{
	vector<int> sorted_truth_muons(3,-1);
	if(!isMC) return sorted_truth_muons;
	
	TMapdi muMapTruAll;
	TMapdi muMapTruAcc;
	unsigned int nparticles = mc_pdgId->size();
	for(unsigned int imc=0 ; imc<nparticles ; imc++)
	{
		int pdgId = mc_pdgId->at(imc);
		int status = mc_status->at(imc);
		
		bool istau = ((pdgId==15 || pdgId==-15) && (status==2 || status==10902));
		if(!istau) continue;
		
		int nchildren = mc_children->at(imc).size();
		if(nchildren!=3) continue;
		
		int imu1 = mc_children->at(imc)[0]; int pdgId1 = mc_pdgId->at(imu1); int status1 = mc_status->at(imu1);
		int imu2 = mc_children->at(imc)[1]; int pdgId2 = mc_pdgId->at(imu2); int status2 = mc_status->at(imu2);
		int imu3 = mc_children->at(imc)[2]; int pdgId3 = mc_pdgId->at(imu3); int status3 = mc_status->at(imu3);
		bool ismu1 = ((pdgId1==13 || pdgId1==-13) && status1==1);
		bool ismu2 = ((pdgId2==13 || pdgId2==-13) && status2==1);
		bool ismu3 = ((pdgId3==13 || pdgId3==-13) && status3==1);
		if((ismu1+ismu2+ismu3)!=3) continue;

		int nparents = mc_parents->at(imc).size();
		if(nparents!=1) continue;
		
		int iparent = mc_parents->at(imc)[0];
		int pdgIdParent = mc_pdgId->at(iparent);
		if(!(pdgIdParent==24 || pdgIdParent==-24)) continue;
		
		// all signal triplets
		muMapTruAll.insert(make_pair(mc_pt->at(imu1),imu1));
		muMapTruAll.insert(make_pair(mc_pt->at(imu2),imu2));
		muMapTruAll.insert(make_pair(mc_pt->at(imu3),imu3));
		
		// in acceptance signal muons
		bool ismu1acc = isInAcceptance(mc_pt->at(imu1),mc_eta->at(imu1));
		bool ismu2acc = isInAcceptance(mc_pt->at(imu2),mc_eta->at(imu2));
		bool ismu3acc = isInAcceptance(mc_pt->at(imu3),mc_eta->at(imu3));
		if(ismu1acc) muMapTruAcc.insert(make_pair(mc_pt->at(imu1),imu1));
		if(ismu2acc) muMapTruAcc.insert(make_pair(mc_pt->at(imu2),imu2));
		if(ismu3acc) muMapTruAcc.insert(make_pair(mc_pt->at(imu3),imu3));
		
		/////////////////////////////////////////////////////
		// allow only one signal W->tau->3mu per event... ///
		break; //////////////////////////////////////////////
		/////////////////////////////////////////////////////
	}
	
	_DEBUG("");
	
	if(type=="All")
	{
		for(TMapdi::reverse_iterator rit=muMapTruAll.rbegin() ; rit!=muMapTruAll.rend() ; ++rit)
		{
			int itru = rit->second;
			int order = getorder(muMapTruAll,itru);
			if(order>=0) sorted_truth_muons[order] = itru;
		}
	}
	else if(type=="Acc")
	{
		for(TMapdi::reverse_iterator rit=muMapTruAcc.rbegin() ; rit!=muMapTruAcc.rend() ; ++rit)
		{
			int itru = rit->second;
			int order = getorder(muMapTruAll,itru);
			if(order>=0) sorted_truth_muons[order] = itru;
		}
	}
	else _FATAL("type="+(string)type+" is not supported");
	
	return sorted_truth_muons;
}
void simulateTripletBuilder(TString channel, TMapTSP2TH1& histos1, TString primaryMuon, TString primaryTP, TString secondaryTP)
{
	int iPrimaryMuon = -1;
	int iPrimaryTP   = -1;
	int iSecondaryTP = -1;
	
	_DEBUG("");
	
	TString histname = channel+"_tripletCategories_noVertexing";
	TString histname1 = channel+"_tripletCategories_norm_noVertexing";
	
	///////////////////////////////////////////////////////////////
	if     (primaryMuon=="Muons")               iPrimaryMuon = MUONS;
	else if(primaryMuon=="MuidMuonCollection")  iPrimaryMuon = MUID;
	else _FATAL("no muon collection named: "+(string)primaryMuon);
	
	if     (primaryMuon=="Muons"               && primaryTP=="CombinedFitMuonParticles" && secondaryTP=="SegmentTagTrackParticles")       { iPrimaryTP = MUONSTPA; iSecondaryTP = MUONSTPB; }
	else if(primaryMuon=="Muons"               && primaryTP=="CombinedFitMuonParticles" && secondaryTP=="MuGirlRefittedTrackParticles")   { iPrimaryTP = MUONSTPA; iSecondaryTP = MUONSTPC; }
	else if(primaryMuon=="MuidMuonCollection"  && primaryTP=="MuidCombTrackParticles"   && secondaryTP=="MuTagIMOTrackParticles")         { iPrimaryTP = MUIDTPA; iSecondaryTP = MUIDTPB; }
	else if(primaryMuon=="MuidMuonCollection"  && primaryTP=="MuidCombTrackParticles"   && secondaryTP=="MuGirlRefittedTrackParticles")   { iPrimaryTP = MUIDTPA; iSecondaryTP = MUIDTPC; }
	else _FATAL("primary TP collection named: "+(string)primaryTP+" or secondary TP collection named: "+(string)secondaryTP+" is inconsistent with this muon collection: "+(string)primaryMuon);
	///////////////////////////////////////////////////////////////
	
	bool isMuid = primaryMuon.Contains("Muid");
	vector<double>* muref_pt     = isMuid ? muid_pt                    : muons_pt;
	vector<double>* muref_eta    = isMuid ? muid_eta                   : muons_eta;
	vector<double>* muref_phi    = isMuid ? muid_phi                   : muons_phi;
	vector<int>*    muref_charge = isMuid ? muid_charge                : muons_charge;
	vector<bool>*   muref_isCB   = isMuid ? muid_isCombined            : muons_isCombined;
	vector<bool>*   muref_hasID  = isMuid ? muid_hasInDetTrackParticle : muons_hasInDetTrackParticle;
	
	// histos1[histname]->AddBinContent(1); // all events
	
	// if(!istrigger) return;
	// histos1[histname]->AddBinContent(2);
	
	// vector<int> sorted_truth_muons_inacc = getSignalMuons(isMC,"Acc");
	// bool isinacceptance = (isMC) ? (sorted_truth_muons_inacc[0]>=0 && sorted_truth_muons_inacc[1]>=0 && sorted_truth_muons_inacc[2]>=0) : true;
	// if(!isinacceptance) return;
	// histos1[histname]->AddBinContent(3);
	
	/////////////////////////////////////////////////////////////////////////
	//// should decide here which Muon collection tu use: Muons/Muid/... ////
	/////////////////////////////////////////////////////////////////////////
	TString muname = getMuonCollection(iPrimaryMuon);
	TString muprefix = muname+"_";
	
	unsigned int nmuons = 0;
	for(unsigned int muRef=0 ; muRef<muref_pt->size() ; ++muRef)
	{
		bool isCB = muref_isCB->at(muRef);
		bool hasID = muref_hasID->at(muRef);
		if(!hasID) continue;
		if(!isCB)  continue;
		nmuons++;
	}

	TString name_tpA = getCollection(iPrimaryTP);
	TString prefix_tpA = name_tpA+"_";
	unsigned int ntpAs = tpmu_vd[prefix_tpA+"pt"]->size();
	
	TString name_tpB = getCollection(iSecondaryTP);
	TString prefix_tpB = name_tpB+"_";
	unsigned int ntpBs = tpmu_vd[prefix_tpB+"pt"]->size();	

	unsigned int ncalos = calo_muons_pt->size();
	
	_DEBUG("");
	

	// remove muon-tpA/tpB overlaps
	vector<int> matched_tpAs;
	vector<int> matched_tpBs;
	vector<int> matched_calo;
	vector<int> nonmatched_tpAs;
	vector<int> nonmatched_tpBs;
	vector<int> nonmatched_calo;
	for(unsigned int muRef=0 ; muRef<muref_pt->size() ; ++muRef)
	{
		bool isCB = muref_isCB->at(muRef);
		bool hasID = muref_hasID->at(muRef);
		if(!hasID) continue;
		if(!isCB)  continue;
		
		TLorentzVector refMuon;
		refMuon.SetPtEtaPhiM(muref_pt->at(muRef),muref_eta->at(muRef),muref_phi->at(muRef),muonMassMeV);
		for(unsigned int muTst=0 ; muTst<ntpAs ; ++muTst)
		{
			TLorentzVector tstMuon;
			tstMuon.SetPtEtaPhiM(tpmu_vd[prefix_tpA+"pt"]->at(muTst),tpmu_vd[prefix_tpA+"eta"]->at(muTst),tpmu_vd[prefix_tpA+"phi"]->at(muTst),muonMassMeV);
			double dr  = refMuon.DeltaR(tstMuon);
			double dq  = muref_charge->at(muRef)-qtrk(tpmu_vd[prefix_tpA+"qOverP"]->at(muTst));
			double dpt = (muref_pt->at(muRef)!=0.) ? fabs(muref_pt->at(muRef)-tpmu_vd[prefix_tpA+"pt"]->at(muTst))/muref_pt->at(muRef) : 1.e20;
			bool matched = isMatchedTP2Muon(dr,dpt,dq);
			if(matched && findInVector(matched_tpAs,muTst)==-1) matched_tpAs.push_back(muTst);
		}
		for(unsigned int muTst=0 ; muTst<ntpBs ; ++muTst)
		{
			TLorentzVector tstMuon;
			tstMuon.SetPtEtaPhiM(tpmu_vd[prefix_tpB+"pt"]->at(muTst),tpmu_vd[prefix_tpB+"eta"]->at(muTst),tpmu_vd[prefix_tpB+"phi"]->at(muTst),muonMassMeV);
			double dr  = refMuon.DeltaR(tstMuon);
			double dq  = muref_charge->at(muRef)-qtrk(tpmu_vd[prefix_tpB+"qOverP"]->at(muTst));
			double dpt = (muref_pt->at(muRef)!=0.) ? fabs(muref_pt->at(muRef)-tpmu_vd[prefix_tpB+"pt"]->at(muTst))/muref_pt->at(muRef) : 1.e20;
			bool matched = isMatchedTP2Muon(dr,dpt,dq);
			if(matched && findInVector(matched_tpBs,muTst)==-1) matched_tpBs.push_back(muTst);
		}
		for(unsigned int muTst=0 ; muTst<ncalos ; ++muTst)
		{
			if(fabs(calo_muons_eta->at(muTst))>0.1) continue;
			
			TLorentzVector tstMuon;
			tstMuon.SetPtEtaPhiM(calo_muons_pt->at(muTst),calo_muons_eta->at(muTst),calo_muons_phi->at(muTst),muonMassMeV);
			double dr  = refMuon.DeltaR(tstMuon);
			double dq  = muref_charge->at(muRef)-calo_muons_charge->at(muTst);
			double dpt = (muref_pt->at(muRef)!=0.) ? fabs(muref_pt->at(muRef)-calo_muons_pt->at(muTst))/muref_pt->at(muRef) : 1.e20;
			bool matched = isMatchedTP2Muon(dr,dpt,dq);
			if(matched && findInVector(matched_calo,muTst)==-1) matched_calo.push_back(muTst);
		}
	}
	for(unsigned int k=0 ; k<ntpAs ; ++k) { if(findInVector(matched_tpAs,k)==-1) nonmatched_tpAs.push_back(k); }
	for(unsigned int k=0 ; k<ntpBs ; ++k) { if(findInVector(matched_tpBs,k)==-1) nonmatched_tpBs.push_back(k); }
	unsigned int ncaloseta015 = 0;
	for(unsigned int k=0 ; k<ncalos ; ++k)
	{
		if(fabs(calo_muons_eta->at(k))>0.1) continue;
		ncaloseta015++;
		if(findInVector(matched_calo,k)==-1) nonmatched_calo.push_back(k);
	}
	unsigned int nmatchedtpAs = matched_tpAs.size();
	unsigned int nmatchedtpBs = matched_tpBs.size();
	unsigned int nmatchedcalo = matched_calo.size();
	unsigned int nnonmatchedtpAs = nonmatched_tpAs.size(); if(nnonmatchedtpAs!=ntpAs-nmatchedtpAs)        _FATAL("nnonmatchedtpAs!=ndtpAs-nmatchedtpAs");
	unsigned int nnonmatchedtpBs = nonmatched_tpBs.size(); if(nnonmatchedtpBs!=ntpBs-nmatchedtpBs)        _FATAL("nnonmatchedtpBs!=ndtpBs-nmatchedtpBs");
	unsigned int nnonmatchedcalo = nonmatched_calo.size(); if(nnonmatchedcalo!=ncaloseta015-nmatchedcalo) _FATAL("nnonmatchedcalo!=ncaloseta015-nmatchedcalo");
	
	_DEBUG("");
	
	// remove tpA-tpB overlaps
	vector<int> matched_tpA_tpBs;
	vector<int> nonmatched_tpA_tpBs;
	for(unsigned int j=0 ; j<nnonmatchedtpAs ; ++j)
	{
		int muRef = nonmatched_tpAs[j];
		TLorentzVector refMuon;
		refMuon.SetPtEtaPhiM(tpmu_vd[prefix_tpA+"pt"]->at(muRef),tpmu_vd[prefix_tpA+"eta"]->at(muRef),tpmu_vd[prefix_tpA+"phi"]->at(muRef),muonMassMeV);
		for(unsigned int k=0 ; k<nnonmatchedtpBs ; ++k)
		{
			int muTst = nonmatched_tpBs[k];
			TLorentzVector tstMuon;
			tstMuon.SetPtEtaPhiM(tpmu_vd[prefix_tpB+"pt"]->at(muTst),tpmu_vd[prefix_tpB+"eta"]->at(muTst),tpmu_vd[prefix_tpB+"phi"]->at(muTst),muonMassMeV);
			double dr  = refMuon.DeltaR(tstMuon);
			double dq  = qtrk(tpmu_vd[prefix_tpA+"qOverP"]->at(muRef))-qtrk(tpmu_vd[prefix_tpB+"qOverP"]->at(muTst));
			double dpt = (tpmu_vd[prefix_tpA+"pt"]->at(muRef)!=0.) ? fabs(tpmu_vd[prefix_tpA+"pt"]->at(muRef)-tpmu_vd[prefix_tpB+"pt"]->at(muTst))/tpmu_vd[prefix_tpA+"pt"]->at(muRef) : 1.e20;
			bool matched = isMatchedTP2Muon(dr,dpt,dq);
			if(matched && findInVector(matched_tpA_tpBs,muTst)==-1) matched_tpA_tpBs.push_back(muTst);
		}
	}
	for(unsigned int k=0 ; k<nnonmatchedtpBs ; ++k)
	{
		int imu = nonmatched_tpBs[k];
		if(findInVector(matched_tpA_tpBs,imu)==-1) nonmatched_tpA_tpBs.push_back(imu);
	}
	// unsigned int nmatchedtpAtpBs    = matched_tpA_tpBs.size();
	unsigned int nnonmatchedtpAtpBs = nonmatched_tpA_tpBs.size();
	
	_DEBUG("");
	
	// remove tpA-calo overlaps
	vector<int> matched_tpA_calo;
	vector<int> nonmatched_tpA_calo;
	for(unsigned int j=0 ; j<nnonmatchedtpAs ; ++j)
	{
		int muRef = nonmatched_tpAs[j];
		TLorentzVector refMuon;
		refMuon.SetPtEtaPhiM(tpmu_vd[prefix_tpA+"pt"]->at(muRef),tpmu_vd[prefix_tpA+"eta"]->at(muRef),tpmu_vd[prefix_tpA+"phi"]->at(muRef),muonMassMeV);
		for(unsigned int k=0 ; k<nnonmatchedcalo ; ++k)
		{
			int muTst = nonmatched_calo[k];
			TLorentzVector tstMuon;
			tstMuon.SetPtEtaPhiM(calo_muons_pt->at(muTst),calo_muons_eta->at(muTst),calo_muons_phi->at(muTst),muonMassMeV);
			double dr  = refMuon.DeltaR(tstMuon);
			double dq  = qtrk(tpmu_vd[prefix_tpA+"qOverP"]->at(muRef))-calo_muons_charge->at(muTst);
			double dpt = (tpmu_vd[prefix_tpA+"pt"]->at(muRef)!=0.) ? fabs(tpmu_vd[prefix_tpA+"pt"]->at(muRef)-calo_muons_pt->at(muTst))/tpmu_vd[prefix_tpA+"pt"]->at(muRef) : 1.e20;
			bool matched = isMatchedTP2Muon(dr,dpt,dq);
			if(matched && findInVector(matched_tpA_calo,muTst)==-1) matched_tpA_calo.push_back(muTst);
		}
	}
	for(unsigned int k=0 ; k<nnonmatchedcalo ; ++k)
	{
		int imu = nonmatched_calo[k];
		if(findInVector(matched_tpA_calo,imu)==-1) nonmatched_tpA_calo.push_back(imu);
	}
	// unsigned int nmatchedtpAcalo    = matched_tpA_calo.size();
	// unsigned int nnonmatchedtpAcalo = nonmatched_tpA_calo.size();
	
	_DEBUG("");
	
	// fill the bins
	if(nmuons>2)                                                 { histos1[histname]->AddBinContent(1);  histos1[histname1]->AddBinContent(1);  return; } // 3mu
	if(nmuons==2 && nnonmatchedtpAs>0)                           { histos1[histname]->AddBinContent(2);  histos1[histname1]->AddBinContent(2);  return; } // 2mu+1tpA
	if(nmuons==2 && nnonmatchedtpAs==0 && nnonmatchedtpBs>0)     { histos1[histname]->AddBinContent(3);  histos1[histname1]->AddBinContent(3);  return; } // 2mu+1tpB | 0tpA
	if(nmuons==1 && nnonmatchedtpAs>1)                           { histos1[histname]->AddBinContent(4);  histos1[histname1]->AddBinContent(4);  return; } // 1mu+2tpA 
	if(nmuons==1 && nnonmatchedtpAs==0 && nnonmatchedtpBs>1)     { histos1[histname]->AddBinContent(5);  histos1[histname1]->AddBinContent(5);  return; } // 1mu+2tpB | 0tpA
	if(nmuons==1 && nnonmatchedtpAs==1 && nnonmatchedtpAtpBs>0)  { histos1[histname]->AddBinContent(6);  histos1[histname1]->AddBinContent(6);  return; } // 1mu+1tpA+1tpB
	if(nmuons==2 && nnonmatchedtpAs==0 && nnonmatchedcalo>0)     { histos1[histname]->AddBinContent(7);  histos1[histname1]->AddBinContent(7);  return; } // 2mu+1calo | 0tpA
	if(nmuons==0 && nnonmatchedtpAs>2)                           { histos1[histname]->AddBinContent(8);  histos1[histname1]->AddBinContent(8);  return; } // 3tpA | 0mu
	if(nmuons==0 && nnonmatchedtpAs==0 && nnonmatchedtpAtpBs>2)  { histos1[histname]->AddBinContent(9);  histos1[histname1]->AddBinContent(9);  return; } // 3tpB | 0mu,0tpA
	if(nmuons==0 && nnonmatchedtpAs==2 && nnonmatchedtpAtpBs>0)  { histos1[histname]->AddBinContent(10); histos1[histname1]->AddBinContent(10); return; } // 2tpA+1tpB | 0mu
	if(nmuons==0 && nnonmatchedtpAs==1 && nnonmatchedtpAtpBs>1)  { histos1[histname]->AddBinContent(11); histos1[histname1]->AddBinContent(11); return; } // 1tpA+2tpB | 0mu
	if(nmuons==0 && nnonmatchedtpAs==2 && nnonmatchedcalo>0)     { histos1[histname]->AddBinContent(12); histos1[histname1]->AddBinContent(12); return; } // 2tpA+1calo | 0mu
}
void fillCategories(vector<unsigned int>& ivtx, TString channel, TString hname, TMapTSP2TH1& histos1)
{
	TMapTSi categories;
	categories.insert(make_pair("3muons", 0));
	categories.insert(make_pair("2muons1tpmuA",0));
	categories.insert(make_pair("2muons1tpmuB",0));
	categories.insert(make_pair("1muons2tpmuA",0));
	categories.insert(make_pair("1muons2tpmuB",0));
	categories.insert(make_pair("1muons1tpmuA1tpmuB",0));
	categories.insert(make_pair("2muons1calomu",0));
	categories.insert(make_pair("0muons3tpmuA",0));
	categories.insert(make_pair("0muons3tpmuB",0));
	categories.insert(make_pair("0muons2tpmuA1tmpmuB",0));
	categories.insert(make_pair("0muons1tpmuA2tmpmuB",0));
	categories.insert(make_pair("0muons2tpmuA1calomu",0));
	
	_DEBUG("");
	
	for(unsigned int i=0 ; i<ivtx.size() ; i++)
	{	
		unsigned int vtx = ivtx[i];
		TString shorttype = classifyTripletShort(vtx);
		if(!validatedVertex(vtx)) continue;
		categories[shorttype]++;
	}
	
	_DEBUG("");
	
	TString histname = channel+"_"+hname;
	if(categories["3muons"]>0)              { histos1[histname]->AddBinContent(1);  return; } // 3+ muons
	if(categories["2muons1tpmuA"]>0)        { histos1[histname]->AddBinContent(2);  return; } // 2mu+1tpA
	if(categories["2muons1tpmuB"]>0)        { histos1[histname]->AddBinContent(3);  return; } // 2mu+1tpB | 0tpA
	if(categories["1muons2tpmuA"]>0)        { histos1[histname]->AddBinContent(4);  return; } // 1mu+2tpA 
	if(categories["1muons2tpmuB"]>0)        { histos1[histname]->AddBinContent(5);  return; } // 1mu+2tpB | 0tpA
	if(categories["1muons1tpmuA1tpmuB"]>0)  { histos1[histname]->AddBinContent(6);  return; } // 1mu+1tpA+1tpB
	if(categories["2muons1calomu"]>0)       { histos1[histname]->AddBinContent(7);  return; } // 2mu+1calo | 0tpA
	if(categories["0muons3tpmuA"]>0)        { histos1[histname]->AddBinContent(8);  return; } // 3tpA | 0mu
	if(categories["0muons3tpmuB"]>0)        { histos1[histname]->AddBinContent(9);  return; } // 3tpB | 0mu,0tpA
	if(categories["0muons2tpmuA1tmpmuB"]>0) { histos1[histname]->AddBinContent(10); return; } // 2tpA+1tpB | 0mu
	if(categories["0muons1tpmuA2tmpmuB"]>0) { histos1[histname]->AddBinContent(11); return; } // 1tpA+2tpB | 0mu
	if(categories["0muons2tpmuA1calomu"]>0) { histos1[histname]->AddBinContent(12); return; } // 2tpA+1calo | 0mu
	
	_DEBUG("");
}
void fillCategories(unsigned int vtx, TString channel, TString hname, TMapTSP2TH1& histos1)
{
	TMapTSi categories;
	categories.insert(make_pair("3muons", 0));
	categories.insert(make_pair("2muons1tpmuA",0));
	categories.insert(make_pair("2muons1tpmuB",0));
	categories.insert(make_pair("1muons2tpmuA",0));
	categories.insert(make_pair("1muons2tpmuB",0));
	categories.insert(make_pair("1muons1tpmuA1tpmuB",0));
	categories.insert(make_pair("2muons1calomu",0));
	categories.insert(make_pair("0muons3tpmuA",0));
	categories.insert(make_pair("0muons3tpmuB",0));
	categories.insert(make_pair("0muons2tpmuA1tmpmuB",0));
	categories.insert(make_pair("0muons1tpmuA2tmpmuB",0));
	categories.insert(make_pair("0muons2tpmuA1calomu",0));
	
	_DEBUG("");
	
	TString shorttype = classifyTripletShort(vtx);
	if(!validatedVertex(vtx)) _FATAL("Not a valid chain");
	categories[shorttype]++;
	
	_DEBUG("");
	
	TString histname = channel+"_"+hname;
	if(categories["3muons"]>0)              { histos1[histname]->AddBinContent(1);  return; } // 3+ muons
	if(categories["2muons1tpmuA"]>0)        { histos1[histname]->AddBinContent(2);  return; } // 2mu+1tpA
	if(categories["2muons1tpmuB"]>0)        { histos1[histname]->AddBinContent(3);  return; } // 2mu+1tpB | 0tpA
	if(categories["1muons2tpmuA"]>0)        { histos1[histname]->AddBinContent(4);  return; } // 1mu+2tpA 
	if(categories["1muons2tpmuB"]>0)        { histos1[histname]->AddBinContent(5);  return; } // 1mu+2tpB | 0tpA
	if(categories["1muons1tpmuA1tpmuB"]>0)  { histos1[histname]->AddBinContent(6);  return; } // 1mu+1tpA+1tpB
	if(categories["2muons1calomu"]>0)       { histos1[histname]->AddBinContent(7);  return; } // 2mu+1calo | 0tpA
	if(categories["0muons3tpmuA"]>0)        { histos1[histname]->AddBinContent(8);  return; } // 3tpA | 0mu
	if(categories["0muons3tpmuB"]>0)        { histos1[histname]->AddBinContent(9);  return; } // 3tpB | 0mu,0tpA
	if(categories["0muons2tpmuA1tmpmuB"]>0) { histos1[histname]->AddBinContent(10); return; } // 2tpA+1tpB | 0mu
	if(categories["0muons1tpmuA2tmpmuB"]>0) { histos1[histname]->AddBinContent(11); return; } // 1tpA+2tpB | 0mu
	if(categories["0muons2tpmuA1calomu"]>0) { histos1[histname]->AddBinContent(12); return; } // 2tpA+1calo | 0mu
	
	_DEBUG("");
}
void fillCategoriesdRmin(vector<unsigned int>& ivtx, TString channel, TMapTSP2TH1& histos1)
{	
	_DEBUG("");
	
	TString hname = channel+"_triplet_dRmin_";
	
	for(unsigned int i=0 ; i<ivtx.size() ; i++)
	{	
		unsigned int vtx = ivtx[i];
		TString shorttype = classifyTripletShort(vtx);
		if(!validatedVertex(vtx)) continue;

		_DEBUG("i="+_s(i));

		sources src;
		getSrc(vtx,src);
		TLorentzVector pSum = src.srcTlv[0]+src.srcTlv[1]+src.srcTlv[2];

		_DEBUG("after sum");

		double dRmin = +1.e20;
		dRmin = (pSum.DeltaR(src.srcTlv[0])<dRmin) ? pSum.DeltaR(src.srcTlv[0]) : dRmin;
		dRmin = (pSum.DeltaR(src.srcTlv[1])<dRmin) ? pSum.DeltaR(src.srcTlv[1]) : dRmin;
		dRmin = (pSum.DeltaR(src.srcTlv[2])<dRmin) ? pSum.DeltaR(src.srcTlv[2]) : dRmin;

		_DEBUG("dRmin="+_s(dRmin));

		if(shorttype=="3muons")              histos1[hname+"3muons"]->Fill(dRmin);      
		if(shorttype=="2muons1tpmuA")        histos1[hname+"2muons_1tpmu"]->Fill(dRmin);
		if(shorttype=="2muons1tpmuB")        histos1[hname+"2muons_1tpmu"]->Fill(dRmin);
		if(shorttype=="1muons2tpmuA")        histos1[hname+"1muons_2tpmu"]->Fill(dRmin);
		if(shorttype=="1muons2tpmuB")        histos1[hname+"1muons_2tpmu"]->Fill(dRmin);
		if(shorttype=="1muons1tpmuA1tpmuB")  histos1[hname+"1muons_2tpmu"]->Fill(dRmin);
		if(shorttype=="2muons1calomu")       histos1[hname+"2muons_1calo"]->Fill(dRmin);
		if(shorttype=="0muons3tpmuA")        histos1[hname+"0muons_3tpmu"]->Fill(dRmin);
		if(shorttype=="0muons3tpmuB")        histos1[hname+"0muons_3tpmu"]->Fill(dRmin);
		if(shorttype=="0muons2tpmuA1tmpmuB") histos1[hname+"0muons_3tpmu"]->Fill(dRmin);
		if(shorttype=="0muons1tpmuA2tmpmuB") histos1[hname+"0muons_3tpmu"]->Fill(dRmin);
		if(shorttype=="0muons2tpmuA1calomu") histos1[hname+"0muons_3tpmu"]->Fill(dRmin);
		
		_DEBUG("");
	}
		
	_DEBUG("");
}



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

int findPVindex(int pvtype=1)
{
	int vtxpvindex = -1;
	for(unsigned int pv=0 ; pv<pv_index->size() ; pv++)
	{
		if(pv_type->at(pv)==pvtype)
		{
			vtxpvindex = pv_index->at(pv);
			break;
		}
	}
	return vtxpvindex;
}
void findPVindex(vector<int>& vpv, int pvtype=1)
{
	vpv.clear();
	for(unsigned int pv=0 ; pv<pv_index->size() ; pv++)
	{
		if(pv_type->at(pv)==pvtype) vpv.push_back( pv_index->at(pv) );
	}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


TMapVL getDoubletsFromTriplet(unsigned int vtx, int& iUniqueCharge)
{
	TMapVL pairsTLV;
	if(fabs(vtx_charge->at(vtx))!=1.) return pairsTLV;
	
	sources src;
	getSrc(vtx,src);
	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
	
	int srcorder1 = src.srcOrder[0];
	int srcorder2 = src.srcOrder[1];
	int srcorder3 = src.srcOrder[2];
	int i1=-1;
	int i2=-1;
	int i3=-1;
	if     (srcorder1==1 && srcorder2==2 && srcorder3==3) { i1=0; i2=1; i3=2; }
	else if(srcorder1==1 && srcorder2==3 && srcorder3==2) { i1=0; i2=2; i3=1; }
	else if(srcorder1==2 && srcorder2==1 && srcorder3==3) { i1=1; i2=0; i3=2; }
	else if(srcorder1==2 && srcorder2==3 && srcorder3==1) { i1=1; i2=2; i3=0; }
	else if(srcorder1==3 && srcorder2==1 && srcorder3==2) { i1=2; i2=0; i3=1; }
	else if(srcorder1==3 && srcorder2==2 && srcorder3==1) { i1=2; i2=1; i3=0; }
	else _FATAL("cannot sort the tracks using the `order` parameter");
	if(i1==-1 || i2==-1 || i3==-1) _FATAL("i1==-1 || i2==-1 || i3==-1");
	
	/////////////////////////////////////
	//// below, everything is sorted ////
	/////////////////////////////////////
	
	double px1 = vtx_reftrks_px->at(vtx)[i1];
	double px2 = vtx_reftrks_px->at(vtx)[i2];
	double px3 = vtx_reftrks_px->at(vtx)[i3];
	double py1 = vtx_reftrks_py->at(vtx)[i1];
	double py2 = vtx_reftrks_py->at(vtx)[i2];
	double py3 = vtx_reftrks_py->at(vtx)[i3];
	double pz1 = vtx_reftrks_pz->at(vtx)[i1];
	double pz2 = vtx_reftrks_pz->at(vtx)[i2];
	double pz3 = vtx_reftrks_pz->at(vtx)[i3];
	TLorentzVector p1,p2,p3, pOS1, pOS2, pSS;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
	
	
	int itrk1 = src.trkIndex[i1];
	int itrk2 = src.trkIndex[i2];
	int itrk3 = src.trkIndex[i3];
	
	// TString src1 = src.srcName[i1];
	// TString src2 = src.srcName[i2];
	// TString src3 = src.srcName[i3];
	// 
	// int isrc1 = src.srcIndex[i1];
	// int isrc2 = src.srcIndex[i2];
	// int isrc3 = src.srcIndex[i3];
	// 
	// bool isMuon1 = src.isMuon[i1];
	// bool isMuon2 = src.isMuon[i2];
	// bool isMuon3 = src.isMuon[i3];
	// 
	// bool isCalo1 = src.isCalo[i1];
	// bool isCalo2 = src.isCalo[i2];
	// bool isCalo3 = src.isCalo[i3];
	// 
	// bool isTPmu1 = src.isTPmu[i1];
	// bool isTPmu2 = src.isTPmu[i2];
	// bool isTPmu3 = src.isTPmu[i3];
	
	// double q1 = (isMuon1 || isCalo1) ? qtrk(trks_qoverp->at(itrk1)) : qtrk(tpmu_vd[src1+"_qOverP"]->at(isrc1));
	// double q2 = (isMuon2 || isCalo2) ? qtrk(trks_qoverp->at(itrk2)) : qtrk(tpmu_vd[src2+"_qOverP"]->at(isrc2));
	// double q3 = (isMuon3 || isCalo3) ? qtrk(trks_qoverp->at(itrk3)) : qtrk(tpmu_vd[src3+"_qOverP"]->at(isrc3));
	
	double q1 = qtrk(trks_qoverp->at(itrk1)); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles)
	double q2 = qtrk(trks_qoverp->at(itrk2)); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles)
	double q3 = qtrk(trks_qoverp->at(itrk3)); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles)
	double q12 = q1+q2;
	double q13 = q1+q3;
	double q23 = q2+q3;
	
	if     (q12==0. && q13==0.) { pOS1=p1+p2; pOS2=p1+p3; pSS=p2+p3; iUniqueCharge = 1; }
	else if(q12==0. && q23==0.) { pOS1=p1+p2; pOS2=p2+p3; pSS=p1+p3; iUniqueCharge = 2; }
	else if(q13==0. && q23==0.) { pOS1=p1+p3; pOS2=p2+p3; pSS=p1+p2; iUniqueCharge = 3; }
	else return pairsTLV;
	pairsTLV.insert(make_pair("OS1", pOS1));
	pairsTLV.insert(make_pair("OS2", pOS2));
	pairsTLV.insert(make_pair("SS",  pSS));
	return pairsTLV;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

bool findTrigger(vector<string>* v, string x)
{
	vector<string>::iterator it = find(v->begin(),v->end(),x);
	if(it!=v->end()) return true;
	return false;
}
void buildTriggerbits()
{
	// should be called only once, at initialization
	triggerorder.clear();
	
	triggerorder.insert(make_pair(1,           "EF_3mu4T"               ));
	triggerorder.insert(make_pair(2,           "EF_3mu6"                ));
	triggerorder.insert(make_pair(3,           "EF_3mu6_MSonly"         ));
	triggerorder.insert(make_pair(4,           "EF_2mu13"               ));
	triggerorder.insert(make_pair(5,           "EF_mu18_tight_mu8_EFFS" ));
	triggerorder.insert(make_pair(6,           "EF_mu18_tight_2mu4_EFFS"));
	triggerorder.insert(make_pair(7,           "EF_2mu8_EFxe30_tclcw"   ));
	if(skim) triggerorder.insert(make_pair(8,  "EF_mu24_tight_EFxe40"   ));
	if(skim) triggerorder.insert(make_pair(9,  "EF_mu24i_tight"         ));
	if(skim) triggerorder.insert(make_pair(10, "EF_mu36_tight"          ));

	/*
	triggerorder.insert(make_pair(1,  "EF_mu18"));
	triggerorder.insert(make_pair(2,  "EF_mu18_MG"));
	triggerorder.insert(make_pair(3,  "EF_mu40_MSonly_barrel"));
	triggerorder.insert(make_pair(4,  "EF_mu15_mu10_EFFS"));
	triggerorder.insert(make_pair(5,  "EF_2mu10_loose"));
	triggerorder.insert(make_pair(6,  "EF_mu15_xe30_noMu"));
	triggerorder.insert(make_pair(7,  "EF_2mu4_DiMu"));
	triggerorder.insert(make_pair(8,  "EF_mu4Tmu6_DiMu"));;
	triggerorder.insert(make_pair(9,  "EF_2mu4T_Bmumu"));
	triggerorder.insert(make_pair(10, "EF_3mu4T"));
	triggerorder.insert(make_pair(11, "EF_3mu6"));
	triggerorder.insert(make_pair(12, "EF_3mu4_MSonly"));
	triggerorder.insert(make_pair(13, "EF_2mu6_MSonly_g10_loose"));
	triggerorder.insert(make_pair(14, "EF_e10_medium_2mu6"));
	triggerorder.insert(make_pair(15, "EF_e10_medium_mu6"));
	triggerorder.insert(make_pair(16, "EF_2mu4T_xe30_noMu"));
	*/
	
	cout << "----------------- Trigger bits -----------------" << endl;
	for(TMapuiTS::iterator it=triggerorder.begin() ; it!=triggerorder.end() ; ++it)
	{
		TString name = it->second;
		cout << "Trigger bit: " << name << endl;
		triggerbits.insert(make_pair(name,0));
	}
	cout << "------------------------------------------------" << endl;
}
void clearTriggerbits()
{
	if(triggerorder.size()==0) _FATAL("Need to call buildTriggerbits() at initialization");
	if(triggerbits.size()==0)  _FATAL("Need to call buildTriggerbits() at initialization");
	for(TMapTSui::iterator it=triggerbits.begin() ; it!=triggerbits.end() ; ++it) it->second = 0;
}
void fillTriggerbits()
{
	if(triggerorder.size()==0) _FATAL("Need to call buildTriggerbits() at initialization");
	if(triggerbits.size()==0)  _FATAL("Need to call buildTriggerbits() at initialization");
	///////////////////////
	clearTriggerbits(); ///
	///////////////////////
	for(TMapTSui::iterator it=triggerbits.begin() ; it!=triggerbits.end() ; ++it)
	{
		string name = (string)it->first;
		bool decision = 0;
		if(it->first.Contains("L1_")) decision = findTrigger(L1_trigger_name,name);
		if(it->first.Contains("EF_")) decision = findTrigger(EF_trigger_name,name);
		it->second = decision;
	}
}
bool getTriggerbit(TString name)
{
	if(triggerorder.size()==0) _FATAL("Need to call buildTriggerbits() at initialization");
	if(triggerbits.size()==0)  _FATAL("Need to call buildTriggerbits() at initialization");	
	return triggerbits[name];
}
bool getUniqueTriggerbit(TString name)
{
	if(triggerorder.size()==0) _FATAL("Need to call buildTriggerbits() at initialization");
	if(triggerbits.size()==0)  _FATAL("Need to call buildTriggerbits() at initialization");
	
	// retrun false if it didn't fire
	bool fired = triggerbits[name];
	if(!fired) return false;
		
	for(TMapuiTS::iterator ii=triggerorder.begin() ; ii!=triggerorder.end() ; ++ii)
	{
		if(ii->second==name)        continue;     // ignore "this" trigger
		if(triggerbits[ii->second]) return false; // it is not unique -> return false
	}
	
	return true; // otherwise, it has fired and it is unique -> return true
}
void fillTriggerbitsHist(TH1* h, TString type="")
{
	for(TMapuiTS::iterator it=triggerorder.begin() ; it!=triggerorder.end() ; ++it)
	{
		unsigned int bin = it->first;
		TString tname = it->second;
		bool decision = (triggerbits[tname]);
		if(decision)
		{
			if     (type=="") h->AddBinContent(bin);
			else if(type=="unique")
			{
				bool isUnique = true; // assume it is unique
				for(TMapuiTS::iterator ii=triggerorder.begin() ; ii!=triggerorder.end() ; ++ii)
				{
					if(ii->first==it->first) continue; // ignore "this" trigger
					if(triggerbits[ii->second]) { isUnique = false; break; } // check if other triggers have fired as well
				}
				if(isUnique) h->AddBinContent(bin);
			}
			else _FATAL("unknown type="+(string)type);
		}
	}
}
bool acceptTrigger()
{
	for(TMapTSui::iterator it=triggerbits.begin() ; it!=triggerbits.end() ; ++it)
	{
		if(it->second) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


void buildTPmus()
{	
	tpmus.insert(make_pair(CALO,     "CaloMuonCollection"));
	tpmus.insert(make_pair(MUID,     "MuidMuonCollection"));
	tpmus.insert(make_pair(MUONS,    "Muons"));
	tpmus.insert(make_pair(MUONSTPA, "CombinedFitMuonParticles"));
	tpmus.insert(make_pair(MUONSTPB, "SegmentTagTrackParticles"));
	tpmus.insert(make_pair(MUONSTPC, "MuGirlRefittedTrackParticles"));
	tpmus.insert(make_pair(MUIDTPA,  "MuidCombTrackParticles"));
	tpmus.insert(make_pair(MUIDTPB,  "MuTagIMOTrackParticles"));
	tpmus.insert(make_pair(MUIDTPC,  "MuGirlRefittedTrackParticles")); // same as for MUONS
}


void initOutTrees(TFile* fout, TMapTSP2TTREE& otrees, TMapTSP2TCHAIN& chains, TMapTSP2TCHAIN& chainfriends)
{
	if(!skim || !fout)
	{
		_ERROR("Calling init out trees with NULL root file");
		return;
	}
	
	_DEBUG("");
	
	fout->cd();
	
	for(TMapTSP2TCHAIN::iterator it=chains.begin() ; it!=chains.end() ; it++)
	{
		TString treename = it->second->GetName();
		otrees.insert(make_pair(treename,it->second->CloneTree(0)));
		cout << "\nMaster tree: " << treename << endl;
	}
	
	// for(TMapTSP2TTREE::iterator it=otrees.begin() ; it!=otrees.end() ; it++)
	// {
	// 	_INFO("otrees -> "+(string)it->first);
	// }
	
	_DEBUG("");
	
	int nfriends = 0;
	TMapTSP2TTREE::iterator oit = otrees.begin();
	// TList* listoffreinds = (TList*)oit->second->GetListOfFriends();	
	TList* listoffreinds;
	if(!oit->second) _FATAL("otrees["+(string)oit->first+"] is NULL ==> nothing to process !");
	if(oit->second->GetListOfFriends()) listoffreinds = (TList*)oit->second->GetListOfFriends();
	else                                _FATAL("list of friends is NULL");
	
	for(TMapTSP2TCHAIN::iterator it=chainfriends.begin() ; it!=chainfriends.end() ; it++)
	{
		TString friendname = it->second->GetName();
		cout << "Friend tree [" << nfriends << "]: " << friendname;
		if(friendname.Contains("PLS1TRK")) continue; // remove the quadruplets tree from the output because it is huge...
		listoffreinds->Remove(listoffreinds->FindObject(friendname)); // remove before clone
		otrees.insert(make_pair(friendname,it->second->CloneTree(0)));
		nfriends++;
		cout << " -> successful." << endl;
	}
		
	olddir->cd();
	
	_DEBUG("");
}
void finitOutTrees(TFile* fout, TMapTSP2TTREE& otrees)
{
	if(!skim || !fout)
	{
		_ERROR("Calling finite out trees with NULL root file");
		return;
	}
	
	_DEBUG("");
	
	fout->cd();
	for(TMapTSP2TTREE::iterator it=otrees.begin() ; it!=otrees.end() ; it++)
	{
		it->second->AutoSave();
		// tout->Write("", TObject::kOverwrite);
	}
	fout->Write();
	fout->Close();
	olddir->cd();
	
	_DEBUG("");
}
void fillOutTrees(TMapTSP2TTREE& otrees)
{
	if(!skim) return;
	
	_DEBUG("");
	
	for(TMapTSP2TTREE::iterator it=otrees.begin() ; it!=otrees.end() ; it++)
	{
		TFile* fout = it->second->GetCurrentFile();
		fout->cd();
		it->second->Fill();

		int nAll = counters["nPassing_evt_all"];
		if(nAll%10000==0 && nAll!=0)
		{
			it->second->FlushBaskets();
			// it->second->Write("", TObject::kOverwrite);
		}
	}
	olddir->cd();
	
	_DEBUG("");
}
void initFlatoutTree(TString name)
{
	if(skim) return;

	_DEBUG("");

	TString flatout_fname = foutname;
	flatout_fname.ReplaceAll("out.","flatout.");
	if(flatout_fname.Contains("muons")) flatout_fname.ReplaceAll("muons.","muons."+smethod+".");
	if(flatout_fname.Contains("muid"))  flatout_fname.ReplaceAll("muid.","muid."+smethod+".");
	flatout_init(flatout_fname,name,olddir);
	flatout_book(olddir);
	
	_DEBUG("");
}
void fillFlatoutTree(vector<vertex>& vertices, int allPassing)
{
	///////////////////////
	// don't produce it ///
	// at the skim level //
	if(skim) return; //////
	///////////////////////
	
	/////////////////////
	// first, clear. ////
	flatout_clear(); ////
	/////////////////////
	
	flatout_ints["vtx_n"] = vertices.size();
	
	for(unsigned int i=0 ; i<vertices.size() ; ++i)
	{
		vertex v = vertices[i];
		unsigned int vtx = v.vtxIndex();
		
		_DEBUG("");
		
		for(unsigned int m=0 ; m<nMaxTracks ; ++m)
		{
			TString x = _s((float)(m+1),0);
			if(m>2 && isFourthTrack[vtx][m]) {  flatout_vfloats["vtx_mQuad"+x]->push_back(v.vtxMquad(m-3)); }
		}
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_order"+x]          ->push_back(v.trkOrder(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_type"+x]           ->push_back(v.trkType(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_pt"+x]           ->push_back(v.trkP(m).Pt()); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_eta"+x]          ->push_back(v.trkP(m).Eta()); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_phi"+x]          ->push_back(v.trkP(m).Phi()); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_sctangsig"+x]    ->push_back(v.trkSctAngMu(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_sctngbsig"+x]    ->push_back(v.trkSctNgbMu(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_pbalsig"+x]      ->push_back(v.trkPbalMu(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_chi2trkfit"+x]   ->push_back(v.trkChi2(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_ndftrkfit"+x]    ->push_back(v.trkNdf(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_chi2ndftrkfit"+x]->push_back(v.trkChi2Ndf(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_pvaltrkfit"+x]   ->push_back(v.trkPval(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_srcqoverp"+x]    ->push_back(v.trkQoverPsrc(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_trkqoverp"+x]    ->push_back(v.trkQoverPtrk(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_ptfrac"+x]       ->push_back(v.trkPtfrac(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vfloats["mu_pixeldEdx"+x]    ->push_back(v.trkPixeldEdx(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_isMedium"+x]       ->push_back(v.trkIsMediumMu(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nPIXhits"+x]       ->push_back(v.trkPIXhits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nDeadPIX"+x]       ->push_back(v.trkDeadPIX(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nPIXholes"+x]      ->push_back(v.trkPIXholes(m)); }	
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nSCThits"+x]       ->push_back(v.trkSCThits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nDeadSCT"+x]       ->push_back(v.trkDeadSCT(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nSCTholes"+x]      ->push_back(v.trkSCTholes(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nTRThits"+x]       ->push_back(v.trkTRThits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nTRToutliers"+x]   ->push_back(v.trkTRToutliers(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_htTRThits"+x]      ->push_back(v.trkHtTRThits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nUsedHitsdEdx"+x]  ->push_back(v.trkUsedHitsdEdx(m)); }
	 	for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nMDThits"+x]->push_back(v.trkMDThits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nTGCPhiHits"+x]->push_back(v.trkTGCPhiHits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nTGCEtaHits"+x]->push_back(v.trkTGCEtaHits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nCSCPhiHits"+x]->push_back(v.trkCSCPhiHits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nCSCEtaHits"+x]->push_back(v.trkCSCEtaHits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nRPCPhiHits"+x]->push_back(v.trkRPCPhiHits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nRPCEtaHits"+x]->push_back(v.trkRPCEtaHits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nCSCEtaHoles"+x]->push_back(v.trkCSCEtaHoles(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nCSCPhiHoles"+x]->push_back(v.trkCSCPhiHoles(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nRPCEtaHoles"+x]->push_back(v.trkRPCEtaHoles(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nRPCPhiHoles"+x]->push_back(v.trkRPCPhiHoles(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nMDTHoles"+x]->push_back(v.trkMDTholes(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nTGCEtaHoles"+x]->push_back(v.trkTGCEtaHoles(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nTGCPhiHoles"+x]->push_back(v.trkTGCPhiHoles(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nOutliersOnTrack"+x]->push_back(v.trkOutliersOnTrack(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_standardDeviationOfChi2OS"+x]->push_back(v.trkStdDevOfChi2OS(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nPrecisionHits"+x]->push_back(v.trkPrecisionHits(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nPhiLayers"+x]->push_back(v.trkPhiLayers(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nEtaPhiLayers"+x]->push_back(v.trkEtaPhiLayers(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nPrecisionHoles"+x]->push_back(v.trkPrecisionHoles(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nEtaTriggerHoleLayers"+x]->push_back(v.trkEtaTriggerHoleLayers(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nPhiHoleLayers"+x]->push_back(v.trkPhiHoleLayers(m)); }
		for(unsigned int m=0 ; m<nMaxTracks ; ++m) { TString x = _s((float)(m+1),0); if(isFourthTrack[vtx][m]) flatout_vints["mu_nPrecisionOutliers"+x]->push_back(v.trkPrecisionOutliers(m)); }
		
		_DEBUG("");

		flatout_vfloats["jet_pt1"]      ->push_back((v.jetN()>0) ? v.jetPM(0).Pt() : 0.);
		flatout_vfloats["jet_pt2"]      ->push_back((v.jetN()>1) ? v.jetPM(1).Pt() : 0.);
		flatout_vfloats["jet_pt3"]      ->push_back((v.jetN()>2) ? v.jetPM(2).Pt() : 0.);
		flatout_vfloats["jet_pt4"]      ->push_back((v.jetN()>3) ? v.jetPM(3).Pt() : 0.);
		flatout_vfloats["jet_eta1"]     ->push_back((v.jetN()>0) ? v.jetPM(0).Eta() : -999.);
		flatout_vfloats["jet_eta2"]     ->push_back((v.jetN()>1) ? v.jetPM(1).Eta() : -999.);
		flatout_vfloats["jet_eta3"]     ->push_back((v.jetN()>2) ? v.jetPM(2).Eta() : -999.);
		flatout_vfloats["jet_eta4"]     ->push_back((v.jetN()>3) ? v.jetPM(3).Eta() : -999.);
		flatout_vfloats["jet_phi1"]     ->push_back((v.jetN()>0) ? v.jetPM(0).Phi() : -999.);
		flatout_vfloats["jet_phi2"]     ->push_back((v.jetN()>1) ? v.jetPM(1).Phi() : -999.);
		flatout_vfloats["jet_phi3"]     ->push_back((v.jetN()>2) ? v.jetPM(2).Phi() : -999.);
		flatout_vfloats["jet_phi4"]     ->push_back((v.jetN()>3) ? v.jetPM(3).Phi() : -999.);
		flatout_vfloats["jet_m1"]       ->push_back((v.jetN()>0) ? v.jetPM(0).M() : 0.);
		flatout_vfloats["jet_m2"]       ->push_back((v.jetN()>1) ? v.jetPM(1).M() : 0.);
		flatout_vfloats["jet_m3"]       ->push_back((v.jetN()>2) ? v.jetPM(2).M() : 0.);
		flatout_vfloats["jet_m4"]       ->push_back((v.jetN()>3) ? v.jetPM(3).M() : 0.);
		flatout_vfloats["jet_E1"]       ->push_back((v.jetN()>0) ? v.jetPM(0).E() : 0.);
		flatout_vfloats["jet_E2"]       ->push_back((v.jetN()>1) ? v.jetPM(1).E() : 0.);
		flatout_vfloats["jet_E3"]       ->push_back((v.jetN()>2) ? v.jetPM(2).E() : 0.);
		flatout_vfloats["jet_E4"]       ->push_back((v.jetN()>3) ? v.jetPM(3).E() : 0.);
		flatout_vfloats["jet_MV1w1"]    ->push_back( v.jetMV1(0));
		flatout_vfloats["jet_MV1w2"]    ->push_back( v.jetMV1(1));
		flatout_vfloats["jet_MV1w3"]    ->push_back( v.jetMV1(2));
		flatout_vfloats["jet_MV1w4"]    ->push_back( v.jetMV1(3));
		flatout_vfloats["jet_vtxf1"]    ->push_back( v.jetVtxFrac(0));
		flatout_vfloats["jet_vtxf2"]    ->push_back( v.jetVtxFrac(1));
		flatout_vfloats["jet_vtxf3"]    ->push_back( v.jetVtxFrac(2));
		flatout_vfloats["jet_vtxf4"]    ->push_back( v.jetVtxFrac(3));
		flatout_vfloats["jet_sumpt12"]  ->push_back( v.jetSumPt());
		flatout_vfloats["jet_dphi3muJ1"]->push_back( v.jetDphi3body());
		flatout_vfloats["jet_dR3muJ1"]  ->push_back( v.jetDR3body());
		flatout_vfloats["jet_dphiJ1J2"] ->push_back( v.jetDphi12());
		flatout_vfloats["jet_dRJ1J2"]   ->push_back( v.jetDR12());
		
		_DEBUG("");
		
		flatout_vstrings["vtx_type"]    ->push_back((string)v.vtxType());
		flatout_vints["vtx_code"]       ->push_back(v.vtxCode());
		flatout_vfloats["vtx_charge"]   ->push_back(v.vtxQ());
		flatout_vfloats["vtx_mass"]     ->push_back(v.vtxM());
		flatout_vfloats["vtx_mOS1"]     ->push_back(v.vtxMOS1());
		flatout_vfloats["vtx_mOS2"]     ->push_back(v.vtxMOS2());
		flatout_vfloats["vtx_mSS"]      ->push_back(v.vtxMSS());
		flatout_vfloats["vtx_pt"]       ->push_back(v.vtxPt());
		flatout_vfloats["vtx_rapidity"] ->push_back(v.vtxP().Rapidity());
		flatout_vfloats["vtx_chi2"]     ->push_back(v.vtxChi2());
		flatout_vfloats["vtx_ndf"]      ->push_back(v.vtxNdf());
		flatout_vfloats["vtx_chi2ndf"]  ->push_back(v.vtxChi2Ndf());
		flatout_vfloats["vtx_pval"]     ->push_back(v.vtxPvalue());
		flatout_vfloats["vtx_lxy"]      ->push_back(v.vtxLxy());
		flatout_vfloats["vtx_lxySig"]   ->push_back((v.vtxLxyErr()!=0.) ? v.vtxLxy()/v.vtxLxyErr() : 0.);
		flatout_vfloats["vtx_a0"]       ->push_back(v.vtxA0());
		flatout_vfloats["vtx_a0xy"]     ->push_back(v.vtxA0xy());
		flatout_vfloats["vtx_cosT"]     ->push_back(fabs(v.vtxCosT()));
		flatout_vfloats["vtx_cosTxy"]   ->push_back(fabs(v.vtxCosTxy()));
		flatout_vfloats["vtx_tau"]      ->push_back(v.vtxTau());
		flatout_vfloats["vtx_ptfrac12"]  ->push_back(v.vtxPtFrac12());
		flatout_vfloats["vtx_ptfrac23"]  ->push_back(v.vtxPtFrac23());
		flatout_vfloats["vtx_ptfrac13"]  ->push_back(v.vtxPtFrac13());
		flatout_vfloats["vtx_dpt12"]     ->push_back(v.vtxDPt12());
		flatout_vfloats["vtx_dpt23"]     ->push_back(v.vtxDPt23());
		flatout_vfloats["vtx_dpt13"]     ->push_back(v.vtxDPt13());
		flatout_vfloats["vtx_dRmax"]     ->push_back(v.vtxDRmax());
		flatout_vfloats["vtx_dRmin"]     ->push_back(v.vtxDRmin());
		
		flatout_vfloats["vtx_isolation000"]->push_back(v.vtxIsolation(0));
		flatout_vfloats["vtx_isolation001"]->push_back(v.vtxIsolation(1));
		flatout_vfloats["vtx_isolation002"]->push_back(v.vtxIsolation(2));
		flatout_vfloats["vtx_isolation003"]->push_back(v.vtxIsolation(3));
		flatout_vfloats["vtx_isolation004"]->push_back(v.vtxIsolation(4));
		flatout_vfloats["vtx_isolation005"]->push_back(v.vtxIsolation(5));
		flatout_vfloats["vtx_isolation006"]->push_back(v.vtxIsolation(6));
		flatout_vfloats["vtx_isolation007"]->push_back(v.vtxIsolation(7));
		flatout_vfloats["vtx_isolation008"]->push_back(v.vtxIsolation(8));
		flatout_vfloats["vtx_isolation009"]->push_back(v.vtxIsolation(9));
		flatout_vfloats["vtx_isolation010"]->push_back(v.vtxIsolation(10));
		flatout_vfloats["vtx_isolation012"]->push_back(v.vtxIsolation(11));
		flatout_vfloats["vtx_isolation014"]->push_back(v.vtxIsolation(12));
		flatout_vfloats["vtx_isolation016"]->push_back(v.vtxIsolation(13));
		flatout_vfloats["vtx_isolation018"]->push_back(v.vtxIsolation(14));
		flatout_vfloats["vtx_isolation020"]->push_back(v.vtxIsolation(15));
		flatout_vfloats["vtx_isolation022"]->push_back(v.vtxIsolation(16));
		flatout_vfloats["vtx_isolation024"]->push_back(v.vtxIsolation(17));
		flatout_vfloats["vtx_isolation026"]->push_back(v.vtxIsolation(18));
		flatout_vfloats["vtx_isolation028"]->push_back(v.vtxIsolation(19));
		flatout_vfloats["vtx_isolation030"]->push_back(v.vtxIsolation(20));
		
		flatout_vfloats["met_reffinal_mT"]     ->push_back(v.metMt());
		flatout_vfloats["met_reffinal_dPhi3mu"]->push_back(v.metDphi3body());
	}
	
	_DEBUG("");
	
	flatout_floats["met_reffinal_et"]      = MET_RefFinal_et;
	flatout_floats["met_reffinal_phi"]     = MET_RefFinal_phi;
	
	_DEBUG("");
	
	flatout_floats["wgt_shapeFONLL"] = weights["wFONLLshape"];
	flatout_floats["wgt_normFONLL"]  = weights["wFONLLflat"];
	flatout_floats["wgt_luminosity"] = weights["wLumi"];
	flatout_floats["wgt_kfactor"]    = weights["wKfac"];
	flatout_floats["wgt_dijets"]     = weights["wDijet"];
	flatout_floats["wgt_pileup"]     = weights["wPileup"];
	flatout_floats["wgt_total"]      = weights["wgt"];
	
	_DEBUG("");
	
	flatout_ints["evt_RunNumber"]         = RunNumber;
	flatout_ints["evt_lbn"]               = lbn;
	flatout_ints["evt_EventNumber"]       = EventNumber;
	flatout_ints["evt_actualIntPerXing"]  = aux_actualIntPerXing;
	flatout_ints["evt_averageIntPerXing"] = aux_averageIntPerXing;
	
	_DEBUG("");
	
	flatout_ints["EF_3mu4T"]                = findTrigger(EF_trigger_name, "EF_3mu4T");
	flatout_ints["EF_3mu6"]                 = findTrigger(EF_trigger_name, "EF_3mu6");
	flatout_ints["EF_3mu6_MSonly"]          = findTrigger(EF_trigger_name, "EF_3mu6_MSonly");
	flatout_ints["EF_2mu13"]                = findTrigger(EF_trigger_name, "EF_2mu13");
	flatout_ints["EF_mu18_tight_mu8_EFFS"]  = findTrigger(EF_trigger_name, "EF_mu18_tight_mu8_EFFS");
	flatout_ints["EF_mu18_tight_2mu4_EFFS"] = findTrigger(EF_trigger_name, "EF_mu18_tight_2mu4_EFFS");
	flatout_ints["EF_2mu8_EFxe30_tclcw"]    = findTrigger(EF_trigger_name, "EF_2mu8_EFxe30_tclcw");
	flatout_ints["EF_mu24_tight_EFxe40"]    = findTrigger(EF_trigger_name, "EF_mu24_tight_EFxe40");
	flatout_ints["EF_mu24i_tight"]          = findTrigger(EF_trigger_name, "EF_mu24i_tight");
	flatout_ints["EF_mu36_tight"]           = findTrigger(EF_trigger_name, "EF_mu36_tight");
	
	_DEBUG("");
	
	////////////////////////////////////////////////
	// then fill the tree. /////////////////////////
	if(!skim) flatout_fill(allPassing,olddir); /////
	////////////////////////////////////////////////
}


// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BtagAnalysis172Winter2013
bool isStandardBtag(float flavor_weight_MV1)
{
	return (flavor_weight_MV1>=0.3900); // 80% nominal efficiency
}
bool isMediumBtag(float flavor_weight_MV1)
{
	return (flavor_weight_MV1>=0.8119); // 70% nominal efficiency
}
bool isTightBtag(float flavor_weight_MV1)
{
	return (flavor_weight_MV1>=0.9867); // 60% nominal efficiency
}
bool isMV1btag(bool isTight=true)
{
	// for(unsigned int i=0 ; i<jets_pt->size() ; ++i)
	// {
	// 	bool isUgly = (jets_isUgly->at(i));
	// 	if(isUgly)  continue;
	// 	
	// 	bool isBtag = (isTight) ? isTightBtag(jets_flwgt_MV1->at(i)) : isStandardBtag(jets_flwgt_MV1->at(i));
	// 	if(isBtag) return true;
	// }
	// return false;

	for(int i=0 ; i<AntiKt4LCTopoLets_n ; i++)
	{
		if(!isGoodJet(i,JetQuality))  continue;
		bool isBtag = (isTight) ? isTightBtag(AntiKt4LCTopoLets_flavor_weight_MV1->at(i)) : isStandardBtag(AntiKt4LCTopoLets_flavor_weight_MV1->at(i));
		if(isBtag) return true;
	}
	return false;
}
int bestOverlapingJet(TLorentzVector lvX)
{
	TLorentzVector lvJet;
	float dRmin = +1.e20;
	int   idRmin = -1;
	
	// for(unsigned int i=0 ; i<jets_pt->size() ; ++i)
	// {
	// 	bool isUgly = (jets_isUgly->at(i));
	// 	if(isUgly)  continue;
	// 	
	// 	lvJet.SetPtEtaPhiE(jets_pt->at(i),jets_eta->at(i),jets_phi->at(i),jets_e->at(i));
	// 	float dR = lvJet.DeltaR(lvX);
	// 	if(dR<dRmin) { dRmin=dR; idRmin=i; }
	// }
	// return idRmin;
	
	for(int i=0 ; i<AntiKt4LCTopoLets_n ; i++)
	{
		if(!isGoodJet(i,JetQuality))  continue;
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(i),AntiKt4LCTopoLets_eta->at(i),AntiKt4LCTopoLets_phi->at(i),AntiKt4LCTopoLets_E->at(i));
		float dR = lvJet.DeltaR(lvX);
		if(dR<dRmin) { dRmin=dR; idRmin=i; }
	}
	return idRmin;
}



void setLumis()
{	
	periodlumi.insert(make_pair("A",0.794017));
	periodlumi.insert(make_pair("B",5.09468));
	periodlumi.insert(make_pair("C",1.40602));
	periodlumi.insert(make_pair("D",3.28839));
	periodlumi.insert(make_pair("E",2.52628));
	periodlumi.insert(make_pair("G",1.27481));
	periodlumi.insert(make_pair("H",1.44493));
	periodlumi.insert(make_pair("I",1.01626));
	periodlumi.insert(make_pair("J",2.59634));
	periodlumi.insert(make_pair("L",0.839765));
	// periodlumi.insert(make_pair("M",0));
}
void enablePeriod(TString period)
{
	if(periodlumi.size()==0) _FATAL("Before you can call enablePeriod() you need to call setLumis()");
	bool isOK        = (periodlumi.find(period)!=periodlumi.end());
	bool isAlreadyIn = (periodenable.find(period)!=periodenable.end());
	if(isAlreadyIn) _FATAL("period name: "+(string)period+" was already enabled");
	if(isOK) periodenable.insert(make_pair(period,true));
	else _FATAL("period name: "+(string)period+" not found in periodlumi map");
}
double getTotalLumi()
{
	if(periodenable.size()==0) _FATAL("Before you can call getTotalLumi() you need to call enablePeriod() at least once to enable at least one period");
	double totallumi = 0;
	for(TMapTSb::iterator it=periodenable.begin() ; it!=periodenable.end() ; ++it)
	{
		TString period = it->first;
		totallumi += periodlumi[period];
	}
	return totallumi;
}
void enableBinnedMC(TString binnedmcchannel)
{
	bool isAlreadyIn = (binnedmcenable.find(binnedmcchannel)!=binnedmcenable.end());
	if(isAlreadyIn) _FATAL("binned mc channel name: "+(string)binnedmcchannel+" was already enabled");
	binnedmcenable.insert(make_pair(binnedmcchannel,true));
}

void properties(TString channel,
				TString   label,  Color_t   color,  Int_t    pattern,  TString   legoption,  Int_t isdraw,
				TMapTSTS& labels, TMapTSTC& colors, TMapTSI& patterns, TMapTSTS& legoptions)
{
	labels.insert(make_pair(channel,     label));
	colors.insert(make_pair(channel,     color));
	patterns.insert(make_pair(channel,   pattern));
	legoptions.insert(make_pair(channel, legoption));
	drawchannels.insert(make_pair(channel, isdraw));
}
bool isWsignal(TString name)
{
	if(name=="Wtaunu_3mu") return true;
	return false;
}
bool isSignal(TString name)
{
	if(name.Contains("3mu")) return true;
	return false;
}
bool isData(TString name)
{
	if(name.Contains("Data") || name.Contains("period")) return true;
	return false;
}
bool isBGsum(TString name)
{
	if(name.Contains("Backgrounds")) return true;
	return false;
}
bool isSIGsum(TString name)
{
	if(name.Contains("Signals")) return true;
	return false;
}
bool isBinned(TString name)
{
	if(name.Contains("uNp")) { _INFO("Binned sample! "+(string)name); return true; }
	if(name.Contains("JZ"))  { _INFO("Binned sample! "+(string)name); return true; }
	// if(name.Contains("period")) { _INFO("Binned sample! "+(string)name); return true; }
	if(isData(name))         { _INFO("Binned sample! "+(string)name); return true; }
	// if(name.Contains("...")) { _INFO("Binned sample! "+(string)name); return true; }
	// if(name.Contains("...")) { _INFO("Binned sample! "+(string)name); return true; }
	return false;
}
TString getUnbinnedName(TString name)
{
	TString bindex = name;	
	if(name.Contains("WmunuNp"))
	{
		bindex.ReplaceAll("WmunuNp","");
		name.ReplaceAll(bindex,"X");
	}
	else if(name.Contains("ZmumuNp"))
	{
		bindex.ReplaceAll("ZmumuNp","");
		name.ReplaceAll(bindex,"X");
	}
	else if(name.Contains("WtaunuNp"))
	{
		bindex.ReplaceAll("WtaunuNp","");
		name.ReplaceAll(bindex,"X");
	}
	else if(name.Contains("JZ"))
	{
		bindex.ReplaceAll("JZ","");
		bindex.ReplaceAll("W", "");
		name.ReplaceAll(bindex,"x");
	}
	else if(name.Contains("period"))
	{
		name = "Data";
	}
	else
	{
		_FATAL("calling getUnbinnedName on an unbinned sample is forbidden.");
	}
	_INFO("running with name="+(string)name);
	return name;
}



void addCounter(TString name, int order)
{
	counters.insert(make_pair(name,0));
	counters_weighted.insert(make_pair(name,0.));
	counters_evtvisited.insert(make_pair(name,0));
	counters_ordered.insert(make_pair(order,name));
}
void initCounters()
{
	//////////////////////////////////////////
	//// cutflow counters (have order>=0) ////
	//////////////////////////////////////////
	
	//// event preselection
	addCounter("nPassing_evt_all",          0);
	addCounter("nPassing_evt_trigger",      1);
	addCounter("nPassing_evt_goodtriplets", 2);
	
	//// 2nd skim level cuts - may be redundant ?
	if(skim) addCounter("nPassing_skim2_m3mu",       10);
	if(skim) addCounter("nPassing_skim2_pT3mu",      11);
	// if(skim) addCounter("nPassing_skim2_mcp",     12);
	if(skim) addCounter("nPassing_skim2_met",        13);
	if(skim) addCounter("nPassing_skim2_mt",         14);
	
	//// single muon
	addCounter("nPassing_mu_mcp",        20);
	addCounter("nPassing_mu_pt",         21);
	addCounter("nPassing_mu_eta",        22);
	addCounter("nPassing_mu_sctangsig",  23);
	addCounter("nPassing_mu_pbalsig",    24);
	// addCounter("nPassing_mu_type",    25);
	addCounter("nPassing_mu_trkquality", 26);
	
	//// triplet cuts
	addCounter("nPassing_triplet_pvalueLoose", 50);
	addCounter("nPassing_triplet_m",           51);
	addCounter("nPassing_triplet_ptLoose",     52);
	addCounter("nPassing_triplet_charge",      53);
	addCounter("nPassing_triplet_mOSmSS",      54);
	// addCounter("nPassing_triplet_ptfrac",   55);
	
	//// Hadronic cleaning (W-oriented) cuts
	addCounter("nPassing_jets_bjetveto_std",    70);
	addCounter("nPassing_jets_coljetvetoLoose", 71);
	addCounter("nPassing_jets_dijetvetoLoose",  72);
	addCounter("nPassing_met_metLoose",         73);
	
	//// one triplet per event (easier for MVA)
	if(!skim && !doMVA) addCounter("nPassing_evt_onetriplet", 80);
	
	//// cuts mode ?
	if(!doMVA) addCounter("nPassing_jets_coljetveto",   100);
	if(!doMVA) addCounter("nPassing_jets_dijetveto",    101);
	if(!doMVA) addCounter("nPassing_triplet_pt",        102);
	if(!doMVA) addCounter("nPassing_triplet_pvalue",    103);
	if(!doMVA) addCounter("nPassing_triplet_vtxclean",  104);
	if(!doMVA) addCounter("nPassing_triplet_isolation", 105);
	if(!doMVA) addCounter("nPassing_W_MET",             106);
	if(!doMVA) addCounter("nPassing_W_dphi3muMET",      107);
	if(!doMVA) addCounter("nPassing_W_mT",              108);

	//// MVA mode ?
	if(doMVA) addCounter("nPassing_MVA_vtxmet",          200);

	addCounter("nPassing_evt_sidebands",                 300);


	/////////////////////////////////////
	/// other counters (have order<0) ///
	/////////////////////////////////////
	addCounter("nPassing_trusigtriplet_matchedmuons",     -102);
	addCounter("nPassing_trusigtriplet_atleast3recmuon",  -103);
	addCounter("nPassing_trusigtriplet_acc",              -104);
	addCounter("nPassing_trusigtriplet_all",              -105);
}
bool isCounter(TString name)
{
	return (counters.find(name)!=counters.end());
}
void clearCounters()
{
	counters.clear();
	counters_weighted.clear();
	counters_evtvisited.clear();
	counters_ordered.clear();
}
void resetCounterFlags()
{
	for(TMapTSui::iterator it=counters_evtvisited.begin() ; it!=counters_evtvisited.end() ; it++)
	{
		it->second = 0;
	}
}
void incrementCounter(TString name, double weight=1.)
{
	if(!isCounter(name)) return;
	
	if(!counters_evtvisited[name])
	{
		counters[name]++;
		counters_weighted[name]+=weight;
		counters_evtvisited[name] = 1;
	}
}
int getCounter(TString name, bool doWeighted=false)
{
	if(!isCounter(name)) return -9999;
	
	if(doWeighted) return counters_weighted[name];
	return counters[name];
}
string printCounters(TString type, TString pname, bool doWeighted=false)
{
	string allLines = "";
	int maxlength = 0;
	int maxcounterlength = 0;
	for(TMapTSui::iterator it = counters.begin() ; it!=counters.end() ; it++)
	{
		TString name = it->first;
		int icounter = counters[name];
		double dcounter = counters_weighted[name];
		TString counter = (doWeighted) ? (TString)_s(dcounter,3) : (TString)_s(icounter);
		maxlength        = (name.Length()>maxlength)           ? name.Length()    : maxlength;
		maxcounterlength = (counter.Length()>maxcounterlength) ? counter.Length() : maxcounterlength;
	}
	
	string header = "";
	if     (type=="cutflow" && !doWeighted) header = "================== CUTFLOW "+ pname +" =================";
	else if(type=="cutflow" && doWeighted)  header = "================== WEIGHTED-CUTFLOW "+ pname +" =================";
	else if(type=="objects" && !doWeighted) header = "================== OBJECTS "+pname+" =================";
	else if(type=="objects" && doWeighted)  header = "================== WEIGHTED-OBJECTS "+pname+" =================";
	else                     _ERROR("Unsupported option: "+(string)type+" for printCounters(TString,TString,bool) method");
	cout << header << endl;
	allLines += header+"\n";
	
	
	TString previousname = "";
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; it++)
	{
		TString name = it->second;
		int order    = it->first;
		
		int icounter = counters[name];
		double dcounter = counters_weighted[name];
		
		TString thiscounter = (doWeighted) ? (TString)_s(dcounter,3) : (TString)_s(icounter);
		int thislength = thiscounter.Length();
		int namelength = it->second.Length();
		if(type=="cutflow" && order<0)  continue;
		if(type=="objects" && order>=0) continue;
		
		string spaces = "";
		for(int j=1 ; j<(maxlength-namelength+1) ; j++) spaces+=" ";
		string digits = "";
		for(int j=1 ; j<(maxcounterlength-thislength+1) ; j++) digits+=" ";
		
		if(previousname=="") previousname = name; // only for the 1st cut
		double thiscut = (doWeighted) ? dcounter : icounter; // counters[name];
		double prevcut = (doWeighted) ? counters_weighted[previousname] : counters[previousname];
		double reldiff = 100*(prevcut-thiscut)/prevcut;
		
		string sthiscut = (doWeighted) ? _s(thiscut,3) : _s(thiscut,0);
		string cutline = (string)name+" "+spaces+digits+sthiscut+" [-"+_s(reldiff,1)+"\%]";
		cout << cutline << endl;
		allLines += cutline+"\n";
		
		previousname = name;
	}
	string footer = "";
	if(!doWeighted)  footer = "================== CUTFLOW "+pname+" =================";
	else             footer = "=========================== "+pname+" ========================";
	cout << footer << endl;
	allLines += footer+"\n";
	
	return allLines;
}
string getTime()
{
	time_t rawtime;
	// struct tm* timeinfo;
	time(&rawtime);
	string timenow = "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n<<<<< Current local time and date: "+(string)ctime(&rawtime);
	return timenow;
}
void writeCoutners(string fname, string lines)
{	
	ofstream ofs;
	ofs.open (fname.c_str(), ofstream::out | ofstream::app);
	ofs << lines;
	ofs.close();
}



void addHist(TMapTSP2TH1& histos, TString channel, TString name, TString titles, TMapTSTS& labels,
			 int nbins, float xmin, float xmax, bool addchannelName=true)
{
	TString hname = channel+"_"+name;
	histos.insert(make_pair(hname, new TH1F(channel+"_"+name,titles,nbins,xmin,xmax)));
	TString title = histos[hname]->GetTitle();
	title = (title!="" && addchannelName) ? labels[channel]+" "+titles : titles;
	histos[hname]->SetTitle(title);
}
void addHist(TMapTSP2TH2& histos, TString channel, TString name, TString titles, TMapTSTS& labels,
			 int nbinsx, float xmin, float xmax,
			 int nbinsy, float ymin, float ymax, bool addchannelName=true)
{
	TString hname = channel+"_"+name;
	histos.insert(make_pair(hname, new TH2F(channel+"_"+name,titles,nbinsx,xmin,xmax,nbinsy,ymin,ymax)));
	TString title = histos[hname]->GetTitle();
	title = (title!="" && addchannelName) ? labels[channel]+" "+titles : titles;
	histos[hname]->SetTitle(title);
}
void makeCutflowHisto(TString channel, TMapTSP2TH1& histos, TMapTSTS& labels, bool addchannelName=true)
{
	vector<TString> cutnames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it)
	{
		if(it->first>=0) cutnames.push_back(it->second);
	}
	Int_t ncuts = cutnames.size();
	
	TString label = (skim) ? "after 1st skim" : "after 2nd skim";
	
	addHist(histos,channel,"cutflow_normalized","Normalized cutflow "+label+";;Normalized",labels,ncuts,0,ncuts,false);
	addHist(histos,channel,"cutflow_absolute","Absolute cutflow "+label+";;Events",labels,ncuts,0,ncuts,addchannelName);
	addHist(histos,channel,"cutflow_weighted","Weighted cutflow "+label+";;Events",labels,ncuts,0,ncuts,addchannelName);
	
	for(Int_t b=1 ; b<=ncuts ; b++)
	{
		TString cutname = cutnames[b-1];
		cutname.ReplaceAll("nPassing_","");
		histos[channel+"_cutflow_normalized"]->GetXaxis()->SetBinLabel(b,cutname);
		histos[channel+"_cutflow_absolute"]->GetXaxis()->SetBinLabel(b,cutname);
		histos[channel+"_cutflow_weighted"]->GetXaxis()->SetBinLabel(b,cutname);
	}
}
void makeCountersHisto(TString channel, TMapTSP2TH1& histos, TMapTSTS& labels, bool addchannelName=true)
{
	vector<TString> counternames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it)
	{
		if(it->first<0) counternames.push_back(it->second);
	}
	Int_t ncounters = counternames.size();
	
	TString label = (skim) ? "after 1st skim" : "after 2nd skim";
	
	addHist(histos,channel,"counters_normalized","Normalized counters "+label+";;Normalized",labels,ncounters,0,ncounters,false);
	addHist(histos,channel,"counters_absolute","Absolute counters "+label+";;Events",labels,ncounters,0,ncounters,addchannelName);
	addHist(histos,channel,"counters_weighted","Weighted counters "+label+";;Events",labels,ncounters,0,ncounters,addchannelName);
	
	for(Int_t b=1 ; b<=ncounters ; b++)
	{
		TString countername = counternames[b-1];
		countername.ReplaceAll("nPassing_","");
		histos[channel+"_counters_normalized"]->GetXaxis()->SetBinLabel(b,countername);
		histos[channel+"_counters_absolute"]->GetXaxis()->SetBinLabel(b,countername);
		histos[channel+"_counters_weighted"]->GetXaxis()->SetBinLabel(b,countername);
	}
}
void makeCategoriesHisto(TString channel, TMapTSP2TH1& histos, TMapTSTS& labels, bool addchannelName=false)
{
	vector<TString> categories;
	// categories.push_back("All");
	// categories.push_back("Triggered");
	// categories.push_back("in acceptance");
	categories.push_back("3mu");
	categories.push_back("2mu+1tpA");
	categories.push_back("2mu+1tpB | 0tpA");
	categories.push_back("1mu+2tpA");
	categories.push_back("1mu+2tpB | 0tpA");
	categories.push_back("1mu+1tpA+1tpB");
	categories.push_back("2mu+1calo(|#eta|<0.1) | 0tpA");
	categories.push_back("3tpA | 0mu");
	categories.push_back("3tpB | 0mu | 0tpA");
	categories.push_back("2tpA+1tpB | 0mu");
	categories.push_back("1tpA+1tpB | 0mu");
	categories.push_back("2tpA+1calo | 0mu");
	unsigned int ncat = categories.size();
	TString hname = "";
	TString fullname = "";
	
	_DEBUG("");
	
	hname = "tripletCategories_noVertexing";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories no vertexing;;Normalized to 1st bin", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	hname = "tripletCategories_norm_noVertexing";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories no vertexing;;Normalized", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	_DEBUG("");
	
	hname = "tripletCategories_afterVertexing";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after vertexing;;Normalized to 1st bin", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	hname = "tripletCategories_norm_afterVertexing";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after vertexing;;Normalized", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);

	_DEBUG("");
	
	hname = "tripletCategories_after_muons";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after muon cuts;;Normalized to 1st bin", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	hname = "tripletCategories_norm_after_muons";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after muon cuts;;Normalized", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	_DEBUG("");
	
	hname = "tripletCategories_after_triplet";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after triplet cuts;;Normalized to 1st bin", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	hname = "tripletCategories_norm_after_triplet";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after triplet cuts;;Normalized", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	_DEBUG("");
	
	hname = "tripletCategories_after_hadclean";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after hadronic cleaning;;Normalized to 1st bin", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	hname = "tripletCategories_norm_after_hadclean";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories after hadronic cleaning;;Normalized", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	_DEBUG("");
	
	hname = "tripletCategories_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories end of selection;;Normalized to 1st bin", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	hname = "tripletCategories_norm_endOfSelection";
	fullname = channel+"_"+hname;
	addHist(histos,channel,hname, "Triplet categories end of selection;;Normalized", labels, ncat,0,ncat, addchannelName);
	for(unsigned int i=1 ; i<=ncat ; ++i)  histos[fullname]->GetXaxis()->SetBinLabel(i,categories[i-1]);
	
	_DEBUG("");
}
void fillCutFlowHisto(TString channel, TMapTSP2TH1& histos)
{	
	vector<TString> cutnames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it) { if(it->first>=0) cutnames.push_back(it->second); }
	Int_t ncuts = cutnames.size();
	for(Int_t b=1 ; b<=ncuts ; b++)
	{
		unsigned int counter  = getCounter(cutnames[b-1]);
		unsigned int wcounter = getCounter(cutnames[b-1],true);
		counter  = (std::isinf(counter)  || std::isnan(counter)  || counter<1)  ? 1.e-3 : counter;
		wcounter = (std::isinf(wcounter) || std::isnan(wcounter) || wcounter<1) ? 1.e-3 : wcounter;
		Double_t bincontent = histos[channel+"_cutflow_absolute"]->GetBinContent(b); // that is nonzero only for channel="Data"
		histos[channel+"_cutflow_normalized"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_cutflow_absolute"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_cutflow_weighted"]->SetBinContent(b,wcounter+bincontent);
	}
}
void fillCountersHisto(TString channel, TMapTSP2TH1& histos)
{	
	vector<TString> counternames;
	for(TMapiTS::iterator it=counters_ordered.begin() ; it!=counters_ordered.end() ; ++it) { if(it->first<0) counternames.push_back(it->second); }
	Int_t ncounters = counternames.size();
	for(Int_t b=1 ; b<=ncounters ; b++)
	{
		unsigned int counter  = getCounter(counternames[b-1]);
		unsigned int wcounter = getCounter(counternames[b-1],true);
		counter  = (std::isinf(counter)  || std::isnan(counter)  || counter<1)  ? 1.e-3 : counter;
		wcounter = (std::isinf(wcounter) || std::isnan(wcounter) || wcounter<1) ? 1.e-3 : wcounter;
		Double_t bincontent = histos[channel+"_counters_absolute"]->GetBinContent(b); // that is nonzero only for channel="Data"
		histos[channel+"_counters_normalized"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_counters_absolute"]->SetBinContent(b,counter+bincontent);
		histos[channel+"_counters_weighted"]->SetBinContent(b,wcounter+bincontent);
	}
}




int vetoHF(TString name)
{
	//////////////////////////////////////
	// only for JZxW /////////////////////
	if(!name.Contains("JZ")) return 0; ///
	//////////////////////////////////////
	
	for(unsigned int i=0 ; i<mc_pdgId->size() ; i++)
	{
		if( isBhadron((unsigned int)mc_pdgId->at(i)) ) return (int)mc_pdgId->at(i);
	}
	return 0;
}
double fFONLLbbTomu15(double pt) // pt in GeV !!!
{
	double a0 = 8.30678e-01;
	double a1 = 4.03935e+02;
	double a2 = 1.19803e-01;
	double a3 = 9.88432e-04;
	double a4 = 1.76653e-02;
	return a0*TMath::Power(pt+a1,a2)-TMath::Power(pt,a3)/TMath::Exp(a4*pt);
}
double fFONLLccTomu15(double pt) // pt in GeV !!!
{
	double a0 = 1.36119e-01;
	double a1 = 4.32113e+03;
	double a2 = 2.36906e-01;
	double a3 = -1.94770e-07;
	double a4 = 5.75564e-02;
	return a0*TMath::Power(pt+a1,a2)-a3*TMath::Power(pt,4)/TMath::Exp(a4*pt);
}
double getLeadingPtGeVTruthMu()
{
	double ptmax = -1.e20;
	for(unsigned int imc=0 ; imc<mc_pdgId->size() ; imc++)
	{				
		if(fabs(mc_pdgId->at(imc))!=PDTMU) continue;
		if(mc_status->at(imc)!=1)          continue;
		if(fabs(mc_eta->at(imc))>2.5)      continue;
		
		double pt = mc_pt->at(imc);
		ptmax = (pt>ptmax) ? pt : ptmax;
	}
	return ptmax;
}
double getFONLLShapeWeight(TString name)
{
	if(!name.Contains("ccTomu15") && !name.Contains("bbTomu15")) _FATAL("FONLL is to be used either with ccTomu15 or bbTomu15 samples !");
	double pt = getLeadingPtGeVTruthMu()*MeV2GeV;
	double weight = -1.;
	if(name.Contains("ccTomu15")) weight = fFONLLccTomu15(pt);
	if(name.Contains("bbTomu15")) weight = fFONLLbbTomu15(pt);	
	return weight;
}
double getFONLLFlatWeight(TString name)
{
	if(!name.Contains("ccTomu15") && !name.Contains("bbTomu15"))
	{
		_FATAL("FONLL is to be used either with ccTomu15 or bbTomu15 samples !");
	}
	return 1.;
}
double getHFnormalization(TMapuiTS& channels, TMapTSP2TH1& histos, double xmin, double xmax, TString varname)
{
	double EWX     = 0.;
	double HF      = 0.;
	double errEWX  = 0.;
	double errHF   = 0.;
	double data    = 0.;
	double errdata = 0.;
	for(TMapuiTS::iterator it=channels.begin() ; it!=channels.end() ; it++)
	{
		TString name = it->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if     (isSignal(name)) continue; // skip the signal
		else if(isBGsum(name))  continue; // skip the bkgs sum
		else if(isSIGsum(name)) continue; // skip the signals sum
		else if(isData(name))
		{
			data    = getEntries(histos[hname], xmin,xmax);
			errdata = getErrors(histos[hname], xmin,xmax);
		}
		else if(name.Contains("bbTomu15") || name.Contains("ccTomu15"))
		{
			HF += getEntries(histos[hname], xmin,xmax);
			errHF += getErrors(histos[hname], xmin,xmax);
		}
		else
		{
			EWX += getEntries(histos[hname], xmin,xmax);
			errEWX += getErrors(histos[hname], xmin,xmax);
		}
	}
	
	// Pythia's HF true normalization is incorrect.
	// Normalize it to the data and therefore get HF'.
	// Work in the 1st few bins of the pT rec histo,
	// where the HF component is dominant over the
	// rest of the contributing channels to the 1mu pT.
	// Or, do the same on the J/psi peak.
	// Reqire: data/(HF'+EWX) = 1 (in the 1st few pT bins or on the J/psi peak)
	// ==>     data = HF'+EWX
	// ==>     HF' = data-EWX
	// ==>     normaliztion = HF'/HF = (data-EWX)/HF
	
	double normaliztion = (data-EWX)/HF;
	double normaliztionErr = normaliztion*TMath::Sqrt(TMath::Power((errdata+errEWX)/(data-EWX),2) + TMath::Power(errHF/HF,2));
	cout << "\n[hist=" << varname << "]: HF Normalization = " << normaliztion << " +- " << normaliztionErr << endl;
	cout << "Now that you have this number, check that it is correct in : getFONLLFlatWeight(TString name)" << endl;
	cout << "and run this process again if it is different than the hardcoded number there." << endl;
	return normaliztion;
}
double getKfactorWeight(TString name)
{
	if(!name.Contains("WmunuNp") && !name.Contains("ZmumuNp"))
	{
		_FATAL("K factor weight can be used either with WmunuNpX or with ZmumuNpX !");
	}
	return (name.Contains("WmunuNp")) ? 1.19 : 1.23;
}
double getDijetWeight(TString name)
{
	if(!name.Contains("JZ"))
	{
		_FATAL("dijet weight can be used only with JZxW, with x=(0),1,2,3");
	}
	return mc_event_weight;
}

double getSampleWeight(TString name)
{
	double sigma,Nevents,genEff;
	double sigmaSM=-1.;
	if     (name.Contains("Wtaunu_3mu"))    { sigma=9.5753E+00*nb2fb; Nevents=100000.;   genEff=1.0000E+00; sigmaSM=9.1081E+00*nb2fb; }
	else if(name.Contains("bbTotau10_3mu")) { sigma=3.5388E+02*nb2fb; Nevents=100000.;   genEff=1.0000E+00; sigmaSM=3.5388E+02*nb2fb; }
	else if(name.Contains("ccTotau10_3mu")) { sigma=2.9674E+02*nb2fb; Nevents=100000.;   genEff=1.0000E+00; sigmaSM=2.9674E+02*nb2fb; }
	else if(name.Contains("bb_mu4mu4"))     { sigma=1.1464E+02*nb2fb; Nevents=19980978.; genEff=1.0000E+00;}
	else if(name.Contains("bbTomu15"))      { sigma=1.9898E+02*nb2fb; Nevents=4998088.;  genEff=1.0000E+00;}
	else if(name.Contains("ccTomu15"))      { sigma=8.0088E+01*nb2fb; Nevents=4998690.;  genEff=1.0000E+00;}
	else if(name.Contains("bb_Jpsimu4mu4")) { sigma=2.0874E+02*nb2fb; Nevents=9969994.;  genEff=1.0000E+00;} // cross section is wrong !
	else if(name.Contains("WmunuNp0"))      { sigma=8.0463E+00*nb2fb; Nevents=3469591.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp1"))      { sigma=1.5812E+00*nb2fb; Nevents=2499893.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp2"))      { sigma=4.7741E-01*nb2fb; Nevents=3769890.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp3"))      { sigma=1.3391E-01*nb2fb; Nevents=1009896.;  genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp4"))      { sigma=3.5982E-02*nb2fb; Nevents=255000.;   genEff=1.0000E+00;}
	else if(name.Contains("WmunuNp5"))      { sigma=1.0409E-02*nb2fb; Nevents=20000.;    genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp0"))      { sigma=7.1211E-01*nb2fb; Nevents=6609982.;  genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp1"))      { sigma=1.5477E-01*nb2fb; Nevents=1334897.;  genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp2"))      { sigma=4.8912E-02*nb2fb; Nevents=404897.;   genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp3"))      { sigma=1.4226E-02*nb2fb; Nevents=110000.;   genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp4"))      { sigma=3.7838E-03*nb2fb; Nevents=29999.;    genEff=1.0000E+00;}
	else if(name.Contains("ZmumuNp5"))      { sigma=1.1148E-03*nb2fb; Nevents=10000.;    genEff=1.0000E+00;}
	else if(name.Contains("JZ0W"))          { sigma=7.2850E+07*nb2fb; Nevents=3998693.;  genEff=4.2413E-04;}
	else if(name.Contains("JZ1W"))          { sigma=4.1440E+06*nb2fb; Nevents=1999694.;  genEff=3.4217E-05;}
	else if(name.Contains("JZ2W"))          { sigma=5.0147E+03*nb2fb; Nevents=2499692.;  genEff=7.0769E-04;}
	else if(name.Contains("JZ3W"))          { sigma=5.4418E+02*nb2fb; Nevents=998688.;   genEff=6.7835E-05;}
	
	else _FATAL("sample name: "+(string)name+" is not supported");
	double Lmc = (sigmaSM>0.) ? (Nevents/genEff)/(sigmaSM*(BRbelle*BRfactor)) : (Nevents/genEff)/sigma;
	// double weight = luminosity/Lmc;
	double weight = totalLumi/Lmc;
	return weight;
}

bool MCP(unsigned int index)
{
	// !expectBLayerHit OR numberOfBLayerHits > 0,  ===> REMOVED !!!
	// number of pixel hits + number of crossed dead pixel sensors > 0,
	// number of SCT hits + number of crossed dead SCT sensors > 4,
	// Number of pixel holes + number of SCT holes < 3,
	// TRT hits > 5 for 0.1<|eta|<1.9, outlier fraction < 0.9.
	// If TRT hits > 5 for |eta|<=0.1 or |eta|>=1.9, outlier fraction < 0.9  ===> REMOVED !!!

	bool passedPIX = ((trks_nPix->at(index)+trks_nDeadPixels->at(index)) > 0);
	if(!passedPIX) return false;

	bool passedSCT = ((trks_nSCT->at(index)+trks_nDeadSCT->at(index)) > 4);
	if(!passedSCT) return false;

	bool passedHOL = ((trks_nPixHoles->at(index)+trks_nSCTHoles->at(index)) < 3);
	if(!passedHOL) return false;

	float TRTratio = 0.9;
	float nTRTmin  = 5;
	float nTRT     = trks_nTRT->at(index)+trks_nTRTOutliers->at(index);
	float nTRTol   = trks_nTRTOutliers->at(index);
	float etamu    = trks_eta->at(index);
	if(fabs(etamu)>0.1 && fabs(etamu)<1.9)
	{
		if(nTRT>nTRTmin)
		{
			if(nTRTol>=TRTratio*nTRT) return false;
		}
	}

	//// otherwise return true
	return true;
}



void setBranches(TString tType, TChain* t)
{
	if(tType=="physics")
	{
		AntiKt4LCTopoLets_E = 0;
		AntiKt4LCTopoLets_pt = 0;
		AntiKt4LCTopoLets_m = 0;
		AntiKt4LCTopoLets_eta = 0;
		AntiKt4LCTopoLets_phi = 0;
		
		AntiKt4LCTopoLets_emfrac = 0;
		AntiKt4LCTopoLets_hecf = 0;
		AntiKt4LCTopoLets_LArQuality = 0;
		AntiKt4LCTopoLets_HECQuality = 0;
		AntiKt4LCTopoLets_Timing = 0;
		AntiKt4LCTopoLets_sumPtTrk_pv0_500MeV = 0;
		AntiKt4LCTopoLets_eta = 0;
		AntiKt4LCTopoLets_pt = 0;
		AntiKt4LCTopoLets_fracSamplingMax = 0;
		AntiKt4LCTopoLets_NegativeE = 0;
		AntiKt4LCTopoLets_AverageLArQF = 0;
		
		AntiKt4LCTopoLets_EtaOrigin = 0;
		AntiKt4LCTopoLets_PhiOrigin = 0;
		AntiKt4LCTopoLets_MOrigin = 0;
		AntiKt4LCTopoLets_WIDTH = 0;
		AntiKt4LCTopoLets_n90 = 0;
		AntiKt4LCTopoLets_nTrk = 0;
		AntiKt4LCTopoLets_sumPtTrk = 0;
		AntiKt4LCTopoLets_OriginIndex = 0;
        
		AntiKt4LCTopoLets_BCH_CORR_CELL = 0;
		AntiKt4LCTopoLets_BCH_CORR_DOTX = 0;
		AntiKt4LCTopoLets_BCH_CORR_JET = 0;
		AntiKt4LCTopoLets_BCH_CORR_JET_FORCELL = 0;
		AntiKt4LCTopoLets_ENG_BAD_CELLS = 0;
		AntiKt4LCTopoLets_N_BAD_CELLS = 0;
		AntiKt4LCTopoLets_N_BAD_CELLS_CORR = 0;
		AntiKt4LCTopoLets_BAD_CELLS_CORR_E = 0;
        
		AntiKt4LCTopoLets_isUgly = 0;
		AntiKt4LCTopoLets_isBadLooseMinus = 0;
		AntiKt4LCTopoLets_isBadLoose = 0;
		AntiKt4LCTopoLets_isBadMedium = 0;
		AntiKt4LCTopoLets_isBadTight = 0;
        
		AntiKt4LCTopoLets_Offset = 0;
		AntiKt4LCTopoLets_EMJES = 0;
		AntiKt4LCTopoLets_EMJES_EtaCorr = 0;
		AntiKt4LCTopoLets_EMJESnooffset = 0;
        
		AntiKt4LCTopoLets_emscale_E = 0;
		AntiKt4LCTopoLets_emscale_pt = 0;
		AntiKt4LCTopoLets_emscale_m = 0;
		AntiKt4LCTopoLets_emscale_eta = 0;
		AntiKt4LCTopoLets_emscale_phi = 0;
		
		AntiKt4LCTopoLets_jvtx_x = 0;
		AntiKt4LCTopoLets_jvtx_y = 0;
		AntiKt4LCTopoLets_jvtx_z = 0;
		AntiKt4LCTopoLets_jvtxf = 0;
		
		AntiKt4LCTopoLets_LikeLihood_0 = 0;
		AntiKt4LCTopoLets_ActiveArea = 0;
		AntiKt4LCTopoLets_ActiveAreaPx = 0;
		AntiKt4LCTopoLets_ActiveAreaPy = 0;
		AntiKt4LCTopoLets_ActiveAreaPz = 0;
		AntiKt4LCTopoLets_ActiveAreaE = 0;
		AntiKt4LCTopoLets_VoronoiArea = 0;
		AntiKt4LCTopoLets_VoronoiAreaPx = 0;
		AntiKt4LCTopoLets_VoronoiAreaPy = 0;
		AntiKt4LCTopoLets_VoronoiAreaPz = 0;
		AntiKt4LCTopoLets_VoronoiAreaE = 0;
		AntiKt4LCTopoLets_LowEtConstituentsFrac = 0;
		
		AntiKt4LCTopoLets_constscale_E = 0;
		AntiKt4LCTopoLets_constscale_pt = 0;
		AntiKt4LCTopoLets_constscale_m = 0;
		AntiKt4LCTopoLets_constscale_eta = 0;
		AntiKt4LCTopoLets_constscale_phi = 0;
		AntiKt4LCTopoLets_LArBadHVEnergy = 0;
		AntiKt4LCTopoLets_LArBadHVRatio = 0;
		
		AntiKt4LCTopoLets_flavor_weight_MV1 = 0;
		
		trk_pt = 0;
		trk_eta = 0;
		trk_phi_wrtPV = 0;
        
		mu_muons_E = 0;
		mu_muons_pt = 0;
		mu_muons_m = 0;
		mu_muons_eta = 0;
		mu_muons_phi = 0;
		mu_muons_charge = 0;
		mu_muons_ms_phi = 0;
		mu_muons_ms_theta = 0;
		mu_muons_ms_qoverp = 0;
        
		mu_staco_E = 0;
		mu_staco_pt = 0;
		mu_staco_m = 0;
		mu_staco_eta = 0;
		mu_staco_phi = 0;
		mu_staco_charge = 0;
		mu_staco_ms_phi = 0;
		mu_staco_ms_theta = 0;
		mu_staco_ms_qoverp = 0;
        
		mu_muid_E = 0;
		mu_muid_pt = 0;
		mu_muid_m = 0;
		mu_muid_eta = 0;
		mu_muid_phi = 0;
		mu_muid_charge = 0;
		mu_muid_ms_phi = 0;
		mu_muid_ms_theta = 0;
		mu_muid_ms_qoverp = 0;
        
		el_E = 0;
		el_Et = 0;
		el_pt = 0;
		el_m = 0;
		el_eta = 0;
		el_phi = 0;
		
		ph_E = 0;
		ph_Et = 0;
		ph_pt = 0;
		ph_m = 0;
		ph_eta = 0;
		ph_phi = 0;
		
		tau_Et = 0;
		tau_pt = 0;
		tau_m = 0;
		tau_eta = 0;
		tau_phi = 0;
        
		el_MET_wpx = 0;
		el_MET_wpy = 0;
		el_MET_wet = 0;
		el_MET_statusWord = 0;
		ph_MET_wpx = 0;
		ph_MET_wpy = 0;
		ph_MET_wet = 0;
		ph_MET_statusWord = 0;
		mu_staco_MET_wpx = 0;
		mu_staco_MET_wpy = 0;
		mu_staco_MET_wet = 0;
		mu_staco_MET_statusWord = 0;
		mu_muid_MET_wpx = 0;
		mu_muid_MET_wpy = 0;
		mu_muid_MET_wet = 0;
		mu_muid_MET_statusWord = 0;
		mu_MET_wpx = 0;
		mu_MET_wpy = 0;
		mu_MET_wet = 0;
		mu_MET_statusWord = 0;
		tau_MET_wpx = 0;
		tau_MET_wpy = 0;
		tau_MET_wet = 0;
		tau_MET_statusWord = 0;
		jet_antikt4LCtopo_MET_wpx = 0;
		jet_antikt4LCtopo_MET_wpy = 0;
		jet_antikt4LCtopo_MET_wet = 0;
		jet_antikt4LCtopo_MET_statusWord = 0;
		cl_MET_wpx = 0;
		cl_MET_wpy = 0;
		cl_MET_wet = 0;
		cl_MET_statusWord = 0;
		trk_MET_wpx = 0;
		trk_MET_wpy = 0;
		trk_MET_wet = 0;
		trk_MET_statusWord = 0;
		
		musp_eta = 0;
		musp_phi = 0;
		musp_innerSegments = 0;
		musp_middleSegments = 0;
		musp_outerSegments = 0;
		
		vxp_nTracks = 0;
		
		
		
		t->SetBranchAddress("RunNumber",&phys_RunNumber);
		t->SetBranchAddress("EventNumber",&phys_EventNumber);
		if(glob_isMC) t->SetBranchAddress("mc_channel_number",&phys_mc_channel_number);
		if(glob_isMC) t->SetBranchAddress("mc_event_number",&phys_mc_event_number);
		if(glob_isMC) t->SetBranchAddress("mc_event_weight",&phys_mc_event_weight);
		t->SetBranchAddress("lbn",&phys_lbn);
		t->SetBranchAddress("actualIntPerXing",&phys_actualIntPerXing);
		t->SetBranchAddress("averageIntPerXing",&phys_averageIntPerXing);

		t->SetBranchAddress("AntiKt4LCTopoLets_n", &AntiKt4LCTopoLets_n);
		t->SetBranchAddress("AntiKt4LCTopoLets_E", &AntiKt4LCTopoLets_E);
		t->SetBranchAddress("AntiKt4LCTopoLets_pt", &AntiKt4LCTopoLets_pt);
		t->SetBranchAddress("AntiKt4LCTopoLets_m", &AntiKt4LCTopoLets_m);
		t->SetBranchAddress("AntiKt4LCTopoLets_eta", &AntiKt4LCTopoLets_eta);
		t->SetBranchAddress("AntiKt4LCTopoLets_phi", &AntiKt4LCTopoLets_phi);
		
		t->SetBranchAddress("AntiKt4LCTopoLets_emfrac", &AntiKt4LCTopoLets_emfrac);
		t->SetBranchAddress("AntiKt4LCTopoLets_hecf", &AntiKt4LCTopoLets_hecf);
		t->SetBranchAddress("AntiKt4LCTopoLets_LArQuality", &AntiKt4LCTopoLets_LArQuality);
		t->SetBranchAddress("AntiKt4LCTopoLets_HECQuality", &AntiKt4LCTopoLets_HECQuality);
		t->SetBranchAddress("AntiKt4LCTopoLets_Timing", &AntiKt4LCTopoLets_Timing);
		t->SetBranchAddress("AntiKt4LCTopoLets_sumPtTrk_pv0_500MeV", &AntiKt4LCTopoLets_sumPtTrk_pv0_500MeV);
		t->SetBranchAddress("AntiKt4LCTopoLets_eta", &AntiKt4LCTopoLets_eta);
		t->SetBranchAddress("AntiKt4LCTopoLets_pt", &AntiKt4LCTopoLets_pt);
		t->SetBranchAddress("AntiKt4LCTopoLets_fracSamplingMax", &AntiKt4LCTopoLets_fracSamplingMax);
		t->SetBranchAddress("AntiKt4LCTopoLets_NegativeE", &AntiKt4LCTopoLets_NegativeE);
		t->SetBranchAddress("AntiKt4LCTopoLets_AverageLArQF", &AntiKt4LCTopoLets_AverageLArQF);
		
		t->SetBranchAddress("AntiKt4LCTopoLets_EtaOrigin", &AntiKt4LCTopoLets_EtaOrigin);
		t->SetBranchAddress("AntiKt4LCTopoLets_PhiOrigin", &AntiKt4LCTopoLets_PhiOrigin);
		t->SetBranchAddress("AntiKt4LCTopoLets_MOrigin", &AntiKt4LCTopoLets_MOrigin);
		t->SetBranchAddress("AntiKt4LCTopoLets_WIDTH", &AntiKt4LCTopoLets_WIDTH);
		t->SetBranchAddress("AntiKt4LCTopoLets_n90", &AntiKt4LCTopoLets_n90);
		t->SetBranchAddress("AntiKt4LCTopoLets_nTrk", &AntiKt4LCTopoLets_nTrk);
		t->SetBranchAddress("AntiKt4LCTopoLets_sumPtTrk", &AntiKt4LCTopoLets_sumPtTrk);
		t->SetBranchAddress("AntiKt4LCTopoLets_OriginIndex", &AntiKt4LCTopoLets_OriginIndex);

		t->SetBranchAddress("AntiKt4LCTopoLets_BCH_CORR_CELL", &AntiKt4LCTopoLets_BCH_CORR_CELL);
		t->SetBranchAddress("AntiKt4LCTopoLets_BCH_CORR_DOTX", &AntiKt4LCTopoLets_BCH_CORR_DOTX);
		t->SetBranchAddress("AntiKt4LCTopoLets_BCH_CORR_JET", &AntiKt4LCTopoLets_BCH_CORR_JET);
		t->SetBranchAddress("AntiKt4LCTopoLets_BCH_CORR_JET_FORCELL", &AntiKt4LCTopoLets_BCH_CORR_JET_FORCELL);
		t->SetBranchAddress("AntiKt4LCTopoLets_ENG_BAD_CELLS", &AntiKt4LCTopoLets_ENG_BAD_CELLS);
		t->SetBranchAddress("AntiKt4LCTopoLets_N_BAD_CELLS", &AntiKt4LCTopoLets_N_BAD_CELLS);
		t->SetBranchAddress("AntiKt4LCTopoLets_N_BAD_CELLS_CORR", &AntiKt4LCTopoLets_N_BAD_CELLS_CORR);
		t->SetBranchAddress("AntiKt4LCTopoLets_BAD_CELLS_CORR_E", &AntiKt4LCTopoLets_BAD_CELLS_CORR_E);

		t->SetBranchAddress("AntiKt4LCTopoLets_isUgly", &AntiKt4LCTopoLets_isUgly);
		t->SetBranchAddress("AntiKt4LCTopoLets_isBadLooseMinus", &AntiKt4LCTopoLets_isBadLooseMinus);
		t->SetBranchAddress("AntiKt4LCTopoLets_isBadLoose", &AntiKt4LCTopoLets_isBadLoose);
		t->SetBranchAddress("AntiKt4LCTopoLets_isBadMedium", &AntiKt4LCTopoLets_isBadMedium);
		t->SetBranchAddress("AntiKt4LCTopoLets_isBadTight", &AntiKt4LCTopoLets_isBadTight);

		t->SetBranchAddress("AntiKt4LCTopoLets_Offset", &AntiKt4LCTopoLets_Offset);
		t->SetBranchAddress("AntiKt4LCTopoLets_EMJES", &AntiKt4LCTopoLets_EMJES);
		t->SetBranchAddress("AntiKt4LCTopoLets_EMJES_EtaCorr", &AntiKt4LCTopoLets_EMJES_EtaCorr);
		t->SetBranchAddress("AntiKt4LCTopoLets_EMJESnooffset", &AntiKt4LCTopoLets_EMJESnooffset);

		t->SetBranchAddress("AntiKt4LCTopoLets_emscale_E", &AntiKt4LCTopoLets_emscale_E);
		t->SetBranchAddress("AntiKt4LCTopoLets_emscale_pt", &AntiKt4LCTopoLets_emscale_pt);
		t->SetBranchAddress("AntiKt4LCTopoLets_emscale_m", &AntiKt4LCTopoLets_emscale_m);
		t->SetBranchAddress("AntiKt4LCTopoLets_emscale_eta", &AntiKt4LCTopoLets_emscale_eta);
		t->SetBranchAddress("AntiKt4LCTopoLets_emscale_phi", &AntiKt4LCTopoLets_emscale_phi);
		
		t->SetBranchAddress("AntiKt4LCTopoLets_jvtx_x", &AntiKt4LCTopoLets_jvtx_x);
		t->SetBranchAddress("AntiKt4LCTopoLets_jvtx_y", &AntiKt4LCTopoLets_jvtx_y);
		t->SetBranchAddress("AntiKt4LCTopoLets_jvtx_z", &AntiKt4LCTopoLets_jvtx_z);
		t->SetBranchAddress("AntiKt4LCTopoLets_jvtxf", &AntiKt4LCTopoLets_jvtxf);
		
		t->SetBranchAddress("AntiKt4LCTopoLets_LikeLihood_0", &AntiKt4LCTopoLets_LikeLihood_0);
		t->SetBranchAddress("AntiKt4LCTopoLets_ActiveArea", &AntiKt4LCTopoLets_ActiveArea);
		t->SetBranchAddress("AntiKt4LCTopoLets_ActiveAreaPx", &AntiKt4LCTopoLets_ActiveAreaPx);
		t->SetBranchAddress("AntiKt4LCTopoLets_ActiveAreaPy", &AntiKt4LCTopoLets_ActiveAreaPy);
		t->SetBranchAddress("AntiKt4LCTopoLets_ActiveAreaPz", &AntiKt4LCTopoLets_ActiveAreaPz);
		t->SetBranchAddress("AntiKt4LCTopoLets_ActiveAreaE", &AntiKt4LCTopoLets_ActiveAreaE);
		t->SetBranchAddress("AntiKt4LCTopoLets_VoronoiArea", &AntiKt4LCTopoLets_VoronoiArea);
		t->SetBranchAddress("AntiKt4LCTopoLets_VoronoiAreaPx", &AntiKt4LCTopoLets_VoronoiAreaPx);
		t->SetBranchAddress("AntiKt4LCTopoLets_VoronoiAreaPy", &AntiKt4LCTopoLets_VoronoiAreaPy);
		t->SetBranchAddress("AntiKt4LCTopoLets_VoronoiAreaPz", &AntiKt4LCTopoLets_VoronoiAreaPz);
		t->SetBranchAddress("AntiKt4LCTopoLets_VoronoiAreaE", &AntiKt4LCTopoLets_VoronoiAreaE);
		t->SetBranchAddress("AntiKt4LCTopoLets_LowEtConstituentsFrac", &AntiKt4LCTopoLets_LowEtConstituentsFrac);
		
		t->SetBranchAddress("AntiKt4LCTopoLets_constscale_E", &AntiKt4LCTopoLets_constscale_E);
		t->SetBranchAddress("AntiKt4LCTopoLets_constscale_pt", &AntiKt4LCTopoLets_constscale_pt);
		t->SetBranchAddress("AntiKt4LCTopoLets_constscale_m", &AntiKt4LCTopoLets_constscale_m);
		t->SetBranchAddress("AntiKt4LCTopoLets_constscale_eta", &AntiKt4LCTopoLets_constscale_eta);
		t->SetBranchAddress("AntiKt4LCTopoLets_constscale_phi", &AntiKt4LCTopoLets_constscale_phi);
		t->SetBranchAddress("AntiKt4LCTopoLets_LArBadHVEnergy", &AntiKt4LCTopoLets_LArBadHVEnergy);
		t->SetBranchAddress("AntiKt4LCTopoLets_LArBadHVRatio", &AntiKt4LCTopoLets_LArBadHVRatio);
		
		t->SetBranchAddress("AntiKt4LCTopoLets_flavor_weight_MV1", &AntiKt4LCTopoLets_flavor_weight_MV1);
		
		t->SetBranchAddress("MET_RefFinal_etx", &MET_RefFinal_etx);
		t->SetBranchAddress("MET_RefFinal_ety", &MET_RefFinal_ety);
		t->SetBranchAddress("MET_RefFinal_phi", &MET_RefFinal_phi);
		t->SetBranchAddress("MET_RefFinal_et", &MET_RefFinal_et);
		t->SetBranchAddress("MET_RefFinal_sumet", &MET_RefFinal_sumet);
		
		
		if(glob_isSig)
		{
			t->SetBranchAddress("MET_RefMuon_etx", &MET_RefMuon_etx);
			t->SetBranchAddress("MET_RefMuon_ety", &MET_RefMuon_ety);
			t->SetBranchAddress("MET_RefMuon_phi", &MET_RefMuon_phi);
			t->SetBranchAddress("MET_RefMuon_et", &MET_RefMuon_et);
			t->SetBranchAddress("MET_RefMuon_sumet", &MET_RefMuon_sumet);
        	
			t->SetBranchAddress("MET_RefTau_etx", &MET_RefTau_etx);
			t->SetBranchAddress("MET_RefTau_ety", &MET_RefTau_ety);
			t->SetBranchAddress("MET_RefTau_phi", &MET_RefTau_phi);
			t->SetBranchAddress("MET_RefTau_et", &MET_RefTau_et);
			t->SetBranchAddress("MET_RefTau_sumet", &MET_RefTau_sumet);
        	
			t->SetBranchAddress("MET_CellOut_Eflow_etx", &MET_CellOut_Eflow_etx);
			t->SetBranchAddress("MET_CellOut_Eflow_ety", &MET_CellOut_Eflow_ety);
			t->SetBranchAddress("MET_CellOut_Eflow_phi", &MET_CellOut_Eflow_phi);
			t->SetBranchAddress("MET_CellOut_Eflow_et", &MET_CellOut_Eflow_et);
			t->SetBranchAddress("MET_CellOut_Eflow_sumet", &MET_CellOut_Eflow_sumet);
        	
			t->SetBranchAddress("MET_Track_etx", &MET_Track_etx);
			t->SetBranchAddress("MET_Track_ety", &MET_Track_ety);
			t->SetBranchAddress("MET_Track_phi", &MET_Track_phi);
			t->SetBranchAddress("MET_Track_et", &MET_Track_et);
			t->SetBranchAddress("MET_Track_sumet", &MET_Track_sumet);
        	
			t->SetBranchAddress("MET_Staco_etx", &MET_Staco_etx);
			t->SetBranchAddress("MET_Staco_ety", &MET_Staco_ety);
			t->SetBranchAddress("MET_Staco_phi", &MET_Staco_phi);
			t->SetBranchAddress("MET_Staco_et", &MET_Staco_et);
			t->SetBranchAddress("MET_Staco_sumet", &MET_Staco_sumet);
        	
			t->SetBranchAddress("MET_Muid_etx", &MET_Muid_etx);
			t->SetBranchAddress("MET_Muid_ety", &MET_Muid_ety);
			t->SetBranchAddress("MET_Muid_phi", &MET_Muid_phi);
			t->SetBranchAddress("MET_Muid_et", &MET_Muid_et);
			t->SetBranchAddress("MET_Muid_sumet", &MET_Muid_sumet);
        	
			t->SetBranchAddress("MET_RefGamma_etx", &MET_RefGamma_etx);
			t->SetBranchAddress("MET_RefGamma_ety", &MET_RefGamma_ety);
			t->SetBranchAddress("MET_RefGamma_phi", &MET_RefGamma_phi);
			t->SetBranchAddress("MET_RefGamma_et", &MET_RefGamma_et);
			t->SetBranchAddress("MET_RefGamma_sumet", &MET_RefGamma_sumet);
        	
			t->SetBranchAddress("MET_RefEle_etx", &MET_RefEle_etx);
			t->SetBranchAddress("MET_RefEle_ety", &MET_RefEle_ety);
			t->SetBranchAddress("MET_RefEle_phi", &MET_RefEle_phi);
			t->SetBranchAddress("MET_RefEle_et", &MET_RefEle_et);
			t->SetBranchAddress("MET_RefEle_sumet", &MET_RefEle_sumet);
        	
			t->SetBranchAddress("MET_RefJet_etx", &MET_RefJet_etx);
			t->SetBranchAddress("MET_RefJet_ety", &MET_RefJet_ety);
			t->SetBranchAddress("MET_RefJet_phi", &MET_RefJet_phi);
			t->SetBranchAddress("MET_RefJet_et", &MET_RefJet_et);
			t->SetBranchAddress("MET_RefJet_sumet", &MET_RefJet_sumet);
        	
			t->SetBranchAddress("MET_SoftJets_etx", &MET_SoftJets_etx);
			t->SetBranchAddress("MET_SoftJets_ety", &MET_SoftJets_ety);
			t->SetBranchAddress("MET_SoftJets_phi", &MET_SoftJets_phi);
			t->SetBranchAddress("MET_SoftJets_et", &MET_SoftJets_et);
			t->SetBranchAddress("MET_SoftJets_sumet", &MET_SoftJets_sumet);
        	
			t->SetBranchAddress("MET_Truth_NonInt_etx", &MET_Truth_NonInt_etx);
			t->SetBranchAddress("MET_Truth_NonInt_ety", &MET_Truth_NonInt_ety);
			t->SetBranchAddress("MET_Truth_NonInt_phi", &MET_Truth_NonInt_phi);
			t->SetBranchAddress("MET_Truth_NonInt_et", &MET_Truth_NonInt_et);
			t->SetBranchAddress("MET_Truth_NonInt_sumet", &MET_Truth_NonInt_sumet);
			
			// t->SetBranchAddress("trk_n", &trk_n);
			// t->SetBranchAddress("trk_pt", &trk_pt);
			// t->SetBranchAddress("trk_eta", &trk_eta);
			// t->SetBranchAddress("trk_phi_wrtPV", &trk_phi_wrtPV);
			//         	
			// t->SetBranchAddress("mu_n", &mu_muons_n);
			// t->SetBranchAddress("mu_E", &mu_muons_E);
			// t->SetBranchAddress("mu_pt", &mu_muons_pt);
			// t->SetBranchAddress("mu_m", &mu_muons_m);
			// t->SetBranchAddress("mu_eta", &mu_muons_eta);
			// t->SetBranchAddress("mu_phi", &mu_muons_phi);
			// t->SetBranchAddress("mu_charge", &mu_muons_charge);
			// t->SetBranchAddress("mu_ms_phi", &mu_muons_ms_phi);
			// t->SetBranchAddress("mu_ms_theta", &mu_muons_ms_theta);
			// t->SetBranchAddress("mu_ms_qoverp", &mu_muons_ms_qoverp);
        	
			t->SetBranchAddress("mu_staco_n", &mu_staco_n);
			t->SetBranchAddress("mu_staco_E", &mu_staco_E);
			t->SetBranchAddress("mu_staco_pt", &mu_staco_pt);
			t->SetBranchAddress("mu_staco_m", &mu_staco_m);
			t->SetBranchAddress("mu_staco_eta", &mu_staco_eta);
			t->SetBranchAddress("mu_staco_phi", &mu_staco_phi);
			t->SetBranchAddress("mu_staco_charge", &mu_staco_charge);
			t->SetBranchAddress("mu_staco_ms_phi", &mu_staco_ms_phi);
			t->SetBranchAddress("mu_staco_ms_theta", &mu_staco_ms_theta);
			t->SetBranchAddress("mu_staco_ms_qoverp", &mu_staco_ms_qoverp);
        	
			// t->SetBranchAddress("mu_muid_n", &mu_muid_n);
			// t->SetBranchAddress("mu_muid_E", &mu_muid_E);
			// t->SetBranchAddress("mu_muid_pt", &mu_muid_pt);
			// t->SetBranchAddress("mu_muid_m", &mu_muid_m);
			// t->SetBranchAddress("mu_muid_eta", &mu_muid_eta);
			// t->SetBranchAddress("mu_muid_phi", &mu_muid_phi);
			// t->SetBranchAddress("mu_muid_charge", &mu_muid_charge);
			// t->SetBranchAddress("mu_muid_ms_phi", &mu_muid_ms_phi);
			// t->SetBranchAddress("mu_muid_ms_theta", &mu_muid_ms_theta);
			// t->SetBranchAddress("mu_muid_ms_qoverp", &mu_muid_ms_qoverp);
			//        	
			t->SetBranchAddress("el_n", &el_n);
			t->SetBranchAddress("el_E", &el_E);
			t->SetBranchAddress("el_Et", &el_Et);
			t->SetBranchAddress("el_pt", &el_pt);
			t->SetBranchAddress("el_m", &el_m);
			t->SetBranchAddress("el_eta", &el_eta);
			t->SetBranchAddress("el_phi", &el_phi);
			         	
			t->SetBranchAddress("ph_n", &ph_n);
			t->SetBranchAddress("ph_E", &ph_E);
			t->SetBranchAddress("ph_Et", &ph_Et);
			t->SetBranchAddress("ph_pt", &ph_pt);
			t->SetBranchAddress("ph_m", &ph_m);
			t->SetBranchAddress("ph_eta", &ph_eta);
			t->SetBranchAddress("ph_phi", &ph_phi);
			         	
			// t->SetBranchAddress("tau_n", &tau_n);
			// t->SetBranchAddress("tau_Et", &tau_Et);
			// t->SetBranchAddress("tau_pt", &tau_pt);
			// t->SetBranchAddress("tau_m", &tau_m);
			// t->SetBranchAddress("tau_eta", &tau_eta);
			// t->SetBranchAddress("tau_phi", &tau_phi);

			t->SetBranchAddress("el_MET_n", &el_MET_n);
			t->SetBranchAddress("el_MET_wpx", &el_MET_wpx);
			t->SetBranchAddress("el_MET_wpy", &el_MET_wpy);
			t->SetBranchAddress("el_MET_wet", &el_MET_wet);
			t->SetBranchAddress("el_MET_statusWord", &el_MET_statusWord);
			t->SetBranchAddress("ph_MET_n", &ph_MET_n);
			t->SetBranchAddress("ph_MET_wpx", &ph_MET_wpx);
			t->SetBranchAddress("ph_MET_wpy", &ph_MET_wpy);
			t->SetBranchAddress("ph_MET_wet", &ph_MET_wet);
			t->SetBranchAddress("ph_MET_statusWord", &ph_MET_statusWord);
			t->SetBranchAddress("mu_staco_MET_n", &mu_staco_MET_n);
			t->SetBranchAddress("mu_staco_MET_wpx", &mu_staco_MET_wpx);
			t->SetBranchAddress("mu_staco_MET_wpy", &mu_staco_MET_wpy);
			t->SetBranchAddress("mu_staco_MET_wet", &mu_staco_MET_wet);
			t->SetBranchAddress("mu_staco_MET_statusWord", &mu_staco_MET_statusWord);
			// t->SetBranchAddress("mu_muid_MET_n", &mu_muid_MET_n);
			// t->SetBranchAddress("mu_muid_MET_wpx", &mu_muid_MET_wpx);
			// t->SetBranchAddress("mu_muid_MET_wpy", &mu_muid_MET_wpy);
			// t->SetBranchAddress("mu_muid_MET_wet", &mu_muid_MET_wet);
			// t->SetBranchAddress("mu_muid_MET_statusWord", &mu_muid_MET_statusWord);
			// t->SetBranchAddress("mu_MET_n", &mu_MET_n);
			// t->SetBranchAddress("mu_MET_wpx", &mu_MET_wpx);
			// t->SetBranchAddress("mu_MET_wpy", &mu_MET_wpy);
			// t->SetBranchAddress("mu_MET_wet", &mu_MET_wet);
			// t->SetBranchAddress("mu_MET_statusWord", &mu_MET_statusWord);
			// t->SetBranchAddress("tau_MET_n", &tau_MET_n);
			// t->SetBranchAddress("tau_MET_wpx", &tau_MET_wpx);
			// t->SetBranchAddress("tau_MET_wpy", &tau_MET_wpy);
			// t->SetBranchAddress("tau_MET_wet", &tau_MET_wet);
			// t->SetBranchAddress("tau_MET_statusWord", &tau_MET_statusWord);
			t->SetBranchAddress("jet_antikt4LCtopo_MET_n", &jet_antikt4LCtopo_MET_n);
			t->SetBranchAddress("jet_antikt4LCtopo_MET_wpx", &jet_antikt4LCtopo_MET_wpx);
			t->SetBranchAddress("jet_antikt4LCtopo_MET_wpy", &jet_antikt4LCtopo_MET_wpy);
			t->SetBranchAddress("jet_antikt4LCtopo_MET_wet", &jet_antikt4LCtopo_MET_wet);
			t->SetBranchAddress("jet_antikt4LCtopo_MET_statusWord", &jet_antikt4LCtopo_MET_statusWord);
			// t->SetBranchAddress("cl_MET_n", &cl_MET_n);
			// t->SetBranchAddress("cl_MET_wpx", &cl_MET_wpx);
			// t->SetBranchAddress("cl_MET_wpy", &cl_MET_wpy);
			// t->SetBranchAddress("cl_MET_wet", &cl_MET_wet);
			// t->SetBranchAddress("cl_MET_statusWord", &cl_MET_statusWord);
			// t->SetBranchAddress("trk_MET_n", &trk_MET_n);
			// t->SetBranchAddress("trk_MET_wpx", &trk_MET_wpx);
			// t->SetBranchAddress("trk_MET_wpy", &trk_MET_wpy);
			// t->SetBranchAddress("trk_MET_wet", &trk_MET_wet);
			// t->SetBranchAddress("trk_MET_statusWord", &trk_MET_statusWord);
		}
		
		t->SetBranchAddress("rhorhoKt3EM", &rhorhoKt3EM);
		t->SetBranchAddress("rhorhoKt4EM", &rhorhoKt4EM);
		t->SetBranchAddress("rhorhoKt3LC", &rhorhoKt3LC);
		t->SetBranchAddress("rhorhoKt4LC", &rhorhoKt4LC);
		
		if(glob_isSig)
		{
			t->SetBranchAddress("musp_eta", &musp_eta);
			t->SetBranchAddress("musp_phi", &musp_phi);
			t->SetBranchAddress("musp_innerSegments", &musp_innerSegments);
			t->SetBranchAddress("musp_middleSegments", &musp_middleSegments);
			t->SetBranchAddress("musp_outerSegments", &musp_outerSegments);
		
			t->SetBranchAddress("vxp_nTracks", &vxp_nTracks);
		}

		if(skim)
		{
		  	t->SetBranchStatus("trig_*_mu18it_*", 0);
		  	t->SetBranchStatus("trig_*_mu4T_j*",  0);
		  	t->SetBranchStatus("trig_*_mu6i*",  0);
		  	t->SetBranchStatus("trig_*_Jpsi*",  0);
		  	t->SetBranchStatus("trig_*_a4tchad*",  0);
		  	t->SetBranchStatus("trig_*2mu8_EFxe40*",  0);
		  	t->SetBranchStatus("trig_*l2muonSA*",  0);
		  	t->SetBranchStatus("trig_*_Trk_*",  0);
		  	t->SetBranchStatus("trig_*Upsi*",  0);
		  	t->SetBranchStatus("trig_*Bmumux*",  0);
		  	t->SetBranchStatus("trig_*_c4cchad_*",  0);
		  	t->SetBranchStatus("AntiKt4LCTopoLets_ptconst_default*",  0);
		
			if(glob_isSig)
			{	
				t->SetBranchStatus("trk_*", 0);
				//t->SetBranchStatus("el_*", 0);
				//t->SetBranchStatus("ph_*", 0);
				t->SetBranchStatus("tau_*", 0);
				t->SetBranchStatus("cl_*", 0);
		    	
				//t->SetBranchStatus("el_MET_*", 0);
				//t->SetBranchStatus("ph_MET_*", 0);
				t->SetBranchStatus("mu_muid_MET_*", 0);
				t->SetBranchStatus("mu_MET_*", 0);
				t->SetBranchStatus("tau_MET_*", 0);
				t->SetBranchStatus("cl_MET_*", 0);
				t->SetBranchStatus("trk_MET_*", 0);
			}
		}
	}
	

	if(tType=="MUONS_TRIPLET" || tType=="MUID_TRIPLET")
	{
		TString prefix = "VTX_";

		vtx_index = 0;
		vtx_charge = 0;
		vtx_chi2 = 0;
		vtx_ndf = 0;
		vtx_mass = 0;
		vtx_pt = 0;
		vtx_rapidity = 0;
		
		vtx_x = 0;
		vtx_y = 0;
		vtx_z = 0;
		vtx_xErr = 0;
		vtx_yErr = 0;
		vtx_zErr = 0;
		
		vtx_refPVx = 0;
		vtx_refPVy = 0;
		vtx_refPVz = 0;
		vtx_refPVcovxx = 0;
		vtx_refPVcovyy = 0;
		vtx_refPVcovzz = 0;
		
		vtx_reftrks_px = 0;
		vtx_reftrks_py = 0;
		vtx_reftrks_pz = 0;
		
		vtx_srcIndex = 0;
		vtx_trkIndex = 0;
		vtx_srcName = 0;
		
		vtx_lxy = 0;
		vtx_lxyErr = 0;
		vtx_tau = 0;
		vtx_tauErr = 0;
		vtx_a0 = 0;
		vtx_a0Err = 0;
		vtx_a0XY = 0;
		vtx_a0XYErr = 0;
		vtx_cosTheta = 0;
		vtx_cosThetaXY = 0;
		
		vtx_trks_index = 0;
		vtx_muons_index = 0;

		t->SetBranchAddress("runNumber",&RunNumber);
		t->SetBranchAddress("lumiBlock",&lbn);
		t->SetBranchAddress("eventNumber",&EventNumber);
		
		t->SetBranchAddress("index",&vtx_index);
		t->SetBranchAddress("charge",&vtx_charge);
		t->SetBranchAddress(prefix+"pt",&vtx_pt);
		t->SetBranchAddress(prefix+"mass",&vtx_mass);
		t->SetBranchAddress(prefix+"rapidity",&vtx_rapidity);
		t->SetBranchAddress(prefix+"chi2",&vtx_chi2);
		t->SetBranchAddress(prefix+"NDF",&vtx_ndf);
		
		t->SetBranchAddress(prefix+"xposition",&vtx_x);
		t->SetBranchAddress(prefix+"yposition",&vtx_y);
		t->SetBranchAddress(prefix+"zposition",&vtx_z);
		t->SetBranchAddress(prefix+"xxPosCov",&vtx_xErr);
		t->SetBranchAddress(prefix+"yyPosCov",&vtx_yErr);
		t->SetBranchAddress(prefix+"zzPosCov",&vtx_zErr);
		
		t->SetBranchAddress(prefix+"refTrksPx",&vtx_reftrks_px);
		t->SetBranchAddress(prefix+"refTrksPy",&vtx_reftrks_py);
		t->SetBranchAddress(prefix+"refTrksPz",&vtx_reftrks_pz);
		
		t->SetBranchAddress(prefix+"srcIndex",&vtx_srcIndex);
		t->SetBranchAddress(prefix+"trkIndex",&vtx_trkIndex);
		t->SetBranchAddress(prefix+"srcName",&vtx_srcName);
		
		t->SetBranchAddress(prefix+"lxy",&vtx_lxy);
		t->SetBranchAddress(prefix+"lxyError",&vtx_lxyErr);
		t->SetBranchAddress(prefix+"tau",&vtx_tau);
		t->SetBranchAddress(prefix+"tauError",&vtx_tauErr);
		t->SetBranchAddress(prefix+"a0",&vtx_a0);
		t->SetBranchAddress(prefix+"a0Error",&vtx_a0Err);
		t->SetBranchAddress(prefix+"a0xy",&vtx_a0XY);
		t->SetBranchAddress(prefix+"a0xyError",&vtx_a0XYErr);
		t->SetBranchAddress(prefix+"cosTheta",&vtx_cosTheta);
		t->SetBranchAddress(prefix+"cosThetaXY",&vtx_cosThetaXY);

		t->SetBranchAddress("TRKS_index",&vtx_trks_index);
		t->SetBranchAddress("MUONS_index",&vtx_muons_index);
		

		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"refPVx", 0);
			t->SetBranchStatus(prefix+"refPVy", 0);
			t->SetBranchStatus(prefix+"refPVz", 0);
			t->SetBranchStatus(prefix+"refPVcovxx", 0);
			t->SetBranchStatus(prefix+"refPVcovyy", 0);
			t->SetBranchStatus(prefix+"refPVcovzz", 0);
			
			t->SetBranchStatus(prefix+"ptError",   0);
			t->SetBranchStatus(prefix+"massError",   0);
			t->SetBranchStatus(prefix+"pxMassHyp",   0);
			t->SetBranchStatus(prefix+"pyMassHyp",   0);
			t->SetBranchStatus(prefix+"pzMassHyp",   0);
			t->SetBranchStatus(prefix+"thetaStar",   0);
			t->SetBranchStatus(prefix+"cosThetaStar",   0);
			t->SetBranchStatus(prefix+"phiStar",   0);
			t->SetBranchStatus(prefix+"fullCovDiagonal",   0);
			t->SetBranchStatus(prefix+"tauMassHyp",   0);
			t->SetBranchStatus(prefix+"tauErrorMassHyp",   0);
			t->SetBranchStatus(prefix+"tauCorrected",   0);
			t->SetBranchStatus(prefix+"tauErrorCorrected",   0);
			t->SetBranchStatus(prefix+"massTauCov",   0);
			
			t->SetBranchStatus("isSignal",         0);
			t->SetBranchStatus("label",            0);
			t->SetBranchStatus("childrenIndices",  0);
			t->SetBranchStatus("TRUE_TRKS_index",  0);
			t->SetBranchStatus("MUON_TRKS_index",  0);
			t->SetBranchStatus("ELECTRONS_index",  0);
			t->SetBranchStatus("ELECTRONS_GSFTRKS_index",  0);
			t->SetBranchStatus("TRUTH_pdgID",  0);
			t->SetBranchStatus("TRUTH_barcode",  0);
		}
	}
	

	
	if(tType.Contains("TRKPRT"))
	{	
		TString prefix = tType+"_TRKS_";
		TString tTypeShort = tType;
		tTypeShort.ReplaceAll("TRKPRT_","");
		
		tpmu_vd.insert(make_pair(tTypeShort+"_chi2"        ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_qOverPErr"   ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_BS_qOverPErr",(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_qOverP"      ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_BS_qOverP"   ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_pt"          ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_eta"         ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_px"          ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_py"          ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_pz"          ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_pixeldEdx"   ,(vector<double>*)0));
		tpmu_vd.insert(make_pair(tTypeShort+"_phi"         ,new vector<double>)); // this is not a branch
		
		tpmu_vi.insert(make_pair(tTypeShort+"_ndf"                       ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nBLayer"                   ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nPix"                      ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nSCT"                      ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nTRT"                      ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nTRTOutliers"              ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nHighThresholdTRTHits"     ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nPixHoles"                 ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nSCTHoles"                 ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nDeadPixels"               ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nDeadSCT"                  ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_expectBLayer"              ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_nUsedHitsdEdx"             ,(vector<int>*)0));

		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfMdtHits"           ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfTgcPhiHits"        ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfTgcEtaHits"        ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfCscPhiHits"        ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfCscEtaHits"        ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfRpcPhiHits"        ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfRpcEtaHits"        ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfCscEtaHoles"       ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfCscPhiHoles"       ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfRpcEtaHoles"       ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfRpcPhiHoles"       ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfMdtHoles"          ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfTgcEtaHoles"       ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfTgcPhiHoles"       ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_numberOfOutliersOnTrack"   ,(vector<int>*)0));
		tpmu_vi.insert(make_pair(tTypeShort+"_standardDeviationOfChi2OS" ,(vector<int>*)0));

		tpmu_vvi.insert(make_pair(tTypeShort+"_sectors"               ,(vector<vector<int> >*)0));
		tpmu_vvi.insert(make_pair(tTypeShort+"_nprecisionHits"        ,(vector<vector<int> >*)0));
		tpmu_vvi.insert(make_pair(tTypeShort+"_nphiLayers"            ,(vector<vector<int> >*)0));
		tpmu_vvi.insert(make_pair(tTypeShort+"_netaPhiLayers"         ,(vector<vector<int> >*)0));
		tpmu_vvi.insert(make_pair(tTypeShort+"_nprecisionHoles"       ,(vector<vector<int> >*)0));
		tpmu_vvi.insert(make_pair(tTypeShort+"_netaTriggerHoleLayers" ,(vector<vector<int> >*)0));
		tpmu_vvi.insert(make_pair(tTypeShort+"_nphiHoleLayers"        ,(vector<vector<int> >*)0));
		tpmu_vvi.insert(make_pair(tTypeShort+"_nprecisionOutliers"    ,(vector<vector<int> >*)0));

		tpmu_vvs.insert(make_pair(tTypeShort+"_stationName"           ,(vector<vector<string> >*)0));
		tpmu_vvs.insert(make_pair(tTypeShort+"_phiName"               ,(vector<vector<string> >*)0));


		
		t->SetBranchAddress(prefix+"chi2",         &tpmu_vd[tTypeShort+"_chi2"        ]);
		t->SetBranchAddress(prefix+"qOverPErr",    &tpmu_vd[tTypeShort+"_qOverPErr"   ]);
		t->SetBranchAddress(prefix+"BS_qOverPErr", &tpmu_vd[tTypeShort+"_BS_qOverPErr"]);
		t->SetBranchAddress(prefix+"qOverP",       &tpmu_vd[tTypeShort+"_qOverP"      ]);
		t->SetBranchAddress(prefix+"BS_qOverP",    &tpmu_vd[tTypeShort+"_BS_qOverP"   ]);
		t->SetBranchAddress(prefix+"pt",           &tpmu_vd[tTypeShort+"_pt"          ]);
		t->SetBranchAddress(prefix+"eta",          &tpmu_vd[tTypeShort+"_eta"         ]);
		t->SetBranchAddress(prefix+"px",           &tpmu_vd[tTypeShort+"_px"          ]);
		t->SetBranchAddress(prefix+"py",           &tpmu_vd[tTypeShort+"_py"          ]);
		t->SetBranchAddress(prefix+"pz",           &tpmu_vd[tTypeShort+"_pz"          ]);
		t->SetBranchAddress(prefix+"pixeldEdx",    &tpmu_vd[tTypeShort+"_pixeldEdx"   ]);
		
		t->SetBranchAddress(prefix+"ndf",                   &tpmu_vi[tTypeShort+"_ndf"                  ]);
		t->SetBranchAddress(prefix+"nBLayer",               &tpmu_vi[tTypeShort+"_nBLayer"              ]);
		t->SetBranchAddress(prefix+"nPix",                  &tpmu_vi[tTypeShort+"_nPix"                 ]);
		t->SetBranchAddress(prefix+"nSCT",                  &tpmu_vi[tTypeShort+"_nSCT"                 ]);
		t->SetBranchAddress(prefix+"nTRT",                  &tpmu_vi[tTypeShort+"_nTRT"                 ]);
		t->SetBranchAddress(prefix+"nTRTOutliers",          &tpmu_vi[tTypeShort+"_nTRTOutliers"         ]);
		t->SetBranchAddress(prefix+"nHighThresholdTRTHits", &tpmu_vi[tTypeShort+"_nHighThresholdTRTHits"]);
		t->SetBranchAddress(prefix+"nPixHoles",             &tpmu_vi[tTypeShort+"_nPixHoles"            ]);
		t->SetBranchAddress(prefix+"nSCTHoles",             &tpmu_vi[tTypeShort+"_nSCTHoles"            ]);
		t->SetBranchAddress(prefix+"nDeadPixels",           &tpmu_vi[tTypeShort+"_nDeadPixels"          ]);
		t->SetBranchAddress(prefix+"nDeadSCT",              &tpmu_vi[tTypeShort+"_nDeadSCT"             ]);
		t->SetBranchAddress(prefix+"expectBLayer",          &tpmu_vi[tTypeShort+"_expectBLayer"         ]);
		t->SetBranchAddress(prefix+"nUsedHitsdEdx",         &tpmu_vi[tTypeShort+"_nUsedHitsdEdx"        ]);

		t->SetBranchAddress(prefix+"numberOfMdtHits",            &tpmu_vi[tTypeShort+"_numberOfMdtHits"]);
		t->SetBranchAddress(prefix+"numberOfTgcPhiHits",         &tpmu_vi[tTypeShort+"_numberOfTgcPhiHits"]);
		t->SetBranchAddress(prefix+"numberOfTgcEtaHits",         &tpmu_vi[tTypeShort+"_numberOfTgcEtaHits"]);
		t->SetBranchAddress(prefix+"numberOfCscPhiHits",         &tpmu_vi[tTypeShort+"_numberOfCscPhiHits"]);
		t->SetBranchAddress(prefix+"numberOfCscEtaHits",         &tpmu_vi[tTypeShort+"_numberOfCscEtaHits"]);
		t->SetBranchAddress(prefix+"numberOfRpcPhiHits",         &tpmu_vi[tTypeShort+"_numberOfRpcPhiHits"]);
		t->SetBranchAddress(prefix+"numberOfRpcEtaHits",         &tpmu_vi[tTypeShort+"_numberOfRpcEtaHits"]);
		t->SetBranchAddress(prefix+"numberOfCscEtaHoles",        &tpmu_vi[tTypeShort+"_numberOfCscEtaHoles"]);
		t->SetBranchAddress(prefix+"numberOfCscPhiHoles",        &tpmu_vi[tTypeShort+"_numberOfCscPhiHoles"]);
		t->SetBranchAddress(prefix+"numberOfRpcEtaHoles",        &tpmu_vi[tTypeShort+"_numberOfRpcEtaHoles"]);
		t->SetBranchAddress(prefix+"numberOfRpcPhiHoles",        &tpmu_vi[tTypeShort+"_numberOfRpcPhiHoles"]);
		t->SetBranchAddress(prefix+"numberOfMdtHoles",           &tpmu_vi[tTypeShort+"_numberOfMdtHoles"]);
		t->SetBranchAddress(prefix+"numberOfTgcEtaHoles",        &tpmu_vi[tTypeShort+"_numberOfTgcEtaHoles"]);
		t->SetBranchAddress(prefix+"numberOfTgcPhiHoles",        &tpmu_vi[tTypeShort+"_numberOfTgcPhiHoles"]);
		t->SetBranchAddress(prefix+"numberOfOutliersOnTrack",    &tpmu_vi[tTypeShort+"_numberOfOutliersOnTrack"]);
		t->SetBranchAddress(prefix+"standardDeviationOfChi2OS",  &tpmu_vi[tTypeShort+"_standardDeviationOfChi2OS"]);
		
		t->SetBranchAddress(prefix+"sectors",                &tpmu_vvi[tTypeShort+"_sectors"]);
		t->SetBranchAddress(prefix+"nprecisionHits",         &tpmu_vvi[tTypeShort+"_nprecisionHits"]);
		t->SetBranchAddress(prefix+"nphiLayers",             &tpmu_vvi[tTypeShort+"_nphiLayers"]);
		t->SetBranchAddress(prefix+"netaPhiLayers",          &tpmu_vvi[tTypeShort+"_netaPhiLayers"]);
		t->SetBranchAddress(prefix+"nprecisionHoles",        &tpmu_vvi[tTypeShort+"_nprecisionHoles"]);
		t->SetBranchAddress(prefix+"netaTriggerHoleLayers",  &tpmu_vvi[tTypeShort+"_netaTriggerHoleLayers"]);
		t->SetBranchAddress(prefix+"nphiHoleLayers",         &tpmu_vvi[tTypeShort+"_nphiHoleLayers"]);
		t->SetBranchAddress(prefix+"nprecisionOutliers",     &tpmu_vvi[tTypeShort+"_nprecisionOutliers"]);
		
		t->SetBranchAddress(prefix+"stationName",            &tpmu_vvs[tTypeShort+"_stationName"]);
		t->SetBranchAddress(prefix+"phiName",                &tpmu_vvs[tTypeShort+"_phiName"]);

		
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"index", 0);
			t->SetBranchStatus(prefix+"vtxIndex", 0);
			t->SetBranchStatus(prefix+"BS_d0", 0);
			t->SetBranchStatus(prefix+"BS_z0", 0);
			t->SetBranchStatus(prefix+"BS_phi0", 0);
			t->SetBranchStatus(prefix+"BS_theta", 0);
			t->SetBranchStatus(prefix+"BS_d0Err", 0);
			t->SetBranchStatus(prefix+"BS_z0Err", 0);
			t->SetBranchStatus(prefix+"BS_phi0Err", 0);
			t->SetBranchStatus(prefix+"BS_thetaErr", 0);
			
			t->SetBranchStatus(prefix+"d0", 0);
			t->SetBranchStatus(prefix+"z0", 0);
			t->SetBranchStatus(prefix+"extrapZ0", 0);
			t->SetBranchStatus(prefix+"extrapD0", 0);
			t->SetBranchStatus(prefix+"phi0", 0);
			t->SetBranchStatus(prefix+"theta", 0);
			t->SetBranchStatus(prefix+"d0Err", 0);
			t->SetBranchStatus(prefix+"z0Err", 0);
			t->SetBranchStatus(prefix+"phi0Err", 0);
			t->SetBranchStatus(prefix+"thetaErr", 0);
		}
	}
	
	
	if(tType=="SEL_TRACKS")
	{
		TString prefix = "SEL_TRACKS_TRKS_";
		
		trks_chi2 = 0;
		trks_qoverp = 0;
		trks_qoverpErr = 0;
		trks_pt = 0;
		trks_eta = 0;
		trks_px = 0;
		trks_py = 0;
		trks_pz = 0;
		
		trks_ndf = 0;
		trks_nBLayer = 0;
		trks_nPix = 0;
		trks_nSCT = 0;
		trks_nTRT = 0;
		trks_nTRTOutliers = 0;
		trks_nHighThresholdTRTHits = 0;
		trks_nPixHoles = 0;
		trks_nSCTHoles = 0;
		trks_nDeadPixels = 0;
		trks_nDeadSCT = 0;
		trks_expectBLayer = 0;
		
		trks_pixeldEdx = 0;
		trks_nUsedHitsdEdx = 0;
		
		t->SetBranchAddress(prefix+"chi2",&trks_chi2);
		t->SetBranchAddress(prefix+"qOverPErr",&trks_qoverpErr);
		t->SetBranchAddress(prefix+"qOverP",&trks_qoverp);
		t->SetBranchAddress(prefix+"pt",&trks_pt);
		t->SetBranchAddress(prefix+"eta",&trks_eta);
		t->SetBranchAddress(prefix+"px",&trks_px);
		t->SetBranchAddress(prefix+"py",&trks_py);
		t->SetBranchAddress(prefix+"pz",&trks_pz);
		
		t->SetBranchAddress(prefix+"ndf",&trks_ndf);
		t->SetBranchAddress(prefix+"nBLayer",&trks_nBLayer);
		t->SetBranchAddress(prefix+"nPix",&trks_nPix);
		t->SetBranchAddress(prefix+"nSCT",&trks_nSCT);
		t->SetBranchAddress(prefix+"nTRT",&trks_nTRT);
		t->SetBranchAddress(prefix+"nTRTOutliers",&trks_nTRTOutliers);
		t->SetBranchAddress(prefix+"nHighThresholdTRTHits",&trks_nHighThresholdTRTHits);
		t->SetBranchAddress(prefix+"nPixHoles",&trks_nPixHoles);
		t->SetBranchAddress(prefix+"nSCTHoles",&trks_nSCTHoles);
		t->SetBranchAddress(prefix+"nDeadPixels",&trks_nDeadPixels);
		t->SetBranchAddress(prefix+"nDeadSCT",&trks_nDeadSCT);
		t->SetBranchAddress(prefix+"expectBLayer",&trks_expectBLayer);
		t->SetBranchAddress(prefix+"pixeldEdx",    &trks_pixeldEdx);
		t->SetBranchAddress(prefix+"nUsedHitsdEdx",&trks_nUsedHitsdEdx);
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"index", 0);
			t->SetBranchStatus(prefix+"vtxIndex", 0);
			t->SetBranchStatus(prefix+"BS_d0", 0);
			t->SetBranchStatus(prefix+"BS_z0", 0);
			t->SetBranchStatus(prefix+"BS_phi0", 0);
			t->SetBranchStatus(prefix+"BS_theta", 0);
			t->SetBranchStatus(prefix+"BS_qOverP", 0);
			t->SetBranchStatus(prefix+"BS_d0Err", 0);
			t->SetBranchStatus(prefix+"BS_z0Err", 0);
			t->SetBranchStatus(prefix+"BS_phi0Err", 0);
			t->SetBranchStatus(prefix+"BS_thetaErr", 0);
			t->SetBranchStatus(prefix+"BS_qOverPErr", 0);
			
			t->SetBranchStatus(prefix+"d0", 0);
			t->SetBranchStatus(prefix+"z0", 0);
			t->SetBranchStatus(prefix+"extrapZ0", 0);
			t->SetBranchStatus(prefix+"extrapD0", 0);
			t->SetBranchStatus(prefix+"phi0", 0);
			t->SetBranchStatus(prefix+"theta", 0);
			t->SetBranchStatus(prefix+"d0Err", 0);
			t->SetBranchStatus(prefix+"z0Err", 0);
			t->SetBranchStatus(prefix+"phi0Err", 0);
			t->SetBranchStatus(prefix+"thetaErr", 0);
		}
	}
	
	if(tType=="MC")
	{
		TString prefix = "mc_";
		
		mc_pdgId = 0;
		mc_status = 0;
		mc_barcode = 0;
		mc_productionvtx_barcode = 0;
		mc_decayvtx_barcode = 0;
		mc_has_productionvtx = 0;
		mc_has_decayvtx = 0;
		mc_pt = 0;
		mc_eta = 0;
		mc_phi = 0;
		mc_m = 0;
		mc_charge = 0;
		mc_productionvtx_x = 0;
		mc_productionvtx_y = 0;
		mc_productionvtx_z = 0;
		mc_decayvtx_x = 0;
		mc_decayvtx_y = 0;
		mc_decayvtx_z = 0;
		mc_parents = 0;
		mc_children = 0;

		t->SetBranchAddress(prefix+"channel_number",&mc_channel_number);
		t->SetBranchAddress(prefix+"event_number",&mc_event_number);
		t->SetBranchAddress(prefix+"event_weight",&mc_event_weight);
		t->SetBranchAddress(prefix+"pdgId",&mc_pdgId);
		t->SetBranchAddress(prefix+"status",&mc_status);
		t->SetBranchAddress(prefix+"barcode",&mc_barcode);		
		t->SetBranchAddress(prefix+"pt",&mc_pt);
		t->SetBranchAddress(prefix+"eta",&mc_eta);
		t->SetBranchAddress(prefix+"phi",&mc_phi);
		t->SetBranchAddress(prefix+"charge",&mc_charge);
		t->SetBranchAddress(prefix+"m",&mc_m);
		
		t->SetBranchAddress(prefix+"has_production_vx",&mc_has_productionvtx);
		t->SetBranchAddress(prefix+"production_vx_barcode",&mc_productionvtx_barcode);
		t->SetBranchAddress(prefix+"production_vx_x",&mc_productionvtx_x);
		t->SetBranchAddress(prefix+"production_vx_y",&mc_productionvtx_y);
		t->SetBranchAddress(prefix+"production_vx_z",&mc_productionvtx_z);
		t->SetBranchAddress(prefix+"has_decay_vx",&mc_has_decayvtx);
		t->SetBranchAddress(prefix+"decay_vx_barcode",&mc_decayvtx_barcode);
		t->SetBranchAddress(prefix+"decay_vx_x",&mc_decayvtx_x);
		t->SetBranchAddress(prefix+"decay_vx_y",&mc_decayvtx_y);
		t->SetBranchAddress(prefix+"decay_vx_z",&mc_decayvtx_z);
		
		t->SetBranchAddress(prefix+"parents",&mc_parents);
		t->SetBranchAddress(prefix+"children",&mc_children);
	}
	
	if(tType=="PRIVX")
	{
		TString prefix = "PRIVX_";
		
		pv_mu = 0;
		pv_x = 0;
		pv_y = 0;
		pv_z = 0;
		pv_xErr = 0;
		pv_yErr = 0;
		pv_zErr = 0;
		pv_chi2 = 0;
		pv_ndf = 0;
		pv_ntrk = 0;
		pv_type = 0;
		pv_index = 0;
		
		t->SetBranchAddress(prefix+"mu",&pv_mu);
		t->SetBranchAddress(prefix+"ntrk",&pv_ntrk);
		t->SetBranchAddress(prefix+"type",&pv_type);
		t->SetBranchAddress(prefix+"index",&pv_index);
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"chi2",  0);
			t->SetBranchStatus(prefix+"ndf",   0);
			t->SetBranchStatus("BS_VTX_X",   0);
			t->SetBranchStatus("BS_VTX_Y",   0);
			t->SetBranchStatus("BS_VTX_Z",   0);
			t->SetBranchStatus("BS_POS_X",   0);
			t->SetBranchStatus("BS_POS_Y",   0);
			t->SetBranchStatus("BS_POS_Z",   0);
			t->SetBranchStatus("BS_SIGMA_X",   0);
			t->SetBranchStatus("BS_SIGMA_Y",   0);
			t->SetBranchStatus("BS_SIGMA_Z",   0);
			t->SetBranchStatus("BS_SIGMA_XY",   0);
			t->SetBranchStatus("BS_TILT_XZ",   0);
			t->SetBranchStatus("BS_TILT_YZ",   0);
			t->SetBranchStatus(prefix+"xposition",   0);
			t->SetBranchStatus(prefix+"yposition",   0);
			t->SetBranchStatus(prefix+"zposition",   0);
			t->SetBranchStatus(prefix+"xxPosCov",   0);
			t->SetBranchStatus(prefix+"yyPosCov",   0);
			t->SetBranchStatus(prefix+"zzPosCov",   0);
			t->SetBranchStatus(prefix+"xyPosCov",   0);
			t->SetBranchStatus(prefix+"yzPosCov",   0);
			t->SetBranchStatus(prefix+"xzPosCov",   0);
		}
	}

	
	if(tType=="AUX")
	{
		TString prefix = "AUX_";
				
		t->SetBranchAddress(prefix+"nTriplets",&aux_nTriplets);
		t->SetBranchAddress(prefix+"nDoublets",&aux_nDoublets);
		t->SetBranchAddress(prefix+"nQuadruplets",&aux_nQuadruplets);
		t->SetBranchAddress(prefix+"averageIntPerXing",&aux_averageIntPerXing);
		t->SetBranchAddress(prefix+"actualIntPerXing",&aux_actualIntPerXing);
	}
	
	if(tType=="TRIG")
	{
		TString prefix = "TRIG_";
		
		L1_trigger_name = 0;
		EF_trigger_name = 0;
		
		t->SetBranchAddress(prefix+"L1_trigger_name",&L1_trigger_name);
		t->SetBranchAddress(prefix+"EF_trigger_name",&EF_trigger_name);
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"offtimeTGC",   0);
			t->SetBranchStatus(prefix+"L1_prescale",   0);
			t->SetBranchStatus(prefix+"L2_prescale",   0);
			t->SetBranchStatus(prefix+"EF_prescale",   0);
			t->SetBranchStatus(prefix+"L2_pass_through",   0);
			t->SetBranchStatus(prefix+"EF_pass_through",   0);
		}
	}
	
	if(tType=="JET")
	{
		TString prefix = "JET_";
		
		jets_nTrks = 0;
		jets_LArQuality = 0;
		jets_isGood = 0;
		jets_isBad = 0;
		jets_isUgly = 0;
		jets_eMfrac = 0;
		jets_px = 0;
		jets_py = 0;
		jets_pz = 0;
		jets_pt = 0;
		jets_et = 0;
		jets_e = 0;
		jets_m = 0;
		jets_eta = 0;
		jets_phi = 0;
		jets_flwgt_MV1 = 0;
		
		t->SetBranchAddress(prefix+"nTrks",&jets_nTrks);
		t->SetBranchAddress(prefix+"LArQuality",&jets_LArQuality);
		t->SetBranchAddress(prefix+"isGood",&jets_isGood);
		t->SetBranchAddress(prefix+"isBad",&jets_isBad);
		t->SetBranchAddress(prefix+"isUgly",&jets_isUgly);
		t->SetBranchAddress(prefix+"EMfrac",&jets_eMfrac);
		t->SetBranchAddress(prefix+"px",&jets_px);
		t->SetBranchAddress(prefix+"py",&jets_py);
		t->SetBranchAddress(prefix+"pz",&jets_pz);
		t->SetBranchAddress(prefix+"pt",&jets_pt);
		t->SetBranchAddress(prefix+"et",&jets_et);
		t->SetBranchAddress(prefix+"e",&jets_e);
		t->SetBranchAddress(prefix+"m",&jets_m);
		t->SetBranchAddress(prefix+"eta",&jets_eta);
		t->SetBranchAddress(prefix+"phi",&jets_phi);
		t->SetBranchAddress(prefix+"flwgt_MV1",&jets_flwgt_MV1);
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"trkIdx",   0);
			t->SetBranchStatus(prefix+"flwgt_default",  0);
			t->SetBranchStatus(prefix+"flwgt_GbbNN", 0);
			t->SetBranchStatus(prefix+"flwgt_IP2D",  0);
			t->SetBranchStatus(prefix+"flwgt_IP3D",  0);
			t->SetBranchStatus(prefix+"flwgt_JetFitterCOMBNN",  0);
			t->SetBranchStatus(prefix+"flwgt_JetFitterCharm",  0);
			t->SetBranchStatus(prefix+"flwgt_JetFitterTagNN",  0);
			t->SetBranchStatus(prefix+"flwgt_SV0",  0);
			t->SetBranchStatus(prefix+"flwgt_SV1",  0);
			t->SetBranchStatus(prefix+"flwgt_SV2",  0);
			t->SetBranchStatus(prefix+"flwgt_SecondSoftMuonTagChi2",  0);
			t->SetBranchStatus(prefix+"flwgt_SoftMuonTagChi2",  0);
		}
	}
	
	if(tType=="MET")
	{
		TString prefix = "MET_";
		
		t->SetBranchAddress(prefix+"RefFinal_source",&met_RefFinal_source);
		t->SetBranchAddress(prefix+"RefFinal_etX",&met_RefFinal_etX);
		t->SetBranchAddress(prefix+"RefFinal_etY",&met_RefFinal_etY);
		t->SetBranchAddress(prefix+"RefFinal_sumET",&met_RefFinal_sumet);
		t->SetBranchAddress(prefix+"RefFinal_et",&met_RefFinal_et);
		t->SetBranchAddress(prefix+"RefFinal_phi",&met_RefFinal_phi);
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{

		}
	}
	
	if(tType=="MUONS")
	{
		TString prefix = "MUONS_MU_";
		
		muons_ptcone10 = 0;
		muons_ptcone20 = 0;
		muons_ptcone30 = 0;
		muons_ptcone40 = 0;
		muons_etcone10 = 0;
		muons_etcone20 = 0;
		muons_etcone30 = 0;
		muons_etcone40 = 0;
		muons_eLoss = 0;
		muons_eLossErr = 0;
		muons_charge = 0;
		muons_px = 0;
		muons_py = 0;
		muons_pz = 0;
		muons_e = 0;
		muons_pt = 0;
		muons_eta = 0;
		muons_phi = 0;
		muons_sctangsig = 0;
		muons_sctngbsig = 0;
		muons_pbalsig = 0;
		muons_isCombined = 0;
		muons_hasInDetTrackParticle = 0;
		muons_inDetTrackIndex = 0;
		muons_index = 0;
		muons_chi2 = 0;
		muons_ndf = 0;
		muons_author = 0;
		muons_matchchi2 = 0;
		muons_matchchi2ndf = 0;
		muons_phi_me = 0;
		muons_phi_ie = 0;
		muons_eta_me = 0;
		muons_eta_ie = 0;
		muons_pt_me = 0;
		muons_pt_ie = 0;		
		muons_px_me = 0;
		muons_px_ie = 0;
		muons_py_me = 0;
		muons_py_ie = 0;
		muons_pz_me = 0;
		muons_pz_ie = 0;
		muons_e_me = 0;
		muons_e_ie = 0;
		muons_isLoose = 0;
		muons_isMedium = 0;
		muons_isTight = 0;
		
		muons_numberOfMdtHits = 0;
		muons_numberOfTgcPhiHits = 0;
		muons_numberOfTgcEtaHits = 0;
		muons_numberOfCscPhiHits = 0;
		muons_numberOfCscEtaHits = 0;
		muons_numberOfRpcPhiHits = 0;
		muons_numberOfRpcEtaHits = 0;
		muons_numberOfCscEtaHoles = 0;
		muons_numberOfCscPhiHoles = 0;
		muons_numberOfRpcEtaHoles = 0;
		muons_numberOfRpcPhiHoles = 0;
		muons_numberOfMdtHoles = 0;
		muons_numberOfTgcEtaHoles = 0;
		muons_numberOfTgcPhiHoles = 0;
		muons_numberOfOutliersOnTrack = 0;
		muons_standardDeviationOfChi2OS = 0;

		muons_sectors = 0;
		muons_nprecisionHits = 0;
		muons_nphiLayers = 0;
		muons_netaPhiLayers = 0;
		muons_nprecisionHoles = 0;
		muons_netaTriggerHoleLayers = 0;
		muons_nphiHoleLayers = 0;
		muons_nprecisionOutliers = 0;

		muons_stationName = 0;
		muons_phiName = 0;
		
			
		t->SetBranchAddress(prefix+"ptcone10", &muons_ptcone10);
		t->SetBranchAddress(prefix+"ptcone20", &muons_ptcone20);
		t->SetBranchAddress(prefix+"ptcone30", &muons_ptcone30);
		t->SetBranchAddress(prefix+"ptcone40", &muons_ptcone40);
		t->SetBranchAddress(prefix+"etcone10", &muons_etcone10);
		t->SetBranchAddress(prefix+"etcone20", &muons_etcone20);
		t->SetBranchAddress(prefix+"etcone30", &muons_etcone30);
		t->SetBranchAddress(prefix+"etcone40", &muons_etcone40);
		t->SetBranchAddress(prefix+"eLoss",    &muons_eLoss);
		t->SetBranchAddress(prefix+"eLossErr", &muons_eLossErr);
		t->SetBranchAddress(prefix+"charge", &muons_charge);
		t->SetBranchAddress(prefix+"px", &muons_px);
		t->SetBranchAddress(prefix+"py", &muons_py);
		t->SetBranchAddress(prefix+"pz", &muons_pz);
		t->SetBranchAddress(prefix+"e", &muons_e);
		t->SetBranchAddress(prefix+"pt", &muons_pt);
		t->SetBranchAddress(prefix+"eta", &muons_eta);
		t->SetBranchAddress(prefix+"phi", &muons_phi);
		t->SetBranchAddress(prefix+"scatteringCurvatureSignificance", &muons_sctangsig);
		t->SetBranchAddress(prefix+"scatteringNeighbourSignificance", &muons_sctngbsig);
		t->SetBranchAddress(prefix+"momentumBalanceSignificance", &muons_pbalsig);
		t->SetBranchAddress(prefix+"isCombinedMuon", &muons_isCombined);
		t->SetBranchAddress(prefix+"hasInDetTrackParticle", &muons_hasInDetTrackParticle);
		t->SetBranchAddress(prefix+"inDetTrackIndex", &muons_inDetTrackIndex);
		t->SetBranchAddress(prefix+"index", &muons_index);
		t->SetBranchAddress(prefix+"trkComb_chi2", &muons_chi2);
		t->SetBranchAddress(prefix+"trkComb_nDoF", &muons_ndf);
		t->SetBranchAddress(prefix+"isLoose",  &muons_isLoose);
		t->SetBranchAddress(prefix+"isMedium", &muons_isMedium);
		t->SetBranchAddress(prefix+"isTight",  &muons_isTight);
		t->SetBranchAddress(prefix+"author", &muons_author);
		t->SetBranchAddress(prefix+"matchChi2", &muons_matchchi2);
		t->SetBranchAddress(prefix+"matchChi2OverDoF", &muons_matchchi2ndf);
		t->SetBranchAddress(prefix+"trkMuonExtr_phi", &muons_phi_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_eta", &muons_eta_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_pt", &muons_pt_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_px", &muons_px_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_py", &muons_py_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_pz", &muons_pz_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_e", &muons_e_me);
		t->SetBranchAddress(prefix+"trkInnerExtr_phi", &muons_phi_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_eta", &muons_eta_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_pt", &muons_pt_ie);		
		t->SetBranchAddress(prefix+"trkInnerExtr_px", &muons_px_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_py", &muons_py_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_pz", &muons_pz_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_e", &muons_e_ie);
		
		t->SetBranchAddress(prefix+"numberOfMDTHits",            &muons_numberOfMdtHits);
		t->SetBranchAddress(prefix+"numberOfTGCPhiHits",         &muons_numberOfTgcPhiHits);
		t->SetBranchAddress(prefix+"numberOfTGCEtaHits",         &muons_numberOfTgcEtaHits);
		t->SetBranchAddress(prefix+"numberOfCSCPhiHits",         &muons_numberOfCscPhiHits);
		t->SetBranchAddress(prefix+"numberOfCSCEtaHits",         &muons_numberOfCscEtaHits);
		t->SetBranchAddress(prefix+"numberOfRPCPhiHits",         &muons_numberOfRpcPhiHits);
		t->SetBranchAddress(prefix+"numberOfRPCEtaHits",         &muons_numberOfRpcEtaHits);
		t->SetBranchAddress(prefix+"numberOfCscEtaHoles",        &muons_numberOfCscEtaHoles);
		t->SetBranchAddress(prefix+"numberOfCscPhiHoles",        &muons_numberOfCscPhiHoles);
		t->SetBranchAddress(prefix+"numberOfRpcEtaHoles",        &muons_numberOfRpcEtaHoles);
		t->SetBranchAddress(prefix+"numberOfRpcPhiHoles",        &muons_numberOfRpcPhiHoles);
		t->SetBranchAddress(prefix+"numberOfMdtHoles",           &muons_numberOfMdtHoles);
		t->SetBranchAddress(prefix+"numberOfTgcEtaHoles",        &muons_numberOfTgcEtaHoles);
		t->SetBranchAddress(prefix+"numberOfTgcPhiHoles",        &muons_numberOfTgcPhiHoles);
		t->SetBranchAddress(prefix+"numberOfOutliersOnTrack",    &muons_numberOfOutliersOnTrack);
		t->SetBranchAddress(prefix+"standardDeviationOfChi2OS",  &muons_standardDeviationOfChi2OS);
		
		t->SetBranchAddress(prefix+"sectors",                &muons_sectors);
		t->SetBranchAddress(prefix+"nprecisionHits",         &muons_nprecisionHits);
		t->SetBranchAddress(prefix+"nphiLayers",             &muons_nphiLayers);
		t->SetBranchAddress(prefix+"netaPhiLayers",          &muons_netaPhiLayers);
		t->SetBranchAddress(prefix+"nprecisionHoles",        &muons_nprecisionHoles);
		t->SetBranchAddress(prefix+"netaTriggerHoleLayers",  &muons_netaTriggerHoleLayers);
		t->SetBranchAddress(prefix+"nphiHoleLayers",         &muons_nphiHoleLayers);
		t->SetBranchAddress(prefix+"nprecisionOutliers",     &muons_nprecisionOutliers);
		
		t->SetBranchAddress(prefix+"stationName",            &muons_stationName);
		t->SetBranchAddress(prefix+"phiName",                &muons_phiName);
		
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"trkMS_px", 0);
			t->SetBranchStatus(prefix+"trkComb_px", 0);
			t->SetBranchStatus(prefix+"trkMS_py", 0);
			t->SetBranchStatus(prefix+"trkComb_py", 0);
			t->SetBranchStatus(prefix+"trkMS_pz", 0);
			t->SetBranchStatus(prefix+"trkComb_pz", 0);
			t->SetBranchStatus(prefix+"trkMS_e", 0);
			t->SetBranchStatus(prefix+"trkComb_e", 0);
			t->SetBranchStatus(prefix+"cotTh", 0);
			t->SetBranchStatus(prefix+"trkMS_cotTh", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_cotTh", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_cotTh", 0);
			t->SetBranchStatus(prefix+"trkComb_cotTh", 0);
			t->SetBranchStatus(prefix+"trkMS_pt", 0);
			t->SetBranchStatus(prefix+"trkComb_pt", 0);
			t->SetBranchStatus(prefix+"trkMS_eta", 0);
			t->SetBranchStatus(prefix+"trkComb_eta", 0);
			t->SetBranchStatus(prefix+"trkMS_phi", 0);
			t->SetBranchStatus(prefix+"trkComb_phi", 0);
			t->SetBranchStatus(prefix+"trkMS_chi2", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_chi2", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_chi2", 0);
			t->SetBranchStatus(prefix+"trkMS_nDoF", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_nDoF", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_nDoF", 0);
			t->SetBranchStatus(prefix+"etaAtMsPivot", 0);
			t->SetBranchStatus(prefix+"phiAtMsPivot", 0);
			t->SetBranchStatus(prefix+"sigmaEtaAtMsPivot", 0);
			t->SetBranchStatus(prefix+"sigmaPhiAtMsPivot", 0);
			t->SetBranchStatus(prefix+"muonTrackIndex", 0);
		}
	}
	
	
	
	
	if(tType=="MUID")
	{
		TString prefix = "MUID_MU_";
		
		muid_ptcone10 = 0;
		muid_ptcone20 = 0;
		muid_ptcone30 = 0;
		muid_ptcone40 = 0;
		muid_etcone10 = 0;
		muid_etcone20 = 0;
		muid_etcone30 = 0;
		muid_etcone40 = 0;
		muid_eLoss = 0;
		muid_eLossErr = 0;
		muid_charge = 0;
		muid_px = 0;
		muid_py = 0;
		muid_pz = 0;
		muid_e = 0;
		muid_pt = 0;
		muid_eta = 0;
		muid_phi = 0;
		muid_sctangsig = 0;
		muid_sctngbsig = 0;
		muid_pbalsig = 0;
		muid_isCombined = 0;
		muid_hasInDetTrackParticle = 0;
		muid_inDetTrackIndex = 0;
		muid_index = 0;
		muid_chi2 = 0;
		muid_ndf = 0;
		muid_author = 0;
		muid_matchchi2 = 0;
		muid_matchchi2ndf = 0;
		muid_phi_me = 0;
		muid_phi_ie = 0;
		muid_eta_me = 0;
		muid_eta_ie = 0;
		muid_pt_me = 0;
		muid_pt_ie = 0;
		muid_px_me = 0;
		muid_px_ie = 0;
		muid_py_me = 0;
		muid_py_ie = 0;
		muid_pz_me = 0;
		muid_pz_ie = 0;
		muid_e_me = 0;
		muid_e_ie = 0;
		muid_isLoose = 0;
		muid_isMedium = 0;
		muid_isTight = 0;
		
		muid_numberOfMdtHits = 0;
		muid_numberOfTgcPhiHits = 0;
		muid_numberOfTgcEtaHits = 0;
		muid_numberOfCscPhiHits = 0;
		muid_numberOfCscEtaHits = 0;
		muid_numberOfRpcPhiHits = 0;
		muid_numberOfRpcEtaHits = 0;
		muid_numberOfCscEtaHoles = 0;
		muid_numberOfCscPhiHoles = 0;
		muid_numberOfRpcEtaHoles = 0;
		muid_numberOfRpcPhiHoles = 0;
		muid_numberOfMdtHoles = 0;
		muid_numberOfTgcEtaHoles = 0;
		muid_numberOfTgcPhiHoles = 0;
		muid_numberOfOutliersOnTrack = 0;
		muid_standardDeviationOfChi2OS = 0;

		muid_sectors = 0;
		muid_nprecisionHits = 0;
		muid_nphiLayers = 0;
		muid_netaPhiLayers = 0;
		muid_nprecisionHoles = 0;
		muid_netaTriggerHoleLayers = 0;
		muid_nphiHoleLayers = 0;
		muid_nprecisionOutliers = 0;

		muid_stationName = 0;
		muid_phiName = 0;
		
		t->SetBranchAddress(prefix+"ptcone10", &muid_ptcone10);
		t->SetBranchAddress(prefix+"ptcone20", &muid_ptcone20);
		t->SetBranchAddress(prefix+"ptcone30", &muid_ptcone30);
		t->SetBranchAddress(prefix+"ptcone40", &muid_ptcone40);
		t->SetBranchAddress(prefix+"etcone10", &muid_etcone10);
		t->SetBranchAddress(prefix+"etcone20", &muid_etcone20);
		t->SetBranchAddress(prefix+"etcone30", &muid_etcone30);
		t->SetBranchAddress(prefix+"etcone40", &muid_etcone40);
		t->SetBranchAddress(prefix+"eLoss",    &muid_eLoss);
		t->SetBranchAddress(prefix+"eLossErr", &muid_eLossErr);
		t->SetBranchAddress(prefix+"charge", &muid_charge);
		t->SetBranchAddress(prefix+"px", &muid_px);
		t->SetBranchAddress(prefix+"py", &muid_py);
		t->SetBranchAddress(prefix+"pz", &muid_pz);
		t->SetBranchAddress(prefix+"e", &muid_e);
		t->SetBranchAddress(prefix+"pt", &muid_pt);
		t->SetBranchAddress(prefix+"eta", &muid_eta);
		t->SetBranchAddress(prefix+"phi", &muid_phi);
		t->SetBranchAddress(prefix+"scatteringCurvatureSignificance", &muid_sctangsig);
		t->SetBranchAddress(prefix+"scatteringNeighbourSignificance", &muid_sctngbsig);
		t->SetBranchAddress(prefix+"momentumBalanceSignificance", &muid_pbalsig);
		t->SetBranchAddress(prefix+"isCombinedMuon", &muid_isCombined);
		t->SetBranchAddress(prefix+"hasInDetTrackParticle", &muid_hasInDetTrackParticle);
		t->SetBranchAddress(prefix+"inDetTrackIndex", &muid_inDetTrackIndex);
		t->SetBranchAddress(prefix+"index", &muid_index);
		t->SetBranchAddress(prefix+"trkComb_chi2", &muid_chi2);
		t->SetBranchAddress(prefix+"trkComb_nDoF", &muid_ndf);
		t->SetBranchAddress(prefix+"isLoose",  &muid_isLoose);
		t->SetBranchAddress(prefix+"isMedium", &muid_isMedium);
		t->SetBranchAddress(prefix+"isTight",  &muid_isTight);
		t->SetBranchAddress(prefix+"author", &muid_author);
		t->SetBranchAddress(prefix+"matchChi2", &muid_matchchi2);
		t->SetBranchAddress(prefix+"matchChi2OverDoF", &muid_matchchi2ndf);
		t->SetBranchAddress(prefix+"trkMuonExtr_phi", &muid_phi_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_eta", &muid_eta_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_pt", &muid_pt_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_px", &muid_px_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_py", &muid_py_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_pz", &muid_pz_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_e", &muid_e_me);
		t->SetBranchAddress(prefix+"trkInnerExtr_phi", &muid_phi_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_eta", &muid_eta_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_pt", &muid_pt_ie);		
		t->SetBranchAddress(prefix+"trkInnerExtr_px", &muid_px_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_py", &muid_py_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_pz", &muid_pz_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_e", &muid_e_ie);
		
		t->SetBranchAddress(prefix+"numberOfMDTHits",            &muid_numberOfMdtHits);
		t->SetBranchAddress(prefix+"numberOfTGCPhiHits",         &muid_numberOfTgcPhiHits);
		t->SetBranchAddress(prefix+"numberOfTGCEtaHits",         &muid_numberOfTgcEtaHits);
		t->SetBranchAddress(prefix+"numberOfCSCPhiHits",         &muid_numberOfCscPhiHits);
		t->SetBranchAddress(prefix+"numberOfCSCEtaHits",         &muid_numberOfCscEtaHits);
		t->SetBranchAddress(prefix+"numberOfRPCPhiHits",         &muid_numberOfRpcPhiHits);
		t->SetBranchAddress(prefix+"numberOfRPCEtaHits",         &muid_numberOfRpcEtaHits);
		t->SetBranchAddress(prefix+"numberOfCscEtaHoles",        &muid_numberOfCscEtaHoles);
		t->SetBranchAddress(prefix+"numberOfCscPhiHoles",        &muid_numberOfCscPhiHoles);
		t->SetBranchAddress(prefix+"numberOfRpcEtaHoles",        &muid_numberOfRpcEtaHoles);
		t->SetBranchAddress(prefix+"numberOfRpcPhiHoles",        &muid_numberOfRpcPhiHoles);
		t->SetBranchAddress(prefix+"numberOfMdtHoles",           &muid_numberOfMdtHoles);
		t->SetBranchAddress(prefix+"numberOfTgcEtaHoles",        &muid_numberOfTgcEtaHoles);
		t->SetBranchAddress(prefix+"numberOfTgcPhiHoles",        &muid_numberOfTgcPhiHoles);
		t->SetBranchAddress(prefix+"numberOfOutliersOnTrack",    &muid_numberOfOutliersOnTrack);
		t->SetBranchAddress(prefix+"standardDeviationOfChi2OS",  &muid_standardDeviationOfChi2OS);
		
		t->SetBranchAddress(prefix+"sectors",                &muid_sectors);
		t->SetBranchAddress(prefix+"nprecisionHits",         &muid_nprecisionHits);
		t->SetBranchAddress(prefix+"nphiLayers",             &muid_nphiLayers);
		t->SetBranchAddress(prefix+"netaPhiLayers",          &muid_netaPhiLayers);
		t->SetBranchAddress(prefix+"nprecisionHoles",        &muid_nprecisionHoles);
		t->SetBranchAddress(prefix+"netaTriggerHoleLayers",  &muid_netaTriggerHoleLayers);
		t->SetBranchAddress(prefix+"nphiHoleLayers",         &muid_nphiHoleLayers);
		t->SetBranchAddress(prefix+"nprecisionOutliers",     &muid_nprecisionOutliers);
		
		t->SetBranchAddress(prefix+"stationName",            &muid_stationName);
		t->SetBranchAddress(prefix+"phiName",                &muid_phiName);
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"trkMS_px", 0);
			t->SetBranchStatus(prefix+"trkComb_px", 0);
			t->SetBranchStatus(prefix+"trkMS_py", 0);
			t->SetBranchStatus(prefix+"trkComb_py", 0);
			t->SetBranchStatus(prefix+"trkMS_pz", 0);
			t->SetBranchStatus(prefix+"trkComb_pz", 0);
			t->SetBranchStatus(prefix+"trkMS_e", 0);
			t->SetBranchStatus(prefix+"trkComb_e", 0);
			t->SetBranchStatus(prefix+"cotTh", 0);
			t->SetBranchStatus(prefix+"trkMS_cotTh", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_cotTh", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_cotTh", 0);
			t->SetBranchStatus(prefix+"trkComb_cotTh", 0);
			t->SetBranchStatus(prefix+"trkMS_pt", 0);
			t->SetBranchStatus(prefix+"trkComb_pt", 0);
			t->SetBranchStatus(prefix+"trkMS_eta", 0);
			t->SetBranchStatus(prefix+"trkComb_eta", 0);
			t->SetBranchStatus(prefix+"trkMS_phi", 0);
			t->SetBranchStatus(prefix+"trkComb_phi", 0);
			t->SetBranchStatus(prefix+"trkMS_chi2", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_chi2", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_chi2", 0);
			t->SetBranchStatus(prefix+"trkMS_nDoF", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_nDoF", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_nDoF", 0);
			t->SetBranchStatus(prefix+"etaAtMsPivot", 0);
			t->SetBranchStatus(prefix+"phiAtMsPivot", 0);
			t->SetBranchStatus(prefix+"sigmaEtaAtMsPivot", 0);
			t->SetBranchStatus(prefix+"sigmaPhiAtMsPivot", 0);
			t->SetBranchStatus(prefix+"muonTrackIndex", 0);
		}
	}
	
	
	
	
	
	if(tType=="CALO_MUONS")
	{
		TString prefix = "CALO_MUONS_MU_";
		
		calo_muons_ptcone10 = 0;
		calo_muons_ptcone20 = 0;
		calo_muons_ptcone30 = 0;
		calo_muons_ptcone40 = 0;
		calo_muons_etcone10 = 0;
		calo_muons_etcone20 = 0;
		calo_muons_etcone30 = 0;
		calo_muons_etcone40 = 0;
		calo_muons_eLoss = 0;
		calo_muons_eLossErr = 0;
		calo_muons_charge = 0;
		calo_muons_px = 0;
		calo_muons_py = 0;
		calo_muons_pz = 0;
		calo_muons_e = 0;
		calo_muons_pt = 0;
		calo_muons_eta = 0;
		calo_muons_phi = 0;
		calo_muons_sctangsig = 0;
		calo_muons_sctngbsig = 0;
		calo_muons_pbalsig = 0;
		calo_muons_isCombined = 0;
		calo_muons_hasInDetTrackParticle = 0;
		calo_muons_inDetTrackIndex = 0;
		calo_muons_index = 0;
		calo_muons_chi2 = 0;
		calo_muons_ndf = 0;
		calo_muons_author = 0;
		calo_muons_matchchi2 = 0;
		calo_muons_matchchi2ndf = 0;
		calo_muons_phi_me = 0;
		calo_muons_phi_ie = 0;
		calo_muons_eta_me = 0;
		calo_muons_eta_ie = 0;
		calo_muons_pt_me = 0;
		calo_muons_pt_ie = 0;
		calo_muons_px_me = 0;
		calo_muons_px_ie = 0;
		calo_muons_py_me = 0;
		calo_muons_py_ie = 0;
		calo_muons_pz_me = 0;
		calo_muons_pz_ie = 0;
		calo_muons_e_me = 0;
		calo_muons_e_ie = 0;
		calo_muons_isLoose = 0;
		calo_muons_isMedium = 0;
		calo_muons_isTight = 0;
		
		t->SetBranchAddress(prefix+"ptcone10", &calo_muons_ptcone10);
		t->SetBranchAddress(prefix+"ptcone20", &calo_muons_ptcone20);
		t->SetBranchAddress(prefix+"ptcone30", &calo_muons_ptcone30);
		t->SetBranchAddress(prefix+"ptcone40", &calo_muons_ptcone40);
		t->SetBranchAddress(prefix+"etcone10", &calo_muons_etcone10);
		t->SetBranchAddress(prefix+"etcone20", &calo_muons_etcone20);
		t->SetBranchAddress(prefix+"etcone30", &calo_muons_etcone30);
		t->SetBranchAddress(prefix+"etcone40", &calo_muons_etcone40);
		t->SetBranchAddress(prefix+"eLoss",    &calo_muons_eLoss);
		t->SetBranchAddress(prefix+"eLossErr", &calo_muons_eLossErr);
		t->SetBranchAddress(prefix+"charge", &calo_muons_charge);
		t->SetBranchAddress(prefix+"px", &calo_muons_px);
		t->SetBranchAddress(prefix+"py", &calo_muons_py);
		t->SetBranchAddress(prefix+"pz", &calo_muons_pz);
		t->SetBranchAddress(prefix+"e", &calo_muons_e);
		t->SetBranchAddress(prefix+"pt", &calo_muons_pt);
		t->SetBranchAddress(prefix+"eta", &calo_muons_eta);
		t->SetBranchAddress(prefix+"phi", &calo_muons_phi);
		t->SetBranchAddress(prefix+"scatteringCurvatureSignificance", &calo_muons_sctangsig);
		t->SetBranchAddress(prefix+"scatteringNeighbourSignificance", &calo_muons_sctngbsig);
		t->SetBranchAddress(prefix+"momentumBalanceSignificance", &calo_muons_pbalsig);
		t->SetBranchAddress(prefix+"isCombinedMuon", &calo_muons_isCombined);
		t->SetBranchAddress(prefix+"hasInDetTrackParticle", &calo_muons_hasInDetTrackParticle);
		t->SetBranchAddress(prefix+"inDetTrackIndex", &calo_muons_inDetTrackIndex);
		t->SetBranchAddress(prefix+"index", &calo_muons_index);
		t->SetBranchAddress(prefix+"trkComb_chi2", &calo_muons_chi2);
		t->SetBranchAddress(prefix+"trkComb_nDoF", &calo_muons_ndf);
		t->SetBranchAddress(prefix+"isLoose",  &calo_muons_isLoose);
		t->SetBranchAddress(prefix+"isMedium", &calo_muons_isMedium);
		t->SetBranchAddress(prefix+"isTight",  &calo_muons_isTight);
		t->SetBranchAddress(prefix+"author", &calo_muons_author);
		t->SetBranchAddress(prefix+"matchChi2", &calo_muons_matchchi2);
		t->SetBranchAddress(prefix+"matchChi2OverDoF", &calo_muons_matchchi2ndf);
		t->SetBranchAddress(prefix+"trkMuonExtr_phi", &calo_muons_phi_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_eta", &calo_muons_eta_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_pt", &calo_muons_pt_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_px", &calo_muons_px_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_py", &calo_muons_py_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_pz", &calo_muons_pz_me);
		t->SetBranchAddress(prefix+"trkMuonExtr_e", &calo_muons_e_me);
		t->SetBranchAddress(prefix+"trkInnerExtr_phi", &calo_muons_phi_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_eta", &calo_muons_eta_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_pt", &calo_muons_pt_ie);		
		t->SetBranchAddress(prefix+"trkInnerExtr_px", &calo_muons_px_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_py", &calo_muons_py_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_pz", &calo_muons_pz_ie);
		t->SetBranchAddress(prefix+"trkInnerExtr_e", &calo_muons_e_ie);
		
		///////////////////////////////////
		// disable unnecessary branches ///
		///////////////////////////////////
		if(skim)
		{
			t->SetBranchStatus(prefix+"trkMS_px", 0);
			t->SetBranchStatus(prefix+"trkComb_px", 0);
			t->SetBranchStatus(prefix+"trkMS_py", 0);
			t->SetBranchStatus(prefix+"trkComb_py", 0);
			t->SetBranchStatus(prefix+"trkMS_pz", 0);
			t->SetBranchStatus(prefix+"trkComb_pz", 0);
			t->SetBranchStatus(prefix+"trkMS_e", 0);
			t->SetBranchStatus(prefix+"trkComb_e", 0);
			t->SetBranchStatus(prefix+"cotTh", 0);
			t->SetBranchStatus(prefix+"trkMS_cotTh", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_cotTh", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_cotTh", 0);
			t->SetBranchStatus(prefix+"trkComb_cotTh", 0);
			t->SetBranchStatus(prefix+"trkMS_pt", 0);
			t->SetBranchStatus(prefix+"trkComb_pt", 0);
			t->SetBranchStatus(prefix+"trkMS_eta", 0);
			t->SetBranchStatus(prefix+"trkComb_eta", 0);
			t->SetBranchStatus(prefix+"trkMS_phi", 0);
			t->SetBranchStatus(prefix+"trkComb_phi", 0);
			t->SetBranchStatus(prefix+"trkMS_chi2", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_chi2", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_chi2", 0);
			t->SetBranchStatus(prefix+"trkMS_nDoF", 0);
			t->SetBranchStatus(prefix+"trkMuonExtr_nDoF", 0);
			t->SetBranchStatus(prefix+"trkInnerExtr_nDoF", 0);
			t->SetBranchStatus(prefix+"etaAtMsPivot", 0);
			t->SetBranchStatus(prefix+"phiAtMsPivot", 0);
			t->SetBranchStatus(prefix+"sigmaEtaAtMsPivot", 0);
			t->SetBranchStatus(prefix+"sigmaPhiAtMsPivot", 0);
			t->SetBranchStatus(prefix+"muonTrackIndex", 0);
		}
	}
}
void addChain(TString name, TMapTSP2TCHAIN& chains, TString mastertree="MUONS_TRIPLET")
{
	chains.insert(make_pair(name, new TChain(mastertree)));
}
void makeChain(TString path, TChain* chain)
{
	void* dir = gSystem->OpenDirectory(path);
	if(!dir) _FATAL("Directory: "+(string)path+" could not be opened");
	
	TString sample = gSystem->GetDirEntry(dir);
	TString fullpath = path+"/"+sample+"/";
	dir = gSystem->OpenDirectory(fullpath);
	if(!dir) _FATAL("Directory: "+(string)fullpath+" could not be opened");
	
	Int_t nFiles = 0;
	TString pattern = ".root";

	char* ent;
	vector<TString> vfiles;
	while((ent=(char*)gSystem->GetDirEntry(dir)))
	{
		TString fn = Form("%s%s",fullpath.Data(),ent);
		if(fn.Contains(pattern))
		{
			// FileStat_t st;
			// if(!gSystem->GetPathInfo(fn,st) && R_ISREG(st.fMode)) chain->Add(fn);
			
			// chain->Add(fn);
			vfiles.push_back(fn);
			nFiles++;
			
		}
	}
	if(!nFiles) _FATAL("Directory: "+(string)fullpath+" contains no files with "+(string)pattern+" pattern in their name");
	
	// add the files in chunks if necessary
	Bool_t isValid = (nChunks<=nFiles);
	Bool_t ignoresplit = (splitstr=="" || splitstr=="0:0");
	if(!isValid) _FATAL("nChunks>nFiles...");
	lChunk     = (!ignoresplit) ? ceil((Float_t)nFiles/(Float_t)nChunks) : -1; // chunk size or "length" (rounded up)
	iFirstFile = (!ignoresplit) ? lChunk*(iChunk-1) : -1;
	iLastFile  = (!ignoresplit) ? lChunk*(iChunk-1)+(lChunk-1) : -1;
	if(!chainWasPrinted)
	{
		cout << "--------------------- splitter --------------------" << endl;
		cout << "ignoresplit = " << ignoresplit << endl;
		cout << "nFiles      = " << nFiles << endl;
		for(unsigned int f=0 ; f<vfiles.size() ; f++) cout << "  iFile[" << f << "]: " << vfiles[f] << endl;
		cout << "lChunk      = " << lChunk << endl;
		cout << "iFirstFile  = " << iFirstFile << endl;
		cout << "iLastFile   = " << iLastFile << endl;
		cout << "---------------------------------------------------" << endl;
		chainWasPrinted = true;
	}
	for(Int_t f=0 ; f<(Int_t)vfiles.size() ; f++)
	{	
		if(!ignoresplit && f<iFirstFile) continue;
		if(!ignoresplit && f>iLastFile)  continue;
		chain->Add(vfiles[f]);
	}
	
	// Example of split alg:
	// [1] 37 files
	// [2] splitstr = "10:3"
	// ====> split 37 files into 10 chuncks (subjobs) and now you are working on the 3rd chunck (subjob).
	// ====> number of elements in one chunk: chunkSize = ceil(nFiles/nChunks) = ceil(37/10) = ceil(3.7) = 4
	// ====> chunks: 0-3, 4-7, 8-11, 12-15, 16-19, 20-23, 24-27, 28-31, 32-35, 36
	// ======> first file of the 3rd chunk has index: 8 = lChunk*(iChunk-1) = 4*(3-1) = 8
	// ======> last file has index: lChunk*(iChunk-1) + lChunk-1 = 8 + 3 = 11
	
	gSystem->FreeDirectory(dir);
}
void setFriends(TString cname, TMapTSP2TCHAIN& cfriends)
{
	TString name = "";
	for(TMapTSP2TCHAIN::iterator it=cfriends.begin() ; it!=cfriends.end() ; it++) delete it->second;
	cfriends.clear(); // no need to keep them all (only per ntuple)
	
	bool isdata = isData(cname);
	TString mcdata = (isdata) ? "data/" : "mc/";
	TString master = mastername+"/";
	TString suffix = "/*.root*";
	
	name="MUONS";            cfriends.insert(make_pair(name, new TChain(name)));
	name="MUID";             cfriends.insert(make_pair(name, new TChain(name)));
	name="CALO_MUONS";       cfriends.insert(make_pair(name, new TChain(name)));
	name="PRIVX";            cfriends.insert(make_pair(name, new TChain(name)));
	name="SEL_TRACKS";       cfriends.insert(make_pair(name, new TChain(name)));
	name="AUX";              cfriends.insert(make_pair(name, new TChain(name)));
	name="JET";              cfriends.insert(make_pair(name, new TChain(name)));
	name="MET";              cfriends.insert(make_pair(name, new TChain(name)));
	name="TRIG";             cfriends.insert(make_pair(name, new TChain(name)));
	name="physics";          cfriends.insert(make_pair(name, new TChain(name)));
	if(!isdata) { name="MC"; cfriends.insert(make_pair(name, new TChain(name))); }
	for(TMapiTS::iterator it=tpmus.begin() ; it!=tpmus.end() ; ++it)
	{
		if(it->first<=MUID) continue; // add only the TPs
		name="TRKPRT_"+it->second;     cfriends.insert(make_pair(name, new TChain(name)));
	}
	
	for(TMapTSP2TCHAIN::iterator it=cfriends.begin() ; it!=cfriends.end() ; it++)
	{
		name = it->first;
		TString path = basepath+mcdata+master+cname+"/";
		if(skim) path.ReplaceAll(master,"");
		makeChain(path,it->second);	
	}
}
void setAllBranches(TChain* c, TMapTSP2TCHAIN& cfriends, TString mastertree="MUONS_TRIPLET")
{
	setBranches(mastertree,c);
	
	for(TMapTSP2TCHAIN::iterator it=cfriends.begin() ; it!=cfriends.end() ; it++)
	{
		bool isNULL = (!it->second);
		if(isNULL) _FATAL("Chain friend "+(string)it->first+" is NULL");
		setBranches(it->first,it->second);
	}
}
void addFriends(TChain* c, TMapTSP2TCHAIN& cfriends)
{	
	for(TMapTSP2TCHAIN::iterator it=cfriends.begin() ; it!=cfriends.end() ; it++)
	{
		bool isNULL = (!it->second);
		if(isNULL) _FATAL("Chain friend "+(string)it->first+" is NULL");
		c->AddFriend(it->second);
	}
}
void prepareChains(TString name, TString mastertree, TMapTSP2TCHAIN& chains, TMapTSP2TCHAIN& chainfriends, TFile* fout, TMapTSP2TTREE& otrees)
{
	TString mcdata = (isData(name)) ? "data" : "mc";
	TString path = basepath+mcdata+"/"+mastername+"/"+name+"/";
	if(skim) path.ReplaceAll(mastername+"/","");
	addChain(name,chains,mastertree);
	chainWasPrinted = false;
	makeChain(path,chains[name]);	
	setFriends(name,chainfriends);
	setAllBranches(chains[name],chainfriends,mastertree);
	addFriends(chains[name],chainfriends);
	if(otrees.size()==0 && skim && fout) initOutTrees(fout,otrees,chains,chainfriends);
}



void reBlindAllMassHists(TMapTSP2TH1& histos, double mMin, double mMax)
{
	for(TMapTSP2TH1::iterator it=histos.begin() ; it!=histos.end() ; ++it)
	{
		TString name = it->first;
		if(!isData(name))          continue;
		if(!name.Contains("m3mu")) continue;
		Int_t bMin = it->second->FindBin(mMin);
		Int_t bMax = it->second->FindBin(mMax);
		for(Int_t b=bMin ; b<=bMax ; b++)
		{
			if(it->second->GetBinContent(b)>0)
			{
				it->second->SetBinContent(b,0);
				it->second->SetBinError(b,0);
			}
		}
	}
}
void reBlindAllMassHists(TMapTSP2TH2& histos, double mMin, double mMax)
{
	for(TMapTSP2TH2::iterator it=histos.begin() ; it!=histos.end() ; ++it)
	{
		TString name = it->first;
		if(!isData(name))          continue;
		if(!name.Contains("m3mu")) continue;
		TString xtitle = it->second->GetXaxis()->GetTitle();
		TString ytitle = it->second->GetYaxis()->GetTitle();
		Bool_t isX = (xtitle.Contains("m_{3body}") || xtitle.Contains("m(3#mu)"));
		if(isX)
		{
			Int_t bMin = it->second->GetXaxis()->FindBin(mMin);
			Int_t bMax = it->second->GetXaxis()->FindBin(mMax);
			Int_t nY = it->second->GetNbinsY();
			for(Int_t b=bMin ; b<=bMax ; b++)
			{
				for(Int_t y=1 ; y<=nY ; y++)
				{
					if(it->second->GetBinContent(b,y)>0)
					{
						it->second->SetBinContent(b,y,0);
						it->second->SetBinError(b,y,0);
					}
				}
			}
		}
		else
		{
			Int_t bMin = it->second->GetYaxis()->FindBin(mMin);
			Int_t bMax = it->second->GetYaxis()->FindBin(mMax);
			Int_t nX = it->second->GetNbinsX();
			for(Int_t b=bMin ; b<=bMax ; b++)
			{
				for(Int_t x=1 ; x<=nX ; x++)
				{
					if(it->second->GetBinContent(x,b)>0)
					{
						it->second->SetBinContent(x,b,0);
						it->second->SetBinError(x,b,0);
					}
				}
			}
		}
	}
}
bool reBlind(TH1* h, double m, double mMin, double mMax)
{
	TString name = h->GetName();
	if(!isData(name)) return false;
	Int_t bMin = h->FindBin(mMin);
	Int_t bMax = h->FindBin(mMax);
	double xMin = h->GetBinLowEdge(bMin);
	double xMax = h->GetBinLowEdge(bMax)+h->GetBinWidth(bMax);
	if(m>=xMin && m<xMax) return true;
	return false;
}



bool isBtaggedEvent(int& btagj1, int& btagj2)
{
	btagj1 = -1;
	btagj2 = -1;
	float btagwgt1 = -999;
	float btagwgt2 = -999;
	// for(unsigned int i=0 ; i<jets_pt->size() ; ++i)
	// {
	// 	bool isUgly = (jets_isUgly->at(i));
	// 	if(isUgly)  continue;
	// 	
	// 	float thisbtagwgt = jets_flwgt_MV1->at(i);
	// 
	// 	// highest and second highest flavor weight MV1
	// 	if(thisbtagwgt>btagwgt1)                         { btagwgt1=thisbtagwgt;  btagj1=i; } // highest flavor weight
	// 	if(thisbtagwgt>btagwgt2 && thisbtagwgt<btagwgt1) { btagwgt2=thisbtagwgt;  btagj2=i; } // second highest flavor weight
	// }
	
	for(int i=0 ; i<AntiKt4LCTopoLets_n ; i++)
	{
		if(!isGoodJet(i,JetQuality))  continue;
		
		float thisbtagwgt = AntiKt4LCTopoLets_flavor_weight_MV1->at(i);

		// highest and second highest flavor weight MV1
		if(thisbtagwgt>btagwgt1)                         { btagwgt1=thisbtagwgt;  btagj1=i; } // highest flavor weight
		if(thisbtagwgt>btagwgt2 && thisbtagwgt<btagwgt1) { btagwgt2=thisbtagwgt;  btagj2=i; } // second highest flavor weight
	}
	
	//////////////////////////////////////////////////////////////
	// make the standard tagging decision ////////////////////////
	// bool isBtag  = (btagj1!=-1 && isStandardBtag(btagwgt1)); //
	bool isBtag  = (btagj1!=-1 && isTightBtag(btagwgt1)); ////////
	if(!isBtag) { btagj1 = -1; btagj2 = -1; } ////////////////////
	//////////////////////////////////////////////////////////////

	return isBtag;
}

void bTagging(TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight, TLorentzVector& pMuSum, unsigned int iBtaggingTriplet)
{	
	/////////////////////////////////////
	// this is not the best requirement,
	// because now the histos are filled
	// only for the 1st good vertex that
	// went up to the b-tagging call so
	// other 3mu systems in the same evt
	// are ignored at the moment.
	// However, this is a very small effect  
	if(iBtaggingTriplet>1) return; //////
	/////////////////////////////////////
	
	float dflt = -1.e20;
	float pTjet1  = dflt;
	float pTjet2  = dflt;
	int   ipTjet1 = -1;
	int   ipTjet2 = -1;
	float btagwgt1 = dflt; int btagj1 = -1;
	float btagwgt2 = dflt; int btagj2 = -1;
	float pTjet1_noOvrlp  = dflt;
	float pTjet2_noOvrlp  = dflt;
	int   ipTjet1_noOvrlp = -1;
	int   ipTjet2_noOvrlp = -1;
	float btagwgt1_noOvrlp = dflt; int btagj1_noOvrlp = -1;
	float btagwgt2_noOvrlp = dflt; int btagj2_noOvrlp = -1;
	float pTjet1_stdMV1  = dflt;
	float pTjet2_stdMV1  = dflt;
	int   ipTjet1_stdMV1 = -1;
	int   ipTjet2_stdMV1 = -1;
	float btagwgt1_stdMV1 = dflt; int btagj1_stdMV1 = -1;
	float btagwgt2_stdMV1 = dflt; int btagj2_stdMV1 = -1;
	int nJet    = 0;
	int nJet10  = 0;
	int nJet20  = 0;
	int nBtag   = 0;
	int nBtag10 = 0;
	int nBtag20 = 0;
	TMapTSf maxflwgt;
	maxflwgt.insert(make_pair("MV1", -1.e20));
	
	int iOvlpJet = bestOverlapingJet(pMuSum);
	
	// for(int i=0 ; i<(int)jets_pt->size() ; ++i)
	// {
	// 	bool isUgly = (jets_isUgly->at(i));
	// 	if(isUgly)  continue;
	for(int i=0 ; i<AntiKt4LCTopoLets_n ; i++)
	{
		if(!isGoodJet(i,JetQuality))  continue;
		
		nJet++;

		bool isBtagS = isStandardBtag(AntiKt4LCTopoLets_flavor_weight_MV1->at(i));
		bool isBtagT = isTightBtag(AntiKt4LCTopoLets_flavor_weight_MV1->at(i));
		if(isBtagT) nBtag++;
		
		// different pT thresholds
		if(AntiKt4LCTopoLets_pt->at(i)>10*GeV2MeV) nJet10++;
		if(AntiKt4LCTopoLets_pt->at(i)>20*GeV2MeV) nJet20++;
		if(isBtagT)
		{
			if(AntiKt4LCTopoLets_pt->at(i)>10*GeV2MeV) nBtag10++;
			if(AntiKt4LCTopoLets_pt->at(i)>20*GeV2MeV) nBtag20++;
		}
		
		// some histos
		histos[name+"_jet_E"]->Fill(AntiKt4LCTopoLets_E->at(i),weight);
		histos[name+"_jet_pt"]->Fill(AntiKt4LCTopoLets_pt->at(i),weight);
		histos[name+"_jet_m"]->Fill(AntiKt4LCTopoLets_m->at(i),weight);
		histos[name+"_jet_eta"]->Fill(AntiKt4LCTopoLets_eta->at(i),weight);
		histos[name+"_jet_phi"]->Fill(AntiKt4LCTopoLets_phi->at(i),weight);
		if(isBtagT)
		{
			histos[name+"_bjet_E"]->Fill(AntiKt4LCTopoLets_E->at(i),weight);
			histos[name+"_bjet_pt"]->Fill(AntiKt4LCTopoLets_pt->at(i),weight);
			histos[name+"_bjet_m"]->Fill(AntiKt4LCTopoLets_m->at(i),weight);
			histos[name+"_bjet_eta"]->Fill(AntiKt4LCTopoLets_eta->at(i),weight);
			histos[name+"_bjet_phi"]->Fill(AntiKt4LCTopoLets_phi->at(i),weight);
		}
		
		// find the different tagger MAX weights
		maxflwgt["MV1"] = (AntiKt4LCTopoLets_flavor_weight_MV1->at(i)>maxflwgt["MV1"]) ? AntiKt4LCTopoLets_flavor_weight_MV1->at(i) : maxflwgt["MV1"];
		
		float thisbtagwgt = AntiKt4LCTopoLets_flavor_weight_MV1->at(i);
		float thisjetpt   = AntiKt4LCTopoLets_pt->at(i);
		
		// leading and subleading jet pT
		if(thisjetpt>=pTjet1)
		{
			if(pTjet1!=dflt) { pTjet2=pTjet1;  ipTjet2=ipTjet1; } // first change subleading jet to previous leading jet
			pTjet1=thisjetpt; ipTjet1=i;                          // now change the leading jet
		}
		if(thisjetpt>=pTjet2 && thisjetpt<pTjet1) { pTjet2=thisjetpt;  ipTjet2=i; } // subleading jet
		
		// highest and second highest flavor weight MV1
		if(thisbtagwgt>btagwgt1)
		{
			if(btagwgt1!=dflt) { btagwgt2=btagwgt1;  btagj2=btagj1; } // first change 2nd highest weight to previous highest weight
			btagwgt1=thisbtagwgt; btagj1=i;                           // now change the highest flavor weight
		}
		if(thisbtagwgt>btagwgt2 && thisbtagwgt<btagwgt1) { btagwgt2=thisbtagwgt;  btagj2=i; } // second highest flavor weight
		
		// do the same but ignore the best overlapping jet
		if(i!=iOvlpJet)
		{
			// leading and subleading jet pT
			if(thisjetpt>=pTjet1_noOvrlp)
			{
				if(pTjet1_noOvrlp!=dflt) { pTjet2_noOvrlp=pTjet1_noOvrlp;  ipTjet2_noOvrlp=ipTjet1_noOvrlp; } // first change subleading jet to previous leading jet
				pTjet1_noOvrlp=thisjetpt;  ipTjet1_noOvrlp=i; // leading jet
			}
			if(thisjetpt>=pTjet2_noOvrlp && thisjetpt<pTjet1_noOvrlp) { pTjet2_noOvrlp=thisjetpt;  ipTjet2_noOvrlp=i; } // subleading jet
			
			// highest and second highest flavor weight MV1
			if(thisbtagwgt>=btagwgt1_noOvrlp)
			{
				if(btagwgt1_noOvrlp!=dflt) { btagwgt2_noOvrlp=btagwgt1_noOvrlp;  btagj2_noOvrlp=btagj1_noOvrlp; } // first change 2nd highest weight to previous highest weight
				btagwgt1_noOvrlp=thisbtagwgt;  btagj1_noOvrlp=i; // highest flavor weight
			}
			if(thisbtagwgt>=btagwgt2_noOvrlp && thisbtagwgt<btagwgt1_noOvrlp) { btagwgt2_noOvrlp=thisbtagwgt;  btagj2_noOvrlp=i; } // second highest flavor weight
		}
		
		// do the same but only for standard MV1 b-tagged jets
		if(isBtagS)
		{
			// leading and subleading jet pT
			if(thisjetpt>=pTjet1_stdMV1)
			{
				if(pTjet1_stdMV1!=dflt) { pTjet2_stdMV1=pTjet1_stdMV1;  ipTjet2_stdMV1=ipTjet1_stdMV1; } // first change subleading jet to previous leading jet
				pTjet1_stdMV1=thisjetpt;  ipTjet1_stdMV1=i; // leading jet
			}
			if(thisjetpt>=pTjet2_stdMV1 && thisjetpt<pTjet1_stdMV1) { pTjet2_stdMV1=thisjetpt;  ipTjet2_stdMV1=i; } // subleading jet
			
			// highest and second highest flavor weight MV1
			if(thisbtagwgt>=btagwgt1_stdMV1)
			{
				if(btagwgt1_stdMV1!=dflt) { btagwgt2_stdMV1=btagwgt1_stdMV1;  btagj2_stdMV1=btagj1_stdMV1; } // first change 2nd highest weight to previous highest weight
				btagwgt1_stdMV1=thisbtagwgt;  btagj1_stdMV1=i; // highest flavor weight
			}
			if(thisbtagwgt>=btagwgt2_stdMV1 && thisbtagwgt<btagwgt1_stdMV1) { btagwgt2_stdMV1=thisbtagwgt;  btagj2_stdMV1=i; } // second highest flavor weight
		}
	}
	
	// counters
	histos[name+"_jet_n"]->Fill(nJet,weight);
	histos[name+"_jet_n10"]->Fill(nJet10,weight);
	histos[name+"_jet_n20"]->Fill(nJet20,weight);
	histos[name+"_bjet_n"]->Fill(nBtag,weight);
	histos[name+"_bjet_n10"]->Fill(nBtag10,weight);
	histos[name+"_bjet_n20"]->Fill(nBtag20,weight);
	histos2[name+"_nBtag_vs_nJet"]->Fill(nJet,nBtag,weight);	

	// the jet with highest flavor weight per tagger
	histos[name+"_jet_maxFlvWgt_MV1"]->Fill(maxflwgt["MV1"],weight);
	
	// the pT leading jet flavor weight per tagger
	if(ipTjet1!=-1)
	{		
		histos[name+"_jet_leadingJetFlvWgt_MV1"]->Fill(AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet1), weight);
	}
	// the pT subleading jet flavor weight per tagger
	if(ipTjet2!=-1)
	{		
		histos[name+"_jet_subleadingJetFlvWgt_MV1"]->Fill(AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet2), weight);
	}
		
	TLorentzVector lvJet, lvJet1;
	float dR   = -1;
	float dphi = -1;
	float wj1  = -1; 
	float wj2  = -1; 
	if(ipTjet1!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet1),
						   AntiKt4LCTopoLets_eta->at(ipTjet1),
						   AntiKt4LCTopoLets_phi->at(ipTjet1),
						   AntiKt4LCTopoLets_E->at(ipTjet1));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muLeadingJet"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muLeadingJet"]->Fill(dphi,weight);
	}
	if(ipTjet2!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet2),
						   AntiKt4LCTopoLets_eta->at(ipTjet2),
						   AntiKt4LCTopoLets_phi->at(ipTjet2),
						   AntiKt4LCTopoLets_E->at(ipTjet2));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muSubleadingJet"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muSubleadingJet"]->Fill(dphi,weight);
		
		wj1 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet1);
		wj2 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet2);
		histos2[name+"_jet_wj1_vs_wj2"]->Fill(wj2,wj1,weight);
		lvJet1.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet1),
						   AntiKt4LCTopoLets_eta->at(ipTjet1),
						   AntiKt4LCTopoLets_phi->at(ipTjet1),
						   AntiKt4LCTopoLets_E->at(ipTjet1));
		dphi = fabs(lvJet.DeltaPhi(lvJet1));
		histos[name+"_bjet_dPhi_jet1_jet2"]->Fill(dphi,weight);
	}
	if(btagj1!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj1),
						   AntiKt4LCTopoLets_eta->at(btagj1),
						   AntiKt4LCTopoLets_phi->at(btagj1),
						   AntiKt4LCTopoLets_E->at(btagj1));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muBjet_btagwgt1"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muBjet_btagwgt1"]->Fill(dphi,weight);
		histos2[name+"_bjet_dR3muBjet_btagwgt1_vs_btagwgt1"]->Fill(btagwgt1,dR,weight);
		histos2[name+"_bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1"]->Fill(btagwgt1,dphi,weight);
	}
	if(btagj2!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj2),
						   AntiKt4LCTopoLets_eta->at(btagj2),
						   AntiKt4LCTopoLets_phi->at(btagj2),
						   AntiKt4LCTopoLets_E->at(btagj2));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muBjet_btagwgt2"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muBjet_btagwgt2"]->Fill(dphi,weight);
		histos2[name+"_bjet_dR3muBjet_btagwgt2_vs_btagwgt2"]->Fill(btagwgt2,dR,weight);
		histos2[name+"_bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2"]->Fill(btagwgt2,dphi,weight);
		lvJet1.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj1),
						   AntiKt4LCTopoLets_eta->at(btagj1),
						   AntiKt4LCTopoLets_phi->at(btagj1),
						   AntiKt4LCTopoLets_E->at(btagj1));
		dphi = fabs(lvJet.DeltaPhi(lvJet1));
		histos[name+"_bjet_dPhi_btagwgt1_btagwgt2"]->Fill(dphi,weight);
		histos2[name+"_bjet_dPhib1b2_vs_b1"]->Fill(btagwgt1,dphi,weight);
		histos2[name+"_bjet_dPhib1b2_vs_b2"]->Fill(btagwgt2,dphi,weight);
		histos2[name+"_bjet_wb1_vs_wb2"]->Fill(btagwgt2,btagwgt1,weight);
	}
		
	if(ipTjet1_noOvrlp!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet1_noOvrlp),
						   AntiKt4LCTopoLets_eta->at(ipTjet1_noOvrlp),
						   AntiKt4LCTopoLets_phi->at(ipTjet1_noOvrlp),
						   AntiKt4LCTopoLets_E->at(ipTjet1_noOvrlp));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muLeadingJet_noOvrlp"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muLeadingJet_noOvrlp"]->Fill(dphi,weight);
	}
	if(ipTjet2_noOvrlp!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet2_noOvrlp),
						   AntiKt4LCTopoLets_eta->at(ipTjet2_noOvrlp),
						   AntiKt4LCTopoLets_phi->at(ipTjet2_noOvrlp),
						   AntiKt4LCTopoLets_E->at(ipTjet2_noOvrlp));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muSubleadingJet_noOvrlp"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muSubleadingJet_noOvrlp"]->Fill(dphi,weight);
		
		wj1 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet1_noOvrlp);
		wj2 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet2_noOvrlp);
		histos2[name+"_jet_wj1_vs_wj2_noOvrlp"]->Fill(wj2,wj1,weight);
		lvJet1.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet1_noOvrlp),
						   AntiKt4LCTopoLets_eta->at(ipTjet1_noOvrlp),
						   AntiKt4LCTopoLets_phi->at(ipTjet1_noOvrlp),
						   AntiKt4LCTopoLets_E->at(ipTjet1_noOvrlp));
		dphi = fabs(lvJet.DeltaPhi(lvJet1));
		histos[name+"_bjet_dPhi_jet1_jet2_noOvrlp"]->Fill(dphi,weight);
	}
	if(btagj1_noOvrlp!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj1_noOvrlp),
						   AntiKt4LCTopoLets_eta->at(btagj1_noOvrlp),
						   AntiKt4LCTopoLets_phi->at(btagj1_noOvrlp),
						   AntiKt4LCTopoLets_E->at(btagj1_noOvrlp));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muBjet_btagwgt1_noOvrlp"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muBjet_btagwgt1_noOvrlp"]->Fill(dphi,weight);
		histos2[name+"_bjet_dR3muBjet_btagwgt1_vs_btagwgt1_noOvrlp"]->Fill(btagwgt1_noOvrlp,dR,weight);
		histos2[name+"_bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_noOvrlp"]->Fill(btagwgt1_noOvrlp,dphi,weight);
	}
	if(btagj2_noOvrlp!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj2_noOvrlp),
						   AntiKt4LCTopoLets_eta->at(btagj2_noOvrlp),
						   AntiKt4LCTopoLets_phi->at(btagj2_noOvrlp),
						   AntiKt4LCTopoLets_E->at(btagj2_noOvrlp));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muBjet_btagwgt2_noOvrlp"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muBjet_btagwgt2_noOvrlp"]->Fill(dphi,weight);
		histos2[name+"_bjet_dR3muBjet_btagwgt2_vs_btagwgt2_noOvrlp"]->Fill(btagwgt2_noOvrlp,dR,weight);
		histos2[name+"_bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_noOvrlp"]->Fill(btagwgt2_noOvrlp,dphi,weight);
		lvJet1.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj1_noOvrlp),
						   AntiKt4LCTopoLets_eta->at(btagj1_noOvrlp),
						   AntiKt4LCTopoLets_phi->at(btagj1_noOvrlp),
						   AntiKt4LCTopoLets_E->at(btagj1_noOvrlp));
		dphi = fabs(lvJet.DeltaPhi(lvJet1));
		histos[name+"_bjet_dPhi_btagwgt1_btagwgt2_noOvrlp"]->Fill(dphi,weight);
		histos2[name+"_bjet_dPhib1b2_vs_b1_noOvrlp"]->Fill(btagwgt1_noOvrlp,dphi,weight);
		histos2[name+"_bjet_dPhib1b2_vs_b2_noOvrlp"]->Fill(btagwgt2_noOvrlp,dphi,weight);
		histos2[name+"_bjet_wb1_vs_wb2_noOvrlp"]->Fill(btagwgt2_noOvrlp,btagwgt1_noOvrlp,weight);
	}
	
	if(ipTjet1_stdMV1!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet1_stdMV1),
						   AntiKt4LCTopoLets_eta->at(ipTjet1_stdMV1),
						   AntiKt4LCTopoLets_phi->at(ipTjet1_stdMV1),
						   AntiKt4LCTopoLets_E->at(ipTjet1_stdMV1));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muLeadingJet_stdMV1"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muLeadingJet_stdMV1"]->Fill(dphi,weight);
	}
	if(ipTjet2_stdMV1!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet2_stdMV1),
						   AntiKt4LCTopoLets_eta->at(ipTjet2_stdMV1),
						   AntiKt4LCTopoLets_phi->at(ipTjet2_stdMV1),
						   AntiKt4LCTopoLets_E->at(ipTjet2_stdMV1));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muSubleadingJet_stdMV1"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muSubleadingJet_stdMV1"]->Fill(dphi,weight);
		
		wj1 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet1_stdMV1);
		wj2 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ipTjet2_stdMV1);
		histos2[name+"_jet_wj1_vs_wj2_stdMV1"]->Fill(wj2,wj1,weight);
		lvJet1.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(ipTjet1_stdMV1),
						   AntiKt4LCTopoLets_eta->at(ipTjet1_stdMV1),
						   AntiKt4LCTopoLets_phi->at(ipTjet1_stdMV1),
						   AntiKt4LCTopoLets_E->at(ipTjet1_stdMV1));
		dphi = fabs(lvJet.DeltaPhi(lvJet1));
		histos[name+"_bjet_dPhi_jet1_jet2_stdMV1"]->Fill(dphi,weight);
	}
	if(btagj1_stdMV1!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj1_stdMV1),
						   AntiKt4LCTopoLets_eta->at(btagj1_stdMV1),
						   AntiKt4LCTopoLets_phi->at(btagj1_stdMV1),
						   AntiKt4LCTopoLets_E->at(btagj1_stdMV1));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muBjet_btagwgt1_stdMV1"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muBjet_btagwgt1_stdMV1"]->Fill(dphi,weight);
		histos2[name+"_bjet_dR3muBjet_btagwgt1_vs_btagwgt1_stdMV1"]->Fill(btagwgt1_stdMV1,dR,weight);
		histos2[name+"_bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_stdMV1"]->Fill(btagwgt1_stdMV1,dphi,weight);
	}
	if(btagj2_stdMV1!=-1)
	{
		lvJet.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj2_stdMV1),
						   AntiKt4LCTopoLets_eta->at(btagj2_stdMV1),
						   AntiKt4LCTopoLets_phi->at(btagj2_stdMV1),
						   AntiKt4LCTopoLets_E->at(btagj2_stdMV1));
		dR   = lvJet.DeltaR(pMuSum);
		dphi = fabs(lvJet.DeltaPhi(pMuSum));
		histos[name+"_bjet_dR3muBjet_btagwgt2_stdMV1"]->Fill(dR,weight);
		histos[name+"_bjet_dPhi3muBjet_btagwgt2_stdMV1"]->Fill(dphi,weight);
		histos2[name+"_bjet_dR3muBjet_btagwgt2_vs_btagwgt2_stdMV1"]->Fill(btagwgt2_stdMV1,dR,weight);
		histos2[name+"_bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_stdMV1"]->Fill(btagwgt2_stdMV1,dphi,weight);
		lvJet1.SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(btagj1_stdMV1),
						   AntiKt4LCTopoLets_eta->at(btagj1_stdMV1),
						   AntiKt4LCTopoLets_phi->at(btagj1_stdMV1),
						   AntiKt4LCTopoLets_E->at(btagj1_stdMV1));
		dphi = fabs(lvJet.DeltaPhi(lvJet1));
		histos[name+"_bjet_dPhi_btagwgt1_btagwgt2_stdMV1"]->Fill(dphi,weight);
		histos2[name+"_bjet_dPhib1b2_vs_b1_stdMV1"]->Fill(btagwgt1_stdMV1,dphi,weight);
		histos2[name+"_bjet_dPhib1b2_vs_b2_stdMV1"]->Fill(btagwgt2_stdMV1,dphi,weight);
		histos2[name+"_bjet_wb1_vs_wb2_stdMV1"]->Fill(btagwgt2_stdMV1,btagwgt1_stdMV1,weight);
	}
	
	// kinematics of the jet with the maximum flavor weight
	if(btagj1!=-1)
	{
		histos[name+"_bjet_maxFlvWgt_E"]->Fill(AntiKt4LCTopoLets_E->at(btagj1),weight);
		histos[name+"_bjet_maxFlvWgt_pt"]->Fill(AntiKt4LCTopoLets_pt->at(btagj1),weight);
		histos[name+"_bjet_maxFlvWgt_m"]->Fill(AntiKt4LCTopoLets_m->at(btagj1),weight);
		histos[name+"_bjet_maxFlvWgt_eta"]->Fill(AntiKt4LCTopoLets_eta->at(btagj1),weight);
		histos[name+"_bjet_maxFlvWgt_phi"]->Fill(AntiKt4LCTopoLets_phi->at(btagj1),weight);
	}
}


double tripletpTconeXX(int trk1, int trk2, int trk3, double conemargins=0.)
{
	if(conemargins<0.) _FATAL("Cone margins cannot be negative");
	
	double ptconeXX = 0.;
	TVector3 p1,p2,p3,pRef,pTst;
	p1.SetXYZ(trks_px->at(trk1),trks_py->at(trk1),trks_pz->at(trk1));
	p2.SetXYZ(trks_px->at(trk2),trks_py->at(trk2),trks_pz->at(trk2));
	p3.SetXYZ(trks_px->at(trk3),trks_py->at(trk3),trks_pz->at(trk3));
	pRef = p1+p2+p3;
	double dRmax = -1.e20;
	dRmax = (p1.DeltaR(p2)>dRmax) ? p1.DeltaR(p2) : dRmax;
	dRmax = (p1.DeltaR(p3)>dRmax) ? p1.DeltaR(p3) : dRmax;
	dRmax = (p2.DeltaR(p3)>dRmax) ? p2.DeltaR(p3) : dRmax;
	// dRmax += (dRmax<0.1) ? conemargins : 0.;
	dRmax += conemargins;
	for(int trk=0 ; trk<(int)trks_pt->size() ; ++trk)
	{
		if(trk==trk1)                    continue; // ignore triplet tracks
		if(trk==trk2)                    continue; // ignore triplet tracks
		if(trk==trk3)                    continue; // ignore triplet tracks
		if(trks_pt->at(trk)<0.5*GeV2MeV) continue; // ignore low pT tracks (<0.5 GeV)
		pTst.SetXYZ(trks_px->at(trk),trks_py->at(trk),trks_pz->at(trk));
		if(pRef.DeltaR(pTst)<dRmax) ptconeXX += trks_pt->at(trk);
	}
	return ptconeXX;
}
double tripletIsolation(unsigned int vtx, double conemargins=0.)
{
	if(conemargins<0.) _FATAL("Cone margins cannot be negative");
	
	sources src;
	getSrc(vtx,src);
	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
	// match the tracks
	int itrk1 = src.trkIndex[0];
	int itrk2 = src.trkIndex[1];
	int itrk3 = src.trkIndex[2];
	double tripletptconeXX = tripletpTconeXX(itrk1,itrk2,itrk3,conemargins);
	double tripletiso = (vtx_pt->at(vtx)>0.) ? tripletptconeXX/vtx_pt->at(vtx) : 999.;
	return tripletiso;
}

double ptconeXX(double cone, int reftrk, int ignoretrk1=-1, int ignoretrk2=-1)
{
	double trkisolation = 0.;
	TVector3 pRef,pTst;
	pRef.SetXYZ(trks_px->at(reftrk),trks_py->at(reftrk),trks_pz->at(reftrk));
	for(int trk=0 ; trk<(int)trks_pt->size() ; ++trk)
	{
		if(trk==reftrk)                  continue; // ignore itself
		if(trk==ignoretrk1)              continue; // ignore the other 2 muons
		if(trk==ignoretrk2)              continue; // ignore the other 2 muons
		if(trks_pt->at(trk)<0.5*GeV2MeV) continue; // ignore low pT tracks (<0.5 GeV)
		pTst.SetXYZ(trks_px->at(trk),trks_py->at(trk),trks_pz->at(trk));
		if(pRef.DeltaR(pTst)<cone) trkisolation += trks_pt->at(trk);
	}
	return trkisolation;
}
vector<double> correctedIsolation(unsigned int vtx, double cone)
{
	sources src;
	getSrc(vtx,src);
	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
	// match the tracks
	int itrk1 = src.trkIndex[0];
	int itrk2 = src.trkIndex[1];
	int itrk3 = src.trkIndex[2];
	
	_DEBUG("");
	
	if(cone!=0.1 && cone!=0.2 && cone!=0.3 && cone!=0.4) _FATAL("no cone="+_s(cone)+" is available");
	double ptconeXX1 = ptconeXX(cone,itrk1,itrk2,itrk3);
	double ptconeXX2 = ptconeXX(cone,itrk2,itrk1,itrk3);
	double ptconeXX3 = ptconeXX(cone,itrk3,itrk1,itrk2);
	double iso1 = (trks_pt->at(itrk1)!=0.) ? ptconeXX1/trks_pt->at(itrk1) : 999.;
	double iso2 = (trks_pt->at(itrk2)!=0.) ? ptconeXX2/trks_pt->at(itrk2) : 999.;
	double iso3 = (trks_pt->at(itrk3)!=0.) ? ptconeXX3/trks_pt->at(itrk3) : 999.;
	vector<double> muiso; muiso.push_back(iso1); muiso.push_back(iso2); muiso.push_back(iso3);
	return muiso;
}


void tauLxy(unsigned int vtx, TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight,
			double& type1_tau, double& type1_lxy, double& type3_tau, double& type3_lxy, double& type1_a0XY, double& type1_cosThetaXY)
{
	double chi2   = vtx_chi2->at(vtx);
	double ndf    = vtx_ndf->at(vtx);
	double pvalue = TMath::Prob(chi2,ndf);
	vector<int> pv1; findPVindex(pv1,1);
	vector<int> pv3; findPVindex(pv3,3);
	if(pv1.size()<1) _ERROR("PV not found");
	if(pv1.size()>1) _ERROR("nPV>1");
	type1_tau        = 1.e20;
	type1_lxy        = 1.e20;
	type3_tau        = 1.e20;
	type3_lxy        = 1.e20;
	type1_a0XY       = 1.e20;
	type1_cosThetaXY = 1.e20;
	histos[name+"_vtx_npv_type1"]->Fill(pv1.size(),weight);
	histos[name+"_vtx_npv_type3"]->Fill(pv3.size(),weight);
	// type=1 PVs
	for(unsigned int ipv1=0 ; ipv1<pv1.size() ; ipv1++)
	{
		int pvindex = pv1[ipv1];
		type1_tau        = vtx_tau->at(vtx)[pvindex];
		type1_lxy        = vtx_lxy->at(vtx)[pvindex];
		type1_a0XY       = vtx_a0XY->at(vtx)[pvindex];
		type1_cosThetaXY = fabs(vtx_cosThetaXY->at(vtx)[pvindex]);
		histos[name+"_type1_lxy"]->Fill(type1_lxy,weight);
		histos[name+"_type1_tau"]->Fill(type1_tau,weight);
		histos[name+"_type1_a0"]->Fill(vtx_a0->at(vtx)[pvindex],weight);
		histos[name+"_type1_a0XY"]->Fill(vtx_a0XY->at(vtx)[pvindex],weight);
		histos[name+"_type1_cosTh"]->Fill(vtx_cosTheta->at(vtx)[pvindex],weight);
		histos[name+"_type1_cosThXY"]->Fill(vtx_cosThetaXY->at(vtx)[pvindex],weight);
		histos2[name+"_type1_tau_vs_pval"]->Fill(pvalue,type1_tau,weight);
		histos2[name+"_type1_lxy_vs_pval"]->Fill(pvalue,type1_lxy,weight);
	}
	// type=3 PVs
	for(unsigned int ipv3=0 ; ipv3<pv3.size() ; ipv3++)
	{
		int pvindex = pv3[ipv3];
		type3_tau = vtx_tau->at(vtx)[pvindex];
		type3_lxy = vtx_lxy->at(vtx)[pvindex];
		histos[name+"_type3_lxy"]->Fill(type3_lxy,weight);
		histos[name+"_type3_tau"]->Fill(type3_tau,weight);
		histos[name+"_type3_a0"]->Fill(vtx_a0->at(vtx)[pvindex],weight);
		histos[name+"_type3_a0XY"]->Fill(vtx_a0XY->at(vtx)[pvindex],weight);
		histos[name+"_type3_cosTh"]->Fill(vtx_cosTheta->at(vtx)[pvindex],weight);
		histos[name+"_type3_cosThXY"]->Fill(vtx_cosThetaXY->at(vtx)[pvindex],weight);
		histos2[name+"_type3_tau_vs_pval"]->Fill(pvalue,type3_tau,weight);
		histos2[name+"_type3_lxy_vs_pval"]->Fill(pvalue,type3_lxy,weight);
	}
}
void fillHistsPVandSV(unsigned int vtx, TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight)
{
	double chi2   = vtx_chi2->at(vtx);
	double ndf    = vtx_ndf->at(vtx);
	double pvalue = TMath::Prob(chi2,ndf);
	vector<int> pv1; findPVindex(pv1,1);
	vector<int> pv3; findPVindex(pv3,3);
	if(pv1.size()<1) _ERROR("PV not found");
	if(pv1.size()>1) _ERROR("nPV>1");
	histos[name+"_vtx_npv_type1"]->Fill(pv1.size(),weight);
	histos[name+"_vtx_npv_type3"]->Fill(pv3.size(),weight);
	// type=1 PVs
	for(unsigned int ipv1=0 ; ipv1<pv1.size() ; ipv1++)
	{
		int pvindex = pv1[ipv1];
		histos[name+"_type1_lxy"]->Fill(vtx_lxy->at(vtx)[pvindex],weight);
		histos[name+"_type1_tau"]->Fill(vtx_tau->at(vtx)[pvindex],weight);
		histos[name+"_type1_a0"]->Fill(vtx_a0->at(vtx)[pvindex],weight);
		histos[name+"_type1_a0XY"]->Fill(vtx_a0XY->at(vtx)[pvindex],weight);
		histos[name+"_type1_cosTh"]->Fill(vtx_cosTheta->at(vtx)[pvindex],weight);
		histos[name+"_type1_cosThXY"]->Fill(vtx_cosThetaXY->at(vtx)[pvindex],weight);
		histos2[name+"_type1_tau_vs_pval"]->Fill(pvalue,vtx_tau->at(vtx)[pvindex],weight);
		histos2[name+"_type1_lxy_vs_pval"]->Fill(pvalue,vtx_lxy->at(vtx)[pvindex],weight);
	}
	// type=3 PVs
	for(unsigned int ipv3=0 ; ipv3<pv3.size() ; ipv3++)
	{
		int pvindex = pv3[ipv3];
		histos[name+"_type3_lxy"]->Fill(vtx_lxy->at(vtx)[pvindex],weight);
		histos[name+"_type3_tau"]->Fill(vtx_tau->at(vtx)[pvindex],weight);
		histos[name+"_type3_a0"]->Fill(vtx_a0->at(vtx)[pvindex],weight);
		histos[name+"_type3_a0XY"]->Fill(vtx_a0XY->at(vtx)[pvindex],weight);
		histos[name+"_type3_cosTh"]->Fill(vtx_cosTheta->at(vtx)[pvindex],weight);
		histos[name+"_type3_cosThXY"]->Fill(vtx_cosThetaXY->at(vtx)[pvindex],weight);
		histos2[name+"_type3_tau_vs_pval"]->Fill(pvalue,vtx_tau->at(vtx)[pvindex],weight);
		histos2[name+"_type3_lxy_vs_pval"]->Fill(pvalue,vtx_lxy->at(vtx)[pvindex],weight);
	}
}
void fillHistsDoublets(unsigned int vtx, TString name, TString suffix, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight,
					   TLorentzVector pOS1, TLorentzVector pOS2, TLorentzVector pSS,
					   bool foundOS1, bool foundOS2, bool foundSS)
{
	double m3mu  = vtx_mass->at(vtx);
	double pT3mu = vtx_pt->at(vtx);
	if(foundSS) histos[name+"_mSS"+suffix]->Fill(pSS.M(),weight);
	if(foundOS1 && foundOS2)
	{
		histos[name+"_mOS"+suffix]->Fill(pOS1.M(),weight/2.);
		histos[name+"_mOS"+suffix]->Fill(pOS2.M(),weight/2.);
		
		histos2[name+"_mOS2_vs_mOS1"+suffix]->Fill(pOS1.M(),pOS2.M(),weight); /// NEED TO REBLIND !!!
		
		histos2[name+"_m3mu_vs_mOS"+suffix]->Fill(pOS1.M(),m3mu,weight/2.); /// NEED TO REBLIND !!!
		histos2[name+"_m3mu_vs_mOS"+suffix]->Fill(pOS2.M(),m3mu,weight/2.); /// NEED TO REBLIND !!!
		
		histos2[name+"_pT3mu_vs_mOS"+suffix]->Fill(pOS1.M(),pT3mu,weight/2.); /// NEED TO REBLIND !!!
		histos2[name+"_pT3mu_vs_mOS"+suffix]->Fill(pOS2.M(),pT3mu,weight/2.); /// NEED TO REBLIND !!!
		if(foundSS)
		{
			histos2[name+"_mSS_vs_mOS"+suffix]->Fill(pOS1.M(),pSS.M(),weight/2.);
			histos2[name+"_mSS_vs_mOS"+suffix]->Fill(pOS2.M(),pSS.M(),weight/2.);
		}
	}
	else if(foundOS1  && !foundOS2)
	{
		histos[name+"_mOS"+suffix]->Fill(pOS1.M(),weight);
		histos2[name+"_m3mu_vs_mOS"+suffix]->Fill(pOS1.M(),m3mu,weight); /// NEED TO REBLIND !!!
		histos2[name+"_pT3mu_vs_mOS"+suffix]->Fill(pOS1.M(),pT3mu,weight/2.); /// NEED TO REBLIND !!!
		if(foundSS) histos2[name+"_mSS_vs_mOS"+suffix]->Fill(pOS1.M(),pSS.M(),weight);
	}
	else if(!foundOS1 && foundOS2)
	{
		histos[name+"_mOS"+suffix]->Fill(pOS2.M(),weight);
		histos2[name+"_m3mu_vs_mOS"+suffix]->Fill(pOS2.M(),m3mu,weight); /// NEED TO REBLIND !!!
		histos2[name+"_pT3mu_vs_mOS"+suffix]->Fill(pOS2.M(),pT3mu,weight/2.); /// NEED TO REBLIND !!!
		if(foundSS) histos2[name+"_mSS_vs_mOS"+suffix]->Fill(pOS2.M(),pSS.M(),weight);
	}
}
void fillHistsMassPt3mu(unsigned int vtx, TString name, TString suffix, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight, double mBlindMin, double mBlindMax)
{
	if(false) cout << histos2.begin()->first << endl; // stupid
	double m3mu  = vtx_mass->at(vtx);
	double pT3mu = vtx_pt->at(vtx);
	if(!reBlind(histos[name+"_m3mu"+suffix],m3mu,mBlindMin,mBlindMax))           histos[name+"_m3mu"+suffix]->Fill(m3mu,weight);
	if(!reBlind(histos[name+"_m3mu_lin"+suffix],m3mu,mBlindMin,mBlindMax))       histos[name+"_m3mu_lin"+suffix]->Fill(m3mu,weight);
	if(!reBlind(histos[name+"_m3mu_lin_zoom"+suffix],m3mu,mBlindMin,mBlindMax))  histos[name+"_m3mu_lin_zoom"+suffix]->Fill(m3mu,weight);
	if(!reBlind(histos[name+"_m3mu_sigregion"+suffix],m3mu,mBlindMin,mBlindMax)) histos[name+"_m3mu_sigregion"+suffix]->Fill(m3mu,weight);
	histos[name+"_pT3mu"+suffix]->Fill(pT3mu,weight);
}
void fillHistsMVAvars(unsigned int vtx, TString name, TString suffix, TMapTSP2TH1& histos, TMapTSP2TH2& histos2,
					  TLorentzVector p3body, TLorentzVector pOS1, TLorentzVector pOS2, TLorentzVector pSS,
					  bool foundOS1, bool foundOS2, bool foundSS,
					  double weight, double mBlindMin, double mBlindMax)
{
	// type=1 PVs
	vector<int> pv1; findPVindex(pv1,1);
	if(pv1.size()<1) _ERROR("PV not found");
	if(pv1.size()>1) _ERROR("nPV>1");
	int pvindex = pv1[0];
	double m3mu = p3body.M();
	
	if(foundOS1) histos[name+"_MVA_vars_mOS1"+suffix]->Fill(pOS1.M(),weight);
	if(foundOS2) histos[name+"_MVA_vars_mOS2"+suffix]->Fill(pOS2.M(),weight);
	if(foundSS)  histos[name+"_MVA_vars_mSS"+suffix]->Fill(pSS.M(),weight);
	histos[name+"_MVA_vars_isolation"+suffix]->Fill(isolation,weight);
	histos[name+"_MVA_vars_pT3body"+suffix]->Fill(p3body.Pt(),weight);
	histos[name+"_MVA_vars_pvalue"+suffix]->Fill(TMath::Prob(vtx_chi2->at(vtx),vtx_ndf->at(vtx)),weight);
	histos[name+"_MVA_vars_Lxy"+suffix]->Fill(vtx_lxy->at(vtx)[pvindex],weight);
	histos[name+"_MVA_vars_a0xy"+suffix]->Fill(vtx_a0XY->at(vtx)[pvindex],weight);
	histos[name+"_MVA_vars_cosTxy"+suffix]->Fill(fabs(vtx_cosThetaXY->at(vtx)[pvindex]),weight);
	histos[name+"_MVA_vars_MET"+suffix]->Fill(MET_RefFinal_et,weight);
	histos[name+"_MVA_vars_dPhi3bodyMET"+suffix]->Fill(fabs(dPhi(MET_RefFinal_phi,p3body.Phi())),weight);
	histos[name+"_MVA_vars_mT3bodyMET"+suffix]->Fill(mT(MET_RefFinal_et,MET_RefFinal_phi,p3body.Pt(),p3body.Phi()),weight);
	histos[name+"_MVA_vars_dPhi3bodyJ1"+suffix]->Fill(JetdPhi3bodyJ1,weight);	
	histos[name+"_MVA_vars_dPhiJ1J2"+suffix]->Fill(JetdPhiJ1J2,weight);
	histos[name+"_MVA_vars_pTJ1"+suffix]->Fill(JetPt1,weight);
	histos[name+"_MVA_vars_pTJ2"+suffix]->Fill(JetPt2,weight);
	histos[name+"_MVA_vars_mupt12Fraction"+suffix]->Fill(ptFraction12,weight);		
	histos[name+"_MVA_vars_mupt23Fraction"+suffix]->Fill(ptFraction23,weight);		
	histos[name+"_MVA_vars_mupt13Fraction"+suffix]->Fill(ptFraction13,weight);
	
	if(!reBlind(histos[name+"_MVA_vars_m3mu"+suffix],m3mu,mBlindMin,mBlindMax))           histos[name+"_MVA_vars_m3mu"+suffix]->Fill(m3mu,weight);
	if(!reBlind(histos[name+"_MVA_vars_m3mu_lin"+suffix],m3mu,mBlindMin,mBlindMax))       histos[name+"_MVA_vars_m3mu_lin"+suffix]->Fill(m3mu,weight);
	if(!reBlind(histos[name+"_MVA_vars_m3mu_lin_zoom"+suffix],m3mu,mBlindMin,mBlindMax))  histos[name+"_MVA_vars_m3mu_lin_zoom"+suffix]->Fill(m3mu,weight);
	if(!reBlind(histos[name+"_MVA_vars_m3mu_sigregion"+suffix],m3mu,mBlindMin,mBlindMax)) histos[name+"_MVA_vars_m3mu_sigregion"+suffix]->Fill(m3mu,weight);
	
	histos2[name+"_MVA_vars_dPhi3muJet1_vs_pTjet1"+suffix]->Fill(JetPt1,JetdPhi3bodyJ1,weight);
	histos2[name+"_MVA_vars_dPhiJet1Jet2_vs_sumpTjet12"+suffix]->Fill(JetSumPt12,JetdPhiJ1J2,weight);
}

double getYline(double x1, double y1, double x2, double y2, double x)
{
	if(x1==x2) _FATAL("the 2 points have the same x --> cannot calculate the slope");
	double slope = (y2-y1)/(x2-x1);
	double y = y1 + slope*(x-x1);
	return y;
}

// void clearZ4mus()
// {
// 	for(TMapuivi::iterator it=fourthMuons.begin() ; it!=fourthMuons.end() ; ++it) it->second.clear(); fourthMuons.clear();
// 	for(TMapuivi::iterator it=fourthTPs.begin()   ; it!=fourthTPs.end()   ; ++it) it->second.clear(); fourthTPs.clear();
// }
// void countZ4mus(unsigned int vtx/*, TString name, TMapTSP2TH1& histos, double weight=1.*/)
// {	
// 	double px1 = vtx_reftrks_px->at(vtx)[0];
// 	double px2 = vtx_reftrks_px->at(vtx)[1];
// 	double px3 = vtx_reftrks_px->at(vtx)[2];	
// 	double py1 = vtx_reftrks_py->at(vtx)[0];
// 	double py2 = vtx_reftrks_py->at(vtx)[1];
// 	double py3 = vtx_reftrks_py->at(vtx)[2];
// 	double pz1 = vtx_reftrks_pz->at(vtx)[0];
// 	double pz2 = vtx_reftrks_pz->at(vtx)[1];
// 	double pz3 = vtx_reftrks_pz->at(vtx)[2];
// 	TLorentzVector p1,p2,p3,psum;
// 	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
// 	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
// 	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
// 	psum = p1+p2+p3; // psum = getTlv3mu(vtx);
// 	
// 	// double mass   = vtx_mass->at(vtx); // psum.M();  in principle
// 	// double pTsum  = vtx_pt->at(vtx);   // psum.Pt(); in principle
// 	// double charge = vtx_charge->at(vtx);
// 	// double pvalue = TMath::Prob(vtx_chi2->at(vtx),vtx_ndf->at(vtx));
// 	
// 	_DEBUG("");
// 	
// 	sources src;
// 	getSrc(vtx,src);
// 	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
// 	TString shortType = classifyTripletShort(vtx);
// 	VtxType = classifyTripletCode(shortType);
// 
// 	// match the tracks
// 	int itrk1 = src.trkIndex[0];
// 	int itrk2 = src.trkIndex[1];
// 	int itrk3 = src.trkIndex[2];
// 	// match the sources        
// 	int isrc1 = src.srcIndex[0];
// 	int isrc2 = src.srcIndex[1];
// 	int isrc3 = src.srcIndex[2];
// 	// match the sources        
// 	TString src1 = src.srcName[0];
// 	TString src2 = src.srcName[1];
// 	TString src3 = src.srcName[2];
// 	// test if muons
// 	bool isMuon1 = src.isMuon[0];
// 	bool isMuon2 = src.isMuon[1];
// 	bool isMuon3 = src.isMuon[2];
// 	// test if TPmuon
// 	bool isTPmu1 = src.isTPmu[0]; bool isTPa1 = src.isTPa[0]; bool isTPb1 = src.isTPb[0];
// 	bool isTPmu2 = src.isTPmu[1]; bool isTPa2 = src.isTPa[1]; bool isTPb2 = src.isTPb[1];
// 	bool isTPmu3 = src.isTPmu[2]; bool isTPa3 = src.isTPa[2]; bool isTPb3 = src.isTPb[2];
// 	
// 	TString allowedTP = "CombinedFitMuonParticles";
// 	unsigned int nMuons = muons_pt->size();
// 	unsigned int nTPas  = tpmu_vd[allowedTP+"_pt"]->size();
// 	if(nMuons+nTPas<4) return;
// 	
// 	double margins = 30.*GeV2MeV;
// 
// 	vector<int> vfourthmuons;
// 	for(unsigned int m=0 ; m<nMuons ; ++m)
// 	{
// 		if(isMuon1 && m==isrc1) continue;
// 		if(isMuon2 && m==isrc2) continue;
// 		if(isMuon3 && m==isrc3) continue;
// 		
// 		int itrk = (mastername!="muid") ? muons_inDetTrackIndex->at(m) : muid_inDetTrackIndex->at(m);
// 		if(itrk<0 || itrk>=trks_pt->size()) continue;
// 		
// 		TLorentzVector vZ, v4th;
// 		v4th.SetPtEtaPhiM(muons_pt->at(m),muons_eta->at(m),muons_phi->at(m),muonMassMeV);
// 		vZ = psum+v4th;
// 		double mQuad = vZ.M();
// 		
// 		if(fabs(vZ.M()-91.*GeV2MeV)>margins) continue;
// 		
// 		vfourthmuons.push_back(m);
// 		
// 		// histos[name+"_triplet_mQuad_muons"]->Fill(mQuad,weight);
// 		// histos[name+"_triplet_mQuad_norm_muons"]->Fill(mQuad,weight);
// 		if(fabs(vZ.M()-91.*GeV2MeV)<margins/3.) cout << "Event: " << EventNumber << "(" << isrc1 << "," << isrc2 << "," << isrc3 << "+" << m << ") -> with 4th muon mQuad=" << mQuad << endl;
// 	}
// 	fourthMuons.insert(make_pair(vtx,vfourthmuons));
// 	
// 	vector<int> vfourthtps;
// 	for(unsigned int m=0 ; m<nTPas ; ++m)
// 	{
// 		if(isTPa1 && m==isrc1) continue;
// 		if(isTPa2 && m==isrc2) continue;
// 		if(isTPa3 && m==isrc3) continue;
// 		
// 		bool isoverlap = false;
// 		for(unsigned int n=0 ; n<vfourthmuons.size() ; ++n)
// 		{
// 			unsigned int m4 = vfourthmuons[n];
// 			double dR  = deltaR(muons_eta->at(m4),muons_phi->at(m4),tpmu_vd[allowedTP+"_eta"]->at(m),tpmu_vd[allowedTP+"_phi"]->at(m));
// 			double dpT = fabs(muons_pt->at(m4)-tpmu_vd[allowedTP+"_pt"]->at(m))/muons_pt->at(m4);
// 			double dQ  = muons_charge->at(m4)-qtrk(tpmu_vd[allowedTP+"_qOverP"]->at(m));
// 			
// 			if(dR<0.01 && dpT<0.2 && dQ==0) { isoverlap=true; break; }
// 		}
// 		if(isoverlap) continue;
// 		
// 		TLorentzVector vZ, v4th;
// 		v4th.SetPtEtaPhiM(tpmu_vd[allowedTP+"_pt"]->at(m),tpmu_vd[allowedTP+"_eta"]->at(m),tpmu_vd[allowedTP+"_phi"]->at(m),muonMassMeV);
// 		vZ = psum+v4th;
// 		double mQuad = vZ.M();
// 		
// 		if(fabs(vZ.M()-91.*GeV2MeV)>margins) continue;
// 		
// 		vfourthtps.push_back(m);
// 		
// 		// histos[name+"_triplet_mQuad_TPa"]->Fill(mQuad,weight);
// 		// histos[name+"_triplet_mQuad_norm_TPa"]->Fill(mQuad,weight);
// 		if(fabs(vZ.M()-91.*GeV2MeV)<margins/3.) cout << "Event: " << EventNumber << "(" << isrc1 << "," << isrc2 << "," << isrc3 << "+" << m << ") -> with 4th TPa mQuad=" << mQuad << endl;
// 	}
// 	fourthTPs.insert(make_pair(vtx,vfourthtps));
// }
void vertex::set(unsigned int vtx)
{
	sources src;
	getSrc(vtx,src);
	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
	m_type  = classifyTripletShort(vtx);
	m_code  = classifyTripletCode(m_type);
	m_index = vtx;
	
	_DEBUG("");
	
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];	
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3,psum;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
	psum = p1+p2+p3; // psum = getTlv3mu(vtx);
	m_p4 = psum;
	m_trkP[0] = p1;
	m_trkP[1] = p2;
	m_trkP[2] = p3;
	
	_DEBUG("");
	
	TLorentzVector pOS1, pOS2, pSS;
	bool foundOS1 = false;
	bool foundOS2 = false;
	bool foundSS  = false;
	int iUniqueCharge = -1;
	TMapVL doubletsTLV = getDoubletsFromTriplet(vtx,iUniqueCharge);
	if(doubletsTLV.find("OS1") != doubletsTLV.end()) { pOS1 = doubletsTLV["OS1"]; foundOS1 = true; }
	if(doubletsTLV.find("OS2") != doubletsTLV.end()) { pOS2 = doubletsTLV["OS2"]; foundOS2 = true; }
	if(doubletsTLV.find("SS")  != doubletsTLV.end()) { pSS  = doubletsTLV["SS"];  foundSS  = true; }
	if(foundOS1+foundOS2<2 && fabs(vtx_charge->at(vtx))==1.) _ERROR("could not find OS doublet"); // should be exactly two OS pairs
	if(!foundSS && fabs(vtx_charge->at(vtx))==1.)            _ERROR("could not find SS doublet"); // should be exactly one SS pair
	// // doublets
	// TLorentzVector pOS1, pOS2, pSS;
	// if(fabs(src.q3body)==1.)
	// {
	// 	if(q2body12==0. && q2body13==0. && q2body23!=0.) { pOS1 = src.p2body12; pOS2 = src.p2body13; pSS = src.p2body23; }
	// 	if(q2body12==0. && q2body13!=0. && q2body23==0.) { pOS1 = src.p2body12; pOS2 = src.p2body23; pSS = src.p2body13; }
	// 	if(q2body12!=0. && q2body13==0. && q2body23==0.) { pOS1 = src.p2body13; pOS2 = src.p2body23; pSS = src.p2body12; }
	// }
	m_pOS1 = pOS1;
	m_pOS2 = pOS2;
	m_pSS  = pSS;
	
	_DEBUG("");
	
	m_charge  = vtx_charge->at(vtx);;
	m_pvalue  = TMath::Prob(vtx_chi2->at(vtx),vtx_ndf->at(vtx));
	m_chi2    = vtx_chi2->at(vtx);
	m_ndf     = vtx_ndf->at(vtx);
	m_chi2ndf = vtx_chi2->at(vtx)/vtx_ndf->at(vtx);
	
	_DEBUG("");
	
	vector<int> pv1; findPVindex(pv1,1);
	if(pv1.size()<1) _ERROR("PV not found OR too many (>1) PVs were found: nPV="+_s((int)pv1.size()));
	int pvindex = pv1[0];
	m_lxy    = vtx_lxy->at(vtx)[pvindex];
	m_lxyErr = vtx_lxyErr->at(vtx)[pvindex];
	m_tau    = vtx_tau->at(vtx)[pvindex];
	m_a0     = vtx_a0->at(vtx)[pvindex];
	m_a0xy   = vtx_a0XY->at(vtx)[pvindex];
	m_cosT   = vtx_cosTheta->at(vtx)[pvindex];
	m_cosTxy = vtx_cosThetaXY->at(vtx)[pvindex];
	
	_DEBUG("");
	
	m_isolation[0]  = tripletIsolation(vtx,0.00);
	m_isolation[1]  = tripletIsolation(vtx,0.01);
	m_isolation[2]  = tripletIsolation(vtx,0.02);
	m_isolation[3]  = tripletIsolation(vtx,0.03);
	m_isolation[4]  = tripletIsolation(vtx,0.04);
	m_isolation[5]  = tripletIsolation(vtx,0.05);
	m_isolation[6]  = tripletIsolation(vtx,0.06);
	m_isolation[7]  = tripletIsolation(vtx,0.07);
	m_isolation[8]  = tripletIsolation(vtx,0.08);
	m_isolation[9]  = tripletIsolation(vtx,0.09);
	m_isolation[10] = tripletIsolation(vtx,0.10);
	m_isolation[11] = tripletIsolation(vtx,0.12);
	m_isolation[12] = tripletIsolation(vtx,0.14);
	m_isolation[13] = tripletIsolation(vtx,0.16);
	m_isolation[14] = tripletIsolation(vtx,0.18);
	m_isolation[15] = tripletIsolation(vtx,0.20);
	m_isolation[16] = tripletIsolation(vtx,0.22);
	m_isolation[17] = tripletIsolation(vtx,0.24);
	m_isolation[18] = tripletIsolation(vtx,0.26);
	m_isolation[19] = tripletIsolation(vtx,0.28);
	m_isolation[20] = tripletIsolation(vtx,0.30);
	
	_DEBUG("");
	
	int itrk1 = src.trkIndex[0];
	int itrk2 = src.trkIndex[1];
	int itrk3 = src.trkIndex[2];
	m_itrk[0] = itrk1;
	m_itrk[1] = itrk2;
	m_itrk[2] = itrk3;
	
	int isrc1 = src.srcIndex[0];
	int isrc2 = src.srcIndex[1];
	int isrc3 = src.srcIndex[2];
	m_isrc[0] = isrc1;
	m_isrc[1] = isrc2;
	m_isrc[2] = isrc3;
	      
	TString src1 = src.srcName[0];
	TString src2 = src.srcName[1];
	TString src3 = src.srcName[2];
	m_src[0] = (string)src1;
	m_src[1] = (string)src2;
	m_src[2] = (string)src3;
	
	int order1 = src.srcOrder[0];
	int order2 = src.srcOrder[1];
	int order3 = src.srcOrder[2];
	m_order[0] = order1;
	m_order[1] = order2;
	m_order[2] = order3;
	
	int code1 = src.srcCode[0];
	int code2 = src.srcCode[1];
	int code3 = src.srcCode[2];
	m_trktype[0] = code1;
	m_trktype[1] = code2;
	m_trktype[2] = code3;
	
	bool isMuon1 = src.isMuon[0];
	bool isMuon2 = src.isMuon[1];
	bool isMuon3 = src.isMuon[2];
	m_ismuon[0] = isMuon1;
	m_ismuon[1] = isMuon2;
	m_ismuon[2] = isMuon3;
	
	bool isCalo1 = src.isCalo[0];
	bool isCalo2 = src.isCalo[1];
	bool isCalo3 = src.isCalo[2];
	m_iscalo[0] = isCalo1;
	m_iscalo[1] = isCalo2;
	m_iscalo[2] = isCalo3;
	
	bool isTPmu1 = src.isTPmu[0]; bool isTPa1 = src.isTPa[0]; bool isTPb1 = src.isTPb[0];
	bool isTPmu2 = src.isTPmu[1]; bool isTPa2 = src.isTPa[1]; bool isTPb2 = src.isTPb[1];
	bool isTPmu3 = src.isTPmu[2]; bool isTPa3 = src.isTPa[2]; bool isTPb3 = src.isTPb[2];
	m_istp[0] = isTPmu1; m_istpa[0] = isTPa1; m_istpb[0] = isTPb1; 
	m_istp[1] = isTPmu2; m_istpa[1] = isTPa2; m_istpb[1] = isTPb2; 
	m_istp[2] = isTPmu3; m_istpa[2] = isTPa3; m_istpb[2] = isTPb3; 
		
	bool isCB1 = 0;
	bool isCB2 = 0;
	bool isCB3 = 0;
	if(isMuon1) isCB1 = (!src1.Contains("Muid")) ? muons_isCombined->at(isrc1) : muid_isCombined->at(isrc1);
	if(isMuon2) isCB2 = (!src2.Contains("Muid")) ? muons_isCombined->at(isrc2) : muid_isCombined->at(isrc2);
	if(isMuon3) isCB3 = (!src3.Contains("Muid")) ? muons_isCombined->at(isrc3) : muid_isCombined->at(isrc3);
	m_iscb[0] = isCB1;
	m_iscb[1] = isCB2;
	m_iscb[2] = isCB3;
	
	bool isTight1 = true;
	bool isTight2 = true;
	bool isTight3 = true;
	if(isMuon1) isTight1 = (!src1.Contains("Muid")) ? muons_isTight->at(isrc1) : muid_isTight->at(isrc1);
	if(isMuon2) isTight2 = (!src2.Contains("Muid")) ? muons_isTight->at(isrc2) : muid_isTight->at(isrc2);
	if(isMuon3) isTight3 = (!src3.Contains("Muid")) ? muons_isTight->at(isrc3) : muid_isTight->at(isrc3);
	m_istight[0] = isTight1;
	m_istight[1] = isTight2;
	m_istight[2] = isTight3;
	
	bool isMedium1 = true;
	bool isMedium2 = true;
	bool isMedium3 = true;
	if(isMuon1) isMedium1 = (!src1.Contains("Muid")) ? muons_isMedium->at(isrc1) : muid_isMedium->at(isrc1);
	if(isMuon2) isMedium2 = (!src2.Contains("Muid")) ? muons_isMedium->at(isrc2) : muid_isMedium->at(isrc2);
	if(isMuon3) isMedium3 = (!src3.Contains("Muid")) ? muons_isMedium->at(isrc3) : muid_isMedium->at(isrc3);
	m_ismedium[0] = isMedium1;
	m_ismedium[1] = isMedium2;
	m_ismedium[2] = isMedium3;
	
	bool isLoose1 = true;
	bool isLoose2 = true;
	bool isLoose3 = true;
	if(isMuon1) isLoose1 = (!src1.Contains("Muid")) ? muons_isLoose->at(isrc1) : muid_isLoose->at(isrc1);
	if(isMuon2) isLoose2 = (!src2.Contains("Muid")) ? muons_isLoose->at(isrc2) : muid_isLoose->at(isrc2);
	if(isMuon3) isLoose3 = (!src3.Contains("Muid")) ? muons_isLoose->at(isrc3) : muid_isLoose->at(isrc3);
	m_isloose[0] = isLoose1;
	m_isloose[1] = isLoose2;
	m_isloose[2] = isLoose3;
	
	double sctangsig1 = 999.;
	double sctangsig2 = 999.;
	double sctangsig3 = 999.;
	if(isMuon1) sctangsig1 = (!src1.Contains("Muid")) ? muons_sctangsig->at(isrc1) : muid_sctangsig->at(isrc1);
	if(isMuon2) sctangsig2 = (!src2.Contains("Muid")) ? muons_sctangsig->at(isrc2) : muid_sctangsig->at(isrc2);
	if(isMuon3) sctangsig3 = (!src3.Contains("Muid")) ? muons_sctangsig->at(isrc3) : muid_sctangsig->at(isrc3);
	m_trksctang[0] = sctangsig1;
	m_trksctang[1] = sctangsig2;
	m_trksctang[2] = sctangsig3;
	
	double sctngbsig1 = 999.;
	double sctngbsig2 = 999.;
	double sctngbsig3 = 999.;
	if(isMuon1) sctngbsig1 = (!src1.Contains("Muid")) ? muons_sctngbsig->at(isrc1) : muid_sctngbsig->at(isrc1);
	if(isMuon2) sctngbsig2 = (!src2.Contains("Muid")) ? muons_sctngbsig->at(isrc2) : muid_sctngbsig->at(isrc2);
	if(isMuon3) sctngbsig3 = (!src3.Contains("Muid")) ? muons_sctngbsig->at(isrc3) : muid_sctngbsig->at(isrc3);
	m_trksctngb[0] = sctngbsig1;
	m_trksctngb[1] = sctngbsig2;
	m_trksctngb[2] = sctngbsig3;
	
	double pbalsig1 = 999.;
	double pbalsig2 = 999.;
	double pbalsig3 = 999.;
	if(isMuon1) pbalsig1 = (!src1.Contains("Muid")) ? muons_pbalsig->at(isrc1) : muid_pbalsig->at(isrc1);
	if(isMuon2) pbalsig2 = (!src2.Contains("Muid")) ? muons_pbalsig->at(isrc2) : muid_pbalsig->at(isrc2);
	if(isMuon3) pbalsig3 = (!src3.Contains("Muid")) ? muons_pbalsig->at(isrc3) : muid_pbalsig->at(isrc3);
	m_trkpbal[0] = pbalsig1;
	m_trkpbal[1] = pbalsig2;
	m_trkpbal[2] = pbalsig3;
	
	double matchchi2ndfmu1 = -999.;
	double matchchi2ndfmu2 = -999.;
	double matchchi2ndfmu3 = -999.;
	if(isMuon1) matchchi2ndfmu1 = (!src1.Contains("Muid")) ? muons_matchchi2ndf->at(isrc1) : muid_matchchi2ndf->at(isrc1);
	if(isMuon2) matchchi2ndfmu2 = (!src2.Contains("Muid")) ? muons_matchchi2ndf->at(isrc2) : muid_matchchi2ndf->at(isrc2);
	if(isMuon3) matchchi2ndfmu3 = (!src3.Contains("Muid")) ? muons_matchchi2ndf->at(isrc3) : muid_matchchi2ndf->at(isrc3);
	m_trkMuMatchChi2Ndf[0] = matchchi2ndfmu1;
	m_trkMuMatchChi2Ndf[1] = matchchi2ndfmu2;
	m_trkMuMatchChi2Ndf[2] = matchchi2ndfmu3;
	
	double fitchi2mu1 = -999.;
	double fitchi2mu2 = -999.;
	double fitchi2mu3 = -999.;
	if(isMuon1) fitchi2mu1 = (!src1.Contains("Muid")) ? muons_chi2->at(isrc1) : muid_chi2->at(isrc1);
	if(isMuon2) fitchi2mu2 = (!src2.Contains("Muid")) ? muons_chi2->at(isrc2) : muid_chi2->at(isrc2);
	if(isMuon3) fitchi2mu3 = (!src3.Contains("Muid")) ? muons_chi2->at(isrc3) : muid_chi2->at(isrc3);
	double fitchi2TP1 = -999.;
	double fitchi2TP2 = -999.;
	double fitchi2TP3 = -999.;
	if(isTPmu1) fitchi2TP1 = tpmu_vd[src1+"_chi2"]->at(isrc1);  
	if(isTPmu2) fitchi2TP2 = tpmu_vd[src2+"_chi2"]->at(isrc2);  
	if(isTPmu3) fitchi2TP3 = tpmu_vd[src3+"_chi2"]->at(isrc3);
	m_trkChi2[0] = (isMuon1) ? fitchi2mu1 : fitchi2TP1;
	m_trkChi2[1] = (isMuon2) ? fitchi2mu2 : fitchi2TP2;
	m_trkChi2[2] = (isMuon3) ? fitchi2mu3 : fitchi2TP3;
	
	int fitndfmu1 = -999;
	int fitndfmu2 = -999;
	int fitndfmu3 = -999;
	if(isMuon1) fitndfmu1 = (!src1.Contains("Muid")) ? muons_ndf->at(isrc1) : muid_ndf->at(isrc1);
	if(isMuon2) fitndfmu2 = (!src2.Contains("Muid")) ? muons_ndf->at(isrc2) : muid_ndf->at(isrc2);
	if(isMuon3) fitndfmu3 = (!src3.Contains("Muid")) ? muons_ndf->at(isrc3) : muid_ndf->at(isrc3);
	double fitndfTP1 = -999.;
	double fitndfTP2 = -999.;
	double fitndfTP3 = -999.;
	if(isTPmu1) fitndfTP1 = tpmu_vi[src1+"_ndf"]->at(isrc1);
	if(isTPmu2) fitndfTP2 = tpmu_vi[src2+"_ndf"]->at(isrc2);
	if(isTPmu3) fitndfTP3 = tpmu_vi[src3+"_ndf"]->at(isrc3);
	m_trkNdf[0] = (isMuon1) ? fitndfmu1 : fitndfTP1;
	m_trkNdf[1] = (isMuon2) ? fitndfmu2 : fitndfTP2;
	m_trkNdf[2] = (isMuon3) ? fitndfmu3 : fitndfTP3;

	double fitpvaluemu1 = -999.;
	double fitpvaluemu2 = -999.;
	double fitpvaluemu3 = -999.;
	if(isMuon1) fitpvaluemu1 = TMath::Prob(fitchi2mu1,fitndfmu1);
	if(isMuon2) fitpvaluemu2 = TMath::Prob(fitchi2mu2,fitndfmu2);
	if(isMuon3) fitpvaluemu3 = TMath::Prob(fitchi2mu3,fitndfmu3);
	double fitpvalueTP1 = -999.;
	double fitpvalueTP2 = -999.;
	double fitpvalueTP3 = -999.;
	if(isTPmu1) fitpvalueTP1 = TMath::Prob(tpmu_vd[src1+"_chi2"]->at(isrc1),tpmu_vi[src1+"_ndf"]->at(isrc1));
	if(isTPmu2) fitpvalueTP2 = TMath::Prob(tpmu_vd[src2+"_chi2"]->at(isrc2),tpmu_vi[src2+"_ndf"]->at(isrc2));
	if(isTPmu3) fitpvalueTP3 = TMath::Prob(tpmu_vd[src3+"_chi2"]->at(isrc3),tpmu_vi[src3+"_ndf"]->at(isrc3));
	m_trkPval[0] = (isMuon1) ? fitpvaluemu1 : fitpvalueTP1;
	m_trkPval[1] = (isMuon2) ? fitpvaluemu2 : fitpvalueTP2;
	m_trkPval[2] = (isMuon3) ? fitpvaluemu3 : fitpvalueTP3;
	
	m_trkChi2Ndf[0] = m_trkChi2[0]/m_trkNdf[0];
	m_trkChi2Ndf[1] = m_trkChi2[1]/m_trkNdf[1];
	m_trkChi2Ndf[2] = m_trkChi2[2]/m_trkNdf[2];
	
	_DEBUG("");
	
	TVector3 p1me, p2me, p3me;
	TVector3 p1ie, p2ie, p3ie;
	double q1=-999.; double q2=-999.; double q3=-999.;
	if(isMuon1)
	{
		if(!src1.Contains("Muid")) { p1me.SetXYZ(muons_px_me->at(isrc1), muons_py_me->at(isrc1), muons_pz_me->at(isrc1)); q1 = muons_charge->at(isrc1); }
		else                       { p1me.SetXYZ(muid_px_me->at(isrc1),  muid_py_me->at(isrc1),  muid_pz_me->at(isrc1));  q1 = muid_charge->at(isrc1);  }
		if(!src1.Contains("Muid")) { p1ie.SetXYZ(muons_px_ie->at(isrc1), muons_py_ie->at(isrc1), muons_pz_ie->at(isrc1)); q1 = muons_charge->at(isrc1); }
		else                       { p1ie.SetXYZ(muid_px_ie->at(isrc1),  muid_py_ie->at(isrc1),  muid_pz_ie->at(isrc1));  q1 = muid_charge->at(isrc1);  }
	}                                                                                                                
	if(isMuon2)                                                                                                      
	{                                                                                                                
		if(!src2.Contains("Muid")) { p2me.SetXYZ(muons_px_me->at(isrc2), muons_py_me->at(isrc2), muons_pz_me->at(isrc2)); q2 = muons_charge->at(isrc2); }
		else                       { p2me.SetXYZ(muid_px_me->at(isrc2),  muid_py_me->at(isrc2),  muid_pz_me->at(isrc2));  q2 = muid_charge->at(isrc2);   }
		if(!src2.Contains("Muid")) { p2ie.SetXYZ(muons_px_ie->at(isrc2), muons_py_ie->at(isrc2), muons_pz_ie->at(isrc2)); q2 = muons_charge->at(isrc2); }
		else                       { p2ie.SetXYZ(muid_px_ie->at(isrc2),  muid_py_ie->at(isrc2),  muid_pz_ie->at(isrc2));  q2 = muid_charge->at(isrc2);  }
	}                                                                                                                
	if(isMuon3)                                                                                                      
	{                                                                                                                
		if(!src3.Contains("Muid")) { p3me.SetXYZ(muons_px_me->at(isrc3), muons_py_me->at(isrc3), muons_pz_me->at(isrc3)); q3 = muons_charge->at(isrc3); }
		else                       { p3me.SetXYZ(muid_px_me->at(isrc3),  muid_py_me->at(isrc3),  muid_pz_me->at(isrc3));  q3 = muid_charge->at(isrc3);  }
		if(!src3.Contains("Muid")) { p3ie.SetXYZ(muons_px_ie->at(isrc3), muons_py_ie->at(isrc3), muons_pz_ie->at(isrc3)); q3 = muons_charge->at(isrc3); }
		else                       { p3ie.SetXYZ(muid_px_ie->at(isrc3),  muid_py_ie->at(isrc3),  muid_pz_ie->at(isrc3));  q3 = muid_charge->at(isrc3);  }
	}
	double qopmemu1 = -999.;
	double qopmemu2 = -999.;
	double qopmemu3 = -999.;
	double qopiemu1 = -999.;
	double qopiemu2 = -999.;
	double qopiemu3 = -999.;
	if(isMuon1) { qopmemu1 = q1/p1me.Mag(); qopiemu1 = q1/p1ie.Mag(); }
	if(isMuon2) { qopmemu2 = q2/p2me.Mag(); qopiemu2 = q2/p2ie.Mag(); }
	if(isMuon3) { qopmemu3 = q3/p3me.Mag(); qopiemu3 = q3/p3ie.Mag(); }
	double qopTP1 = -999.;
	double qopTP2 = -999.;
	double qopTP3 = -999.;
	if(isTPmu1) qopTP1 = tpmu_vd[src1+"_qOverP"]->at(isrc1);
	if(isTPmu2) qopTP2 = tpmu_vd[src2+"_qOverP"]->at(isrc2);
	if(isTPmu3) qopTP3 = tpmu_vd[src3+"_qOverP"]->at(isrc3);
	m_srcQoverP[0] = (isMuon1) ? qopmemu1 : qopTP1;
	m_srcQoverP[1] = (isMuon2) ? qopmemu2 : qopTP2;
	m_srcQoverP[2] = (isMuon3) ? qopmemu3 : qopTP3;
	if(0) cout << "qopiemu1=" << qopiemu1 << endl;
	if(0) cout << "qopiemu2=" << qopiemu2 << endl;
	if(0) cout << "qopiemu3=" << qopiemu3 << endl;
	
	_DEBUG("");
	
	m_trkQoverP[0] = trks_qoverp->at(itrk1);
	m_trkQoverP[1] = trks_qoverp->at(itrk2);
	m_trkQoverP[2] = trks_qoverp->at(itrk3);
	
	_DEBUG("");
	
	m_trkPixeldEdx[0] = trks_pixeldEdx->at(itrk1);
	m_trkPixeldEdx[1] = trks_pixeldEdx->at(itrk2);
	m_trkPixeldEdx[2] = trks_pixeldEdx->at(itrk3);
	
	_DEBUG("");
	
	m_trkUsedHitsdEdx[0] = trks_nUsedHitsdEdx->at(itrk1);
	m_trkUsedHitsdEdx[1] = trks_nUsedHitsdEdx->at(itrk2);
	m_trkUsedHitsdEdx[2] = trks_nUsedHitsdEdx->at(itrk3);
	
	_DEBUG("");
	
	// int nPIXHits1 = (isTPmu1) ?  tpmu_vi[src1+"_nPix"]->at(isrc1) : trks_nPix->at(itrk1);        
	// int nPIXHits2 = (isTPmu2) ?  tpmu_vi[src2+"_nPix"]->at(isrc2) : trks_nPix->at(itrk2);
	// int nPIXHits3 = (isTPmu3) ?  tpmu_vi[src3+"_nPix"]->at(isrc3) : trks_nPix->at(itrk3);
	int nPIXHits1 = trks_nPix->at(itrk1);        
	int nPIXHits2 = trks_nPix->at(itrk2);
	int nPIXHits3 = trks_nPix->at(itrk3);
	m_trkPIXhits[0] = nPIXHits1;
	m_trkPIXhits[1] = nPIXHits2;
	m_trkPIXhits[2] = nPIXHits3;
	// int nDeadPIX1 = (isTPmu1) ?  tpmu_vi[src1+"_nDeadPixels"]->at(isrc1) : trks_nDeadPixels->at(itrk1);        
	// int nDeadPIX2 = (isTPmu2) ?  tpmu_vi[src2+"_nDeadPixels"]->at(isrc2) : trks_nDeadPixels->at(itrk2);
	// int nDeadPIX3 = (isTPmu3) ?  tpmu_vi[src3+"_nDeadPixels"]->at(isrc3) : trks_nDeadPixels->at(itrk3);
	int nDeadPIX1 = trks_nDeadPixels->at(itrk1);        
	int nDeadPIX2 = trks_nDeadPixels->at(itrk2);
	int nDeadPIX3 = trks_nDeadPixels->at(itrk3);
	m_trkDeadPIX[0] = nDeadPIX1;
	m_trkDeadPIX[1] = nDeadPIX2;
	m_trkDeadPIX[2] = nDeadPIX3;
	// int nPIXHoles1 = (isTPmu1) ?  tpmu_vi[src1+"_nPixHoles"]->at(isrc1) : trks_nPixHoles->at(itrk1);        
	// int nPIXHoles2 = (isTPmu2) ?  tpmu_vi[src2+"_nPixHoles"]->at(isrc2) : trks_nPixHoles->at(itrk2);
	// int nPIXHoles3 = (isTPmu3) ?  tpmu_vi[src3+"_nPixHoles"]->at(isrc3) : trks_nPixHoles->at(itrk3);
	int nPIXHoles1 = trks_nPixHoles->at(itrk1);        
	int nPIXHoles2 = trks_nPixHoles->at(itrk2);
	int nPIXHoles3 = trks_nPixHoles->at(itrk3);
	m_trkPIXholes[0] = nPIXHoles1;
	m_trkPIXholes[1] = nPIXHoles2;
	m_trkPIXholes[2] = nPIXHoles3;
	
	// int nSCTHits1 = (isTPmu1) ?  tpmu_vi[src1+"_nSCT"]->at(isrc1) : trks_nSCT->at(itrk1);        
	// int nSCTHits2 = (isTPmu2) ?  tpmu_vi[src2+"_nSCT"]->at(isrc2) : trks_nSCT->at(itrk2);
	// int nSCTHits3 = (isTPmu3) ?  tpmu_vi[src3+"_nSCT"]->at(isrc3) : trks_nSCT->at(itrk3);
	int nSCTHits1 = trks_nSCT->at(itrk1);        
	int nSCTHits2 = trks_nSCT->at(itrk2);
	int nSCTHits3 = trks_nSCT->at(itrk3);
	m_trkSCThits[0] = nSCTHits1;
	m_trkSCThits[1] = nSCTHits2;
	m_trkSCThits[2] = nSCTHits3;
	// int nDeadSCT1 = (isTPmu1) ?  tpmu_vi[src1+"_nDeadSCT"]->at(isrc1) : trks_nDeadSCT->at(itrk1);        
	// int nDeadSCT2 = (isTPmu2) ?  tpmu_vi[src2+"_nDeadSCT"]->at(isrc2) : trks_nDeadSCT->at(itrk2);
	// int nDeadSCT3 = (isTPmu3) ?  tpmu_vi[src3+"_nDeadSCT"]->at(isrc3) : trks_nDeadSCT->at(itrk3);
	int nDeadSCT1 = trks_nDeadSCT->at(itrk1);        
	int nDeadSCT2 = trks_nDeadSCT->at(itrk2);
	int nDeadSCT3 = trks_nDeadSCT->at(itrk3);
	m_trkDeadSCT[0] = nDeadSCT1;
	m_trkDeadSCT[1] = nDeadSCT2;
	m_trkDeadSCT[2] = nDeadSCT3;
	// int nSCTHoles1 = (isTPmu1) ?  tpmu_vi[src1+"_nSCTHoles"]->at(isrc1) : trks_nSCTHoles->at(itrk1);        
	// int nSCTHoles2 = (isTPmu2) ?  tpmu_vi[src2+"_nSCTHoles"]->at(isrc2) : trks_nSCTHoles->at(itrk2);
	// int nSCTHoles3 = (isTPmu3) ?  tpmu_vi[src3+"_nSCTHoles"]->at(isrc3) : trks_nSCTHoles->at(itrk3);
	int nSCTHoles1 = trks_nSCTHoles->at(itrk1);        
	int nSCTHoles2 = trks_nSCTHoles->at(itrk2);
	int nSCTHoles3 = trks_nSCTHoles->at(itrk3);
	m_trkSCTholes[0] = nSCTHoles1;
	m_trkSCTholes[1] = nSCTHoles2;
	m_trkSCTholes[2] = nSCTHoles3;
	
	// int nTRTHits1 = (isTPmu1) ?  tpmu_vi[src1+"_nTRT"]->at(isrc1) : trks_nTRT->at(itrk1);        
	// int nTRTHits2 = (isTPmu2) ?  tpmu_vi[src2+"_nTRT"]->at(isrc2) : trks_nTRT->at(itrk2);
	// int nTRTHits3 = (isTPmu3) ?  tpmu_vi[src3+"_nTRT"]->at(isrc3) : trks_nTRT->at(itrk3);
	int nTRTHits1 = trks_nTRT->at(itrk1);        
	int nTRTHits2 = trks_nTRT->at(itrk2);
	int nTRTHits3 = trks_nTRT->at(itrk3);
	m_trkTRThits[0] = nTRTHits1;
	m_trkTRThits[1] = nTRTHits2;
	m_trkTRThits[2] = nTRTHits3;
	// int nTRTOutliers1 = (isTPmu1) ?  tpmu_vi[src1+"_nTRTOutliers"]->at(isrc1) : trks_nTRTOutliers->at(itrk1);        
	// int nTRTOutliers2 = (isTPmu2) ?  tpmu_vi[src2+"_nTRTOutliers"]->at(isrc2) : trks_nTRTOutliers->at(itrk2);
	// int nTRTOutliers3 = (isTPmu3) ?  tpmu_vi[src3+"_nTRTOutliers"]->at(isrc3) : trks_nTRTOutliers->at(itrk3);
	int nTRTOutliers1 = trks_nTRTOutliers->at(itrk1);        
	int nTRTOutliers2 = trks_nTRTOutliers->at(itrk2);
	int nTRTOutliers3 = trks_nTRTOutliers->at(itrk3);
	m_trkTRToutliers[0] = nTRTOutliers1;
	m_trkTRToutliers[1] = nTRTOutliers2;
	m_trkTRToutliers[2] = nTRTOutliers3;
	// int nHighThresholdTRTHits1 = (isTPmu1) ?  tpmu_vi[src1+"_nHighThresholdTRTHits"]->at(isrc1) : trks_nHighThresholdTRTHits->at(itrk1);        
	// int nHighThresholdTRTHits2 = (isTPmu2) ?  tpmu_vi[src2+"_nHighThresholdTRTHits"]->at(isrc2) : trks_nHighThresholdTRTHits->at(itrk2);
	// int nHighThresholdTRTHits3 = (isTPmu3) ?  tpmu_vi[src3+"_nHighThresholdTRTHits"]->at(isrc3) : trks_nHighThresholdTRTHits->at(itrk3);
	int nHighThresholdTRTHits1 = trks_nHighThresholdTRTHits->at(itrk1);        
	int nHighThresholdTRTHits2 = trks_nHighThresholdTRTHits->at(itrk2);
	int nHighThresholdTRTHits3 = trks_nHighThresholdTRTHits->at(itrk3);
	m_trkHtTRThits[0] = nHighThresholdTRTHits1;
	m_trkHtTRThits[1] = nHighThresholdTRTHits2;
	m_trkHtTRThits[2] = nHighThresholdTRTHits3;
	
	_DEBUG("");
	
	m_trkMDThits[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfMdtHits"]->at(isrc1) : -1;
	m_trkMDThits[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfMdtHits"]->at(isrc2) : -1;
	m_trkMDThits[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfMdtHits"]->at(isrc3) : -1;
	if(isMuon1) m_trkMDThits[0] = (!src1.Contains("Muid")) ?  muons_numberOfMdtHits->at(isrc1) : muid_numberOfMdtHits->at(isrc1);
	if(isMuon2) m_trkMDThits[1] = (!src2.Contains("Muid")) ?  muons_numberOfMdtHits->at(isrc2) : muid_numberOfMdtHits->at(isrc2);
	if(isMuon3) m_trkMDThits[2] = (!src3.Contains("Muid")) ?  muons_numberOfMdtHits->at(isrc3) : muid_numberOfMdtHits->at(isrc3);
	
	_DEBUG("");
	
	m_trkTGCPhiHits[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfTgcPhiHits"]->at(isrc1) : -1;
	m_trkTGCPhiHits[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfTgcPhiHits"]->at(isrc2) : -1;
	m_trkTGCPhiHits[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfTgcPhiHits"]->at(isrc3) : -1;
	if(isMuon1) m_trkTGCPhiHits[0] = (!src1.Contains("Muid")) ?  muons_numberOfTgcPhiHits->at(isrc1) : muid_numberOfTgcPhiHits->at(isrc1);
	if(isMuon2) m_trkTGCPhiHits[1] = (!src2.Contains("Muid")) ?  muons_numberOfTgcPhiHits->at(isrc2) : muid_numberOfTgcPhiHits->at(isrc2);
	if(isMuon3) m_trkTGCPhiHits[2] = (!src3.Contains("Muid")) ?  muons_numberOfTgcPhiHits->at(isrc3) : muid_numberOfTgcPhiHits->at(isrc3);
	
	_DEBUG("");
	
	m_trkTGCEtaHits[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfTgcEtaHits"]->at(isrc1) : -1;
	m_trkTGCEtaHits[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfTgcEtaHits"]->at(isrc2) : -1;
	m_trkTGCEtaHits[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfTgcEtaHits"]->at(isrc3) : -1;
	if(isMuon1) m_trkTGCEtaHits[0] = (!src1.Contains("Muid")) ?  muons_numberOfTgcEtaHits->at(isrc1) : muid_numberOfTgcEtaHits->at(isrc1);
	if(isMuon2) m_trkTGCEtaHits[1] = (!src2.Contains("Muid")) ?  muons_numberOfTgcEtaHits->at(isrc2) : muid_numberOfTgcEtaHits->at(isrc2);
	if(isMuon3) m_trkTGCEtaHits[2] = (!src3.Contains("Muid")) ?  muons_numberOfTgcEtaHits->at(isrc3) : muid_numberOfTgcEtaHits->at(isrc3);
	
	_DEBUG("");
	
	m_trkCSCPhiHits[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfCscPhiHits"]->at(isrc1) : -1;
	m_trkCSCPhiHits[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfCscPhiHits"]->at(isrc2) : -1;
	m_trkCSCPhiHits[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfCscPhiHits"]->at(isrc3) : -1;
	if(isMuon1) m_trkCSCPhiHits[0] = (!src1.Contains("Muid")) ?  muons_numberOfCscPhiHits->at(isrc1) : muid_numberOfCscPhiHits->at(isrc1);
	if(isMuon2) m_trkCSCPhiHits[1] = (!src2.Contains("Muid")) ?  muons_numberOfCscPhiHits->at(isrc2) : muid_numberOfCscPhiHits->at(isrc2);
	if(isMuon3) m_trkCSCPhiHits[2] = (!src3.Contains("Muid")) ?  muons_numberOfCscPhiHits->at(isrc3) : muid_numberOfCscPhiHits->at(isrc3);
	
	
	_DEBUG("");
	
	m_trkCSCEtaHits[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfCscEtaHits"]->at(isrc1) : -1;
	m_trkCSCEtaHits[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfCscEtaHits"]->at(isrc2) : -1;
	m_trkCSCEtaHits[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfCscEtaHits"]->at(isrc3) : -1;
	if(isMuon1) m_trkCSCEtaHits[0] = (!src1.Contains("Muid")) ?  muons_numberOfCscEtaHits->at(isrc1) : muid_numberOfCscEtaHits->at(isrc1);
	if(isMuon2) m_trkCSCEtaHits[1] = (!src2.Contains("Muid")) ?  muons_numberOfCscEtaHits->at(isrc2) : muid_numberOfCscEtaHits->at(isrc2);
	if(isMuon3) m_trkCSCEtaHits[2] = (!src3.Contains("Muid")) ?  muons_numberOfCscEtaHits->at(isrc3) : muid_numberOfCscEtaHits->at(isrc3);
	
	_DEBUG("");
	
	m_trkRPCPhiHits[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfRpcPhiHits"]->at(isrc1) : -1;
	m_trkRPCPhiHits[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfRpcPhiHits"]->at(isrc2) : -1;
	m_trkRPCPhiHits[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfRpcPhiHits"]->at(isrc3) : -1;
	if(isMuon1) m_trkRPCPhiHits[0] = (!src1.Contains("Muid")) ?  muons_numberOfRpcPhiHits->at(isrc1) : muid_numberOfRpcPhiHits->at(isrc1);
	if(isMuon2) m_trkRPCPhiHits[1] = (!src2.Contains("Muid")) ?  muons_numberOfRpcPhiHits->at(isrc2) : muid_numberOfRpcPhiHits->at(isrc2);
	if(isMuon3) m_trkRPCPhiHits[2] = (!src3.Contains("Muid")) ?  muons_numberOfRpcPhiHits->at(isrc3) : muid_numberOfRpcPhiHits->at(isrc3);
	
	_DEBUG("");
	
	m_trkRPCEtaHits[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfRpcEtaHits"]->at(isrc1) : -1;
	m_trkRPCEtaHits[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfRpcEtaHits"]->at(isrc2) : -1;
	m_trkRPCEtaHits[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfRpcEtaHits"]->at(isrc3) : -1;
	if(isMuon1) m_trkRPCEtaHits[0] = (!src1.Contains("Muid")) ?  muons_numberOfRpcEtaHits->at(isrc1) : muid_numberOfRpcEtaHits->at(isrc1);
	if(isMuon2) m_trkRPCEtaHits[1] = (!src2.Contains("Muid")) ?  muons_numberOfRpcEtaHits->at(isrc2) : muid_numberOfRpcEtaHits->at(isrc2);
	if(isMuon3) m_trkRPCEtaHits[2] = (!src3.Contains("Muid")) ?  muons_numberOfRpcEtaHits->at(isrc3) : muid_numberOfRpcEtaHits->at(isrc3);
	
	_DEBUG("");
	
	m_trkCSCEtaHoles[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfCscEtaHoles"]->at(isrc1) : -1;
	m_trkCSCEtaHoles[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfCscEtaHoles"]->at(isrc2) : -1;
	m_trkCSCEtaHoles[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfCscEtaHoles"]->at(isrc3) : -1;
	if(isMuon1) m_trkCSCEtaHoles[0] = (!src1.Contains("Muid")) ?  muons_numberOfCscEtaHoles->at(isrc1) : muid_numberOfCscEtaHoles->at(isrc1);
	if(isMuon2) m_trkCSCEtaHoles[1] = (!src2.Contains("Muid")) ?  muons_numberOfCscEtaHoles->at(isrc2) : muid_numberOfCscEtaHoles->at(isrc2);
	if(isMuon3) m_trkCSCEtaHoles[2] = (!src3.Contains("Muid")) ?  muons_numberOfCscEtaHoles->at(isrc3) : muid_numberOfCscEtaHoles->at(isrc3);
	
	_DEBUG("");
	
	m_trkCSCPhiHoles[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfCscPhiHoles"]->at(isrc1) : -1;
	m_trkCSCPhiHoles[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfCscPhiHoles"]->at(isrc2) : -1;
	m_trkCSCPhiHoles[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfCscPhiHoles"]->at(isrc3) : -1;
	if(isMuon1) m_trkCSCPhiHoles[0] = (!src1.Contains("Muid")) ?  muons_numberOfCscPhiHoles->at(isrc1) : muid_numberOfCscPhiHoles->at(isrc1);
	if(isMuon2) m_trkCSCPhiHoles[1] = (!src2.Contains("Muid")) ?  muons_numberOfCscPhiHoles->at(isrc2) : muid_numberOfCscPhiHoles->at(isrc2);
	if(isMuon3) m_trkCSCPhiHoles[2] = (!src3.Contains("Muid")) ?  muons_numberOfCscPhiHoles->at(isrc3) : muid_numberOfCscPhiHoles->at(isrc3);
	
	_DEBUG("");
	
	m_trkRPCEtaHoles[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfRpcEtaHoles"]->at(isrc1) : -1;
	m_trkRPCEtaHoles[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfRpcEtaHoles"]->at(isrc2) : -1;
	m_trkRPCEtaHoles[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfRpcEtaHoles"]->at(isrc3) : -1;
	if(isMuon1) m_trkRPCEtaHoles[0] = (!src1.Contains("Muid")) ?  muons_numberOfRpcEtaHoles->at(isrc1) : muid_numberOfRpcEtaHoles->at(isrc1);
	if(isMuon2) m_trkRPCEtaHoles[1] = (!src2.Contains("Muid")) ?  muons_numberOfRpcEtaHoles->at(isrc2) : muid_numberOfRpcEtaHoles->at(isrc2);
	if(isMuon3) m_trkRPCEtaHoles[2] = (!src3.Contains("Muid")) ?  muons_numberOfRpcEtaHoles->at(isrc3) : muid_numberOfRpcEtaHoles->at(isrc3);
	
	_DEBUG("");
	
	m_trkRPCPhiHoles[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfRpcPhiHoles"]->at(isrc1) : -1;
	m_trkRPCPhiHoles[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfRpcPhiHoles"]->at(isrc2) : -1;
	m_trkRPCPhiHoles[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfRpcPhiHoles"]->at(isrc3) : -1;
	if(isMuon1) m_trkRPCPhiHoles[0] = (!src1.Contains("Muid")) ?  muons_numberOfRpcPhiHoles->at(isrc1) : muid_numberOfRpcPhiHoles->at(isrc1);
	if(isMuon2) m_trkRPCPhiHoles[1] = (!src2.Contains("Muid")) ?  muons_numberOfRpcPhiHoles->at(isrc2) : muid_numberOfRpcPhiHoles->at(isrc2);
	if(isMuon3) m_trkRPCPhiHoles[2] = (!src3.Contains("Muid")) ?  muons_numberOfRpcPhiHoles->at(isrc3) : muid_numberOfRpcPhiHoles->at(isrc3);
	
	_DEBUG("");
	
	m_trkTGCEtaHoles[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfTgcEtaHoles"]->at(isrc1) : -1;
	m_trkTGCEtaHoles[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfTgcEtaHoles"]->at(isrc2) : -1;
	m_trkTGCEtaHoles[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfTgcEtaHoles"]->at(isrc3) : -1;
	if(isMuon1) m_trkTGCEtaHoles[0] = (!src1.Contains("Muid")) ?  muons_numberOfTgcEtaHoles->at(isrc1) : muid_numberOfTgcEtaHoles->at(isrc1);
	if(isMuon2) m_trkTGCEtaHoles[1] = (!src2.Contains("Muid")) ?  muons_numberOfTgcEtaHoles->at(isrc2) : muid_numberOfTgcEtaHoles->at(isrc2);
	if(isMuon3) m_trkTGCEtaHoles[2] = (!src3.Contains("Muid")) ?  muons_numberOfTgcEtaHoles->at(isrc3) : muid_numberOfTgcEtaHoles->at(isrc3);
	
	_DEBUG("");
	
	m_trkTGCPhiHoles[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfTgcPhiHoles"]->at(isrc1) : -1;
	m_trkTGCPhiHoles[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfTgcPhiHoles"]->at(isrc2) : -1;
	m_trkTGCPhiHoles[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfTgcPhiHoles"]->at(isrc3) : -1;
	if(isMuon1) m_trkTGCPhiHoles[0] = (!src1.Contains("Muid")) ?  muons_numberOfTgcPhiHoles->at(isrc1) : muid_numberOfTgcPhiHoles->at(isrc1);
	if(isMuon2) m_trkTGCPhiHoles[1] = (!src2.Contains("Muid")) ?  muons_numberOfTgcPhiHoles->at(isrc2) : muid_numberOfTgcPhiHoles->at(isrc2);
	if(isMuon3) m_trkTGCPhiHoles[2] = (!src3.Contains("Muid")) ?  muons_numberOfTgcPhiHoles->at(isrc3) : muid_numberOfTgcPhiHoles->at(isrc3);
	
	_DEBUG("");
	
	m_trkMDTholes[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfMdtHoles"]->at(isrc1) : -1;
	m_trkMDTholes[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfMdtHoles"]->at(isrc2) : -1;
	m_trkMDTholes[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfMdtHoles"]->at(isrc3) : -1;
	if(isMuon1) m_trkMDTholes[0] = (!src1.Contains("Muid")) ?  muons_numberOfMdtHoles->at(isrc1) : muid_numberOfMdtHoles->at(isrc1);
	if(isMuon2) m_trkMDTholes[1] = (!src2.Contains("Muid")) ?  muons_numberOfMdtHoles->at(isrc2) : muid_numberOfMdtHoles->at(isrc2);
	if(isMuon3) m_trkMDTholes[2] = (!src3.Contains("Muid")) ?  muons_numberOfMdtHoles->at(isrc3) : muid_numberOfMdtHoles->at(isrc3);
	
	_DEBUG("");
	
	m_trkOutliersOnTrack[0] = (isTPmu1) ?  tpmu_vi[src1+"_numberOfOutliersOnTrack"]->at(isrc1) : -1;
	m_trkOutliersOnTrack[1] = (isTPmu2) ?  tpmu_vi[src2+"_numberOfOutliersOnTrack"]->at(isrc2) : -1;
	m_trkOutliersOnTrack[2] = (isTPmu3) ?  tpmu_vi[src3+"_numberOfOutliersOnTrack"]->at(isrc3) : -1;
	if(isMuon1) m_trkOutliersOnTrack[0] = (!src1.Contains("Muid")) ?  muons_numberOfOutliersOnTrack->at(isrc1) : muid_numberOfOutliersOnTrack->at(isrc1);
	if(isMuon2) m_trkOutliersOnTrack[1] = (!src2.Contains("Muid")) ?  muons_numberOfOutliersOnTrack->at(isrc2) : muid_numberOfOutliersOnTrack->at(isrc2);
	if(isMuon3) m_trkOutliersOnTrack[2] = (!src3.Contains("Muid")) ?  muons_numberOfOutliersOnTrack->at(isrc3) : muid_numberOfOutliersOnTrack->at(isrc3);
	
	_DEBUG("");
	
	m_trkStdDevOfChi2OS[0] = (isTPmu1) ?  tpmu_vi[src1+"_standardDeviationOfChi2OS"]->at(isrc1) : -1;
	m_trkStdDevOfChi2OS[1] = (isTPmu2) ?  tpmu_vi[src2+"_standardDeviationOfChi2OS"]->at(isrc2) : -1;
	m_trkStdDevOfChi2OS[2] = (isTPmu3) ?  tpmu_vi[src3+"_standardDeviationOfChi2OS"]->at(isrc3) : -1;
	if(isMuon1) m_trkStdDevOfChi2OS[0] = (!src1.Contains("Muid")) ?  muons_standardDeviationOfChi2OS->at(isrc1) : muid_standardDeviationOfChi2OS->at(isrc1);
	if(isMuon2) m_trkStdDevOfChi2OS[1] = (!src2.Contains("Muid")) ?  muons_standardDeviationOfChi2OS->at(isrc2) : muid_standardDeviationOfChi2OS->at(isrc2);
	if(isMuon3) m_trkStdDevOfChi2OS[2] = (!src3.Contains("Muid")) ?  muons_standardDeviationOfChi2OS->at(isrc3) : muid_standardDeviationOfChi2OS->at(isrc3);
	
	_DEBUG("");
	
	
	
	//////////////////////////////////////////////////////////////////////
	/////////////////
	//// BUG fix ////
	/////////////////
	TMapII bugMU;
	for(unsigned int i=0 ; i<muons_numberOfCscEtaHoles->size() ;  ++i)
	{
		if(muons_numberOfCscEtaHoles->at(i)<0) continue;
		else
		{
			unsigned int index = bugMU.size(); // 1st iteration size=0, 2nd iteration size=1,...
			bugMU.insert(make_pair(i,index));
		}
	}
	bool bugMU_doExsit1 = (isMuon1 && bugMU.find(isrc1)!=bugMU.end());
	bool bugMU_doExsit2 = (isMuon2 && bugMU.find(isrc2)!=bugMU.end());
	bool bugMU_doExsit3 = (isMuon3 && bugMU.find(isrc3)!=bugMU.end());
	_DEBUG("");
	// if(isTPmu1 && tpmu_vi[src1+"_numberOfCscEtaHoles"]->size()!=tpmu_vvi[src1+"_nprecisionHits"]->size()) cout << "sizeCsc=" << tpmu_vi[src1+"_numberOfCscEtaHoles"]->size() << ", sizePcn=" << tpmu_vvi[src1+"_nprecisionHits"]->size() << endl;
	// if(isTPmu2 && tpmu_vi[src2+"_numberOfCscEtaHoles"]->size()!=tpmu_vvi[src2+"_nprecisionHits"]->size()) cout << "sizeCsc=" << tpmu_vi[src2+"_numberOfCscEtaHoles"]->size() << ", sizePcn=" << tpmu_vvi[src2+"_nprecisionHits"]->size() << endl;
	// if(isTPmu3 && tpmu_vi[src3+"_numberOfCscEtaHoles"]->size()!=tpmu_vvi[src3+"_nprecisionHits"]->size()) cout << "sizeCsc=" << tpmu_vi[src3+"_numberOfCscEtaHoles"]->size() << ", sizePcn=" << tpmu_vvi[src3+"_nprecisionHits"]->size() << endl;	
	// cout << "muons_pt->size()=" << muons_pt->size() << endl;
	// cout << "muons_numberOfCscEtaHoles->size()=" << muons_numberOfCscEtaHoles->size() << endl;
	// for(unsigned int i=0 ; i<muons_numberOfCscEtaHoles->size() ;  ++i)
	// {
	// 	int x = (bugMU.find(i)!=bugMU.end()) ? bugMU[i] : -1;
	// 	cout << "["<<i<<"]: " << "bugMU["<< x <<"]" << muons_numberOfCscEtaHoles->at(i) << endl;
	// }
	// cout << "muons_nprecisionHits->size()=" << muons_nprecisionHits->size() << endl;
	// cout << "muons_nphiLayers->size()=" << muons_nphiLayers->size() << endl;
	// cout << "muons_netaPhiLayers->size()=" << muons_netaPhiLayers->size() << endl;
	// cout << "src1=" << src1 << ", isrc1=" << isrc1 << endl;
	// cout << "src2=" << src2 << ", isrc2=" << isrc2 << endl;
	// cout << "src3=" << src3 << ", isrc3=" << isrc3 << endl;
	//////////////////////////////////////////////////////////////////////
	
	
	_DEBUG("");
	unsigned int N1 = 0;
	unsigned int N2 = 0;
	unsigned int N3 = 0;
	N1 = (isTPmu1) ? tpmu_vvi[src1+"_nprecisionHits"]->at(isrc1).size() : 0;
	N2 = (isTPmu2) ? tpmu_vvi[src2+"_nprecisionHits"]->at(isrc2).size() : 0;
	N3 = (isTPmu3) ? tpmu_vvi[src3+"_nprecisionHits"]->at(isrc3).size() : 0;
	_DEBUG("");
	// cout << "muons_nprecisionHits->size()=" << muons_nprecisionHits->size() << ", muons_pt->size()=" << muons_pt->size() << endl;
	// for(TMapII::iterator it=bugMU.begin() ; it!=bugMU.end() ; ++it) cout << "bug["<<it->first<<"]=" << it->second << endl;
	if(isMuon1 && bugMU_doExsit1) N1 = (!src1.Contains("Muid")) ? muons_nprecisionHits->at(bugMU[isrc1]).size() : muid_nprecisionHits->at(bugMU[isrc1]).size();
	_DEBUG("");
	if(isMuon2 && bugMU_doExsit2) N2 = (!src2.Contains("Muid")) ? muons_nprecisionHits->at(bugMU[isrc2]).size() : muid_nprecisionHits->at(bugMU[isrc2]).size();
	_DEBUG("");
	if(isMuon3 && bugMU_doExsit3) N3 = (!src3.Contains("Muid")) ? muons_nprecisionHits->at(bugMU[isrc3]).size() : muid_nprecisionHits->at(bugMU[isrc3]).size();
	_DEBUG("");
	m_trkPrecisionHits[0] = 0;
	m_trkPrecisionHits[1] = 0;
	m_trkPrecisionHits[2] = 0;
	_DEBUG("");
	if(isTPmu1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPrecisionHits[0] += tpmu_vvi[src1+"_nprecisionHits"]->at(isrc1)[i]; }
	if(isTPmu2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPrecisionHits[1] += tpmu_vvi[src2+"_nprecisionHits"]->at(isrc2)[i]; }
	if(isTPmu3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPrecisionHits[2] += tpmu_vvi[src3+"_nprecisionHits"]->at(isrc3)[i]; }
	_DEBUG("");
	if(isMuon1 && bugMU_doExsit1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPrecisionHits[0] += (!src1.Contains("Muid")) ? muons_nprecisionHits->at(bugMU[isrc1])[i] : muid_nprecisionHits->at(bugMU[isrc1])[i]; }
	if(isMuon2 && bugMU_doExsit2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPrecisionHits[1] += (!src2.Contains("Muid")) ? muons_nprecisionHits->at(bugMU[isrc2])[i] : muid_nprecisionHits->at(bugMU[isrc2])[i]; }
	if(isMuon3 && bugMU_doExsit3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPrecisionHits[2] += (!src3.Contains("Muid")) ? muons_nprecisionHits->at(bugMU[isrc3])[i] : muid_nprecisionHits->at(bugMU[isrc3])[i]; }

	_DEBUG("");

	N1 = (isTPmu1) ? tpmu_vvi[src1+"_nphiLayers"]->at(isrc1).size() : 0;
	N2 = (isTPmu2) ? tpmu_vvi[src2+"_nphiLayers"]->at(isrc2).size() : 0;
	N3 = (isTPmu3) ? tpmu_vvi[src3+"_nphiLayers"]->at(isrc3).size() : 0;
	if(isMuon1 && bugMU_doExsit1) N1 = (!src1.Contains("Muid")) ? muons_nphiLayers->at(bugMU[isrc1]).size() : muid_nphiLayers->at(bugMU[isrc1]).size();
	if(isMuon2 && bugMU_doExsit2) N2 = (!src2.Contains("Muid")) ? muons_nphiLayers->at(bugMU[isrc2]).size() : muid_nphiLayers->at(bugMU[isrc2]).size();
	if(isMuon3 && bugMU_doExsit3) N3 = (!src3.Contains("Muid")) ? muons_nphiLayers->at(bugMU[isrc3]).size() : muid_nphiLayers->at(bugMU[isrc3]).size();
	m_trkPhiLayers[0] = 0;
	m_trkPhiLayers[1] = 0;
	m_trkPhiLayers[2] = 0;
	if(isTPmu1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPhiLayers[0] += tpmu_vvi[src1+"_nphiLayers"]->at(isrc1)[i]; }
	if(isTPmu2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPhiLayers[1] += tpmu_vvi[src2+"_nphiLayers"]->at(isrc2)[i]; }
	if(isTPmu3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPhiLayers[2] += tpmu_vvi[src3+"_nphiLayers"]->at(isrc3)[i]; }
	if(isMuon1 && bugMU_doExsit1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPhiLayers[0] += (!src1.Contains("Muid")) ? muons_nphiLayers->at(bugMU[isrc1])[i] : muid_nphiLayers->at(bugMU[isrc1])[i]; }
	if(isMuon2 && bugMU_doExsit2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPhiLayers[1] += (!src2.Contains("Muid")) ? muons_nphiLayers->at(bugMU[isrc2])[i] : muid_nphiLayers->at(bugMU[isrc2])[i]; }
	if(isMuon3 && bugMU_doExsit3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPhiLayers[2] += (!src3.Contains("Muid")) ? muons_nphiLayers->at(bugMU[isrc3])[i] : muid_nphiLayers->at(bugMU[isrc3])[i]; }
	
	_DEBUG("");
	
	N1 = (isTPmu1) ? tpmu_vvi[src1+"_netaPhiLayers"]->at(isrc1).size() : 0;
	N2 = (isTPmu2) ? tpmu_vvi[src2+"_netaPhiLayers"]->at(isrc2).size() : 0;
	N3 = (isTPmu3) ? tpmu_vvi[src3+"_netaPhiLayers"]->at(isrc3).size() : 0;
	if(isMuon1 && bugMU_doExsit1) N1 = (!src1.Contains("Muid")) ? muons_netaPhiLayers->at(bugMU[isrc1]).size() : muid_netaPhiLayers->at(bugMU[isrc1]).size();
	if(isMuon2 && bugMU_doExsit2) N2 = (!src2.Contains("Muid")) ? muons_netaPhiLayers->at(bugMU[isrc2]).size() : muid_netaPhiLayers->at(bugMU[isrc2]).size();
	if(isMuon3 && bugMU_doExsit3) N3 = (!src3.Contains("Muid")) ? muons_netaPhiLayers->at(bugMU[isrc3]).size() : muid_netaPhiLayers->at(bugMU[isrc3]).size();
	m_trkEtaPhiLayers[0] = 0;
	m_trkEtaPhiLayers[1] = 0;
	m_trkEtaPhiLayers[2] = 0;
	if(isTPmu1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkEtaPhiLayers[0] += tpmu_vvi[src1+"_netaPhiLayers"]->at(isrc1)[i]; }
	if(isTPmu2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkEtaPhiLayers[1] += tpmu_vvi[src2+"_netaPhiLayers"]->at(isrc2)[i]; }
	if(isTPmu3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkEtaPhiLayers[2] += tpmu_vvi[src3+"_netaPhiLayers"]->at(isrc3)[i]; }
	if(isMuon1 && bugMU_doExsit1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkEtaPhiLayers[0] += (!src1.Contains("Muid")) ? muons_netaPhiLayers->at(bugMU[isrc1])[i] : muid_netaPhiLayers->at(bugMU[isrc1])[i]; }
	if(isMuon2 && bugMU_doExsit2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkEtaPhiLayers[1] += (!src2.Contains("Muid")) ? muons_netaPhiLayers->at(bugMU[isrc2])[i] : muid_netaPhiLayers->at(bugMU[isrc2])[i]; }
	if(isMuon3 && bugMU_doExsit3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkEtaPhiLayers[2] += (!src3.Contains("Muid")) ? muons_netaPhiLayers->at(bugMU[isrc3])[i] : muid_netaPhiLayers->at(bugMU[isrc3])[i]; }
	
	
	_DEBUG("");
	
	N1 = (isTPmu1) ? tpmu_vvi[src1+"_nprecisionHoles"]->at(isrc1).size() : 0;
	N2 = (isTPmu2) ? tpmu_vvi[src2+"_nprecisionHoles"]->at(isrc2).size() : 0;
	N3 = (isTPmu3) ? tpmu_vvi[src3+"_nprecisionHoles"]->at(isrc3).size() : 0;
	if(isMuon1 && bugMU_doExsit1) N1 = (!src1.Contains("Muid")) ? muons_nprecisionHoles->at(bugMU[isrc1]).size() : muid_nprecisionHoles->at(bugMU[isrc1]).size();
	if(isMuon2 && bugMU_doExsit2) N2 = (!src2.Contains("Muid")) ? muons_nprecisionHoles->at(bugMU[isrc2]).size() : muid_nprecisionHoles->at(bugMU[isrc2]).size();
	if(isMuon3 && bugMU_doExsit3) N3 = (!src3.Contains("Muid")) ? muons_nprecisionHoles->at(bugMU[isrc3]).size() : muid_nprecisionHoles->at(bugMU[isrc3]).size();
	m_trkPrecisionHoles[0] = 0;
	m_trkPrecisionHoles[1] = 0;
	m_trkPrecisionHoles[2] = 0;
	if(isTPmu1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPrecisionHoles[0] += tpmu_vvi[src1+"_nprecisionHoles"]->at(isrc1)[i]; }
	if(isTPmu2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPrecisionHoles[1] += tpmu_vvi[src2+"_nprecisionHoles"]->at(isrc2)[i]; }
	if(isTPmu3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPrecisionHoles[2] += tpmu_vvi[src3+"_nprecisionHoles"]->at(isrc3)[i]; }
	if(isMuon1 && bugMU_doExsit1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPrecisionHoles[0] += (!src1.Contains("Muid")) ? muons_nprecisionHoles->at(bugMU[isrc1])[i] : muid_nprecisionHoles->at(bugMU[isrc1])[i]; }
	if(isMuon2 && bugMU_doExsit2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPrecisionHoles[1] += (!src2.Contains("Muid")) ? muons_nprecisionHoles->at(bugMU[isrc2])[i] : muid_nprecisionHoles->at(bugMU[isrc2])[i]; }
	if(isMuon3 && bugMU_doExsit3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPrecisionHoles[2] += (!src3.Contains("Muid")) ? muons_nprecisionHoles->at(bugMU[isrc3])[i] : muid_nprecisionHoles->at(bugMU[isrc3])[i]; }
	
	_DEBUG("");
	
	N1 = (isTPmu1) ? tpmu_vvi[src1+"_netaTriggerHoleLayers"]->at(isrc1).size() : 0;
	N2 = (isTPmu2) ? tpmu_vvi[src2+"_netaTriggerHoleLayers"]->at(isrc2).size() : 0;
	N3 = (isTPmu3) ? tpmu_vvi[src3+"_netaTriggerHoleLayers"]->at(isrc3).size() : 0;
	if(isMuon1 && bugMU_doExsit1) N1 = (!src1.Contains("Muid")) ? muons_netaTriggerHoleLayers->at(bugMU[isrc1]).size() : muid_netaTriggerHoleLayers->at(bugMU[isrc1]).size();
	if(isMuon2 && bugMU_doExsit2) N2 = (!src2.Contains("Muid")) ? muons_netaTriggerHoleLayers->at(bugMU[isrc2]).size() : muid_netaTriggerHoleLayers->at(bugMU[isrc2]).size();
	if(isMuon3 && bugMU_doExsit3) N3 = (!src3.Contains("Muid")) ? muons_netaTriggerHoleLayers->at(bugMU[isrc3]).size() : muid_netaTriggerHoleLayers->at(bugMU[isrc3]).size();
	m_trkEtaTriggerHoleLayers[0] = 0;
	m_trkEtaTriggerHoleLayers[1] = 0;
	m_trkEtaTriggerHoleLayers[2] = 0;
	if(isTPmu1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkEtaTriggerHoleLayers[0] += tpmu_vvi[src1+"_netaTriggerHoleLayers"]->at(isrc1)[i]; }
	if(isTPmu2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkEtaTriggerHoleLayers[1] += tpmu_vvi[src2+"_netaTriggerHoleLayers"]->at(isrc2)[i]; }
	if(isTPmu3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkEtaTriggerHoleLayers[2] += tpmu_vvi[src3+"_netaTriggerHoleLayers"]->at(isrc3)[i]; }
	if(isMuon1 && bugMU_doExsit1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkEtaTriggerHoleLayers[0] += (!src1.Contains("Muid")) ? muons_netaTriggerHoleLayers->at(bugMU[isrc1])[i] : muid_netaTriggerHoleLayers->at(bugMU[isrc1])[i]; }
	if(isMuon2 && bugMU_doExsit2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkEtaTriggerHoleLayers[1] += (!src2.Contains("Muid")) ? muons_netaTriggerHoleLayers->at(bugMU[isrc2])[i] : muid_netaTriggerHoleLayers->at(bugMU[isrc2])[i]; }
	if(isMuon3 && bugMU_doExsit3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkEtaTriggerHoleLayers[2] += (!src3.Contains("Muid")) ? muons_netaTriggerHoleLayers->at(bugMU[isrc3])[i] : muid_netaTriggerHoleLayers->at(bugMU[isrc3])[i]; }
	
	_DEBUG("");
	
	N1 = (isTPmu1) ? tpmu_vvi[src1+"_nphiHoleLayers"]->at(isrc1).size() : 0;
	N2 = (isTPmu2) ? tpmu_vvi[src2+"_nphiHoleLayers"]->at(isrc2).size() : 0;
	N3 = (isTPmu3) ? tpmu_vvi[src3+"_nphiHoleLayers"]->at(isrc3).size() : 0;
	if(isMuon1 && bugMU_doExsit1) N1 = (!src1.Contains("Muid")) ? muons_nphiHoleLayers->at(bugMU[isrc1]).size() : muid_nphiHoleLayers->at(bugMU[isrc1]).size();
	if(isMuon2 && bugMU_doExsit2) N2 = (!src2.Contains("Muid")) ? muons_nphiHoleLayers->at(bugMU[isrc2]).size() : muid_nphiHoleLayers->at(bugMU[isrc2]).size();
	if(isMuon3 && bugMU_doExsit3) N3 = (!src3.Contains("Muid")) ? muons_nphiHoleLayers->at(bugMU[isrc3]).size() : muid_nphiHoleLayers->at(bugMU[isrc3]).size();
	m_trkPhiHoleLayers[0] = 0;
	m_trkPhiHoleLayers[1] = 0;
	m_trkPhiHoleLayers[2] = 0;
	if(isTPmu1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPhiHoleLayers[0] += tpmu_vvi[src1+"_nphiHoleLayers"]->at(isrc1)[i]; }
	if(isTPmu2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPhiHoleLayers[1] += tpmu_vvi[src2+"_nphiHoleLayers"]->at(isrc2)[i]; }
	if(isTPmu3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPhiHoleLayers[2] += tpmu_vvi[src3+"_nphiHoleLayers"]->at(isrc3)[i]; }
	if(isMuon1 && bugMU_doExsit1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPhiHoleLayers[0] += (!src1.Contains("Muid")) ? muons_nphiHoleLayers->at(bugMU[isrc1])[i] : muid_nphiHoleLayers->at(bugMU[isrc1])[i]; }
	if(isMuon2 && bugMU_doExsit2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPhiHoleLayers[1] += (!src2.Contains("Muid")) ? muons_nphiHoleLayers->at(bugMU[isrc2])[i] : muid_nphiHoleLayers->at(bugMU[isrc2])[i]; }
	if(isMuon3 && bugMU_doExsit3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPhiHoleLayers[2] += (!src3.Contains("Muid")) ? muons_nphiHoleLayers->at(bugMU[isrc3])[i] : muid_nphiHoleLayers->at(bugMU[isrc3])[i]; }
	
	_DEBUG("");
	
	N1 = (isTPmu1) ? tpmu_vvi[src1+"_nprecisionOutliers"]->at(isrc1).size() : 0;
	N2 = (isTPmu2) ? tpmu_vvi[src2+"_nprecisionOutliers"]->at(isrc2).size() : 0;
	N3 = (isTPmu3) ? tpmu_vvi[src3+"_nprecisionOutliers"]->at(isrc3).size() : 0;
	if(isMuon1 && bugMU_doExsit1) N1 = (!src1.Contains("Muid")) ? muons_nprecisionOutliers->at(bugMU[isrc1]).size() : muid_nprecisionOutliers->at(bugMU[isrc1]).size();
	if(isMuon2 && bugMU_doExsit2) N2 = (!src2.Contains("Muid")) ? muons_nprecisionOutliers->at(bugMU[isrc2]).size() : muid_nprecisionOutliers->at(bugMU[isrc2]).size();
	if(isMuon3 && bugMU_doExsit3) N3 = (!src3.Contains("Muid")) ? muons_nprecisionOutliers->at(bugMU[isrc3]).size() : muid_nprecisionOutliers->at(bugMU[isrc3]).size();
	m_trkPrecisionOutliers[0] = 0;
	m_trkPrecisionOutliers[1] = 0;
	m_trkPrecisionOutliers[2] = 0;
	if(isTPmu1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPrecisionOutliers[0] += tpmu_vvi[src1+"_nprecisionOutliers"]->at(isrc1)[i]; }
	if(isTPmu2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPrecisionOutliers[1] += tpmu_vvi[src2+"_nprecisionOutliers"]->at(isrc2)[i]; }
	if(isTPmu3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPrecisionOutliers[2] += tpmu_vvi[src3+"_nprecisionOutliers"]->at(isrc3)[i]; }
	if(isMuon1 && bugMU_doExsit1) { for(unsigned int i=0 ; i<N1 ; ++i) m_trkPrecisionOutliers[0] += (!src1.Contains("Muid")) ? muons_nprecisionOutliers->at(bugMU[isrc1])[i] : muid_nprecisionOutliers->at(bugMU[isrc1])[i]; }
	if(isMuon2 && bugMU_doExsit2) { for(unsigned int i=0 ; i<N2 ; ++i) m_trkPrecisionOutliers[1] += (!src2.Contains("Muid")) ? muons_nprecisionOutliers->at(bugMU[isrc2])[i] : muid_nprecisionOutliers->at(bugMU[isrc2])[i]; }
	if(isMuon3 && bugMU_doExsit3) { for(unsigned int i=0 ; i<N3 ; ++i) m_trkPrecisionOutliers[2] += (!src3.Contains("Muid")) ? muons_nprecisionOutliers->at(bugMU[isrc3])[i] : muid_nprecisionOutliers->at(bugMU[isrc3])[i]; }
	
	_DEBUG("");


if(!skim)
{	
	TString master = mastername;
	
	//// Fill the 4th tracks, allow only 3 per triplet and fill first the muons
	// cout << "\n------------- Event: " << EventNumber << " -------------" << endl;
	/////////////////////////
	vector<int> vtmp;
	isFourthTrack.insert(make_pair(vtx,vtmp));
	for(unsigned int t=1 ; t<=nMaxTracks ; ++t)
	{
		if(t<=3) isFourthTrack[vtx].push_back(1);
		else     isFourthTrack[vtx].push_back(0);
	}
	/////////////////////////
	double margins = 30.*GeV2MeV;
	TString srctp = "CombinedFitMuonParticles";
	int nMuons = (int)muons_pt->size();
	int nTPas  = (int)tpmu_vd[srctp+"_pt"]->size();
	vector<int> vfourthmuons;
	vector<int> vfourthtps;
	if(nMuons>3 || nTPas>3)
	{
		for(int m=0 ; m<nMuons ; ++m)
		{
			if(isMuon1 && m==isrc1) continue;
			if(isMuon2 && m==isrc2) continue;
			if(isMuon3 && m==isrc3) continue;
			
			int itrk = (mastername!="muid") ? muons_inDetTrackIndex->at(m) : muid_inDetTrackIndex->at(m);
			if(itrk<0 || itrk>=(int)trks_pt->size()) continue;
			
			TLorentzVector vZ, v4th;
			v4th.SetPtEtaPhiM(muons_pt->at(m),muons_eta->at(m),muons_phi->at(m),muonMassMeV);
			vZ = psum+v4th;
			// double mQuad = vZ.M();
			if(fabs(vZ.M()-91.*GeV2MeV)>margins) continue;
			
			vfourthmuons.push_back(m);
			// cout << "\tEvent: " << EventNumber << "(" << isrc1 << "," << isrc2 << "," << isrc3 << "+" << m << ") -> with 4th muon mQuad=" << mQuad << endl;
		}
		
		for(int m=0 ; m<nTPas ; ++m)
		{
			if(isTPa1 && m==isrc1) continue;
			if(isTPa2 && m==isrc2) continue;
			if(isTPa3 && m==isrc3) continue;
			
			bool isoverlap = false;
			for(unsigned int n=0 ; n<vfourthmuons.size() ; ++n)
			{
				unsigned int m4 = vfourthmuons[n];
				double dR  = deltaR(muons_eta->at(m4),muons_phi->at(m4),tpmu_vd[srctp+"_eta"]->at(m),tpmu_vd[srctp+"_phi"]->at(m));
				double dpT = fabs(muons_pt->at(m4)-tpmu_vd[srctp+"_pt"]->at(m))/muons_pt->at(m4);
				double dQ  = muons_charge->at(m4)-qtrk(tpmu_vd[srctp+"_qOverP"]->at(m));
				
				if(dR<0.01 && dpT<0.2 && dQ==0) { isoverlap=true; break; }
			}
			if(isoverlap) continue;
			
			TLorentzVector vZ, v4th;
			v4th.SetPtEtaPhiM(tpmu_vd[srctp+"_pt"]->at(m),tpmu_vd[srctp+"_eta"]->at(m),tpmu_vd[srctp+"_phi"]->at(m),muonMassMeV);
			vZ = psum+v4th;
			// double mQuad = vZ.M();
			if(fabs(vZ.M()-91.*GeV2MeV)>margins) continue;
			
			vfourthtps.push_back(m);
			// cout << "\tEvent: " << EventNumber << "(" << isrc1 << "," << isrc2 << "," << isrc3 << "+" << m << ") -> with 4th TPa mQuad=" << mQuad << endl;
		}
	}
	int n4thMuons = (int)vfourthmuons.size();
	int n4thTPs   = (int)vfourthtps.size();
	int nMax4thTracks = nMaxTracks-3;
	int nFilled = 0;
	// if(n4thMuons+n4thTPs>0) cout << "\tn4thMuons=" << n4thMuons << ", n4thTPs=" << n4thTPs << endl;
	for(int i=0 ; (i<n4thMuons && (2+i+1)<(3+nMax4thTracks)) ; ++i) 
	{
		unsigned int ithis = 2+(i+1);
		
		unsigned int m4 = vfourthmuons[i];
		
		int itrk = (master!="muid") ? muons_inDetTrackIndex->at(m4) : muid_inDetTrackIndex->at(m4);
		// if(itrk<0 || itrk>=(int)trks_pt->size()) continue;
		
		isFourthTrack[vtx][ithis] = 1;
		
		TLorentzVector p4;
		p4.SetPtEtaPhiM(muons_pt->at(m4),muons_eta->at(m4),muons_phi->at(m4),muonMassMeV);
		m_trkP[ithis] = p4;

		m_pQuad[ithis-3]           = p4+psum;
		m_isrc[ithis]              = m4;
		m_src[ithis]               = (string)master;
		m_order[ithis]             = -1;
		m_trktype[ithis]           = MUON;
		m_ismuon[ithis]            = 1;
		m_iscalo[ithis]            = 0;
		m_istp[ithis]              = 0;  
		m_istpa[ithis]             = 0;  
		m_istpb[ithis]             = 0;  
		m_iscb[ithis]              = (master!="muid") ? muons_isCombined->at(m4)   : muid_isCombined->at(m4);
		m_istight[ithis]           = (master!="muid") ? muons_isTight->at(m4)      : muid_isTight->at(m4);
		m_ismedium[ithis]          = (master!="muid") ? muons_isMedium->at(m4)     : muid_isMedium->at(m4);
		m_isloose[ithis]           = (master!="muid") ? muons_isLoose->at(m4)      : muid_isLoose->at(m4);
		m_trksctang[ithis]         = (master!="muid") ? muons_sctangsig->at(m4)    : muid_sctangsig->at(m4);
		m_trksctngb[ithis]         = (master!="muid") ? muons_sctngbsig->at(m4)    : muid_sctngbsig->at(m4);
		m_trkpbal[ithis]           = (master!="muid") ? muons_pbalsig->at(m4)      : muid_pbalsig->at(m4);
		m_trkMuMatchChi2Ndf[ithis] = (master!="muid") ? muons_matchchi2ndf->at(m4) : muid_matchchi2ndf->at(m4);
		m_trkChi2[ithis]           = (master!="muid") ? muons_chi2->at(m4)         : muid_chi2->at(m4); // tpmu_vd[src+"_chi2"]->at(m4);  
		m_trkNdf[ithis]            = (master!="muid") ? muons_ndf->at(m4)          : muid_ndf->at(m4); // tpmu_vi[src+"_ndf"]->at(m4);  
		m_trkPval[ithis]           = TMath::Prob(m_trkChi2[ithis],m_trkNdf[ithis]);
		m_trkChi2Ndf[ithis]        = m_trkChi2[ithis]/m_trkNdf[ithis];

		TVector3 pme;
		double q=-999.;
		if(master!="muid") { pme.SetXYZ(muons_px_me->at(m4),muons_py_me->at(m4),muons_pz_me->at(m4)); q = muons_charge->at(m4); }
		else               { pme.SetXYZ(muid_px_me->at(m4), muid_py_me->at(m4), muid_pz_me->at(m4));  q = muid_charge->at(m4);  }
		double qopmemu = q/pme.Mag();
		m_srcQoverP[ithis] = qopmemu;

		m_itrk[ithis]            = itrk;
		m_trkQoverP[ithis]       = trks_qoverp->at(itrk);
		m_trkPixeldEdx[ithis]    = trks_pixeldEdx->at(itrk);
		m_trkUsedHitsdEdx[ithis] = trks_nUsedHitsdEdx->at(itrk);
		m_trkPIXhits[ithis]      = trks_nPix->at(itrk);
		m_trkDeadPIX[ithis]      = trks_nDeadPixels->at(itrk);
		m_trkPIXholes[ithis]     = trks_nPixHoles->at(itrk);
		m_trkSCThits[ithis]      = trks_nSCT->at(itrk);
		m_trkDeadSCT[ithis]      = trks_nDeadSCT->at(itrk);
		m_trkSCTholes[ithis]     = trks_nSCTHoles->at(itrk);
		m_trkTRThits[ithis]      = trks_nTRT->at(itrk);
		m_trkTRToutliers[ithis]  = trks_nTRTOutliers->at(itrk);
		m_trkHtTRThits[ithis]    = trks_nHighThresholdTRTHits->at(itrk);

		m_trkMDThits[ithis]         = -1; // tpmu_vi[src+"_numberOfMdtHits"]->at(m4);
		m_trkTGCPhiHits[ithis]      = -1; // tpmu_vi[src+"_numberOfTgcPhiHits"]->at(m4);
		m_trkTGCEtaHits[ithis]      = -1; // tpmu_vi[src+"_numberOfTgcEtaHits"]->at(m4);
		m_trkCSCPhiHits[ithis]      = -1; // tpmu_vi[src+"_numberOfCscPhiHits"]->at(m4);
		m_trkCSCEtaHits[ithis]      = -1; // tpmu_vi[src+"_numberOfCscEtaHits"]->at(m4);
		m_trkRPCPhiHits[ithis]      = -1; // tpmu_vi[src+"_numberOfRpcPhiHits"]->at(m4);
		m_trkRPCEtaHits[ithis]      = -1; // tpmu_vi[src+"_numberOfRpcEtaHits"]->at(m4);
		m_trkCSCEtaHoles[ithis]     = -1; // tpmu_vi[src+"_numberOfCscEtaHoles"]->at(m4);
		m_trkCSCPhiHoles[ithis]     = -1; // tpmu_vi[src+"_numberOfCscPhiHoles"]->at(m4);
		m_trkRPCEtaHoles[ithis]     = -1; // tpmu_vi[src+"_numberOfRpcEtaHoles"]->at(m4);
		m_trkRPCPhiHoles[ithis]     = -1; // tpmu_vi[src+"_numberOfRpcPhiHoles"]->at(m4);
		m_trkTGCEtaHoles[ithis]     = -1; // tpmu_vi[src+"_numberOfTgcEtaHoles"]->at(m4);
		m_trkTGCPhiHoles[ithis]     = -1; // tpmu_vi[src+"_numberOfTgcPhiHoles"]->at(m4);
		m_trkMDTholes[ithis]        = -1; // tpmu_vi[src+"_numberOfMdtHoles"]->at(m4);
		m_trkOutliersOnTrack[ithis] = -1; // tpmu_vi[src+"_numberOfOutliersOnTrack"]->at(m4);
		m_trkStdDevOfChi2OS[ithis]  = -1; // tpmu_vi[src+"_standardDeviationOfChi2OS"]->at(m4);
		
		m_trkPrecisionHits[ithis] = 0;
		// unsigned int N = tpmu_vvi[src+"_nprecisionHits"]->at(mr).size();
		// for(unsigned int j=0 ; j<N ; ++j) m_trkPrecisionHits[ithis] += tpmu_vvi[src+"_nprecisionHits"]->at(m4)[j];
		
		m_trkPhiLayers[ithis] = 0;
		// N = tpmu_vvi[src+"_nphiLayers"]->at(m4).size();
		// for(unsigned int j=0 ; j<N ; ++j) m_trkPhiLayers[ithis] += tpmu_vvi[src+"_nphiLayers"]->at(m4)[j];

		m_trkEtaPhiLayers[ithis] = 0;
		// N = tpmu_vvi[src+"_netaPhiLayers"]->at(m4).size();
		// for(unsigned int j=0 ; j<N ; ++j) m_trkEtaPhiLayers[ithis] += tpmu_vvi[src+"_netaPhiLayers"]->at(m4)[j];
		
		m_trkPrecisionHoles[ithis] = 0;
		// N = tpmu_vvi[src+"_nprecisionHoles"]->at(m4).size();
		// for(unsigned int j=0 ; j<N ; ++j) m_trkPrecisionHoles[ithis] += tpmu_vvi[src+"_nprecisionHoles"]->at(m4)[j];
	
		m_trkEtaTriggerHoleLayers[ithis] = 0;
		// N = tpmu_vvi[src+"_netaTriggerHoleLayers"]->at(m4).size();
		// for(unsigned int j=0 ; j<N ; ++j) m_trkEtaTriggerHoleLayers[ithis] += tpmu_vvi[src+"_netaTriggerHoleLayers"]->at(m4)[j];

		m_trkPhiHoleLayers[ithis] = 0;
		// N = tpmu_vvi[src+"_nphiHoleLayers"]->at(m4).size();
		// for(unsigned int j=0 ; j<N ; ++j) m_trkPhiHoleLayers[ithis] += tpmu_vvi[src+"_nphiHoleLayers"]->at(m4)[j];

		m_trkPrecisionOutliers[ithis] = 0;
		// N = tpmu_vvi[src+"_nprecisionOutliers"]->at(m4).size();
		// for(unsigned int j=0 ; j<N ; ++j) m_trkPrecisionOutliers[ithis] += tpmu_vvi[src+"_nprecisionOutliers"]->at(m4)[j];
		
		//////////////
		nFilled++; ///
		//////////////
		
		// cout << "\t Adding muon: ithis=" << ithis << ", m4=" << m4 << ", itrk=" << itrk << ", mQuad=" << (p4+psum).M() << ", nFilled=" << nFilled << endl;
	}
	// cout << "\tAdded " << nFilled << " muons" << endl;
	
	int nTPaFilled = 0;
	if(nFilled<nMax4thTracks)
	{
		// cout << "\tGoing to add " << nMax4thTracks-nFilled << " TPas" << endl;
		
		for(int i=0 ; (i<n4thTPs && (2+nFilled+i+1)<(3+nMax4thTracks)) ; ++i) 
		{
			unsigned int ithis = 2+nFilled+(i+1);
			
			unsigned int m4 = vfourthtps[i];
			TLorentzVector p4;
			p4.SetPtEtaPhiM(tpmu_vd[srctp+"_pt"]->at(m4),tpmu_vd[srctp+"_eta"]->at(m4),tpmu_vd[srctp+"_phi"]->at(m4),muonMassMeV);
			m_trkP[ithis] = p4;
			
			isFourthTrack[vtx][ithis] = 1;
			
			m_pQuad[ithis-3]           = p4+psum;
			m_isrc[ithis]              = m4;
			m_src[ithis]               = "";
			m_order[ithis]             = -1;
			m_trktype[ithis]           = TPA;
			m_ismuon[ithis]            = 0;
			m_iscalo[ithis]            = 0;
			m_istp[ithis]              = 1;  
			m_istpa[ithis]             = 1;  
			m_istpb[ithis]             = 0;  
			m_iscb[ithis]              = -1;
			m_istight[ithis]           = -1;
			m_ismedium[ithis]          = -1;
			m_isloose[ithis]           = -1;
			m_trksctang[ithis]         = -1;
			m_trksctngb[ithis]         = -1;
			m_trkpbal[ithis]           = -1;
			m_trkMuMatchChi2Ndf[ithis] = -1;
			m_trkChi2[ithis]           = tpmu_vd[srctp+"_chi2"]->at(m4);  
			m_trkNdf[ithis]            = tpmu_vi[srctp+"_ndf"]->at(m4);  
			m_trkPval[ithis]           = TMath::Prob(m_trkChi2[ithis],m_trkNdf[ithis]);
			m_trkChi2Ndf[ithis]        = m_trkChi2[ithis]/m_trkNdf[ithis];
			m_srcQoverP[ithis]         = tpmu_vd[srctp+"_qOverP"]->at(m4);
			
			m_itrk[ithis]            = -1; // itrk;
			m_trkQoverP[ithis]       = tpmu_vd[srctp+"_qOverP"]->at(m4); // (itrk<0) ? -1 : trks_qoverp->at(itrk);
			m_trkPixeldEdx[ithis]    = tpmu_vd[srctp+"_pixeldEdx"]->at(m4); // (itrk<0) ? -1 : trks_pixeldEdx->at(itrk);
			m_trkUsedHitsdEdx[ithis] = tpmu_vi[srctp+"_nUsedHitsdEdx"]->at(m4); // (itrk<0) ? -1 : trks_nUsedHitsdEdx->at(itrk);
			m_trkPIXhits[ithis]      = tpmu_vi[srctp+"_nPix"]->at(m4); // (itrk<0) ? -1 : trks_nPix->at(itrk);
			m_trkDeadPIX[ithis]      = tpmu_vi[srctp+"_nDeadPixels"]->at(m4); // (itrk<0) ? -1 : trks_nDeadPixels->at(itrk);
			m_trkPIXholes[ithis]     = tpmu_vi[srctp+"_nPixHoles"]->at(m4); // (itrk<0) ? -1 : trks_nPixHoles->at(itrk);
			m_trkSCThits[ithis]      = tpmu_vi[srctp+"_nSCT"]->at(m4); // (itrk<0) ? -1 : trks_nSCT->at(itrk);
			m_trkDeadSCT[ithis]      = tpmu_vi[srctp+"_nDeadSCT"]->at(m4); // (itrk<0) ? -1 : trks_nDeadSCT->at(itrk);
			m_trkSCTholes[ithis]     = tpmu_vi[srctp+"_nSCTHoles"]->at(m4); // (itrk<0) ? -1 : trks_nSCTHoles->at(itrk);
			m_trkTRThits[ithis]      = tpmu_vi[srctp+"_nTRT"]->at(m4); // (itrk<0) ? -1 : trks_nTRT->at(itrk);
			m_trkTRToutliers[ithis]  = tpmu_vi[srctp+"_nTRTOutliers"]->at(m4); // (itrk<0) ? -1 : trks_nTRTOutliers->at(itrk);
			m_trkHtTRThits[ithis]    = tpmu_vi[srctp+"_nHighThresholdTRTHits"]->at(m4); // (itrk<0) ? -1 : trks_nHighThresholdTRTHits->at(itrk);

			m_trkMDThits[ithis]         = tpmu_vi[srctp+"_numberOfMdtHits"]->at(m4);
			m_trkTGCPhiHits[ithis]      = tpmu_vi[srctp+"_numberOfTgcPhiHits"]->at(m4);
			m_trkTGCEtaHits[ithis]      = tpmu_vi[srctp+"_numberOfTgcEtaHits"]->at(m4);
			m_trkCSCPhiHits[ithis]      = tpmu_vi[srctp+"_numberOfCscPhiHits"]->at(m4);
			m_trkCSCEtaHits[ithis]      = tpmu_vi[srctp+"_numberOfCscEtaHits"]->at(m4);
			m_trkRPCPhiHits[ithis]      = tpmu_vi[srctp+"_numberOfRpcPhiHits"]->at(m4);
			m_trkRPCEtaHits[ithis]      = tpmu_vi[srctp+"_numberOfRpcEtaHits"]->at(m4);
			m_trkCSCEtaHoles[ithis]     = tpmu_vi[srctp+"_numberOfCscEtaHoles"]->at(m4);
			m_trkCSCPhiHoles[ithis]     = tpmu_vi[srctp+"_numberOfCscPhiHoles"]->at(m4);
			m_trkRPCEtaHoles[ithis]     = tpmu_vi[srctp+"_numberOfRpcEtaHoles"]->at(m4);
			m_trkRPCPhiHoles[ithis]     = tpmu_vi[srctp+"_numberOfRpcPhiHoles"]->at(m4);
			m_trkTGCEtaHoles[ithis]     = tpmu_vi[srctp+"_numberOfTgcEtaHoles"]->at(m4);
			m_trkTGCPhiHoles[ithis]     = tpmu_vi[srctp+"_numberOfTgcPhiHoles"]->at(m4);
			m_trkMDTholes[ithis]        = tpmu_vi[srctp+"_numberOfMdtHoles"]->at(m4);
			m_trkOutliersOnTrack[ithis] = tpmu_vi[srctp+"_numberOfOutliersOnTrack"]->at(m4);
			m_trkStdDevOfChi2OS[ithis]  = tpmu_vi[srctp+"_standardDeviationOfChi2OS"]->at(m4);

			_DEBUG("");

			m_trkPrecisionHits[ithis] = 0;
			unsigned int N = tpmu_vvi[srctp+"_nprecisionHits"]->at(m4).size();
			for(unsigned int j=0 ; j<N ; ++j) m_trkPrecisionHits[ithis] += tpmu_vvi[srctp+"_nprecisionHits"]->at(m4)[j];

			m_trkPhiLayers[ithis] = 0;
			N = tpmu_vvi[srctp+"_nphiLayers"]->at(m4).size();
			for(unsigned int j=0 ; j<N ; ++j) m_trkPhiLayers[ithis] += tpmu_vvi[srctp+"_nphiLayers"]->at(m4)[j];

			m_trkEtaPhiLayers[ithis] = 0;
			N = tpmu_vvi[srctp+"_netaPhiLayers"]->at(m4).size();
			for(unsigned int j=0 ; j<N ; ++j) m_trkEtaPhiLayers[ithis] += tpmu_vvi[srctp+"_netaPhiLayers"]->at(m4)[j];

			m_trkPrecisionHoles[ithis] = 0;
			N = tpmu_vvi[srctp+"_nprecisionHoles"]->at(m4).size();
			for(unsigned int j=0 ; j<N ; ++j) m_trkPrecisionHoles[ithis] += tpmu_vvi[srctp+"_nprecisionHoles"]->at(m4)[j];

			m_trkEtaTriggerHoleLayers[ithis] = 0;
			N = tpmu_vvi[srctp+"_netaTriggerHoleLayers"]->at(m4).size();
			for(unsigned int j=0 ; j<N ; ++j) m_trkEtaTriggerHoleLayers[ithis] += tpmu_vvi[srctp+"_netaTriggerHoleLayers"]->at(m4)[j];

			m_trkPhiHoleLayers[ithis] = 0;
			N = tpmu_vvi[srctp+"_nphiHoleLayers"]->at(m4).size();
			for(unsigned int j=0 ; j<N ; ++j) m_trkPhiHoleLayers[ithis] += tpmu_vvi[srctp+"_nphiHoleLayers"]->at(m4)[j];

			m_trkPrecisionOutliers[ithis] = 0;
			N = tpmu_vvi[srctp+"_nprecisionOutliers"]->at(m4).size();
			for(unsigned int j=0 ; j<N ; ++j) m_trkPrecisionOutliers[ithis] += tpmu_vvi[srctp+"_nprecisionOutliers"]->at(m4)[j];
			
			_DEBUG("");
			
			/////////////////
			nTPaFilled++; ///
			/////////////////
			
			// cout << "\t Adding TPa: ithis=" << ithis << ", m4=" << m4 << ", mQuad=" << (p4+psum).M() << ", nTPaFilled=" << nTPaFilled << endl;
		}
		
		// cout << "\tAdded " << nTPaFilled << " TPas" << endl;
	}
	// cout << "\tAdded " << nFilled+nTPaFilled << " Muons+TPas" << endl;
	// cout << "--------------------------------------------" << endl;
}
	_DEBUG("");
	
	
	m_trkptfrac[0] = p1.Pt()/psum.Pt();
	m_trkptfrac[1] = p2.Pt()/psum.Pt();
	m_trkptfrac[2] = p3.Pt()/psum.Pt();
	
	_DEBUG("");
	
	double pt1 = 0;
	double pt2 = 0;
	double pt3 = 0;
	float dpt12, dpt13, dpt23, ptfraction12, ptfraction13, ptfraction23;
	if     (order1==1) { pt1=src.srcTlv[0].Pt(); }
	else if(order1==2) { pt2=src.srcTlv[0].Pt(); }
	else if(order1==3) { pt3=src.srcTlv[0].Pt(); }
	if     (order2==1) { pt1=src.srcTlv[1].Pt(); }
	else if(order2==2) { pt2=src.srcTlv[1].Pt(); }
	else if(order2==3) { pt3=src.srcTlv[1].Pt(); }
	if     (order3==1) { pt1=src.srcTlv[2].Pt(); }
	else if(order3==2) { pt2=src.srcTlv[2].Pt(); }
	else if(order3==3) { pt3=src.srcTlv[2].Pt(); }
	dpt12 = fabs(pt1-pt2);
	dpt23 = fabs(pt2-pt3);
	dpt13 = fabs(pt1-pt3);
	ptfraction12 = dpt12/psum.Pt(); 
	ptfraction23 = dpt23/psum.Pt(); 
	ptfraction13 = dpt13/psum.Pt();
	m_dpt12 = dpt12;
	m_dpt13 = dpt23;
	m_dpt23 = dpt13;
	m_ptFrac12 = ptfraction12;
	m_ptFrac13 = ptfraction13;
	m_ptFrac23 = ptfraction23;
	
	_DEBUG("");
	
	double dRmin = +1.e20;
	dRmin = (psum.DeltaR(p1)<dRmin) ? psum.DeltaR(p1) : dRmin;
	dRmin = (psum.DeltaR(p2)<dRmin) ? psum.DeltaR(p2) : dRmin;
	dRmin = (psum.DeltaR(p3)<dRmin) ? psum.DeltaR(p3) : dRmin;
	double dRmax = -1.e20;
	dRmax = (psum.DeltaR(p1)>dRmax) ? psum.DeltaR(p1) : dRmax;
	dRmax = (psum.DeltaR(p2)>dRmax) ? psum.DeltaR(p2) : dRmax;
	dRmax = (psum.DeltaR(p3)>dRmax) ? psum.DeltaR(p3) : dRmax;
	m_drmax = dRmax;
	m_drmin = dRmin;
	
	_DEBUG("");
	
	m_met          = MET_RefFinal_et;
	m_metPhi       = MET_RefFinal_phi;
	m_metDphi3body = fabs(dPhi(MET_RefFinal_phi,psum.Phi()));
	m_metMt        = mT(MET_RefFinal_et,MET_RefFinal_phi,psum.Pt(),psum.Phi());
	
	_DEBUG("");

	multimap<double,int> pt2i;
	vector<int> ijet;
	// for(int i=0 ; i<(int)jets_pt->size() ; i++)
	// {
	// 	bool isUgly = (jets_isUgly->at(i));
	// 	if(isUgly)  continue;
	for(int i=0 ; i<AntiKt4LCTopoLets_n ; i++)
	{
		if(!isGoodJet(i,JetQuality))  continue;
		pt2i.insert(make_pair(AntiKt4LCTopoLets_pt->at(i),i));
	}
	for(multimap<double,int>::reverse_iterator rit=pt2i.rbegin() ; rit!=pt2i.rend() ; ++rit) ijet.push_back(rit->second);
	unsigned int njet = ijet.size();

	double sumptj12     = -999.;
	double dPhiJet1Jet2 = -999.; double dRJet1Jet2 = -999.;
	double dPhi3muJet1  = -999.; double dR3muJet1  = -999.;
	if(njet>0)
	{
		int j1 = ijet[0];
		double phi1 = AntiKt4LCTopoLets_phi->at(j1);	
		dPhi3muJet1 = fabs(dPhi(psum.Phi(),phi1));
		dR3muJet1   = deltaR(psum.Eta(),psum.Phi(),AntiKt4LCTopoLets_eta->at(j1),AntiKt4LCTopoLets_phi->at(j1));
		if(njet>1)
		{
			int j2 = ijet[1];
			sumptj12 = AntiKt4LCTopoLets_pt->at(j1)+AntiKt4LCTopoLets_pt->at(j2);
			double phi2  = AntiKt4LCTopoLets_phi->at(j2);
			dPhiJet1Jet2 = fabs(dPhi(phi1,phi2));
			dRJet1Jet2   = deltaR(AntiKt4LCTopoLets_eta->at(j1),AntiKt4LCTopoLets_phi->at(j1),AntiKt4LCTopoLets_eta->at(j2),AntiKt4LCTopoLets_phi->at(j2));
		}
	}
	
	_DEBUG("");

	m_njets = (njet<5) ? njet : 4;
	for(int i=0 ; i<m_njets ; ++i)
	{
		int j = ijet[i];
		m_jetPE[i].SetPtEtaPhiE(AntiKt4LCTopoLets_pt->at(j),AntiKt4LCTopoLets_eta->at(j),AntiKt4LCTopoLets_phi->at(j),AntiKt4LCTopoLets_E->at(j));
		m_jetPM[i].SetPtEtaPhiM(AntiKt4LCTopoLets_pt->at(j),AntiKt4LCTopoLets_eta->at(j),AntiKt4LCTopoLets_phi->at(j),AntiKt4LCTopoLets_m->at(j));
		m_jetMV1[i] = AntiKt4LCTopoLets_flavor_weight_MV1->at(j);
		m_jetVtxFrac[i] = AntiKt4LCTopoLets_jvtxf->at(j);
	}
	if(m_njets>1) { m_jetDphi12 = dPhiJet1Jet2; m_jetDR12 = dRJet1Jet2; m_jetSumpt12 = sumptj12; }
	if(m_njets>0) { m_jetDphi3body = dPhi3muJet1; m_jetDR3body = dR3muJet1; }
	
	_DEBUG("");
}


bool acceptMuons(unsigned int vtx, TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight=1., double mBlindMin=1500., double mBlindMax=2000.)
{
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];	
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3,psum;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
	psum = p1+p2+p3; // psum = getTlv3mu(vtx);
	
	_DEBUG("");
	
	double mass   = vtx_mass->at(vtx); // psum.M();  in principle
	double pTsum  = vtx_pt->at(vtx);   // psum.Pt(); in principle
	double charge = vtx_charge->at(vtx);
	double pvalue = TMath::Prob(vtx_chi2->at(vtx),vtx_ndf->at(vtx));
	
	_DEBUG("");
	
	sources src;
	getSrc(vtx,src);
	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
	TString shortType = classifyTripletShort(vtx);
	VtxType = classifyTripletCode(shortType);
	
	// // doublets
	// TLorentzVector pOS1, pOS2, pSS;
	// if(fabs(src.q3body)==1.)
	// {
	// 	if(q2body12==0. && q2body13==0. && q2body23!=0.) { pOS1 = src.p2body12; pOS2 = src.p2body13; pSS = src.p2body23; }
	// 	if(q2body12==0. && q2body13!=0. && q2body23==0.) { pOS1 = src.p2body12; pOS2 = src.p2body23; pSS = src.p2body13; }
	// 	if(q2body12!=0. && q2body13==0. && q2body23==0.) { pOS1 = src.p2body13; pOS2 = src.p2body23; pSS = src.p2body12; }
	// }

	// match the tracks
	int itrk1 = src.trkIndex[0];
	int itrk2 = src.trkIndex[1];
	int itrk3 = src.trkIndex[2];
	// match the sources        
	int isrc1 = src.srcIndex[0];
	int isrc2 = src.srcIndex[1];
	int isrc3 = src.srcIndex[2];
	// match the sources        
	TString src1 = src.srcName[0];
	TString src2 = src.srcName[1];
	TString src3 = src.srcName[2];
	// test if muons
	bool isMuon1 = src.isMuon[0];
	bool isMuon2 = src.isMuon[1];
	bool isMuon3 = src.isMuon[2];
	// // test if calo muons
	// bool isCalo1 = src.isCalo[0];
	// bool isCalo2 = src.isCalo[1];
	// bool isCalo3 = src.isCalo[2];
	// test if TPmuon
	bool isTPmu1 = src.isTPmu[0]; bool isTPa1 = src.isTPa[0]; bool isTPb1 = src.isTPb[0];
	bool isTPmu2 = src.isTPmu[1]; bool isTPa2 = src.isTPa[1]; bool isTPb2 = src.isTPb[1];
	bool isTPmu3 = src.isTPmu[2]; bool isTPa3 = src.isTPa[2]; bool isTPb3 = src.isTPb[2];
	
	_DEBUG("");
		
	bool isCB1 = 0;
	bool isCB2 = 0;
	bool isCB3 = 0;

	//cout << "isMuon1=" << isMuon1 << ", isMuon2=" << isMuon2 << ", isMuon3=" << isMuon3 << endl;
	//if(muons_isCombined) cout << "N=" << muons_isCombined->size() << ", isrc1=" << isrc1 << ", isrc2=" << isrc2 << ", isrc3=" << isrc3 << endl;
	//else cout << "muons_isCombined is NULL" << endl;

	if(isMuon1) isCB1 = (!src1.Contains("Muid")) ? muons_isCombined->at(isrc1) : muid_isCombined->at(isrc1);
	if(isMuon2) isCB2 = (!src2.Contains("Muid")) ? muons_isCombined->at(isrc2) : muid_isCombined->at(isrc2);
	if(isMuon3) isCB3 = (!src3.Contains("Muid")) ? muons_isCombined->at(isrc3) : muid_isCombined->at(isrc3);

	_DEBUG("");

	bool isTight1 = true;
	bool isTight2 = true;
	bool isTight3 = true;
	if(isMuon1) isTight1 = (!src1.Contains("Muid")) ? muons_isTight->at(isrc1) : muid_isTight->at(isrc1);
	if(isMuon2) isTight2 = (!src2.Contains("Muid")) ? muons_isTight->at(isrc2) : muid_isTight->at(isrc2);
	if(isMuon3) isTight3 = (!src3.Contains("Muid")) ? muons_isTight->at(isrc3) : muid_isTight->at(isrc3);
	bool isMedium1 = true;
	bool isMedium2 = true;
	bool isMedium3 = true;
	if(isMuon1) isMedium1 = (!src1.Contains("Muid")) ? muons_isMedium->at(isrc1) : muid_isMedium->at(isrc1);
	if(isMuon2) isMedium2 = (!src2.Contains("Muid")) ? muons_isMedium->at(isrc2) : muid_isMedium->at(isrc2);
	if(isMuon3) isMedium3 = (!src3.Contains("Muid")) ? muons_isMedium->at(isrc3) : muid_isMedium->at(isrc3);
	bool isLoose1 = true;
	bool isLoose2 = true;
	bool isLoose3 = true;
	if(isMuon1) isLoose1 = (!src1.Contains("Muid")) ? muons_isLoose->at(isrc1) : muid_isLoose->at(isrc1);
	if(isMuon2) isLoose2 = (!src2.Contains("Muid")) ? muons_isLoose->at(isrc2) : muid_isLoose->at(isrc2);
	if(isMuon3) isLoose3 = (!src3.Contains("Muid")) ? muons_isLoose->at(isrc3) : muid_isLoose->at(isrc3);

	_DEBUG("");

	double sctangsig1 = 999.;
	double sctangsig2 = 999.;
	double sctangsig3 = 999.;
	if(isMuon1) sctangsig1 = (!src1.Contains("Muid")) ? muons_sctangsig->at(isrc1) : muid_sctangsig->at(isrc1);
	if(isMuon2) sctangsig2 = (!src2.Contains("Muid")) ? muons_sctangsig->at(isrc2) : muid_sctangsig->at(isrc2);
	if(isMuon3) sctangsig3 = (!src3.Contains("Muid")) ? muons_sctangsig->at(isrc3) : muid_sctangsig->at(isrc3);
	double sctngbsig1 = 999.;
	double sctngbsig2 = 999.;
	double sctngbsig3 = 999.;
	if(isMuon1) sctngbsig1 = (!src1.Contains("Muid")) ? muons_sctngbsig->at(isrc1) : muid_sctngbsig->at(isrc1);
	if(isMuon2) sctngbsig2 = (!src2.Contains("Muid")) ? muons_sctngbsig->at(isrc2) : muid_sctngbsig->at(isrc2);
	if(isMuon3) sctngbsig3 = (!src3.Contains("Muid")) ? muons_sctngbsig->at(isrc3) : muid_sctngbsig->at(isrc3);
	double pbalsig1 = 999.;
	double pbalsig2 = 999.;
	double pbalsig3 = 999.;
	if(isMuon1) pbalsig1 = (!src1.Contains("Muid")) ? muons_pbalsig->at(isrc1) : muid_pbalsig->at(isrc1);
	if(isMuon2) pbalsig2 = (!src2.Contains("Muid")) ? muons_pbalsig->at(isrc2) : muid_pbalsig->at(isrc2);
	if(isMuon3) pbalsig3 = (!src3.Contains("Muid")) ? muons_pbalsig->at(isrc3) : muid_pbalsig->at(isrc3);

	_DEBUG("");

	double fitchi2mu1 = -999.;
	double fitchi2mu2 = -999.;
	double fitchi2mu3 = -999.;
	if(isMuon1) fitchi2mu1 = (!src1.Contains("Muid")) ? muons_chi2->at(isrc1) : muid_chi2->at(isrc1);
	if(isMuon2) fitchi2mu2 = (!src2.Contains("Muid")) ? muons_chi2->at(isrc2) : muid_chi2->at(isrc2);
	if(isMuon3) fitchi2mu3 = (!src3.Contains("Muid")) ? muons_chi2->at(isrc3) : muid_chi2->at(isrc3);
	int fitndfmu1 = -999;
	int fitndfmu2 = -999;
	int fitndfmu3 = -999;
	if(isMuon1) fitndfmu1 = (!src1.Contains("Muid")) ? muons_ndf->at(isrc1) : muid_ndf->at(isrc1);
	if(isMuon2) fitndfmu2 = (!src2.Contains("Muid")) ? muons_ndf->at(isrc2) : muid_ndf->at(isrc2);
	if(isMuon3) fitndfmu3 = (!src3.Contains("Muid")) ? muons_ndf->at(isrc3) : muid_ndf->at(isrc3);
	double fitpvaluemu1 = -999.;
	double fitpvaluemu2 = -999.;
	double fitpvaluemu3 = -999.;
	if(isMuon1) fitpvaluemu1 = TMath::Prob(fitchi2mu1,fitndfmu1);
	if(isMuon2) fitpvaluemu2 = TMath::Prob(fitchi2mu2,fitndfmu2);
	if(isMuon3) fitpvaluemu3 = TMath::Prob(fitchi2mu3,fitndfmu3);

	_DEBUG("");
	
	double pbalsigTP1 = 999.;
	double pbalsigTP2 = 999.;
	double pbalsigTP3 = 999.;
	//if(isTPmu1) pbalsigTP1 = momentumBalanceSig(itrk1,isrc1,src1);
	//if(isTPmu2) pbalsigTP2 = momentumBalanceSig(itrk2,isrc2,src2);
	//if(isTPmu3) pbalsigTP3 = momentumBalanceSig(itrk3,isrc3,src3);

	_DEBUG("");

	double fitchi2TP1 = -999.;
	double fitchi2TP2 = -999.;
	double fitchi2TP3 = -999.;
	if(isTPmu1) fitchi2TP1 = tpmu_vd[src1+"_chi2"]->at(isrc1);  
	if(isTPmu2) fitchi2TP2 = tpmu_vd[src2+"_chi2"]->at(isrc2);  
	if(isTPmu3) fitchi2TP3 = tpmu_vd[src3+"_chi2"]->at(isrc3);  
	double fitndfTP1 = -999.;
	double fitndfTP2 = -999.;
	double fitndfTP3 = -999.;
	if(isTPmu1) fitndfTP1 = tpmu_vi[src1+"_ndf"]->at(isrc1);
	if(isTPmu2) fitndfTP2 = tpmu_vi[src2+"_ndf"]->at(isrc2);
	if(isTPmu3) fitndfTP3 = tpmu_vi[src3+"_ndf"]->at(isrc3);
	double fitpvalueTP1 = -999.;
	double fitpvalueTP2 = -999.;
	double fitpvalueTP3 = -999.;
	if(isTPmu1) fitpvalueTP1 = TMath::Prob(tpmu_vd[src1+"_chi2"]->at(isrc1),tpmu_vi[src1+"_ndf"]->at(isrc1));
	if(isTPmu2) fitpvalueTP2 = TMath::Prob(tpmu_vd[src2+"_chi2"]->at(isrc2),tpmu_vi[src2+"_ndf"]->at(isrc2));
	if(isTPmu3) fitpvalueTP3 = TMath::Prob(tpmu_vd[src3+"_chi2"]->at(isrc3),tpmu_vi[src3+"_ndf"]->at(isrc3));
	
	_DEBUG("");	

	int order1 = src.srcOrder[0];
	int order2 = src.srcOrder[1];
	int order3 = src.srcOrder[2];
	
	_DEBUG("");
	
	double matchchi2ndfmu1 = -999.;
	double matchchi2ndfmu2 = -999.;
	double matchchi2ndfmu3 = -999.;
	if(isMuon1) matchchi2ndfmu1 = (!src1.Contains("Muid")) ? muons_matchchi2ndf->at(isrc1) : muid_matchchi2ndf->at(isrc1);
	if(isMuon2) matchchi2ndfmu2 = (!src2.Contains("Muid")) ? muons_matchchi2ndf->at(isrc2) : muid_matchchi2ndf->at(isrc2);
	if(isMuon3) matchchi2ndfmu3 = (!src3.Contains("Muid")) ? muons_matchchi2ndf->at(isrc3) : muid_matchchi2ndf->at(isrc3);
	matchchi2ndfmu1 = (matchchi2ndfmu1>20.) ? 19.999 : matchchi2ndfmu1;
	matchchi2ndfmu2 = (matchchi2ndfmu2>20.) ? 19.999 : matchchi2ndfmu2;
	matchchi2ndfmu3 = (matchchi2ndfmu3>20.) ? 19.999 : matchchi2ndfmu3;
	
	_DEBUG("");
	
	TVector3 p1me, p2me, p3me;
	TVector3 p1ie, p2ie, p3ie;
	double q1=-999.; double q2=-999.; double q3=-999.;
	if(isMuon1)
	{
		if(!src1.Contains("Muid")) { p1me.SetXYZ(muons_px_me->at(isrc1), muons_py_me->at(isrc1), muons_pz_me->at(isrc1)); q1 = muons_charge->at(isrc1); }
		else                       { p1me.SetXYZ(muid_px_me->at(isrc1),  muid_py_me->at(isrc1),  muid_pz_me->at(isrc1));  q1 = muid_charge->at(isrc1);  }
		if(!src1.Contains("Muid")) { p1ie.SetXYZ(muons_px_ie->at(isrc1), muons_py_ie->at(isrc1), muons_pz_ie->at(isrc1)); q1 = muons_charge->at(isrc1); }
		else                       { p1ie.SetXYZ(muid_px_ie->at(isrc1),  muid_py_ie->at(isrc1),  muid_pz_ie->at(isrc1));  q1 = muid_charge->at(isrc1);  }
	}                                                                                                                
	if(isMuon2)                                                                                                      
	{                                                                                                                
		if(!src2.Contains("Muid")) { p2me.SetXYZ(muons_px_me->at(isrc2), muons_py_me->at(isrc2), muons_pz_me->at(isrc2)); q2 = muons_charge->at(isrc2); }
		else                       { p2me.SetXYZ(muid_px_me->at(isrc2),  muid_py_me->at(isrc2),  muid_pz_me->at(isrc2));  q2 = muid_charge->at(isrc2);   }
		if(!src2.Contains("Muid")) { p2ie.SetXYZ(muons_px_ie->at(isrc2), muons_py_ie->at(isrc2), muons_pz_ie->at(isrc2)); q2 = muons_charge->at(isrc2); }
		else                       { p2ie.SetXYZ(muid_px_ie->at(isrc2),  muid_py_ie->at(isrc2),  muid_pz_ie->at(isrc2));  q2 = muid_charge->at(isrc2);  }
	}                                                                                                                
	if(isMuon3)                                                                                                      
	{                                                                                                                
		if(!src3.Contains("Muid")) { p3me.SetXYZ(muons_px_me->at(isrc3), muons_py_me->at(isrc3), muons_pz_me->at(isrc3)); q3 = muons_charge->at(isrc3); }
		else                       { p3me.SetXYZ(muid_px_me->at(isrc3),  muid_py_me->at(isrc3),  muid_pz_me->at(isrc3));  q3 = muid_charge->at(isrc3);  }
		if(!src3.Contains("Muid")) { p3ie.SetXYZ(muons_px_ie->at(isrc3), muons_py_ie->at(isrc3), muons_pz_ie->at(isrc3)); q3 = muons_charge->at(isrc3); }
		else                       { p3ie.SetXYZ(muid_px_ie->at(isrc3),  muid_py_ie->at(isrc3),  muid_pz_ie->at(isrc3));  q3 = muid_charge->at(isrc3);  }
	}
	double qopmemu1 = -999.;
	double qopmemu2 = -999.;
	double qopmemu3 = -999.;
	double qopiemu1 = -999.;
	double qopiemu2 = -999.;
	double qopiemu3 = -999.;
	if(isMuon1) { qopmemu1 = q1/p1me.Mag(); qopiemu1 = q1/p1ie.Mag(); }
	if(isMuon2) { qopmemu2 = q2/p2me.Mag(); qopiemu2 = q2/p2ie.Mag(); }
	if(isMuon3) { qopmemu3 = q3/p3me.Mag(); qopiemu3 = q3/p3ie.Mag(); }
	if(fabs(qopmemu1)>7.e-4) qopmemu1 = 6.999e-4;
	if(fabs(qopmemu2)>7.e-4) qopmemu2 = 6.999e-4;
	if(fabs(qopmemu3)>7.e-4) qopmemu3 = 6.999e-4;
	if(fabs(qopiemu1)>7.e-4) qopiemu1 = 6.999e-4;
	if(fabs(qopiemu2)>7.e-4) qopiemu2 = 6.999e-4;
	if(fabs(qopiemu3)>7.e-4) qopiemu3 = 6.999e-4;
	if(0) cout << "qopiemu1=" << qopiemu1 << endl;
	if(0) cout << "qopiemu2=" << qopiemu2 << endl;
	if(0) cout << "qopiemu3=" << qopiemu3 << endl;
	
	_DEBUG("");
	
	double qopTP1 = -999.;
	double qopTP2 = -999.;
	double qopTP3 = -999.;
	if(isTPmu1) qopTP1 = tpmu_vd[src1+"_qOverP"]->at(isrc1);
	if(isTPmu2) qopTP2 = tpmu_vd[src2+"_qOverP"]->at(isrc2);
	if(isTPmu3) qopTP3 = tpmu_vd[src3+"_qOverP"]->at(isrc3);
	if(fabs(qopTP1)>7.e-4) qopTP1 = 6.999e-4;
	if(fabs(qopTP2)>7.e-4) qopTP2 = 6.999e-4;
	if(fabs(qopTP3)>7.e-4) qopTP3 = 6.999e-4;
	
	_DEBUG("");
	
	int nHighThresholdTRTHits1 = (isTPmu1) ?  tpmu_vi[src1+"_nHighThresholdTRTHits"]->at(isrc1) : trks_nHighThresholdTRTHits->at(itrk1);        
	int nHighThresholdTRTHits2 = (isTPmu2) ?  tpmu_vi[src2+"_nHighThresholdTRTHits"]->at(isrc2) : trks_nHighThresholdTRTHits->at(itrk2);
	int nHighThresholdTRTHits3 = (isTPmu3) ?  tpmu_vi[src3+"_nHighThresholdTRTHits"]->at(isrc3) : trks_nHighThresholdTRTHits->at(itrk3);
		
	_DEBUG("");
	
	//////////////////////////
	//// fill some histos ////
	//////////////////////////
	// OS and SS pairs
	TLorentzVector pOS1, pOS2, pSS;
	bool foundOS1 = false;
	bool foundOS2 = false;
	bool foundSS  = false;
	int iUniqueCharge = -1;
	TMapVL doubletsTLV = getDoubletsFromTriplet(vtx,iUniqueCharge);
	if(doubletsTLV.find("OS1") != doubletsTLV.end()) { pOS1 = doubletsTLV["OS1"]; foundOS1 = true; }
	if(doubletsTLV.find("OS2") != doubletsTLV.end()) { pOS2 = doubletsTLV["OS2"]; foundOS2 = true; }
	if(doubletsTLV.find("SS")  != doubletsTLV.end()) { pSS  = doubletsTLV["SS"];  foundSS  = true; }
	if(foundOS1+foundOS2<2 && fabs(vtx_charge->at(vtx))==1.) _ERROR("could not find OS doublet"); // should be exactly two OS pairs
	if(!foundSS && fabs(vtx_charge->at(vtx))==1.)            _ERROR("could not find SS doublet"); // should be exactly one SS pair
	double mOS1 = (foundOS1) ? pOS1.M() : -1.;
	// double mOS2 = (foundOS2) ? pOS2.M() : -1.;
	double mSS  = (foundSS)  ? pSS.M()  : -1.;
	// double weightOS = (foundOS1+foundOS2==2) ? weight/2. : weight;
	fillHistsPVandSV(vtx,name,histos,histos2,weight); // some vertex-PV-SV histos
	fillHistsMassPt3mu(vtx,name,"",histos,histos2,weight,mBlindMin,mBlindMax);
	fillHistsDoublets(vtx,name,"",histos,histos2,weight,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS); 
	
	_DEBUG("");
	
	// before the pT cut
	bool Tight1 = 0; bool Medium1 = 0; bool Loose1 = 0; bool Muon1 = 0; bool TPa1 = 0; bool TPb1 = 0;
	bool Tight2 = 0; bool Medium2 = 0; bool Loose2 = 0; bool Muon2 = 0; bool TPa2 = 0; bool TPb2 = 0;
	bool Tight3 = 0; bool Medium3 = 0; bool Loose3 = 0; bool Muon3 = 0; bool TPa3 = 0; bool TPb3 = 0;
	Pt1 = 0.; Eta1 = +3.; Phi1 = TMath::Pi()*1.2; Code1 = -1;
	Pt2 = 0.; Eta2 = +3.; Phi2 = TMath::Pi()*1.2; Code2 = -1;
	Pt3 = 0.; Eta3 = +3.; Phi3 = TMath::Pi()*1.2; Code3 = -1;
	SctAngSig1 = +10.; SctNgbSig1 = +10.; PbalSig1 = +10.; Chi2TrkFit1 = -1.; NdfTrkFit1 = -1.; PvalTrkFit1 = 1.2; QoverP1 = -10.; NhtTRThits1 = -1;
	SctAngSig2 = +10.; SctNgbSig2 = +10.; PbalSig2 = +10.; Chi2TrkFit2 = -1.; NdfTrkFit2 = -1.; PvalTrkFit2 = 1.2; QoverP2 = -10.; NhtTRThits2 = -1;
	SctAngSig3 = +10.; SctNgbSig3 = +10.; PbalSig3 = +10.; Chi2TrkFit3 = -1.; NdfTrkFit3 = -1.; PvalTrkFit3 = 1.2; QoverP3 = -10.; NhtTRThits3 = -1;
	
	
	if     (order1==1) { Code1=src.srcCode[0]; Muon1=isMuon1; TPa1=isTPa1; TPb1=isTPb1; Tight1=isTight1; Medium1=isMedium1; Loose1=isLoose1; Pt1=src.srcTlv[0].Pt(); Eta1=src.srcTlv[0].Eta(); Phi1=src.srcTlv[0].Phi(); SctAngSig1=sctangsig1; SctNgbSig1=sctngbsig1; PbalSig1=pbalsig1; Chi2TrkFit1 = (isMuon1) ? fitchi2mu1 : fitchi2TP1; NdfTrkFit1 = (isMuon1) ? fitndfmu1 : fitndfTP1; PvalTrkFit1 = (isMuon1) ? fitpvaluemu1 : fitpvalueTP1;  QoverP1 = (isTPmu1) ? qopTP1 : qopmemu1; NhtTRThits1 = nHighThresholdTRTHits1; }
	else if(order1==2) { Code2=src.srcCode[0]; Muon2=isMuon1; TPa2=isTPa1; TPb2=isTPb1; Tight2=isTight1; Medium2=isMedium1; Loose2=isLoose1; Pt2=src.srcTlv[0].Pt(); Eta2=src.srcTlv[0].Eta(); Phi2=src.srcTlv[0].Phi(); SctAngSig2=sctangsig1; SctNgbSig2=sctngbsig1; PbalSig2=pbalsig1; Chi2TrkFit2 = (isMuon1) ? fitchi2mu1 : fitchi2TP1; NdfTrkFit2 = (isMuon1) ? fitndfmu1 : fitndfTP1; PvalTrkFit2 = (isMuon1) ? fitpvaluemu1 : fitpvalueTP1;  QoverP2 = (isTPmu1) ? qopTP1 : qopmemu1; NhtTRThits2 = nHighThresholdTRTHits1; }
	else if(order1==3) { Code3=src.srcCode[0]; Muon3=isMuon1; TPa3=isTPa1; TPb3=isTPb1; Tight3=isTight1; Medium3=isMedium1; Loose3=isLoose1; Pt3=src.srcTlv[0].Pt(); Eta3=src.srcTlv[0].Eta(); Phi3=src.srcTlv[0].Phi(); SctAngSig3=sctangsig1; SctNgbSig3=sctngbsig1; PbalSig3=pbalsig1; Chi2TrkFit3 = (isMuon1) ? fitchi2mu1 : fitchi2TP1; NdfTrkFit3 = (isMuon1) ? fitndfmu1 : fitndfTP1; PvalTrkFit3 = (isMuon1) ? fitpvaluemu1 : fitpvalueTP1;  QoverP3 = (isTPmu1) ? qopTP1 : qopmemu1; NhtTRThits3 = nHighThresholdTRTHits1; }
	if     (order2==1) { Code1=src.srcCode[1]; Muon1=isMuon2; TPa1=isTPa2; TPb1=isTPb2; Tight1=isTight2; Medium1=isMedium2; Loose1=isLoose2; Pt1=src.srcTlv[1].Pt(); Eta1=src.srcTlv[1].Eta(); Phi1=src.srcTlv[1].Phi(); SctAngSig1=sctangsig2; SctNgbSig1=sctngbsig2; PbalSig1=pbalsig2; Chi2TrkFit1 = (isMuon2) ? fitchi2mu2 : fitchi2TP2; NdfTrkFit1 = (isMuon2) ? fitndfmu2 : fitndfTP2; PvalTrkFit1 = (isMuon2) ? fitpvaluemu2 : fitpvalueTP2;  QoverP1 = (isTPmu2) ? qopTP2 : qopmemu2; NhtTRThits1 = nHighThresholdTRTHits2; }
	else if(order2==2) { Code2=src.srcCode[1]; Muon2=isMuon2; TPa2=isTPa2; TPb2=isTPb2; Tight2=isTight2; Medium2=isMedium2; Loose2=isLoose2; Pt2=src.srcTlv[1].Pt(); Eta2=src.srcTlv[1].Eta(); Phi2=src.srcTlv[1].Phi(); SctAngSig2=sctangsig2; SctNgbSig2=sctngbsig2; PbalSig2=pbalsig2; Chi2TrkFit2 = (isMuon2) ? fitchi2mu2 : fitchi2TP2; NdfTrkFit2 = (isMuon2) ? fitndfmu2 : fitndfTP2; PvalTrkFit2 = (isMuon2) ? fitpvaluemu2 : fitpvalueTP2;  QoverP2 = (isTPmu2) ? qopTP2 : qopmemu2; NhtTRThits2 = nHighThresholdTRTHits2; }
	else if(order2==3) { Code3=src.srcCode[1]; Muon3=isMuon2; TPa3=isTPa2; TPb3=isTPb2; Tight3=isTight2; Medium3=isMedium2; Loose3=isLoose2; Pt3=src.srcTlv[1].Pt(); Eta3=src.srcTlv[1].Eta(); Phi3=src.srcTlv[1].Phi(); SctAngSig2=sctangsig2; SctNgbSig3=sctngbsig2; PbalSig3=pbalsig2; Chi2TrkFit3 = (isMuon2) ? fitchi2mu2 : fitchi2TP2; NdfTrkFit3 = (isMuon2) ? fitndfmu2 : fitndfTP2; PvalTrkFit3 = (isMuon2) ? fitpvaluemu2 : fitpvalueTP2;  QoverP3 = (isTPmu2) ? qopTP2 : qopmemu2; NhtTRThits3 = nHighThresholdTRTHits2; }
	if     (order3==1) { Code1=src.srcCode[2]; Muon1=isMuon3; TPa1=isTPa3; TPb1=isTPb3; Tight1=isTight3; Medium1=isMedium3; Loose1=isLoose3; Pt1=src.srcTlv[2].Pt(); Eta1=src.srcTlv[2].Eta(); Phi1=src.srcTlv[2].Phi(); SctAngSig1=sctangsig3; SctNgbSig1=sctngbsig3; PbalSig1=pbalsig3; Chi2TrkFit1 = (isMuon3) ? fitchi2mu3 : fitchi2TP3; NdfTrkFit1 = (isMuon3) ? fitndfmu3 : fitndfTP3; PvalTrkFit1 = (isMuon3) ? fitpvaluemu3 : fitpvalueTP3;  QoverP1 = (isTPmu3) ? qopTP3 : qopmemu3; NhtTRThits1 = nHighThresholdTRTHits3; }
	else if(order3==2) { Code2=src.srcCode[2]; Muon2=isMuon3; TPa2=isTPa3; TPb2=isTPb3; Tight2=isTight3; Medium2=isMedium3; Loose2=isLoose3; Pt2=src.srcTlv[2].Pt(); Eta2=src.srcTlv[2].Eta(); Phi2=src.srcTlv[2].Phi(); SctAngSig2=sctangsig3; SctNgbSig2=sctngbsig3; PbalSig2=pbalsig3; Chi2TrkFit2 = (isMuon3) ? fitchi2mu3 : fitchi2TP3; NdfTrkFit2 = (isMuon3) ? fitndfmu3 : fitndfTP3; PvalTrkFit2 = (isMuon3) ? fitpvaluemu3 : fitpvalueTP3;  QoverP2 = (isTPmu3) ? qopTP3 : qopmemu3; NhtTRThits2 = nHighThresholdTRTHits3; }
	else if(order3==3) { Code3=src.srcCode[2]; Muon3=isMuon3; TPa3=isTPa3; TPb3=isTPb3; Tight3=isTight3; Medium3=isMedium3; Loose3=isLoose3; Pt3=src.srcTlv[2].Pt(); Eta3=src.srcTlv[2].Eta(); Phi3=src.srcTlv[2].Phi(); SctAngSig3=sctangsig3; SctNgbSig3=sctngbsig3; PbalSig3=pbalsig3; Chi2TrkFit3 = (isMuon3) ? fitchi2mu3 : fitchi2TP3; NdfTrkFit3 = (isMuon3) ? fitndfmu3 : fitndfTP3; PvalTrkFit3 = (isMuon3) ? fitpvaluemu3 : fitpvalueTP3;  QoverP3 = (isTPmu3) ? qopTP3 : qopmemu3; NhtTRThits3 = nHighThresholdTRTHits3; }

	

	dPt12 = fabs(Pt1-Pt2);
	dPt23 = fabs(Pt2-Pt3);
	dPt13 = fabs(Pt1-Pt3);
	ptFraction1 = Pt1/psum.Pt(); 
	ptFraction2 = Pt2/psum.Pt(); 
	ptFraction3 = Pt3/psum.Pt();
	ptFraction12 = dPt12/psum.Pt(); 
	ptFraction23 = dPt23/psum.Pt(); 
	ptFraction13 = dPt13/psum.Pt();
	
	
	////////////////////////
	// single muons cuts ///
	////////////////////////
	
	_DEBUG("");
	
	histos[name+"_mu_isTight1"]->Fill(Tight1,weight); histos[name+"_mu_isMedium1"]->Fill(Medium1,weight); histos[name+"_mu_isLoose1"]->Fill(Loose1,weight);
	histos[name+"_mu_isTight2"]->Fill(Tight2,weight); histos[name+"_mu_isMedium2"]->Fill(Medium2,weight); histos[name+"_mu_isLoose2"]->Fill(Loose2,weight);
	histos[name+"_mu_isTight3"]->Fill(Tight3,weight); histos[name+"_mu_isMedium3"]->Fill(Medium3,weight); histos[name+"_mu_isLoose3"]->Fill(Loose3,weight);
	histos[name+"_mu_isMedium"]->Fill(Medium1,weight);
	histos[name+"_mu_isMedium"]->Fill(Medium2,weight);
	histos[name+"_mu_isMedium"]->Fill(Medium3,weight);
	
	_DEBUG("");
	
	histos[name+"_mu_pt"]->Fill(Pt1,weight); histos[name+"_mu_eta1"]->Fill(Eta1,weight); histos[name+"_mu_phi1"]->Fill(Phi1,weight);
	histos[name+"_mu_pt"]->Fill(Pt2,weight); histos[name+"_mu_eta2"]->Fill(Eta2,weight); histos[name+"_mu_phi2"]->Fill(Phi2,weight);
	histos[name+"_mu_pt"]->Fill(Pt3,weight); histos[name+"_mu_eta3"]->Fill(Eta3,weight); histos[name+"_mu_phi3"]->Fill(Phi3,weight);
	histos[name+"_mu_type"]->Fill(Code1,weight);
	histos[name+"_mu_type"]->Fill(Code2,weight);
	histos[name+"_mu_type"]->Fill(Code3,weight);
	histos2[name+"_mupt1_type_vs_mupt1"]->Fill(Pt1,Code1,weight);
	histos2[name+"_mupt2_type_vs_mupt2"]->Fill(Pt2,Code2,weight);
	histos2[name+"_mupt3_type_vs_mupt3"]->Fill(Pt3,Code3,weight);
	
	_DEBUG("");
	

	if(Muon1) histos[name+"_muOnly_pt"]->Fill(Pt1,weight);
	if(Muon2) histos[name+"_muOnly_pt"]->Fill(Pt2,weight);
	if(Muon3) histos[name+"_muOnly_pt"]->Fill(Pt3,weight);
	if(TPa1) histos[name+"_tpA_pt"]->Fill(Pt1,weight);
	if(TPa2) histos[name+"_tpA_pt"]->Fill(Pt2,weight);
	if(TPa3) histos[name+"_tpA_pt"]->Fill(Pt3,weight);
	if(TPb1) histos[name+"_tpB_pt"]->Fill(Pt1,weight);
	if(TPb2) histos[name+"_tpB_pt"]->Fill(Pt2,weight);
	if(TPb3) histos[name+"_tpB_pt"]->Fill(Pt3,weight);

	
	_DEBUG("");
	
	histos[name+"_mu_pt1"]->Fill(Pt1,weight);
	histos[name+"_mu_pt2"]->Fill(Pt2,weight);
	histos[name+"_mu_pt3"]->Fill(Pt3,weight);
	if(getTriggerbit("EF_3mu4T"))
	{
		histos[name+"_mu_pt1_3mu4T"]->Fill(Pt1,weight);
		histos[name+"_mu_pt2_3mu4T"]->Fill(Pt2,weight);
		histos[name+"_mu_pt3_3mu4T"]->Fill(Pt3,weight);
		if(getUniqueTriggerbit("EF_3mu4T"))
		{
			histos[name+"_mu_pt1_unique_3mu4T"]->Fill(Pt1,weight);
			histos[name+"_mu_pt2_unique_3mu4T"]->Fill(Pt2,weight);
			histos[name+"_mu_pt3_unique_3mu4T"]->Fill(Pt3,weight);
		}
	}
	if(getTriggerbit("EF_2mu8_EFxe30_tclcw"))
	{
		histos[name+"_mu_pt1_2mu8_EFxe30_tclcw"]->Fill(Pt1,weight);
		histos[name+"_mu_pt2_2mu8_EFxe30_tclcw"]->Fill(Pt2,weight);
		histos[name+"_mu_pt3_2mu8_EFxe30_tclcw"]->Fill(Pt3,weight);
		if(getUniqueTriggerbit("EF_2mu8_EFxe30_tclcw"))
		{
			histos[name+"_mu_pt1_unique_2mu8_EFxe30_tclcw"]->Fill(Pt1,weight);
			histos[name+"_mu_pt2_unique_2mu8_EFxe30_tclcw"]->Fill(Pt2,weight);
			histos[name+"_mu_pt3_unique_2mu8_EFxe30_tclcw"]->Fill(Pt3,weight);
		}
	}
	// if(getTriggerbit("EF_mu24i_tight"))
	// {
	// 	histos[name+"_mu_pt1_mu24i_tight"]->Fill(Pt1,weight);
	// 	histos[name+"_mu_pt2_mu24i_tight"]->Fill(Pt2,weight);
	// 	histos[name+"_mu_pt3_mu24i_tight"]->Fill(Pt3,weight);
	// 	if(getUniqueTriggerbit("EF_mu24i_tight"))
	// 	{
	// 		histos[name+"_mu_pt1_unique_mu24i_tight"]->Fill(Pt1,weight);
	// 		histos[name+"_mu_pt2_unique_mu24i_tight"]->Fill(Pt2,weight);
	// 		histos[name+"_mu_pt3_unique_mu24i_tight"]->Fill(Pt3,weight);
	// 	}
	// }
	// if(getTriggerbit("EF_mu36_tight"))
	// {
	// 	histos[name+"_mu_pt1_mu36_tight"]->Fill(Pt1,weight);
	// 	histos[name+"_mu_pt2_mu36_tight"]->Fill(Pt2,weight);
	// 	histos[name+"_mu_pt3_mu36_tight"]->Fill(Pt3,weight);
	// 	if(getUniqueTriggerbit("EF_mu36_tight"))
	// 	{
	// 		histos[name+"_mu_pt1_unique_mu36_tight"]->Fill(Pt1,weight);
	// 		histos[name+"_mu_pt2_unique_mu36_tight"]->Fill(Pt2,weight);
	// 		histos[name+"_mu_pt3_unique_mu36_tight"]->Fill(Pt3,weight);
	// 	}
	// }
	
	
	
	
	
	
	
	
	
	
	_DEBUG("");
	
	// MCP ID hits cut - !!! NOT SORTED !!!
	bool mcp1 = MCP(itrk1);
	bool mcp2 = MCP(itrk2);
	bool mcp3 = MCP(itrk3);
	mcp1 = (isTPb1) ?  (!skim && mcp1 && nHighThresholdTRTHits1<5) : mcp1;        
	mcp2 = (isTPb2) ?  (!skim && mcp2 && nHighThresholdTRTHits2<5) : mcp2;
	mcp3 = (isTPb3) ?  (!skim && mcp3 && nHighThresholdTRTHits3<5) : mcp3;
	mcp1 = (isMuon1 || isTPa1) ?  (!skim && mcp1 && nHighThresholdTRTHits1<10) : mcp1;        
	mcp2 = (isMuon2 || isTPa2) ?  (!skim && mcp2 && nHighThresholdTRTHits2<10) : mcp2;
	mcp3 = (isMuon3 || isTPa3) ?  (!skim && mcp3 && nHighThresholdTRTHits3<10) : mcp3;
	histos[name+"_mu_MCP"]->Fill(mcp1,weight);
	histos[name+"_mu_MCP"]->Fill(mcp2,weight);
	histos[name+"_mu_MCP"]->Fill(mcp3,weight);
	if(isTPb1) histos[name+"_tpB_nHighThresholdTRTHits"]->Fill(tpmu_vi[src1+"_nHighThresholdTRTHits"]->at(isrc1),weight);
	if(isTPb2) histos[name+"_tpB_nHighThresholdTRTHits"]->Fill(tpmu_vi[src2+"_nHighThresholdTRTHits"]->at(isrc2),weight);
	if(isTPb3) histos[name+"_tpB_nHighThresholdTRTHits"]->Fill(tpmu_vi[src3+"_nHighThresholdTRTHits"]->at(isrc3),weight);
	if(isTPa1) histos[name+"_tpA_nHighThresholdTRTHits"]->Fill(tpmu_vi[src1+"_nHighThresholdTRTHits"]->at(isrc1),weight);
	if(isTPa2) histos[name+"_tpA_nHighThresholdTRTHits"]->Fill(tpmu_vi[src2+"_nHighThresholdTRTHits"]->at(isrc2),weight);
	if(isTPa3) histos[name+"_tpA_nHighThresholdTRTHits"]->Fill(tpmu_vi[src3+"_nHighThresholdTRTHits"]->at(isrc3),weight);
	if(isMuon1) histos[name+"_mu_nHighThresholdTRTHits"]->Fill(trks_nHighThresholdTRTHits->at(itrk1),weight);
	if(isMuon2) histos[name+"_mu_nHighThresholdTRTHits"]->Fill(trks_nHighThresholdTRTHits->at(itrk2),weight);
	if(isMuon3) histos[name+"_mu_nHighThresholdTRTHits"]->Fill(trks_nHighThresholdTRTHits->at(itrk3),weight);
	if(isCounter("nPassing_mu_mcp") && (mcp1+mcp2+mcp3)!=3) return false;
	incrementCounter("nPassing_mu_mcp",weight);
	
	_DEBUG("");
	
	// trk pT cut !!! SORTED !!!
	double pTmuCut1 = (skim) ? 2.0*GeV2MeV : 5.5*GeV2MeV;
	double pTmuCut2 = (skim) ? 2.0*GeV2MeV : 3.5*GeV2MeV;
	double pTmuCut3 = (skim) ? 2.0*GeV2MeV : 2.5*GeV2MeV;
	bool pt1 = (fabs(Pt1)>pTmuCut1);
	bool pt2 = (fabs(Pt2)>pTmuCut2);
	bool pt3 = (fabs(Pt3)>pTmuCut3);
	pt1 = (TPb1) ? (!skim && fabs(Pt1)>3.0*GeV2MeV) : pt1;
	pt2 = (TPb2) ? (!skim && fabs(Pt2)>3.0*GeV2MeV) : pt2;
	pt3 = (TPb3) ? (!skim && fabs(Pt3)>3.0*GeV2MeV) : pt3;
	if(isCounter("nPassing_mu_pt") && (pt1+pt2+pt3)!=3) return false;
	incrementCounter("nPassing_mu_pt",weight);
	
	_DEBUG("");
	



	// trk eta cut  !!! SORTED !!!
	histos[name+"_mu_eta"]->Fill(Eta1,weight);
	histos[name+"_mu_eta"]->Fill(Eta2,weight);
	histos[name+"_mu_eta"]->Fill(Eta3,weight);
	if(Muon1) histos[name+"_muOnly_eta"]->Fill(Eta1,weight);
	if(Muon2) histos[name+"_muOnly_eta"]->Fill(Eta2,weight);
	if(Muon3) histos[name+"_muOnly_eta"]->Fill(Eta3,weight);
	if(TPa1) histos[name+"_tpA_eta"]->Fill(Eta1,weight);
	if(TPa2) histos[name+"_tpA_eta"]->Fill(Eta2,weight);
	if(TPa3) histos[name+"_tpA_eta"]->Fill(Eta3,weight);
	if(TPb1) histos[name+"_tpB_eta"]->Fill(Eta1,weight);
	if(TPb2) histos[name+"_tpB_eta"]->Fill(Eta2,weight);
	if(TPb3) histos[name+"_tpB_eta"]->Fill(Eta3,weight);
	bool eta1 = (fabs(Eta1)<2.5);
	bool eta2 = (fabs(Eta2)<2.5);
	bool eta3 = (fabs(Eta3)<2.5);
	eta1 = (TPa1) ? (!skim && eta1 && !(fabs(Eta1)>1.05 && fabs(Eta1)<1.6)) : eta1;
	eta2 = (TPa2) ? (!skim && eta2 && !(fabs(Eta2)>1.05 && fabs(Eta2)<1.6)) : eta2;
	eta3 = (TPa3) ? (!skim && eta3 && !(fabs(Eta3)>1.05 && fabs(Eta3)<1.6)) : eta3;
	if(isCounter("nPassing_mu_eta") && (eta1+eta2+eta3)!=3) return false;
	incrementCounter("nPassing_mu_eta",weight);


	
	_DEBUG("");



	// !!! NOT SORTED !!!
	if(isMuon1) { histos[name+"_mu_sctangsig"]->Fill(sctangsig1,weight); histos[name+"_mu_sctngbsig"]->Fill(sctngbsig1,weight); histos[name+"_mu_chi2"]->Fill(fitchi2mu1,weight); histos[name+"_mu_chi2_zoom"]->Fill(fitchi2mu1,weight); histos[name+"_mu_chi2ndf"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_matchchi2ndf"]->Fill(matchchi2ndfmu1,weight); }
	if(isMuon2) { histos[name+"_mu_sctangsig"]->Fill(sctangsig2,weight); histos[name+"_mu_sctngbsig"]->Fill(sctngbsig2,weight); histos[name+"_mu_chi2"]->Fill(fitchi2mu2,weight); histos[name+"_mu_chi2_zoom"]->Fill(fitchi2mu2,weight); histos[name+"_mu_chi2ndf"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_matchchi2ndf"]->Fill(matchchi2ndfmu2,weight); }
	if(isMuon3) { histos[name+"_mu_sctangsig"]->Fill(sctangsig3,weight); histos[name+"_mu_sctngbsig"]->Fill(sctngbsig3,weight); histos[name+"_mu_chi2"]->Fill(fitchi2mu3,weight); histos[name+"_mu_chi2_zoom"]->Fill(fitchi2mu3,weight); histos[name+"_mu_chi2ndf"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_matchchi2ndf"]->Fill(matchchi2ndfmu3,weight); }
	if(isMuon1) { histos2[name+"_mu_chi2_vs_sctangsig"]->Fill(fitchi2mu1,sctangsig1,weight); histos2[name+"_mu_pvalue_vs_sctangsig"]->Fill(fitpvaluemu1,sctangsig1,weight); }
	if(isMuon2) { histos2[name+"_mu_chi2_vs_sctangsig"]->Fill(fitchi2mu2,sctangsig2,weight); histos2[name+"_mu_pvalue_vs_sctangsig"]->Fill(fitpvaluemu2,sctangsig2,weight); }
	if(isMuon3) { histos2[name+"_mu_chi2_vs_sctangsig"]->Fill(fitchi2mu3,sctangsig3,weight); histos2[name+"_mu_pvalue_vs_sctangsig"]->Fill(fitpvaluemu3,sctangsig3,weight); }
	// Scattering angle significance cut !!! NOT SORTED !!!
	bool psig1 = (isMuon1) ? (fabs(sctangsig1)<5.) : true;
	bool psig2 = (isMuon2) ? (fabs(sctangsig2)<5.) : true;
	bool psig3 = (isMuon3) ? (fabs(sctangsig3)<5.) : true;
	if(isCounter("nPassing_mu_sctangsig") && (psig1+psig2+psig3)!=3) return false;
	incrementCounter("nPassing_mu_sctangsig",weight);
	
	_DEBUG("");
	
	// !!! NOT SORTED !!!
	if(isMuon1) histos[name+"_mu_pbalsig"]->Fill(pbalsig1,weight);	
	if(isMuon2) histos[name+"_mu_pbalsig"]->Fill(pbalsig2,weight);	
	if(isMuon3) histos[name+"_mu_pbalsig"]->Fill(pbalsig3,weight);
	if(isTPmu1) histos[name+"_tpmu_pbalsig"]->Fill(pbalsigTP1,weight);	
	if(isTPmu2) histos[name+"_tpmu_pbalsig"]->Fill(pbalsigTP2,weight);	
	if(isTPmu3) histos[name+"_tpmu_pbalsig"]->Fill(pbalsigTP3,weight);
	// momentum balance significance cut !!! NOT SORTED !!!
	bool pbal1 = (isMuon1 && isCB1) ? (fabs(pbalsig1)<3.) : true;
	bool pbal2 = (isMuon2 && isCB2) ? (fabs(pbalsig2)<3.) : true;
	bool pbal3 = (isMuon3 && isCB3) ? (fabs(pbalsig3)<3.) : true;
	if(isCounter("nPassing_mu_pbalsig") && (pbal1+pbal2+pbal3)!=3) return false;
	incrementCounter("nPassing_mu_pbalsig",weight);
	
	_DEBUG("");
	
	// !!! NOT SORTED !!!
	if(isMuon1) { histos[name+"_mu_qoverp_me"]->Fill(qopmemu1,weight); histos[name+"_mu_qoverp_ie"]->Fill(qopiemu1,weight); }
	if(isMuon2) { histos[name+"_mu_qoverp_me"]->Fill(qopmemu2,weight); histos[name+"_mu_qoverp_ie"]->Fill(qopiemu2,weight); }
	if(isMuon3) { histos[name+"_mu_qoverp_me"]->Fill(qopmemu3,weight); histos[name+"_mu_qoverp_ie"]->Fill(qopiemu3,weight); }
	
	_DEBUG("");
	
	// !!! NOT SORTED !!!
	if(isTPmu1) { histos[name+"_tpmu_ndf"]->Fill(fitndfTP1,weight);	histos2[name+"_tpmu_pvalue_vs_chi2ndf"]->Fill(fitchi2TP1/fitndfTP1,fitpvalueTP1,weight); }
	if(isTPmu2) { histos[name+"_tpmu_ndf"]->Fill(fitndfTP2,weight);	histos2[name+"_tpmu_pvalue_vs_chi2ndf"]->Fill(fitchi2TP2/fitndfTP2,fitpvalueTP2,weight); }
	if(isTPmu3) { histos[name+"_tpmu_ndf"]->Fill(fitndfTP3,weight); histos2[name+"_tpmu_pvalue_vs_chi2ndf"]->Fill(fitchi2TP3/fitndfTP3,fitpvalueTP3,weight); }
	if(isTPmu1) { histos[name+"_tpmu_chi2"]->Fill(fitchi2TP1,weight); histos[name+"_tpmu_chi2_zoom"]->Fill(fitchi2TP1,weight); histos[name+"_tpmu_chi2ndf"]->Fill(fitchi2TP1/fitndfTP1,weight); }
	if(isTPmu2) { histos[name+"_tpmu_chi2"]->Fill(fitchi2TP2,weight); histos[name+"_tpmu_chi2_zoom"]->Fill(fitchi2TP2,weight); histos[name+"_tpmu_chi2ndf"]->Fill(fitchi2TP2/fitndfTP2,weight); }  
	if(isTPmu3) { histos[name+"_tpmu_chi2"]->Fill(fitchi2TP3,weight); histos[name+"_tpmu_chi2_zoom"]->Fill(fitchi2TP3,weight); histos[name+"_tpmu_chi2ndf"]->Fill(fitchi2TP3/fitndfTP3,weight); }
	if(isTPmu1) { histos[name+"_tpmu_pvalue"]->Fill(fitpvalueTP1,weight); histos[name+"_tpmu_pvalue_zoom"]->Fill(fitpvalueTP1,weight); }	
	if(isTPmu2) { histos[name+"_tpmu_pvalue"]->Fill(fitpvalueTP2,weight); histos[name+"_tpmu_pvalue_zoom"]->Fill(fitpvalueTP2,weight); }	
	if(isTPmu3) { histos[name+"_tpmu_pvalue"]->Fill(fitpvalueTP3,weight); histos[name+"_tpmu_pvalue_zoom"]->Fill(fitpvalueTP3,weight); }
	if(isTPa1) { histos2[name+"_tpA_pvalue_vs_chi2ndf"]->Fill(fitchi2TP1/fitndfTP1,fitpvalueTP1,weight); histos[name+"_tpA_qoverp"]->Fill(qopTP1,weight); histos[name+"_tpA_chi2"]->Fill(fitchi2TP1,weight); histos[name+"_tpA_chi2_zoom"]->Fill(fitchi2TP1,weight); histos[name+"_tpA_chi2ndf"]->Fill(fitchi2TP1/fitndfTP1,weight); }
	if(isTPa2) { histos2[name+"_tpA_pvalue_vs_chi2ndf"]->Fill(fitchi2TP2/fitndfTP2,fitpvalueTP2,weight); histos[name+"_tpA_qoverp"]->Fill(qopTP2,weight); histos[name+"_tpA_chi2"]->Fill(fitchi2TP2,weight); histos[name+"_tpA_chi2_zoom"]->Fill(fitchi2TP2,weight); histos[name+"_tpA_chi2ndf"]->Fill(fitchi2TP2/fitndfTP2,weight); }  
	if(isTPa3) { histos2[name+"_tpA_pvalue_vs_chi2ndf"]->Fill(fitchi2TP3/fitndfTP3,fitpvalueTP3,weight); histos[name+"_tpA_qoverp"]->Fill(qopTP3,weight); histos[name+"_tpA_chi2"]->Fill(fitchi2TP3,weight); histos[name+"_tpA_chi2_zoom"]->Fill(fitchi2TP3,weight); histos[name+"_tpA_chi2ndf"]->Fill(fitchi2TP3/fitndfTP3,weight); }
	if(isTPb1) { histos2[name+"_tpB_pvalue_vs_chi2ndf"]->Fill(fitchi2TP1/fitndfTP1,fitpvalueTP1,weight); histos[name+"_tpB_qoverp"]->Fill(qopTP1,weight); histos[name+"_tpB_chi2"]->Fill(fitchi2TP1,weight); histos[name+"_tpB_chi2_zoom"]->Fill(fitchi2TP1,weight); histos[name+"_tpB_chi2ndf"]->Fill(fitchi2TP1/fitndfTP1,weight); }
	if(isTPb2) { histos2[name+"_tpB_pvalue_vs_chi2ndf"]->Fill(fitchi2TP2/fitndfTP2,fitpvalueTP2,weight); histos[name+"_tpB_qoverp"]->Fill(qopTP2,weight); histos[name+"_tpB_chi2"]->Fill(fitchi2TP2,weight); histos[name+"_tpB_chi2_zoom"]->Fill(fitchi2TP2,weight); histos[name+"_tpB_chi2ndf"]->Fill(fitchi2TP2/fitndfTP2,weight); }  
	if(isTPb3) { histos2[name+"_tpB_pvalue_vs_chi2ndf"]->Fill(fitchi2TP3/fitndfTP3,fitpvalueTP3,weight); histos[name+"_tpB_qoverp"]->Fill(qopTP3,weight); histos[name+"_tpB_chi2"]->Fill(fitchi2TP3,weight); histos[name+"_tpB_chi2_zoom"]->Fill(fitchi2TP3,weight); histos[name+"_tpB_chi2ndf"]->Fill(fitchi2TP3/fitndfTP3,weight); }
	if(isTPa1) { histos[name+"_tpA_pvalue"]->Fill(fitpvalueTP1,weight); histos[name+"_tpA_pvalue_zoom"]->Fill(fitpvalueTP1,weight); }	
	if(isTPa2) { histos[name+"_tpA_pvalue"]->Fill(fitpvalueTP2,weight); histos[name+"_tpA_pvalue_zoom"]->Fill(fitpvalueTP2,weight); }	
	if(isTPa3) { histos[name+"_tpA_pvalue"]->Fill(fitpvalueTP3,weight); histos[name+"_tpA_pvalue_zoom"]->Fill(fitpvalueTP3,weight); }
	if(isTPb1) { histos[name+"_tpB_pvalue"]->Fill(fitpvalueTP1,weight); histos[name+"_tpB_pvalue_zoom"]->Fill(fitpvalueTP1,weight); }	
	if(isTPb2) { histos[name+"_tpB_pvalue"]->Fill(fitpvalueTP2,weight); histos[name+"_tpB_pvalue_zoom"]->Fill(fitpvalueTP2,weight); }	
	if(isTPb3) { histos[name+"_tpB_pvalue"]->Fill(fitpvalueTP3,weight); histos[name+"_tpB_pvalue_zoom"]->Fill(fitpvalueTP3,weight); }
	
	_DEBUG("");


	// !!! NOT SORTED !!!
	if(fabs(charge)==1. && foundOS1 && mOS1>680. && mOS1<1100.)
	{	
		histos[name+"_mOS_rhoomegaphi"]->Fill(mOS1,weight);
		histos2[name+"_mSS_vs_mOS_rhoomegaphi"]->Fill(mOS1,mSS,weight);
		histos2[name+"_m3mu_vs_mOS_rhoomegaphi"]->Fill(mOS1,mass,weight);
		histos2[name+"_pT3mu_vs_mOS_rhoomegaphi"]->Fill(mOS1,pTsum,weight);
		
		histos[name+"_mSS_rhoomegaphi"]->Fill(mSS,weight);
		histos[name+"_m3mu_rhoomegaphi"]->Fill(mass,weight);
		histos[name+"_pT3mu_rhoomegaphi"]->Fill(pTsum,weight);
		histos[name+"_pvalue_rhoomegaphi"]->Fill(pvalue,weight);
		
		if(isMuon1 && !isMedium1) { histos[name+"_mu_chi2ndf_failMedium_rhoomegaphi"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue_failMedium_rhoomegaphi"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom_failMedium_rhoomegaphi"]->Fill(fitpvaluemu1,weight); }
		if(isMuon2 && !isMedium2) { histos[name+"_mu_chi2ndf_failMedium_rhoomegaphi"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue_failMedium_rhoomegaphi"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom_failMedium_rhoomegaphi"]->Fill(fitpvaluemu2,weight); }
		if(isMuon3 && !isMedium3) { histos[name+"_mu_chi2ndf_failMedium_rhoomegaphi"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue_failMedium_rhoomegaphi"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom_failMedium_rhoomegaphi"]->Fill(fitpvaluemu3,weight); }
		
		if(isMuon1 && isMedium1) { histos[name+"_mu_chi2ndf_passMedium_rhoomegaphi"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue_passMedium_rhoomegaphi"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom_passMedium_rhoomegaphi"]->Fill(fitpvaluemu1,weight); }
		if(isMuon2 && isMedium2) { histos[name+"_mu_chi2ndf_passMedium_rhoomegaphi"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue_passMedium_rhoomegaphi"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom_passMedium_rhoomegaphi"]->Fill(fitpvaluemu2,weight); }
		if(isMuon3 && isMedium3) { histos[name+"_mu_chi2ndf_passMedium_rhoomegaphi"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue_passMedium_rhoomegaphi"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom_passMedium_rhoomegaphi"]->Fill(fitpvaluemu3,weight); }
		
		if(isTPmu1) { histos[name+"_tpmu_chi2ndf_rhoomegaphi"]->Fill(fitchi2TP1/fitndfTP1,weight); histos[name+"_tpmu_pvalue_rhoomegaphi"]->Fill(fitpvalueTP1,weight); histos[name+"_tpmu_pvalue_zoom_rhoomegaphi"]->Fill(fitpvalueTP1,weight); }
		if(isTPmu2) { histos[name+"_tpmu_chi2ndf_rhoomegaphi"]->Fill(fitchi2TP2/fitndfTP2,weight); histos[name+"_tpmu_pvalue_rhoomegaphi"]->Fill(fitpvalueTP2,weight); histos[name+"_tpmu_pvalue_zoom_rhoomegaphi"]->Fill(fitpvalueTP2,weight); }
		if(isTPmu3) { histos[name+"_tpmu_chi2ndf_rhoomegaphi"]->Fill(fitchi2TP3/fitndfTP3,weight); histos[name+"_tpmu_pvalue_rhoomegaphi"]->Fill(fitpvalueTP3,weight); histos[name+"_tpmu_pvalue_zoom_rhoomegaphi"]->Fill(fitpvalueTP3,weight); }	   
	}
	if(fabs(charge)==1. && foundOS1 && mOS1>2900. && mOS1<3300.)
	{	
		histos[name+"_mOS_Jpsi"]->Fill(mOS1,weight);
		histos2[name+"_mSS_vs_mOS_Jpsi"]->Fill(mOS1,mSS,weight);
		histos2[name+"_m3mu_vs_mOS_Jpsi"]->Fill(mOS1,mass,weight);
		histos2[name+"_pT3mu_vs_mOS_Jpsi"]->Fill(mOS1,pTsum,weight);
		
		histos[name+"_mSS_Jpsi"]->Fill(mSS,weight);
		histos[name+"_m3mu_Jpsi"]->Fill(mass,weight);
		histos[name+"_pT3mu_Jpsi"]->Fill(pTsum,weight);
		histos[name+"_pvalue_Jpsi"]->Fill(pvalue,weight);
		
		if(isMuon1 && !isMedium1) { histos[name+"_mu_chi2ndf_failMedium_Jpsi"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue_failMedium_Jpsi"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom_failMedium_Jpsi"]->Fill(fitpvaluemu1,weight); }
		if(isMuon2 && !isMedium2) { histos[name+"_mu_chi2ndf_failMedium_Jpsi"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue_failMedium_Jpsi"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom_failMedium_Jpsi"]->Fill(fitpvaluemu2,weight); }
		if(isMuon3 && !isMedium3) { histos[name+"_mu_chi2ndf_failMedium_Jpsi"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue_failMedium_Jpsi"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom_failMedium_Jpsi"]->Fill(fitpvaluemu3,weight); }
		
		if(isMuon1 && isMedium1) { histos[name+"_mu_chi2ndf_passMedium_Jpsi"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue_passMedium_Jpsi"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom_passMedium_Jpsi"]->Fill(fitpvaluemu1,weight); }
		if(isMuon2 && isMedium2) { histos[name+"_mu_chi2ndf_passMedium_Jpsi"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue_passMedium_Jpsi"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom_passMedium_Jpsi"]->Fill(fitpvaluemu2,weight); }
		if(isMuon3 && isMedium3) { histos[name+"_mu_chi2ndf_passMedium_Jpsi"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue_passMedium_Jpsi"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom_passMedium_Jpsi"]->Fill(fitpvaluemu3,weight); }
		
		if(isTPmu1) { histos[name+"_tpmu_chi2ndf_Jpsi"]->Fill(fitchi2TP1/fitndfTP1,weight); histos[name+"_tpmu_pvalue_Jpsi"]->Fill(fitpvalueTP1,weight); histos[name+"_tpmu_pvalue_zoom_Jpsi"]->Fill(fitpvalueTP1,weight); }
		if(isTPmu2) { histos[name+"_tpmu_chi2ndf_Jpsi"]->Fill(fitchi2TP2/fitndfTP2,weight); histos[name+"_tpmu_pvalue_Jpsi"]->Fill(fitpvalueTP2,weight); histos[name+"_tpmu_pvalue_zoom_Jpsi"]->Fill(fitpvalueTP2,weight); }
		if(isTPmu3) { histos[name+"_tpmu_chi2ndf_Jpsi"]->Fill(fitchi2TP3/fitndfTP3,weight); histos[name+"_tpmu_pvalue_Jpsi"]->Fill(fitpvalueTP3,weight); histos[name+"_tpmu_pvalue_zoom_Jpsi"]->Fill(fitpvalueTP3,weight); }
	}
	
	_DEBUG("");
	
	// !!! NOT SORTED !!!
	if(isMuon1) histos2[name+"_mu_pvalue_vs_chi2ndf"]->Fill(fitchi2mu1/fitndfmu1,fitpvaluemu1,weight);
	if(isMuon2) histos2[name+"_mu_pvalue_vs_chi2ndf"]->Fill(fitchi2mu2/fitndfmu2,fitpvaluemu2,weight);
	if(isMuon3) histos2[name+"_mu_pvalue_vs_chi2ndf"]->Fill(fitchi2mu3/fitndfmu3,fitpvaluemu3,weight);
	
	// !!! NOT SORTED !!!
	if(isMuon1 && !isMedium1) { histos[name+"_mu_chi2ndf_failMedium"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue_failMedium"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom_failMedium"]->Fill(fitpvaluemu1,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_failMedium"]->Fill(fitchi2mu1/fitndfmu1,fitpvaluemu1,weight); }
	if(isMuon2 && !isMedium2) { histos[name+"_mu_chi2ndf_failMedium"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue_failMedium"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom_failMedium"]->Fill(fitpvaluemu2,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_failMedium"]->Fill(fitchi2mu2/fitndfmu2,fitpvaluemu2,weight); }
	if(isMuon3 && !isMedium3) { histos[name+"_mu_chi2ndf_failMedium"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue_failMedium"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom_failMedium"]->Fill(fitpvaluemu3,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_failMedium"]->Fill(fitchi2mu3/fitndfmu3,fitpvaluemu3,weight); }
	
	if(isMuon1 && isMedium1) { histos[name+"_mu_chi2ndf_passMedium"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue_passMedium"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom_passMedium"]->Fill(fitpvaluemu1,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_passMedium"]->Fill(fitchi2mu1/fitndfmu1,fitpvaluemu1,weight); }
	if(isMuon2 && isMedium2) { histos[name+"_mu_chi2ndf_passMedium"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue_passMedium"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom_passMedium"]->Fill(fitpvaluemu2,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_passMedium"]->Fill(fitchi2mu2/fitndfmu2,fitpvaluemu2,weight); }
	if(isMuon3 && isMedium3) { histos[name+"_mu_chi2ndf_passMedium"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue_passMedium"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom_passMedium"]->Fill(fitpvaluemu3,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_passMedium"]->Fill(fitchi2mu3/fitndfmu3,fitpvaluemu3,weight); }
	
	_DEBUG("");
	
	// !!! NOT SORTED !!!
	if(isMuon1 && isMedium1) { histos[name+"_mu_chi2_ismedium"]->Fill(fitchi2mu1,weight); histos[name+"_mu_chi2_zoom_ismedium"]->Fill(fitchi2mu1,weight); histos[name+"_mu_chi2ndf_ismedium"]->Fill(fitchi2mu1/fitndfmu1,weight); histos[name+"_mu_pvalue_ismedium"]->Fill(fitpvaluemu1,weight); histos[name+"_mu_pvalue_zoom_ismedium"]->Fill(fitpvaluemu1,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_ismedium"]->Fill(fitchi2mu1/fitndfmu1,fitpvaluemu1,weight); }
	if(isMuon2 && isMedium2) { histos[name+"_mu_chi2_ismedium"]->Fill(fitchi2mu2,weight); histos[name+"_mu_chi2_zoom_ismedium"]->Fill(fitchi2mu2,weight); histos[name+"_mu_chi2ndf_ismedium"]->Fill(fitchi2mu2/fitndfmu2,weight); histos[name+"_mu_pvalue_ismedium"]->Fill(fitpvaluemu2,weight); histos[name+"_mu_pvalue_zoom_ismedium"]->Fill(fitpvaluemu2,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_ismedium"]->Fill(fitchi2mu2/fitndfmu2,fitpvaluemu2,weight); }
	if(isMuon3 && isMedium3) { histos[name+"_mu_chi2_ismedium"]->Fill(fitchi2mu3,weight); histos[name+"_mu_chi2_zoom_ismedium"]->Fill(fitchi2mu3,weight); histos[name+"_mu_chi2ndf_ismedium"]->Fill(fitchi2mu3/fitndfmu3,weight); histos[name+"_mu_pvalue_ismedium"]->Fill(fitpvaluemu3,weight); histos[name+"_mu_pvalue_zoom_ismedium"]->Fill(fitpvaluemu3,weight); histos2[name+"_mu_pvalue_vs_chi2ndf_ismedium"]->Fill(fitchi2mu3/fitndfmu3,fitpvaluemu3,weight); }
	
	_DEBUG("");
	
	// the muon/tp track-fit quality cut !!! NOT SORTED !!!
	bool trkquality1 = (isTPmu1) ? ((isTPa1 && fitpvalueTP1>0.01) || (isTPb1 && fitpvalueTP1>0.1)) : isMedium1;
	bool trkquality2 = (isTPmu2) ? ((isTPa2 && fitpvalueTP2>0.01) || (isTPb2 && fitpvalueTP2>0.1)) : isMedium2;
	bool trkquality3 = (isTPmu3) ? ((isTPa3 && fitpvalueTP3>0.01) || (isTPb3 && fitpvalueTP3>0.1)) : isMedium3;
	if(isCounter("nPassing_mu_trkquality") && (trkquality1+trkquality2+trkquality3)!=3) return false;
	incrementCounter("nPassing_mu_trkquality",weight);
	
	_DEBUG("");

	///////////////////////////////////////////////
	// fill some more histos after all muon cuts //
	///////////////////////////////////////////////
	fillHistsMassPt3mu(vtx,name,"_after_muons",histos,histos2,weight,mBlindMin,mBlindMax);
	fillHistsDoublets(vtx,name,"_after_muons",histos,histos2,weight,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS);
	fillCategories(vtx,name,"tripletCategories_after_muons",histos);
	fillCategories(vtx,name,"tripletCategories_norm_after_muons",histos);

	return true;
}

bool acceptVtxMET(TString method, unsigned int vtx, vector<vertex>& vertices, TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight=1., double mBlindMin=1500., double mBlindMax=2000.)
{
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];	
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
	TLorentzVector pMuSum = getTlv3mu(vtx);	
	vector<int> pv1; findPVindex(pv1,1);
	if(pv1.size()<1) _ERROR("PV not found OR too many (>1) PVs were found: nPV="+_s((int)pv1.size()));
	int pvindex = pv1[0];
	unsigned int nvtxpv = vtx_lxy->size();
	double m3mu  = vtx_mass->at(vtx); // pMuSum.M();
	double pT3mu = vtx_pt->at(vtx);   // pMuSum.Pt();
	
	_DEBUG("");
	
	// OS and SS pairs
	TLorentzVector pOS1, pOS2, pSS;
	bool foundOS1 = false;
	bool foundOS2 = false;
	bool foundSS  = false;
	int iUniqueCharge = -1;
	TMapVL doubletsTLV = getDoubletsFromTriplet(vtx,iUniqueCharge);
	if(doubletsTLV.find("OS1") != doubletsTLV.end()) { pOS1 = doubletsTLV["OS1"]; foundOS1 = true; }
	if(doubletsTLV.find("OS2") != doubletsTLV.end()) { pOS2 = doubletsTLV["OS2"]; foundOS2 = true; }
	if(doubletsTLV.find("SS")  != doubletsTLV.end()) { pSS  = doubletsTLV["SS"];  foundSS  = true; }
	if(foundOS1+foundOS2<2 && fabs(vtx_charge->at(vtx))==1.) _ERROR("could not find OS doublet"); // should be exactly two OS pairs
	if(!foundSS && fabs(vtx_charge->at(vtx))==1.)            _ERROR("could not find SS doublet"); // should be exactly one SS pair
	double mOS1 = (foundOS1) ? pOS1.M() : -1.;
	double mOS2 = (foundOS2) ? pOS2.M() : -1.;
	double mSS  = (foundSS)  ? pSS.M()  : -1.;
	double weightOS = (foundOS1+foundOS2==2) ? weight/2. : weight;
	
	_DEBUG("");
	
	double dRmin = +1.e20;
	dRmin = (pMuSum.DeltaR(p1)<dRmin) ? pMuSum.DeltaR(p1) : dRmin;
	dRmin = (pMuSum.DeltaR(p2)<dRmin) ? pMuSum.DeltaR(p2) : dRmin;
	dRmin = (pMuSum.DeltaR(p3)<dRmin) ? pMuSum.DeltaR(p3) : dRmin;
	double dRmax = -1.e20;
	dRmax = (pMuSum.DeltaR(p1)>dRmax) ? pMuSum.DeltaR(p1) : dRmax;
	dRmax = (pMuSum.DeltaR(p2)>dRmax) ? pMuSum.DeltaR(p2) : dRmax;
	dRmax = (pMuSum.DeltaR(p3)>dRmax) ? pMuSum.DeltaR(p3) : dRmax;
	
	_DEBUG("");
	
	// for the MVA
	float vertex_pt        = vtx_pt->at(vtx);
	float vertex_mass      = vtx_mass->at(vtx);
	float vertex_mOS1      = mOS1;
	float vertex_mOS2      = mOS2;
	float vertex_mSS       = mSS;
	// float jets_pt1         = JetPt1;
	// float jets_pt2         = JetPt2;
	// float jets_dphi3muJ1   = JetdPhi3bodyJ1;
	// float jets_dphiJ1J2    = JetdPhiJ1J2;
	// float vtx_ptfrac12     = ptFraction12;
	// float vtx_ptfrac13     = ptFraction13;
	// float vtx_ptfrac23     = ptFraction23;
	// float vertex_isolation = isolation;
	// float vertex_rapidity  = vtx_rapidity->at(vtx);
	// float vertex_charge    = vtx_charge->at(vtx);
	float vertex_pval      = TMath::Prob(vtx_chi2->at(vtx),vtx_ndf->at(vtx));
	float vertex_lxy       = vtx_lxy->at(vtx)[pvindex];
	// float vertex_lxyErr = vtx_lxyErr->at(vtx)[pvindex];
	// float vertex_lxySig = (vertex_lxyErr!=0.) ? vertex_lxy/vertex_lxyErr : 0.;
	// float vertex_a0     = vtx_a0->at(vtx)[pvindex];
	float vertex_a0xy      = vtx_a0XY->at(vtx)[pvindex];
	// float vertex_cosT   = fabs(vtx_cosTheta->at(vtx)[pvindex]);
	float vertex_cosTxy    = fabs(vtx_cosThetaXY->at(vtx)[pvindex]);
	
	// float met_reffinal_sumet   = MET_RefFinal_sumet;
	// float met_reffinal_et      = MET_RefFinal_et;
	// float met_reffinal_phi     = MET_RefFinal_phi;
	float met_reffinal_mT      = mT(MET_RefFinal_et,MET_RefFinal_phi,pMuSum.Pt(),pMuSum.Phi());
	float met_reffinal_dPhi3mu = fabs(dPhi(MET_RefFinal_phi,pMuSum.Phi()));
	
	_DEBUG("");
	
	/////////////////////////////////////////////////////////////
	//// fill some histos before cutting or applying the MVA ////
	/////////////////////////////////////////////////////////////
	histos[name+"_met_sumet"]->Fill(MET_RefFinal_sumet,weight);
	histos[name+"_met_et"]->Fill(MET_RefFinal_et,weight);
	histos[name+"_met_mt_et3mu"]->Fill(met_reffinal_mT,weight);
	histos[name+"_met_dphi3mu"]->Fill(met_reffinal_dPhi3mu,weight);
	histos[name+"_pvalue"]->Fill(vertex_pval,weight);
	histos[name+"_chi2"]->Fill(vtx_chi2->at(vtx),weight);
	histos[name+"_vtx_npv"]->Fill(nvtxpv,weight);
	
	_DEBUG("");
	
	fillHistsMVAvars(vtx,name,"_before_vtxmet",histos,histos2,pMuSum,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS,weight,mBlindMin,mBlindMax);
	
	_DEBUG("");
	
	////////////////////////////
	//// now apply the cuts ////
	////////////////////////////
	if(!doMVA)
	{
		_DEBUG("");

		// multimap<double,int> pt2i;
		// vector<int> ijet;
		// for(int i=0 ; i<(int)jets_pt->size() ; i++)
		// {
		// 	bool isUgly = (jets_isUgly->at(i));
		// 	// cout << "Jet["<<i<<"] isUgly(pt="<< jets_pt->at(i) <<")=" << jets_isUgly->at(i) << ", passCleaning("<< AntiKt4LCTopoLets_pt->at(i) << ")=" << isGoodJet(i,"TightBad") << "|" << isGoodJet(i,"MediumBad") << "|" << isGoodJet(i,JetQuality) << "|" << isGoodJet(i,"VeryLooseBad") << endl;
		// 	if(isUgly)  continue;
		// 	pt2i.insert(make_pair(jets_pt->at(i),i));
		// }
		// for(multimap<double,int>::reverse_iterator rit=pt2i.rbegin() ; rit!=pt2i.rend() ; ++rit) ijet.push_back(rit->second);
		// unsigned int njet = ijet.size();
		
		multimap<double,int> pt2i;
		vector<int> ijet;
		for(int i=0 ; i<AntiKt4LCTopoLets_n ; i++)
		{
			if(!isGoodJet(i,JetQuality))  continue;
			pt2i.insert(make_pair(AntiKt4LCTopoLets_pt->at(i),i));
		}
		for(multimap<double,int>::reverse_iterator rit=pt2i.rbegin() ; rit!=pt2i.rend() ; ++rit) ijet.push_back(rit->second);
		unsigned int njet = ijet.size();
		

		_DEBUG("");
		
		// Fill some histos before doing full had-cleaning		
		if(njet>0 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV)                                                  { histos2[name+"_dPhi3muJet1_vs_pTjet1_before_tightveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0]),fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); histos[name+"_jets_dPhi3bodyJ1_before_tightveto"]->Fill(fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); }
		if(njet>1 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV && AntiKt4LCTopoLets_pt->at(ijet[1])>20.*GeV2MeV) { histos2[name+"_dPhiJet1Jet2_vs_sumpTjet12_before_tightveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0])+AntiKt4LCTopoLets_pt->at(ijet[1]),fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); histos[name+"_jets_dPhiJ1J2_before_tightveto"]->Fill(fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); }
		
		_DEBUG("");
	

		if(isWsignal(name))
		{
			(*ofstr1) << EventNumber << " " << lbn << endl;
			if(njet>0) (*ofstr1) << "  pt1=" << AntiKt4LCTopoLets_pt->at(ijet[0]) << ", dphi3muj1=" << fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))) << endl;
			if(njet>1) (*ofstr1) << "  pt1=" << AntiKt4LCTopoLets_pt->at(ijet[0]) << ", pt2=" << AntiKt4LCTopoLets_pt->at(ijet[1]) << ", dphij1j2=" << fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))) << endl;
			(*ofstr1) << "-------------------------------------------" << endl;
		}
	
		// colinear jet veto
		/////////////////////////////////////////////////////////////////////////////////////
		bool isJet3muOverlap = (njet>0 && AntiKt4LCTopoLets_pt->at(ijet[0])>25.*GeV2MeV) ? (fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0])))<0.2 || fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0])))>2.5) : false;
		if(isCounter("nPassing_jets_coljetveto") && isJet3muOverlap) return false; /////
		incrementCounter("nPassing_jets_coljetveto",weight); /////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////


		if(isWsignal(name)) (*ofstr) << EventNumber << " " << lbn << endl;

		_DEBUG("");

		// di-jet veto
		//////////////////////////////////////////////////////////////////////////
		bool isDijets = (njet>1 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV && AntiKt4LCTopoLets_pt->at(ijet[1])>20.*GeV2MeV) ? (AntiKt4LCTopoLets_pt->at(ijet[0])+AntiKt4LCTopoLets_pt->at(ijet[1])>60.*GeV2MeV && fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1])))>2.5) : false;
		if(isCounter("nPassing_jets_dijetveto") && isDijets) return false; //////
		incrementCounter("nPassing_jets_dijetveto",weight); //////////////////////
		//////////////////////////////////////////////////////////////////////////
		
		_DEBUG("");
		
		// Fill some histos before doing full had-cleaning		
		if(njet>0 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV)                                                  { histos2[name+"_dPhi3muJet1_vs_pTjet1_after_tightveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0]),fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); histos[name+"_jets_dPhi3bodyJ1_after_tightveto"]->Fill(fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); }
		if(njet>1 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV && AntiKt4LCTopoLets_pt->at(ijet[1])>20.*GeV2MeV) { histos2[name+"_dPhiJet1Jet2_vs_sumpTjet12_after_tightveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0])+AntiKt4LCTopoLets_pt->at(ijet[1]),fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); histos[name+"_jets_dPhiJ1J2_after_tightveto"]->Fill(fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); }
		
		_DEBUG("");
		
		// 3body vertex pT
		if(isCounter("nPassing_triplet_pt") && pT3mu<25.*GeV2MeV) return false;
		incrementCounter("nPassing_triplet_pt",weight);
		
		_DEBUG("");
		
		// triplet p-value
		////////////////////////////////////////////////////////////////////////////////
		if(isCounter("nPassing_triplet_pvalue") && vertex_pval<0.1) return false; //////
		incrementCounter("nPassing_triplet_pvalue",weight); ////////////////////////////
		////////////////////////////////////////////////////////////////////////////////
		histos[name+"_mOS_after_pval"]->Fill(mOS1,weight/2.);
		histos[name+"_mOS_after_pval"]->Fill(mOS2,weight/2.);
		histos[name+"_mSS_after_pval"]->Fill(mSS,weight);
		fillHistsMassPt3mu(vtx,name,"_after_pval",histos,histos2,weight,mBlindMin,mBlindMax);
		fillHistsDoublets(vtx,name,"_after_pval",histos,histos2,weight,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS);

		_DEBUG("");

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(isCounter("nPassing_triplet_vtxclean") && (vertex_lxy<0. || vertex_a0xy>0.05 || vertex_cosTxy<0.96)) return false; /////
		incrementCounter("nPassing_triplet_vtxclean",weight); /////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		histos[name+"_triplet_dRmin_after_vtxclean"]->Fill(dRmin,weight);
		histos[name+"_triplet_dRmax_after_vtxclean"]->Fill(dRmax,weight);
		histos2[name+"_m3mu_vs_dRmin_after_vtxclean"]->Fill(dRmin,m3mu,weight); // need to blind Y axis
		histos2[name+"_m3mu_vs_dRmax_after_vtxclean"]->Fill(dRmax,m3mu,weight); // need to blind Y axis
		if(foundOS1) histos[name+"_mOS_after_vtxclean"]->Fill(mOS1,weightOS);
		if(foundOS2) histos[name+"_mOS_after_vtxclean"]->Fill(mOS2,weightOS);
		if(foundSS)  histos[name+"_mSS_after_vtxclean"]->Fill(mSS,weight);
		fillHistsMassPt3mu(vtx,name,"_after_vtxclean",histos,histos2,weight,mBlindMin,mBlindMax);
		fillHistsDoublets(vtx,name,"_after_vtxclean",histos,histos2,weight,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS);
		
		_DEBUG("");
		
		histos[name+"_triplet_isolation_before_isolation"]->Fill(isolation,weight);
		histos[name+"_triplet_isolation_zoom_before_isolation"]->Fill(isolation,weight);
		double vtxisolation = tripletIsolation(vtx,0.03); // used for the cut !!!
		//////////////////////////////////////////////////////////////////////////////////////	
		if(isCounter("nPassing_triplet_isolation") && vtxisolation>0.08) return false; //////////
		incrementCounter("nPassing_triplet_isolation",weight); ///////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////
		if(!reBlind(histos[name+"_m3mu_after_isolation"],m3mu,mBlindMin,mBlindMax))           histos[name+"_m3mu_after_isolation"]->Fill(m3mu,weight);
		if(!reBlind(histos[name+"_m3mu_lin_after_isolation"],m3mu,mBlindMin,mBlindMax))       histos[name+"_m3mu_lin_after_isolation"]->Fill(m3mu,weight);
		if(!reBlind(histos[name+"_m3mu_lin_zoom_after_isolation"],m3mu,mBlindMin,mBlindMax))  histos[name+"_m3mu_lin_zoom_after_isolation"]->Fill(m3mu,weight);
		if(!reBlind(histos[name+"_m3mu_sigregion_after_isolation"],m3mu,mBlindMin,mBlindMax)) histos[name+"_m3mu_sigregion_after_isolation"]->Fill(m3mu,weight);
		histos[name+"_pT3mu_after_isolation"]->Fill(pT3mu,weight);
		
		
		_DEBUG("");
		
		// before any W cuts are applied
		fillHistsMassPt3mu(vtx,name,"_before_W",histos,histos2,weight,mBlindMin,mBlindMax);
		///////////////////////////////////////////////////////////////////////////////////
		if(isCounter("nPassing_W_MET") && MET_RefFinal_et<15.*GeV2MeV) return false; //////
		incrementCounter("nPassing_W_MET",weight); ////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////
		fillHistsMassPt3mu(vtx,name,"_after_MET",histos,histos2,weight,mBlindMin,mBlindMax);
		
		_DEBUG("");
		
		//////////////////////////////////////////////////////////////////////////////////
		if(isCounter("nPassing_W_dphi3muMET") && met_reffinal_dPhi3mu<2) return false; ///
		incrementCounter("nPassing_W_dphi3muMET",weight); ////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////
		fillHistsMassPt3mu(vtx,name,"_after_dphi3muMET",histos,histos2,weight,mBlindMin,mBlindMax);
		
		_DEBUG("");
		
		///////////////////////////////////////////////////////////////////////////////
		if(isCounter("nPassing_W_mT") && met_reffinal_mT<60.*GeV2MeV) return false; ///
		incrementCounter("nPassing_W_mT",weight); /////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////
		fillHistsMassPt3mu(vtx,name,"_after_mT",histos,histos2,weight,mBlindMin,mBlindMax);
		
		_DEBUG("");
	}
	else if(doMVA)
	{
		_DEBUG("");
		
		/////////////////////////////////
		method.ReplaceAll("MVA:",""); ///
		/////////////////////////////////
		
		_DEBUG("");
		
		fillHistsMassPt3mu(vtx,name,"_before_MVA",histos,histos2,weight,mBlindMin,mBlindMax);
		/////////////////////////////////////////////////////////////////////////////////////
		unsigned int ivtx = 0;
		for(unsigned int i=0 ; i<vertices.size() ; ++i) { ivtx=vertices[i].vtxIndex(); if(ivtx==vtx) { ivtx=i; break; } }
		setMVAvars(ivtx, flatout_ints, flatout_floats, flatout_vints, flatout_vfloats); /////
		setMVAspect(ivtx, flatout_ints, flatout_floats, flatout_vints, flatout_vfloats); ////
		Double_t MVAscore = getMVAscore(); //////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		histos[name+"_MVA_score_all"]->Fill(MVAscore,weight);
		if(vertex_mass<mBlindMin) histos[name+"_MVA_score_left_sideband"]->Fill(MVAscore,weight);
		if(vertex_mass>mBlindMax) histos[name+"_MVA_score_right_sideband"]->Fill(MVAscore,weight);
		histos2[name+"_m3mu_vs_MVA_score"]->Fill(MVAscore,vertex_mass,weight);
		histos2[name+"_pT3mu_vs_MVA_score"]->Fill(MVAscore,vertex_pt,weight);
		histos2[name+"_mOS1_vs_MVA_score"]->Fill(MVAscore,vertex_mOS1,weight);
		histos2[name+"_mOS2_vs_MVA_score"]->Fill(MVAscore,vertex_mOS2,weight);
		histos2[name+"_mSS_vs_MVA_score"]->Fill(MVAscore,vertex_mSS,weight);
		/////////////////////////////////////////////////////////////////////////////////////
		// the MVA cut itself ///////////////////////////////////////////////////////////////
		if(isCounter("nPassing_MVA_vtxmet") && MVAscore<0.985) return false; ////////////////
		incrementCounter("nPassing_MVA_vtxmet",weight); /////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		fillHistsMassPt3mu(vtx,name,"_after_MVA",histos,histos2,weight,mBlindMin,mBlindMax);
		
		_DEBUG("");
	}
	
	_DEBUG("");
	
	bool isSidebands = (vertex_mass>1300. && vertex_mass<2300.);
	if(isCounter("nPassing_evt_sidebands") && isSidebands) incrementCounter("nPassing_evt_sidebands",weight); // do not return false
	
	_DEBUG("");
	
	//////////////////////////
	//// fill some histos ////
	//////////////////////////
	fillTriggerbitsHist(histos[name+"_triggers_after_vtxmet"]);                 // after vertex cleaning cut
	fillTriggerbitsHist(histos[name+"_unique_triggers_after_vtxmet"],"unique"); // after vertex cleaning cut
	fillHistsMassPt3mu(vtx,name,"_after_vtxmet",histos,histos2,weight,mBlindMin,mBlindMax);
	fillHistsDoublets(vtx,name,"_after_vtxmet",histos,histos2,weight,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS);
	fillHistsMVAvars(vtx,name,"_after_vtxmet",histos,histos2,pMuSum,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS,weight,mBlindMin,mBlindMax);
	
	_DEBUG("");
	
	return true;
}



bool acceptTriplet(unsigned int vtx, TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight=1., double mBlindMin=1500., double mBlindMax=2000.)
{	
	_DEBUG("");
	
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];	
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3,psum;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
	psum = p1+p2+p3; // psum = getTlv3mu(vtx); // could be done also like that
	
	_DEBUG("");
	
	// OS and SS pairs
	TLorentzVector pOS1, pOS2, pSS;
	bool foundOS1 = false;
	bool foundOS2 = false;
	bool foundSS  = false;
	int iUniqueCharge = -1;
	TMapVL doubletsTLV = getDoubletsFromTriplet(vtx,iUniqueCharge);
	if(doubletsTLV.find("OS1") != doubletsTLV.end()) { pOS1 = doubletsTLV["OS1"]; foundOS1 = true; }
	if(doubletsTLV.find("OS2") != doubletsTLV.end()) { pOS2 = doubletsTLV["OS2"]; foundOS2 = true; }
	if(doubletsTLV.find("SS")  != doubletsTLV.end()) { pSS  = doubletsTLV["SS"];  foundSS  = true; }
	if(foundOS1+foundOS2<2 && fabs(vtx_charge->at(vtx))==1.) _ERROR("could not find OS doublet"); // should be exactly two OS pairs
	if(!foundSS && fabs(vtx_charge->at(vtx))==1.)            _ERROR("could not find SS doublet"); // should be exactly one OS pair
	double mOS1 = (foundOS1) ? pOS1.M() : -1.;
	double mOS2 = (foundOS2) ? pOS2.M() : -1.;
	double mSS  = (foundSS)  ? pSS.M()  : -1.;
	double weightOS = (foundOS1+foundOS2==2) ? weight/2. : weight;
	
	_DEBUG("");
	
	double mass   = vtx_mass->at(vtx); // psum.M();  in principle
	double pTsum  = vtx_pt->at(vtx);   // psum.Pt(); in principle
	double charge = vtx_charge->at(vtx);
	double pvalue = TMath::Prob(vtx_chi2->at(vtx),vtx_ndf->at(vtx));

	_DEBUG("");
	
	
	// some histos before the triplet cuts
	double dRmin = +1.e20;
	dRmin = (psum.DeltaR(p1)<dRmin) ? psum.DeltaR(p1) : dRmin;
	dRmin = (psum.DeltaR(p2)<dRmin) ? psum.DeltaR(p2) : dRmin;
	dRmin = (psum.DeltaR(p3)<dRmin) ? psum.DeltaR(p3) : dRmin;
	double dRmax = -1.e20;
	dRmax = (psum.DeltaR(p1)>dRmax) ? psum.DeltaR(p1) : dRmax;
	dRmax = (psum.DeltaR(p2)>dRmax) ? psum.DeltaR(p2) : dRmax;
	dRmax = (psum.DeltaR(p3)>dRmax) ? psum.DeltaR(p3) : dRmax;
	histos[name+"_triplet_dRmin"]->Fill(dRmin,weight);
	histos2[name+"_m3mu_vs_dRmin"]->Fill(dRmin,mass,weight);
	histos[name+"_triplet_dRmax"]->Fill(dRmax,weight);
	histos2[name+"_m3mu_vs_dRmax"]->Fill(dRmax,mass,weight);
	histos[name+"_pT3mu_before_triplet"]->Fill(pTsum,weight); // NEED TO BLIND THIS !!!
	
	_DEBUG("");
	
	// corrected muon isolation
	iso10 = correctedIsolation(vtx,0.1);
	iso20 = correctedIsolation(vtx,0.2);
	iso30 = correctedIsolation(vtx,0.3);
	iso40 = correctedIsolation(vtx,0.4);
	isolation = tripletIsolation(vtx,0.03); // used for the cut !!!

	histos[name+"_triplet_isolation"]->Fill(isolation,weight);
	histos[name+"_triplet_isolation_zoom"]->Fill(isolation,weight);
	if(iso10.size()>0)
	{
		histos[name+"_mu_ptcone10.1"]->Fill(iso10[0],weight);
		histos[name+"_mu_ptcone10.2"]->Fill(iso10[1],weight);
		histos[name+"_mu_ptcone10.3"]->Fill(iso10[2],weight);
		histos[name+"_mu_ptcone10.1_zoom"]->Fill(iso10[0],weight);
		histos[name+"_mu_ptcone10.2_zoom"]->Fill(iso10[1],weight);	
		histos[name+"_mu_ptcone10.3_zoom"]->Fill(iso10[2],weight);
	}
	if(iso20.size()>0)
	{
		histos[name+"_mu_ptcone20.1"]->Fill(iso20[0],weight);
		histos[name+"_mu_ptcone20.2"]->Fill(iso20[1],weight);
		histos[name+"_mu_ptcone20.3"]->Fill(iso20[2],weight);
		histos[name+"_mu_ptcone20.1_zoom"]->Fill(iso20[0],weight);	
		histos[name+"_mu_ptcone20.2_zoom"]->Fill(iso20[1],weight);	
		histos[name+"_mu_ptcone20.3_zoom"]->Fill(iso20[2],weight);
	}
	if(iso30.size()>0)
	{
		histos[name+"_mu_ptcone30.1"]->Fill(iso30[0],weight);
		histos[name+"_mu_ptcone30.2"]->Fill(iso30[1],weight);
		histos[name+"_mu_ptcone30.3"]->Fill(iso30[2],weight);
		histos[name+"_mu_ptcone30.1_zoom"]->Fill(iso30[0],weight);	
		histos[name+"_mu_ptcone30.2_zoom"]->Fill(iso30[1],weight);	
		histos[name+"_mu_ptcone30.3_zoom"]->Fill(iso30[2],weight);
	}
	if(iso40.size()>0)
	{
		histos[name+"_mu_ptcone40.1"]->Fill(iso40[0],weight);
		histos[name+"_mu_ptcone40.2"]->Fill(iso40[1],weight);
		histos[name+"_mu_ptcone40.3"]->Fill(iso40[2],weight);
		histos[name+"_mu_ptcone40.1_zoom"]->Fill(iso40[0],weight);	
		histos[name+"_mu_ptcone40.2_zoom"]->Fill(iso40[1],weight);	
		histos[name+"_mu_ptcone40.3_zoom"]->Fill(iso40[2],weight);
	}
	// can calculate the "3mu" isolation given the dRmax cone per event
	// for that, I need the entire track container (==> do it on the AOD level)
	
	_DEBUG("");
	
	
	//////////////////////////////////////////
	//////////////  triplet cuts /////////////
	//////////////////////////////////////////
	
	// remove very bad-quality vertices
	histos[name+"_pvalue_all"]->Fill(pvalue,weight);
	if(isCounter("nPassing_triplet_pvalueLoose") && pvalue<1e-4) return false;
	incrementCounter("nPassing_triplet_pvalueLoose",weight);
	
	_DEBUG("");
	
	// 3mu vertex mass
	if(isCounter("nPassing_triplet_m") && mass>4000.) return false;
	incrementCounter("nPassing_triplet_m",weight);
	
	_DEBUG("");
	
	// 3body vertex pT loose cut
	if(isCounter("nPassing_triplet_ptLoose") && pTsum<15.*GeV2MeV) return false;
	incrementCounter("nPassing_triplet_ptLoose",weight);
	
	_DEBUG("");
	
	// 3mu vertex sum charge
	if(isCounter("nPassing_triplet_charge") && fabs(charge)!=1.) return false;
	incrementCounter("nPassing_triplet_charge",weight);

	_DEBUG("");

	// fill the cosThetaBoost histo
	TLorentzVector pUniqueQ;
	if     (iUniqueCharge==1) pUniqueQ = p1;
	else if(iUniqueCharge==2) pUniqueQ = p2;
	else if(iUniqueCharge==3) pUniqueQ = p3;
	else _ERROR("vertex |charge|=1 but no unique charge track was found");
	TVector3 boostVector = psum.BoostVector(); // this is the 3vector of the 3mu (tau)
	pUniqueQ.Boost( -boostVector ); // boost p to the 3mu (tau) CM rest frame
	Double_t cosTh = pUniqueQ.Vect()*psum.Vect()/(pUniqueQ.P()*psum.P());
	histos[name+"_triplet_cosTheta"]->Fill(cosTh);

	_DEBUG("");	

	// remove rho/omega and phi
	if(foundOS1 && foundSS) histos2[name+"_mSS_vs_mOS_beforeCut"]->Fill(mOS1,mSS,weightOS);
	if(foundOS2 && foundSS) histos2[name+"_mSS_vs_mOS_beforeCut"]->Fill(mOS2,mSS,weightOS);
	bool mHigh     = ((mOS1>2400. || mOS2>2400.) || mSS>2400.);
	// bool mLow      = (mSS<500. && (mOS1<1100. || (mOS2<1100. && mOS1<1100.)));
	bool mDiagLow  = (mSS<(500.-mOS1)  || mSS<(500.-mOS2));
	bool mDiagHigh = (mSS>(3300.-mOS1) || mSS>(3300.-mOS2));
	bool fail = (foundOS1 && foundOS2 && foundSS && (mHigh || /*mLow ||*/ mDiagHigh || mDiagLow));
	if(isCounter("nPassing_triplet_mOSmSS") && fail) return false;
	incrementCounter("nPassing_triplet_mOSmSS",weight);
	if(foundOS1 && foundSS) histos2[name+"_mSS_vs_mOS_afterCut"]->Fill(mOS1,mSS,weightOS);
	if(foundOS2 && foundSS) histos2[name+"_mSS_vs_mOS_afterCut"]->Fill(mOS2,mSS,weightOS);

	_DEBUG("");
	
	// pT fraction cut
	histos[name+"_mupt1Fraction"]->Fill(ptFraction1,weight);   histos2[name+"_mupt1Fraction_vs_mupt1"]->Fill(Pt1,ptFraction1,weight);
	histos[name+"_mupt2Fraction"]->Fill(ptFraction2,weight);   histos2[name+"_mupt2Fraction_vs_mupt2"]->Fill(Pt2,ptFraction2,weight);
	histos[name+"_mupt3Fraction"]->Fill(ptFraction3,weight);   histos2[name+"_mupt3Fraction_vs_mupt3"]->Fill(Pt3,ptFraction3,weight);
	histos[name+"_mupt12Fraction"]->Fill(ptFraction12,weight); histos2[name+"_mupt12Fraction_vs_mupt12"]->Fill(dPt12,ptFraction12,weight);
	histos[name+"_mupt23Fraction"]->Fill(ptFraction23,weight); histos2[name+"_mupt23Fraction_vs_mupt23"]->Fill(dPt23,ptFraction23,weight);
	histos[name+"_mupt13Fraction"]->Fill(ptFraction13,weight); histos2[name+"_mupt13Fraction_vs_mupt13"]->Fill(dPt13,ptFraction13,weight);
	histos2[name+"_mupt23Fraction_vs_mupt12Fraction"]->Fill(ptFraction12,ptFraction23,weight);
	histos2[name+"_mupt13Fraction_vs_mupt12Fraction"]->Fill(ptFraction12,ptFraction13,weight);
	double slope1 = (1.00-0.3)/(30000.-5000.); double cutoff1 = 0.3 - slope1*5000.; TString formula1 = "x*"+_s(slope1)+"+"+_s(cutoff1);
	double slope2 = (0.50-0.0)/(12500.-0.000); double cutoff2 = 0.0 - slope2*0.000; TString formula2 = "x*"+_s(slope2)+"+"+_s(cutoff2);
	double slope3 = (0.35-0.1)/(10000.-0.000); double cutoff3 = 0.1 - slope3*0.000; TString formula3 = "x*"+_s(slope3)+"+"+_s(cutoff3);
	TF1 f1("f1",formula1); bool ptfrac1 = (ptFraction1>f1.Eval(Pt1));
	TF1 f2("f2",formula2); bool ptfrac2 = (ptFraction2>f2.Eval(Pt2));
	TF1 f3("f3",formula3); bool ptfrac3 = (ptFraction3>f3.Eval(Pt3));
	if(isCounter("nPassing_triplet_ptfrac") && (ptfrac1 || ptfrac2 || ptfrac3)) return false;
	incrementCounter("nPassing_triplet_ptfrac",weight);
	if(pTsum>25.*GeV2MeV)
	{
		histos[name+"_mupt1Fraction_afterCut"]->Fill(ptFraction1,weight); histos2[name+"_mupt1Fraction_vs_mupt1_afterCut"]->Fill(Pt1,ptFraction1,weight);
		histos[name+"_mupt2Fraction_afterCut"]->Fill(ptFraction2,weight); histos2[name+"_mupt2Fraction_vs_mupt2_afterCut"]->Fill(Pt2,ptFraction2,weight);
		histos[name+"_mupt3Fraction_afterCut"]->Fill(ptFraction3,weight); histos2[name+"_mupt3Fraction_vs_mupt3_afterCut"]->Fill(Pt3,ptFraction3,weight);
	}
	_DEBUG("");

	
	sources src;
	getSrc(vtx,src);
	if(!validateSrcChain(vtx)) _FATAL("wrong chain");
	double mRho    = 770.;  // MeV
	double mOmega  = 882.;  // MeV
	double mPhi    = 1020.; // MeV
	double margins = 25.;   // +/- (MeV)
	
	
	TLorentzVector pOSa, pOSb, p1a,p2a,p3a, p1b,p2b,p3b;
	// double q1a,q2a,q3a, q1b,q2b,q3b;
	double q1a,q1b;
	bool okOS = true;
	if(src.q2body12==0. && src.q2body13==0.)
	{
		pOSa = src.p2body12; p1a = src.srcTlv[0]; p2a = src.srcTlv[1]; p3a = src.srcTlv[2];
		pOSb = src.p2body13; p1b = src.srcTlv[0]; p2a = src.srcTlv[2]; p3a = src.srcTlv[1];
		q1a = qtrk(trks_qoverp->at(src.trkIndex[0])); q1b = qtrk(trks_qoverp->at(src.trkIndex[0]));
		// q2a = qtrk(trks_qoverp->at(src.trkIndex[1])); q2b = qtrk(trks_qoverp->at(src.trkIndex[2]));
		// q3a = qtrk(trks_qoverp->at(src.trkIndex[2])); q3b = qtrk(trks_qoverp->at(src.trkIndex[1]));
	}
	else if(src.q2body12==0. && src.q2body23==0.)
	{
		pOSa = src.p2body12; p1a = src.srcTlv[0]; p2a = src.srcTlv[1]; p3a = src.srcTlv[2];
		pOSb = src.p2body23; p1b = src.srcTlv[1]; p2a = src.srcTlv[2]; p3a = src.srcTlv[0];
		q1a = qtrk(trks_qoverp->at(src.trkIndex[0])); q1b = qtrk(trks_qoverp->at(src.trkIndex[1]));
		// q2a = qtrk(trks_qoverp->at(src.trkIndex[1])); q2b = qtrk(trks_qoverp->at(src.trkIndex[2]));
		// q3a = qtrk(trks_qoverp->at(src.trkIndex[2])); q3b = qtrk(trks_qoverp->at(src.trkIndex[0]));
	}
	else if(src.q2body13==0. && src.q2body23==0.)
	{
		pOSa = src.p2body13; p1a = src.srcTlv[0]; p2a = src.srcTlv[2]; p3a = src.srcTlv[1];
		pOSb = src.p2body23; p1b = src.srcTlv[1]; p2a = src.srcTlv[2]; p3a = src.srcTlv[0];
		q1a = qtrk(trks_qoverp->at(src.trkIndex[0])); q1b = qtrk(trks_qoverp->at(src.trkIndex[1]));
		// q2a = qtrk(trks_qoverp->at(src.trkIndex[2])); q2b = qtrk(trks_qoverp->at(src.trkIndex[2]));
		// q3a = qtrk(trks_qoverp->at(src.trkIndex[1])); q3b = qtrk(trks_qoverp->at(src.trkIndex[0]));
	}
	else { _ERROR("couldn't find 2 OS pairs from this src"); okOS = false; }
	if(okOS)
	{
		bool okOSa = (fabs(pOSa.M()-mRho)<margins || fabs(pOSa.M()-mOmega)<margins || fabs(pOSa.M()-mPhi)<margins);
		bool okOSb = (fabs(pOSb.M()-mRho)<margins || fabs(pOSb.M()-mOmega)<margins || fabs(pOSb.M()-mPhi)<margins);
		double weightOSab = -1.;
		if     (okOSa && okOSb)  weightOSab = weight/2.;
		else if(okOSa && !okOSb) weightOSab = weight;
		else if(!okOSa && okOSb) weightOSab = weight;
		
		TVector3 boostVector1 = pOSa.BoostVector(); // this is the 3vector of the 1st OS pair
		p1a.Boost( -boostVector1 ); // boost p1a to the 1st OS pair CM rest frame
		p2a.Boost( -boostVector1 ); // boost p2a to the 1st OS pair CM rest frame
		p3a.Boost( -boostVector1 ); // boost p3a to the 1st OS pair CM rest frame
		Double_t cosThOS1     = p1a.Vect()*p2a.Vect()/(p1a.P()*p2a.P());
		Double_t cosThOSQ1    = (q1a<0.) ? p1a.Vect()*pOSa.Vect()/(p1a.P()*pOSa.P()) : p2a.Vect()*pOSa.Vect()/(p2a.P()*pOSa.P());
		Double_t cosThOS1trk3 = p3a.Vect()*pOSa.Vect()/(p3a.P()*pOSa.P());
		double dPhiOS1        = fabs(dPhi(p1a.Phi(),p2a.Phi()));
		
		TVector3 boostVector2 = pOSb.BoostVector(); // this is the 3vector of the 2nd OS pair
		p1b.Boost( -boostVector2 ); // boost p1b to the 2nd OS pair CM rest frame
		p2b.Boost( -boostVector2 ); // boost p2b to the 2nd OS pair CM rest frame
		p3b.Boost( -boostVector2 ); // boost p3b to the 2nd OS pair CM rest frame
		Double_t cosThOS2     = p1b.Vect()*p2b.Vect()/(p1b.P()*p2b.P());
		Double_t cosThOSQ2    = (q1b<0.) ? p1b.Vect()*pOSb.Vect()/(p1b.P()*pOSb.P()) : p2b.Vect()*pOSb.Vect()/(p2b.P()*pOSb.P());
		Double_t cosThOS2trk3 = p3b.Vect()*pOSb.Vect()/(p3b.P()*pOSb.P());
		double dPhiOS2        = fabs(dPhi(p1b.Phi(),p2b.Phi()));
		
		double dmOS = fabs(pOSa.M()-pOSb.M());
		
		histos2[name+"_cosThOSQ_vs_mOS"]->Fill(pOSa.M(),cosThOSQ1,weight/2.);
		histos2[name+"_cosThOSQ_vs_mOS"]->Fill(pOSb.M(),cosThOSQ2,weight/2.);
		if(okOSa) histos[name+"_cosThOSQ_lightMesons"]->Fill(cosThOSQ1,weightOSab);
		if(okOSb) histos[name+"_cosThOSQ_lightMesons"]->Fill(cosThOSQ2,weightOSab);
		
		histos2[name+"_cosThOS_vs_mOS"]->Fill(pOSa.M(),cosThOS1,weight/2.);
		histos2[name+"_cosThOS_vs_mOS"]->Fill(pOSb.M(),cosThOS2,weight/2.);
		if(okOSa) histos[name+"_cosThOS_lightMesons"]->Fill(cosThOS1,weightOSab);
		if(okOSb) histos[name+"_cosThOS_lightMesons"]->Fill(cosThOS2,weightOSab);
		
		histos2[name+"_dPhiOS_vs_mOS"]->Fill(pOSa.M(),dPhiOS1,weight/2.);
		histos2[name+"_dPhiOS_vs_mOS"]->Fill(pOSb.M(),dPhiOS2,weight/2.);
		if(okOSa) histos[name+"_dPhiOS_lightMesons"]->Fill(dPhiOS1,weightOSab);
		if(okOSb) histos[name+"_dPhiOS_lightMesons"]->Fill(dPhiOS2,weightOSab);
		
		histos2[name+"_cosThOStrk3_vs_mOS"]->Fill(pOSa.M(),cosThOS1trk3,weight/2.);
		histos2[name+"_cosThOStrk3_vs_mOS"]->Fill(pOSb.M(),cosThOS2trk3,weight/2.);
		if(okOSa) histos[name+"_cosThOStrk3_lightMesons"]->Fill(cosThOS1trk3,weightOSab);
		if(okOSb) histos[name+"_cosThOStrk3_lightMesons"]->Fill(cosThOS2trk3,weightOSab);
		
		histos2[name+"_dmOS_vs_mOS"]->Fill(pOSa.M(),dmOS,weight);
		histos2[name+"_dmOS_vs_mOS"]->Fill(pOSb.M(),dmOS,weight);
		histos[name+"_dmOS_lightMesons"]->Fill(dmOS,weight);
	}
	
	
	_DEBUG("");
	
	
	//////////////////////////
	//// fill some histos ////
	//////////////////////////
	fillHistsMassPt3mu(vtx,name,"_after_triplet",histos,histos2,weight,mBlindMin,mBlindMax);
	fillHistsDoublets(vtx,name,"_after_triplet",histos,histos2,weight,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS);
	fillCategories(vtx,name,"tripletCategories_after_triplet",histos);
	fillCategories(vtx,name,"tripletCategories_norm_after_triplet",histos);
	
	return true;
}
bool acceptHadClean(unsigned int vtx, TString name, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, double weight=1., double mBlindMin=1500., double mBlindMax=2000.)
{
	// the triplet...
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];	
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV);
	p2.SetXYZM(px2,py2,pz2,muonMassMeV);
	p3.SetXYZM(px3,py3,pz3,muonMassMeV);
	TLorentzVector pMuSum = getTlv3mu(vtx);
	double m3mu           = vtx_mass->at(vtx); // pMuSum.M();
	// double pT3mu       = vtx_pt->at(vtx);   // pMuSum.Pt();
	
	// OS and SS pairs
	TLorentzVector pOS1, pOS2, pSS;
	bool foundOS1 = false;
	bool foundOS2 = false;
	bool foundSS  = false;
	int iUniqueCharge = -1;
	TMapVL doubletsTLV = getDoubletsFromTriplet(vtx,iUniqueCharge);
	if(doubletsTLV.find("OS1") != doubletsTLV.end()) { pOS1 = doubletsTLV["OS1"]; foundOS1 = true; }
	if(doubletsTLV.find("OS2") != doubletsTLV.end()) { pOS2 = doubletsTLV["OS2"]; foundOS2 = true; }
	if(doubletsTLV.find("SS")  != doubletsTLV.end()) { pSS  = doubletsTLV["SS"];  foundSS  = true; }
	if(foundOS1+foundOS2<2 && fabs(vtx_charge->at(vtx))==1.) _ERROR("could not find OS doublet"); // should be exactly two OS pairs
	if(!foundSS && fabs(vtx_charge->at(vtx))==1.)            _ERROR("could not find SS doublet"); // should be exactly one SS pair

	
	// jet stuff
	int nStdBtaggedJets30GeV = 0;
	int nStdBtaggedJets20GeV = 0;
	int nMedBtaggedJets30GeV = 0;
	int nMedBtaggedJets20GeV = 0;
	int nTgtBtaggedJets30GeV = 0;
	int nTgtBtaggedJets20GeV = 0;
	nDijets = 0;
	nDijetsLoose = 0;
	nJet3muOverlaps = 0;
	nJet3muOverlapsLoose = 0;

	// multimap<double,int> pt2i;
	// vector<int> ijet;
	// for(int i=0 ; i<(int)jets_pt->size() ; i++)
	// {
	// 	bool isUgly = (jets_isUgly->at(i));
	// 	if(isUgly)  continue;
	// 	pt2i.insert(make_pair(jets_pt->at(i),i));
	// }
	// for(multimap<double,int>::reverse_iterator rit=pt2i.rbegin() ; rit!=pt2i.rend() ; ++rit) ijet.push_back(rit->second);
	// unsigned int njet = ijet.size();
	
	multimap<double,int> pt2i;
	vector<int> ijet;
	for(int i=0 ; i<AntiKt4LCTopoLets_n ; i++)
	{
		if(!isGoodJet(i,JetQuality))  continue;
		pt2i.insert(make_pair(AntiKt4LCTopoLets_pt->at(i),i));
	}
	for(multimap<double,int>::reverse_iterator rit=pt2i.rbegin() ; rit!=pt2i.rend() ; ++rit) ijet.push_back(rit->second);
	unsigned int njet = ijet.size();
	

	// for b-jet veto
	for(unsigned int i=0 ; i<njet && i<4 ; ++i)
	{
		int j = ijet[i];
		double mj  = AntiKt4LCTopoLets_m->at(j);
		double ptj = AntiKt4LCTopoLets_pt->at(j);
		double MV1 = AntiKt4LCTopoLets_flavor_weight_MV1->at(j);
		histos2[name+"_flwMV1_vs_mjet"]->Fill(mj,MV1,weight);
		histos2[name+"_flwMV1_vs_pTjet"]->Fill(ptj,MV1,weight);
		if(ptj>20.*GeV2MeV && mj>4.*GeV2MeV && isStandardBtag(MV1)) nStdBtaggedJets20GeV++;
		if(ptj>30.*GeV2MeV && mj>4.*GeV2MeV && isStandardBtag(MV1)) nStdBtaggedJets30GeV++;
		if(ptj>20.*GeV2MeV && mj>4.*GeV2MeV && isMediumBtag(MV1))   nMedBtaggedJets20GeV++;
		if(ptj>30.*GeV2MeV && mj>4.*GeV2MeV && isMediumBtag(MV1))   nMedBtaggedJets30GeV++;
		if(ptj>20.*GeV2MeV && mj>4.*GeV2MeV && isTightBtag(MV1))    nTgtBtaggedJets20GeV++;
		if(ptj>30.*GeV2MeV && mj>4.*GeV2MeV && isTightBtag(MV1))    nTgtBtaggedJets30GeV++;
	}

	double sumptj12     = -999.;
	double dPhiJet1Jet2 = -999.;
	double dPhi3muJet1  = -999.;
	if(njet>0)
	{
		int j1 = ijet[0];
		double phi1 = AntiKt4LCTopoLets_phi->at(j1);
		dPhi3muJet1 = fabs(dPhi(pMuSum.Phi(),phi1));
		double pt1  = AntiKt4LCTopoLets_pt->at(j1);

		// for leadingjet-3mu(tau) overlap veto
		if(pt1>25.*GeV2MeV && (dPhi3muJet1<0.2 || dPhi3muJet1>2.5))  nJet3muOverlaps++;
		if(pt1>35.*GeV2MeV && (dPhi3muJet1<0.2 || dPhi3muJet1>2.95)) nJet3muOverlapsLoose++;

		if(njet>1)
		{
			int j2 = ijet[1];
			sumptj12 = AntiKt4LCTopoLets_pt->at(j1)+AntiKt4LCTopoLets_pt->at(j2);
			double phi2  = AntiKt4LCTopoLets_phi->at(j2);
			double pt2   = AntiKt4LCTopoLets_phi->at(j2);
			dPhiJet1Jet2 = fabs(dPhi(phi1,phi2));

			// for dijet veto
			if(pt1>20.*GeV2MeV && pt2>20.*GeV2MeV && sumptj12>60.*GeV2MeV && dPhiJet1Jet2>2.5) nDijets++;
			if(pt1>20.*GeV2MeV && pt2>20.*GeV2MeV && sumptj12>75.*GeV2MeV && dPhiJet1Jet2>2.5) nDijetsLoose++;
		}
	}
	
	_DEBUG("");
	
	//// for flatout
	JetPt1 = 0.; JetM1 = 0.; JetE1 = 0.; JetEta1 = +5.; JetPhi1 = TMath::Pi()*1.2; JetMV1w1 = 0.;
	JetPt2 = 0.; JetM2 = 0.; JetE2 = 0.; JetEta2 = +5.; JetPhi2 = TMath::Pi()*1.2; JetMV1w2 = 0.;
	JetSumPt12 = 0.; JetdPhiJ1J2 = TMath::Pi()*1.2; JetdPhi3bodyJ1 = TMath::Pi()*1.2;
	if(njet>0 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV) { JetPt1 = AntiKt4LCTopoLets_pt->at(ijet[0]); JetEta1 = AntiKt4LCTopoLets_eta->at(ijet[0]); JetPhi1 = AntiKt4LCTopoLets_phi->at(ijet[0]); JetM1 = AntiKt4LCTopoLets_m->at(ijet[0]); JetE1 = AntiKt4LCTopoLets_E->at(ijet[0]); JetMV1w1 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ijet[0]); }
	if(njet>1 && AntiKt4LCTopoLets_pt->at(ijet[1])>20.*GeV2MeV) { JetPt2 = AntiKt4LCTopoLets_pt->at(ijet[1]); JetEta2 = AntiKt4LCTopoLets_eta->at(ijet[1]); JetPhi2 = AntiKt4LCTopoLets_phi->at(ijet[1]); JetM2 = AntiKt4LCTopoLets_m->at(ijet[1]); JetE2 = AntiKt4LCTopoLets_E->at(ijet[1]); JetMV1w2 = AntiKt4LCTopoLets_flavor_weight_MV1->at(ijet[1]); }
	if(njet>1 && JetPt1>20.*GeV2MeV && JetPt2>20.*GeV2MeV) { JetSumPt12 = sumptj12; JetdPhiJ1J2 = dPhiJet1Jet2; }
	if(njet>0 && JetPt1>20.*GeV2MeV) { JetdPhi3bodyJ1 = dPhi3muJet1; }
	
	_DEBUG("");
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// if(isCounter("nPassing_jets_bjetveto_std") && nMedBtaggedJets30GeV>0) return false; //////
	if(isCounter("nPassing_jets_bjetveto_std") && nStdBtaggedJets30GeV>0) return false; /////////
	incrementCounter("nPassing_jets_bjetveto_std",weight); //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	
	_DEBUG("");	
	
	// Fill some histos before doing loose had-cleaning
	if(njet>0 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV)                                                  { histos2[name+"_dPhi3muJet1_vs_pTjet1_before_looseveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0]),fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); histos[name+"_jets_dPhi3bodyJ1_before_looseveto"]->Fill(fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); }
	if(njet>1 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV && AntiKt4LCTopoLets_pt->at(ijet[1])>20.*GeV2MeV) { histos2[name+"_dPhiJet1Jet2_vs_sumpTjet12_before_looseveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0])+AntiKt4LCTopoLets_pt->at(ijet[1]),fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); histos[name+"_jets_dPhiJ1J2_before_looseveto"]->Fill(fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); }
	
	_DEBUG("");	

	// if(isWsignal(name) && (EventNumber==39989 || EventNumber==51445 || EventNumber==57442 || EventNumber==63271 || EventNumber==69066 || EventNumber==91929 || EventNumber==95447))
	// {
	// 	(*ofstr1) << EventNumber << " " << lbn << endl;
	// 	if(njet>0) (*ofstr1) << "  Jet1: pt=" << jets_pt->at(ijet[0]) << ", eta=" << jets_eta->at(ijet[0]) << ", phi=" << jets_phi->at(ijet[0]) << endl;
	// 	if(njet>1) (*ofstr1) << "  Jet2: pt=" << jets_pt->at(ijet[1]) << ", eta=" << jets_eta->at(ijet[1]) << ", phi=" << jets_phi->at(ijet[1]) << endl;
	// 	if(njet>0) (*ofstr1) << "  dPhi3muJet1=" << fabs(dPhi(pMuSum.Phi(),jets_phi->at(ijet[0]))) << endl;
	// 	if(njet>1) (*ofstr1) << "  dPhiJet1Jet2=" << fabs(dPhi(jets_phi->at(ijet[0]),jets_phi->at(ijet[1]))) << ", sumptj12=" << jets_pt->at(ijet[0])+jets_pt->at(ijet[1]) << endl;
	// 	(*ofstr1) << "-------------------------------------------" << endl;
	// }

	
	///////////////////////////////////////////////////////////////////////////////////////////
	bool isJet3muOverlapLoose = (njet>0 && AntiKt4LCTopoLets_pt->at(ijet[0])>35.*GeV2MeV) ? (fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0])))<0.2 || fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0])))>2.95) : false;
	if(isCounter("nPassing_jets_coljetvetoLoose") && isJet3muOverlapLoose) return false; //////
	incrementCounter("nPassing_jets_coljetvetoLoose",weight); /////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	
	// if(isWsignal(name)) (*ofstr) << EventNumber << " " << lbn << endl;
		
	_DEBUG("");
	
	////////////////////////////////////////////////////////////////////////////////////
	bool isDijetsLoose = (njet>1 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV && AntiKt4LCTopoLets_pt->at(ijet[1])>20.*GeV2MeV) ? (AntiKt4LCTopoLets_pt->at(ijet[0])+AntiKt4LCTopoLets_pt->at(ijet[1])>75.*GeV2MeV && fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1])))>2.5) : false;
	if(isCounter("nPassing_jets_dijetvetoLoose") && isDijetsLoose) return false; //////
	incrementCounter("nPassing_jets_dijetvetoLoose",weight); ///////////////////////////
	////////////////////////////////////////////////////////////////////////////////////
	
	_DEBUG("");
	
	if(njet>0 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV)                                                  { histos2[name+"_dPhi3muJet1_vs_pTjet1_after_looseveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0]),fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); histos[name+"_jets_dPhi3bodyJ1_after_looseveto"]->Fill(fabs(dPhi(pMuSum.Phi(),AntiKt4LCTopoLets_phi->at(ijet[0]))),weight); }
	if(njet>1 && AntiKt4LCTopoLets_pt->at(ijet[0])>20.*GeV2MeV && AntiKt4LCTopoLets_pt->at(ijet[1])>20.*GeV2MeV) { histos2[name+"_dPhiJet1Jet2_vs_sumpTjet12_after_looseveto"]->Fill(AntiKt4LCTopoLets_pt->at(ijet[0])+AntiKt4LCTopoLets_pt->at(ijet[1]),fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); histos[name+"_jets_dPhiJ1J2_after_looseveto"]->Fill(fabs(dPhi(AntiKt4LCTopoLets_phi->at(ijet[0]),AntiKt4LCTopoLets_phi->at(ijet[1]))),weight); }
	
	_DEBUG("");
	
	//////////////////////////////////////////////////////////////////////////////////////////
	if(isCounter("nPassing_met_metLoose") && MET_RefFinal_et<10.*GeV2MeV) return false; //////
	incrementCounter("nPassing_met_metLoose",weight); ////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	_DEBUG("");
	
	/////////////////////////////////////
	// fill some histos after all cuts //
	/////////////////////////////////////
	
	double dRmin = +1.e20;
	dRmin = (pMuSum.DeltaR(p1)<dRmin) ? pMuSum.DeltaR(p1) : dRmin;
	dRmin = (pMuSum.DeltaR(p2)<dRmin) ? pMuSum.DeltaR(p2) : dRmin;
	dRmin = (pMuSum.DeltaR(p3)<dRmin) ? pMuSum.DeltaR(p3) : dRmin;
	double dRmax = -1.e20;
	dRmax = (pMuSum.DeltaR(p1)>dRmax) ? pMuSum.DeltaR(p1) : dRmax;
	dRmax = (pMuSum.DeltaR(p2)>dRmax) ? pMuSum.DeltaR(p2) : dRmax;
	dRmax = (pMuSum.DeltaR(p3)>dRmax) ? pMuSum.DeltaR(p3) : dRmax;
	histos[name+"_triplet_dRmin_after_hadclean"]->Fill(dRmin,weight);
	histos[name+"_triplet_dRmax_after_hadclean"]->Fill(dRmax,weight);
	histos2[name+"_m3mu_vs_dRmin_after_hadclean"]->Fill(dRmin,m3mu,weight); // need to blind Y axis
	histos2[name+"_m3mu_vs_dRmax_after_hadclean"]->Fill(dRmax,m3mu,weight); // need to blind Y axis
	
	_DEBUG("");
	
	fillTriggerbitsHist(histos[name+"_triggers_after_hadclean"]);
	fillTriggerbitsHist(histos[name+"_unique_triggers_after_hadclean"],"unique");
	fillHistsMassPt3mu(vtx,name,"_after_hadclean",histos,histos2,weight,mBlindMin,mBlindMax);
	fillHistsDoublets(vtx,name,"_after_hadclean",histos,histos2,weight,pOS1,pOS2,pSS,foundOS1,foundOS2,foundSS);
	fillCategories(vtx,name,"tripletCategories_after_hadclean",histos);
	fillCategories(vtx,name,"tripletCategories_norm_after_hadclean",histos);
	
	_DEBUG("");
	
	// primary vertex and pileup histos
	unsigned int npv  = pv_index->size();
	int npv_type1or3 = 0;
	int npv_type1    = 0;
	int npv_type3    = 0;
	int ntrks_type1  = 0;
	for(unsigned int pv=0 ; pv<npv ; pv++)
	{			
		int pvtype = pv_type->at(pv);
		histos[name+"_pv_type"]->Fill(pvtype,weight);
		if(pvtype==1 || pvtype==3) npv_type1or3++;
		if(pvtype==1) { npv_type1++; ntrks_type1 = pv_ntrk->at(pv); }
		if(pvtype==3) npv_type3++;
	}
	histos[name+"_pv_ntrks"]->Fill(ntrks_type1,weight);
	histos[name+"_pv_npv"]->Fill(npv_type1or3,weight);
	histos[name+"_pv_npv_type1"]->Fill(npv_type1,weight);
	histos[name+"_pv_npv_type3"]->Fill(npv_type3,weight);
	histos[name+"_actualInteractionsPerXing"]->Fill(aux_actualIntPerXing,weight);
	
	_DEBUG("");
	
	return true;
}




Bool_t Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.
	return kTRUE;
}
Long64_t LoadTree(Long64_t entry, TChain* chain)
{
	// Set the environment to read one entry
	Int_t fCurrent = -1; //!current Tree number in a TChain
	if (!chain) return -5;
	Long64_t centry = chain->LoadTree(entry);
	if(centry<0) return centry;
	if(chain->GetTreeNumber()!=fCurrent)
	{
		fCurrent = chain->GetTreeNumber();
		Notify();
	}
	return centry;
}
void analysis(TString name, TMapTSP2TCHAIN& chains, TMapTSP2TTREE& otrees, TMapTSP2TH1& histos, TMapTSP2TH2& histos2, Long64_t nentriesMax)
{
	TString bname = name; // can be e.g. WmunuNp0-5 for binned sample or bbTomu15 for unbinned sample
	if(isBinned(name)) name = getUnbinnedName(name);
	// Now for e.g. WmunuNp0-5 we have: name=WmunuNpX always and bname=WmunuNp0-5
	// where for e.g. bbTomu15 we have: name=bname=bbTomu15
	TString pname = (name==bname) ? name : name+"."+bname; // for the printCounters call
	_INFO("name="+(string)name+", bname="+(string)bname);
	
	
	//// TMVA initialization
	if(doMVA)
	{
		TString mvawgtdir = "/afs/cern.ch/user/h/hod/svntest/ThreeMu/D3PDanalysis/Run/Tau3Mu/simple/weights/";
		TString MVAmethod = selectionMethod;
		MVAmethod.ReplaceAll("MVA:","");
		bookMVAvars();
		initTMVA(mvawgtdir+"TMVAClassification."+MVAmethod+".weights.xml",MVAmethod);
	}
	
	//// counters
	clearCounters();
	initCounters();

	bool isdata   = isData(name);
	bool issignal = isSignal(name);

	
	double mBlindMinInitial = 1670.;
	double mBlindMaxInitial = 1890.;
	double mBlindMin = mBlindMinInitial;
	double mBlindMax = mBlindMaxInitial;
	mBlindMinGlob = mBlindMinInitial;
	mBlindMaxGlob = mBlindMaxInitial;
	
	
	if(name=="Wtaunu_3mu" || (isdata && !skim)) initializePileup(); // PU tool is needed for BCH cleaning !
	if(!skim)
	{
		initBCH(isdata); // BCH cleaning should be used either for data or signal
		if(name=="Wtaunu_3mu")
		{
			initJES(false);  // Jet Energy Scale
			initJER();       // Jet Energy Resolution
			initJUN();       // Jet Energy Scale uncertainty
			initMET();       // MET Utility
		}
	}

	
	
	/////////////////////////////////////////////////////////////////////////////////
	// The only place where "bname" is used inside analysis() is here: //////////////
	/////////////////////////////////////////////////////////////////////////////////
	Long64_t chainentires = chains[bname]->GetEntriesFast();
	Long64_t nentries     = (nentriesMax>0 && nentriesMax<chainentires) ? nentriesMax : chainentires;
	Long64_t nbytes=0, nb=0;
	Int_t nEventsAccessed = 0;	
	for(Long64_t jentry=0 ; jentry<nentries ; jentry++)
	{
		Long64_t ientry = LoadTree(jentry,chains[bname]);
		if(ientry<0) break;
		nb = chains[bname]->GetEntry(jentry);
		nbytes += nb;
		nEventsAccessed++;
	
		// TrackParticles do not have builtin "phi"
		// so I have to set it by hand
		setTPmusPhi();
			
		// //////////////////////////////////////////////////
		// // Veto HF in JZxW *only* ////////////////////////
		// // before doing anything else ////////////////////
		// int Bhad = vetoHF(name); /////////////////////////
		// // if(Bhad>0) _INFO("HF event: pdgId="+_s(Bhad));
		// if(Bhad>0) continue; /////////////////////////////
		// //////////////////////////////////////////////////
		
		
		// cout << "Accessing event: " << nEventsAccessed << "\r" << flush;
		int allPassing = getCounter("nPassing_evt_all");
		if(allPassing%5000==0 && allPassing!=0) printCounters("cutflow",pname);


		// weights
		double wgt         = 1.;
		double wFONLLshape = (name.Contains("bbTomu15") || name.Contains("ccTomu15")) ? getFONLLShapeWeight(name) : 1.;
		double wFONLLflat  = (name.Contains("bbTomu15") || name.Contains("ccTomu15")) ? getFONLLFlatWeight(name)  : 1.;
		double wLumi       = (!isdata)                                                ? getSampleWeight(bname)    : 1.;
		double wKfac       = (name.Contains("NpX"))                                   ? getKfactorWeight(name)    : 1.;
		double wDijet      = (name.Contains("JZ"))                                    ? getDijetWeight(name)      : 1.;
		double wPileup     = (name=="Wtaunu_3mu")                                     ? getPileupWeight()         : 1.;
		wgt *= wFONLLshape;
		wgt *= wFONLLflat;
		wgt *= wLumi;
		wgt *= wKfac;
		wgt *= wDijet;
		wgt *= wPileup;
	
		// remember the weights
		weights.clear();
		weights.insert(make_pair("wFONLLshape",wFONLLshape));
		weights.insert(make_pair("wFONLLflat",wFONLLflat));
		weights.insert(make_pair("wLumi",wLumi));
		weights.insert(make_pair("wKfac",wKfac));
		weights.insert(make_pair("wDijet",wDijet));
		weights.insert(make_pair("wPileup",wPileup));
		weights.insert(make_pair("wgt",wgt));



		// counters
		resetCounterFlags();
		incrementCounter("nPassing_evt_all",wgt);
		
		// some local sizes
		unsigned int nvtx = vtx_chi2->size();
		

		
		_DEBUG("");

		// if(!skim && (isdata || name=="Wtaunu_3mu"))
		if(!skim && name=="Wtaunu_3mu")
		{
			seedJER();
			
			calibJets.clear();
			smearedJets.clear();
			badBCHJets.clear();
			calibJetsUnc.clear();
			calibJetsIndex.clear();
			clearCalibratedJets();
			
			unsigned int nJets = AntiKt4LCTopoLets_pt->size();
			TMapVL JetShiftsUp, JetShiftsDwn;
			
			// cout << "\n\n--------- New Event, nJets=" << nJets << " ---------" << endl;
			for(unsigned int jet=0 ; jet<nJets ; jet++)
			{
				TLorentzVector Jet = getJES(jet);
				TLorentzVector smearedJet = Jet;
				calibJets.push_back(Jet);
				addCalibratedJet(Jet);
				
				bool badBCHjet = isBadCalibJetBCH(jet);
				badBCHJets.push_back(badBCHjet);
				
				double quad  = getJUNC(jet, Jet, JetShiftsUp, JetShiftsDwn);
				calibJetsUnc.push_back(quad);
				calibJetsIndex.push_back(jet);
				
				smearedJet = getJER(Jet);
				smearedJets.push_back(smearedJet);
					
				//cout << "  Jet index " << jet << "(pT:" << AntiKt4LCTopoLets_pt->at(jet) << " -> " << Jet.Pt() << " -> " << smearedJet.Pt() << ") -> quadrature=" << quad*100 << "\% " << ", badBCHjet=" << badBCHjet << endl;
				
				fillJetCalibrationHistos(jet, Jet, smearedJet, name, histos, histos2, wgt);
			}
		
			METUtil::METObject calibMET = getMETU();
			// cout << " MET_RefFinal_et : " << MET_RefFinal_et << " -> " << calibMET.et() << endl;
			// cout << " MET_RefFinal_phi: " << MET_RefFinal_phi << " -> " << calibMET.phi() << endl;
			fillMETCalibrationHistos(calibMET, name, histos/*, histos2*/, wgt);
		}



		
		_DEBUG("");
		
		// truth matching
		if(!isdata && issignal)
		{
			_DEBUG("");
			vector<double> tmp;
			vector<vector<double> > vtx_position_mc;
			vector<vector<double> > vtx_position_fit;
			vector<vector<double> > vtx_positionErr_fit;
			for(unsigned int imc=0 ; imc<mc_pdgId->size() ; imc++)
			{				
				if(fabs(mc_pdgId->at(imc))!=15) continue;
				if(mc_status->at(imc)!=2)       continue;
				// if(!mc_has_decayvtx->at(i))  continue;
				vtx_position_mc.push_back(tmp);
				unsigned int index = vtx_position_mc.size()-1;
				vtx_position_mc[index].push_back(mc_decayvtx_x->at(imc));
				vtx_position_mc[index].push_back(mc_decayvtx_y->at(imc));
				vtx_position_mc[index].push_back(mc_decayvtx_z->at(imc));			
			}
			_DEBUG("");
			for(unsigned int vtx=0 ; vtx<nvtx ; vtx++)
			{
				double mass = vtx_mass->at(vtx);
				TLorentzVector pMuSum = getTlv3mu(vtx);
        	
				if(fabs(mass-tauMassMeV)<500.) // |mVtx-mTau|< 500 MeV
				{
					vtx_position_fit.push_back(tmp);
					vtx_positionErr_fit.push_back(tmp);
					unsigned int index = vtx_position_fit.size()-1;
					vtx_position_fit[index].push_back(vtx_x->at(vtx));
					vtx_position_fit[index].push_back(vtx_y->at(vtx));
					vtx_position_fit[index].push_back(vtx_z->at(vtx));
					vtx_positionErr_fit[index].push_back(vtx_xErr->at(vtx));
					vtx_positionErr_fit[index].push_back(vtx_yErr->at(vtx));
					vtx_positionErr_fit[index].push_back(vtx_zErr->at(vtx));
				}
			}
			_DEBUG("");
			for(unsigned int i=0 ; i<vtx_position_mc.size() ; i++)
			{
				double xmc = vtx_position_mc[i][0];
				double ymc = vtx_position_mc[i][1];
				double zmc = vtx_position_mc[i][2];
				for(unsigned int j=0 ; j<vtx_position_fit.size() ; j++)
				{
					double xfit = vtx_position_fit[j][0];
					double yfit = vtx_position_fit[j][1];
					double zfit = vtx_position_fit[j][2];
					double xErrfit = vtx_positionErr_fit[j][0];
					double yErrfit = vtx_positionErr_fit[j][1];
					double zErrfit = vtx_positionErr_fit[j][2];
					double reldifx = (xmc-xfit)/sqrt(xErrfit);
					double reldify = (ymc-yfit)/sqrt(yErrfit);
					double reldifz = (zmc-zfit)/sqrt(zErrfit);
					histos[name+"_vtx_trumatch_x"]->Fill(reldifx,wgt);
					histos[name+"_vtx_trumatch_y"]->Fill(reldify,wgt);
					histos[name+"_vtx_trumatch_z"]->Fill(reldifz,wgt);
				}
			}
		}
		
	
		
		_DEBUG("");
		
		
		
		
		// truth matching - signal only
		if(!skim && name=="Wtaunu_3mu")
		{
			clearTruth();
			getMC(mc_pdgId,mc_status,mc_barcode,mc_children,mc_parents,mc_pt,mc_eta,mc_phi);
			getMU(muons_pt,muons_eta,muons_phi);
			getCAL(calo_muons_pt,calo_muons_eta,calo_muons_phi);
			matchTruth();
			
			if(n_trumu_all==3)  incrementCounter("nPassing_trusigtriplet_all",wgt);
			if(n_trumu_acc==3)  incrementCounter("nPassing_trusigtriplet_acc",wgt);
			if(n_recmu_mat==3)  incrementCounter("nPassing_trusigtriplet_matchedmuons",wgt);
			if(muons_pt->size()>=3) incrementCounter("nPassing_trusigtriplet_atleast3recmuon",wgt);
		}
		
		_DEBUG("");
		
		// Trigger decision (don't cut yet)
		fillTriggerbits(); // only once per event...
		fillTriggerbitsHist(histos[name+"_triggers"]); // before any cut
		fillTriggerbitsHist(histos[name+"_unique_triggers"],"unique"); // before any cut
		
		_DEBUG("");
		
		//// clear the flags telling which 4th track was filled per vertex
		for(TMapuivi::iterator it=isFourthTrack.begin() ; it!=isFourthTrack.end() ; ++it) it->second.clear();
		isFourthTrack.clear();
		
		_DEBUG("");
		
		//// histos of quadruplets
	if(!skim)
	{
		unsigned int nMuons = muons_pt->size();
		if(nMuons>3)
		{
			TLorentzVector m1,m2,m3,m4,p3,p4;
			m1.SetPtEtaPhiM(muons_pt->at(0),muons_eta->at(0),muons_phi->at(0),muonMassMeV);
			m2.SetPtEtaPhiM(muons_pt->at(1),muons_eta->at(1),muons_phi->at(1),muonMassMeV);
			m3.SetPtEtaPhiM(muons_pt->at(2),muons_eta->at(2),muons_phi->at(2),muonMassMeV);
			m4.SetPtEtaPhiM(muons_pt->at(3),muons_eta->at(3),muons_phi->at(3),muonMassMeV);
			p3 = m2+m3+m4;
			p4 = m1+m2+m3+m4;
		
			int nTracks = trks_pt->size();
			int itrk1 = muons_inDetTrackIndex->at(0);
			int itrk2 = muons_inDetTrackIndex->at(1);
			int itrk3 = muons_inDetTrackIndex->at(2);
			int itrk4 = muons_inDetTrackIndex->at(3);
		
			bool ptmu  = (muons_pt->at(0)>4000. && muons_pt->at(1)>4000. && muons_pt->at(2)>4000. && muons_pt->at(3)>4000.);
			bool zmass = (fabs(p4.M()-91.*GeV2MeV)<30.*GeV2MeV);
			bool etamu = (fabs(muons_eta->at(0))<2.5 && fabs(muons_eta->at(1))<2.5 && fabs(muons_eta->at(2))<2.5 && fabs(muons_eta->at(3))<2.5);
			bool itrk  = ((itrk1>=0 && itrk1<nTracks) && (itrk2>=0 && itrk2<nTracks) && (itrk3>=0 && itrk3<nTracks) && (itrk4>=0 && itrk4<nTracks));
			bool mcp   = (itrk) ? (MCP(itrk1) && MCP(itrk2) && MCP(itrk3) && MCP(itrk4)) : false;
			
			if(ptmu && zmass && etamu && mcp)
			{
				//cout << "m4=" << p4.M() << ", m3=" << p3.M() << endl;
				
				histos[name+"_triplet_mQuad_muons"]->Fill(p4.M(),wgt);
				histos[name+"_triplet_mQuad_norm_muons"]->Fill(p4.M(),wgt);
				histos[name+"_triplet_mTrip_muons"]->Fill(p3.M(),wgt);
				histos[name+"_triplet_mTrip_norm_muons"]->Fill(p3.M(),wgt);
				histos[name+"_triplet_pTQuad_muons"]->Fill(p4.Pt(),wgt);
				histos[name+"_triplet_pTQuad1_muons"]->Fill(m1.Pt(),wgt);
				histos[name+"_triplet_pTQuad2_muons"]->Fill(m2.Pt(),wgt);
				histos[name+"_triplet_pTQuad3_muons"]->Fill(m3.Pt(),wgt);
				histos[name+"_triplet_pTQuad4_muons"]->Fill(m4.Pt(),wgt);
			}
		}
	}	
		_DEBUG("");
		
		vector<vertex> vertices;
		for(unsigned int vtx=0 ; vtx<nvtx ; vtx++)
		{
			//// validate the chain
			TString shorttype = classifyTripletShort(vtx);
			if(!validatedVertexChain(vtx)) continue;
			
			//// calculate the 3body mass and decide if to blind it
			double m3body = vtx_mass->at(vtx);
			bool blind = (isdata && (m3body>=mBlindMin && m3body<=mBlindMax));
			
			//// blind the signal region for the data
			if(doBlind && blind) continue;
			
			//// ignore vertices with calo triplets
			if(shorttype.Contains("calo")) continue;
			
			//// ignore vertices with segment-tagged-TPs triplets
			if(shorttype.Contains("tpmuB")) continue;
        
			//// ignore vertices with m3body>4 GeV
			// if(m3body>4.*GeV2MeV) continue;
			
			/////////////////////
			//// write the vertex
			vertex v;
			v.set(vtx);
			vertices.push_back(v);
		}
		
		////////////////////////////////////////////////////
		//// write out only in analysis mode ///////////////
		if(!skim) fillFlatoutTree(vertices,allPassing); ////
		////////////////////////////////////////////////////
		
		////////////////
		//// Trigger cut
		if(isCounter("nPassing_evt_trigger") && !acceptTrigger()) continue;
		incrementCounter("nPassing_evt_trigger",wgt);
		simulateTripletBuilder(name,histos,"Muons","CombinedFitMuonParticles","SegmentTagTrackParticles");
	
		/*	
		_DEBUG("");

		if( (triggerbits["EF_mu24i_tight"] || triggerbits["EF_mu36_tight"] || triggerbits["EF_mu24_tight_EFxe40"]) &&
		   !(triggerbits["EF_3mu4T"] || triggerbits["EF_3mu6"] || triggerbits["EF_3mu6_MSonly"] || triggerbits["EF_2mu13"] || triggerbits["EF_mu18_tight_mu8_EFFS"] || triggerbits["EF_mu18_tight_2mu4_EFFS"]))
		{
			cout << "mu24i("<< triggerbits["EF_mu24i_tight"] <<") / mu36("<< triggerbits["EF_mu36_tight"] <<") / mu24xe("<< triggerbits["EF_mu24_tight_EFxe40"] <<")" << endl;
			cout << "3mu4: " << triggerbits["EF_3mu4T"];
			cout << ", 3mu6: " << triggerbits["EF_3mu6"];
			cout << ", 3mu6_MSonly: " << triggerbits["EF_3mu6_MSonly"];
			cout << ", 2mu13: " << triggerbits["EF_2mu13"];
			cout << ", mu18_mu8: " << triggerbits["EF_mu18_tight_mu8_EFFS"];
			cout << ", mu18_2mu4: " << triggerbits["EF_mu18_tight_2mu4_EFFS"] << endl;
		}
		*/

		_DEBUG("");
		
		//////////////////////////////////////////////////////////////
		// actual analysis goes here /////////////////////////////////
		vector<unsigned int> acceptedtriplets; ///////////////////////
		vector<unsigned int> acceptedvtx; ////////////////////////////
		vector<unsigned int> selectedvtx; ////////////////////////////
		bool acceptSkim2 = false; ////////////////////////////////////
		//////////////////////////////////////////////////////////////
		
		
		unsigned int ngoodtriplets = 0;
		for(unsigned int vtx=0 ; vtx<nvtx ; vtx++)
		{
			if(!validatedVertex(vtx)) continue;
			ngoodtriplets++;
		}
		if(isCounter("nPassing_evt_goodtriplets") && ngoodtriplets<1) continue;
		incrementCounter("nPassing_evt_goodtriplets",wgt);
		
		
		vector<unsigned int> ivtx;
		for(unsigned int i=0 ; i<nvtx ; i++) ivtx.push_back(i);
		fillCategories(ivtx,name,"tripletCategories_afterVertexing",histos);
		fillCategories(ivtx,name,"tripletCategories_norm_afterVertexing",histos);
		fillCategoriesdRmin(ivtx,name,histos);
		
		
		for(unsigned int vtx=0 ; vtx<nvtx ; vtx++)
		{
			///////////////////////////////////////////
			//// validate the vertex chain ////////////
			//// here's the place to select which /////
			//// triplet/vertex categories to keep ////
			if(!validatedVertex(vtx)) continue; ///////
			///////////////////////////////////////////
			
			
			/////////////////////////////////////////////////////////////////////
			// blind the signal region for the data and only then proceed to ////
			// accept the event and fill the histograms /////////////////////////
			double mass = vtx_mass->at(vtx); ////////////////////////////////////
			bool blind = (isdata && (mass>=mBlindMin && mass<=mBlindMax)); //////
			if(doBlind && blind) continue; //////////////////////////////////////
			/////////////////////////////////////////////////////////////////////
			
			if(skim)
			{
				//sources src;
			        //getSrc(vtx,src);
				//int itrk1 = src.trkIndex[0];
			        //int itrk2 = src.trkIndex[1];
			        //int itrk3 = src.trkIndex[2];

				TLorentzVector pMuSum = getTlv3mu(vtx);
				double mt = mT(MET_RefFinal_et,MET_RefFinal_phi,pMuSum.Pt(),pMuSum.Phi());
				
				if(isCounter("nPassing_skim2_m3mu") && mass>4.5*GeV2MeV) continue;
				incrementCounter("nPassing_skim2_m3mu",wgt);

				if(isCounter("nPassing_skim2_pT3mu") && vtx_pt->at(vtx)<7.*GeV2MeV) continue;
				incrementCounter("nPassing_skim2_pT3mu",wgt);
			
				//if(isCounter("nPassing_skim2_mcp") && MCP(itrk1)+MCP(itrk2)+MCP(itrk3)!=3) continue;
				//incrementCounter("nPassing_skim2_mcp",wgt);
				
				if(isCounter("nPassing_skim2_met") && MET_RefFinal_et<7.*GeV2MeV) continue;
				if(skim) incrementCounter("nPassing_skim2_met",wgt);
				
				if(isCounter("nPassing_skim2_mt") && mt<7.*GeV2MeV) continue;
				if(skim) incrementCounter("nPassing_skim2_mt",wgt);
				
				/////////////////////////////////////////////////////
				// cuts of 2nd skim ends here !!!                ////
				// if you got here at least once in skim state,  ////
				// then you can do the skim. otherwise, don't... ////
				acceptSkim2 = true; /////////////////////////////////
				/////////////////////////////////////////////////////
			}
			
			//////////////////////////////////////////////////////////////////////////////////
			// accept triplet constituets ////////////////////////////////////////////////////
			if(!acceptMuons(vtx,name,histos,histos2,wgt,mBlindMin,mBlindMax)) continue; //////
			//////////////////////////////////////////////////////////////////////////////////
			
			//////////////////////////////////////////////////////////////////////////////////	
			// accept triplets ///////////////////////////////////////////////////////////////
			if(!acceptTriplet(vtx,name,histos,histos2,wgt,mBlindMin,mBlindMax)) continue; ////
			//////////////////////////////////////////////////////////////////////////////////
			
			//////////////////////////////////////////////////////////////////////////////////
			// accept hadronic clean /////////////////////////////////////////////////////////
			if(!acceptHadClean(vtx,name,histos,histos2,wgt,mBlindMin,mBlindMax)) continue; ///
			//////////////////////////////////////////////////////////////////////////////////
			
			//////////////////////////////////////
			acceptedtriplets.push_back(vtx); /////
			//////////////////////////////////////
		}
		
		_DEBUG("");
		
		if(!skim && !doMVA)
		{
			///////////////////////////////////////////////////////////////////////////////////////
			if(isCounter("nPassing_evt_onetriplet") && acceptedtriplets.size()!=1) continue; //////
			incrementCounter("nPassing_evt_onetriplet",wgt); //////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////// 
		}

		_DEBUG("");

		
		////////////////////////////////////////////////////
		//// now write out /////////////////////////////////
		if(skim && acceptSkim2) fillOutTrees(otrees); //////
		////////////////////////////////////////////////////
		
		_DEBUG("");
	
		////////////////////////////////////////////////////////////////
		/// count and store the number of accepted muon triplets ///////
		unsigned int nacceptedtriplets = acceptedtriplets.size(); //////
		////////////////////////////////////////////////////////////////
		
		_DEBUG("");

		
		////////////////////////////////////////////////////
		/// now go on to check the vertex properties ///////
		for(unsigned int muvtx=0 ; muvtx<nacceptedtriplets ; muvtx++)
		{
			// get the vertex index
			unsigned int vtx = acceptedtriplets[muvtx];
			
			_DEBUG("");
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			// accept (vertex-MET cut-based or MVA-based) ////////////////////////////////////////////////////////////
			if(!acceptVtxMET(selectionMethod,vtx,vertices,name,histos,histos2,wgt,mBlindMin,mBlindMax)) continue; ////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////

			_DEBUG("");

			////////////////////////////////
			// book-keeping ////////////////
			acceptedvtx.push_back(vtx); ////
			////////////////////////////////
		}
		
		_DEBUG("");
		
		////////////////////////////////////////////////////////////////////////////////////////
		// move to the next event if there's no acepted vertex /////////////////////////////////
		if(acceptedvtx.size()<1) continue; /////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////

		
		_DEBUG("");

	
		// loop on the accepted triplets
		unsigned int nacceptedvtx = acceptedvtx.size();
		for(unsigned int avtx=0 ; avtx<nacceptedvtx ; avtx++)
		{
			/////////////////////////////////////////
			unsigned int vtx = acceptedvtx[avtx]; ///
			/////////////////////////////////////////
			
			// double mass = vtx_mass->at(vtx);
			TLorentzVector pMuSum = getTlv3mu(vtx);
			
			_DEBUG("");
			
			//// b-tagging analysis (before b-tagging cuts !!!)
			unsigned int iBtaggingTriplet = selectedvtx.size()+1;
			/////////////////////////////////////////////////////////////////////
			// b-tag histos should be filled only once, when iBtaggingTriplet = 1.
			// this is not the best requirement, because now the histos are filled
			// only for the 1st good vertex that went up to the b-tagging call so
			// other 3mu systems in the same evt are ignored at the moment.
			// However, this is a very small effect... 
			/////////////////////////////////////////////////////////////////////
			bTagging(name,histos,histos2,wgt,pMuSum,iBtaggingTriplet); //////////
			/////////////////////////////////////////////////////////////////////
			histos[name+"_mu_m3mu_beforeBtag"]->Fill(pMuSum.M(),wgt);
			histos[name+"_mu_m3mu_beforeBtag_norm"]->Fill(pMuSum.M(),wgt);
			int btagj1=-1;
			int btagj2=-1;
			bool isBtag = isBtaggedEvent(btagj1,btagj2);
			if(isBtag)
			{
				histos[name+"_mu_m3mu_afterBtag"]->Fill(pMuSum.M(),wgt);
				histos[name+"_mu_m3mu_afterBtag_norm"]->Fill(pMuSum.M(),wgt);
			}
			
			_DEBUG("");
			
			//////////////////////////////////////////////////
			/////// MORE SELECTIONS SHOULD COME HERE /////////
			//////////////////////////////////////////////////
					
			_DEBUG("");	
			
			//////////////////////////////////////////////
			///// this is the bottom of the cutflow //////
			///// everything below is after all cuts /////
			//////////////////////////////////////////////	
			selectedvtx.push_back(vtx); // important !!!!!
			//////////////////////////////////////////////
		}
		
		_DEBUG("");
		
		vector<unsigned int> iselectedvtx;
		for(unsigned int i=0 ; i<selectedvtx.size() ; i++) iselectedvtx.push_back(selectedvtx[i]);
		fillCategories(iselectedvtx,name,"tripletCategories_endOfSelection",histos);
		fillCategories(iselectedvtx,name,"tripletCategories_norm_endOfSelection",histos);
		
		_DEBUG("");
		
		// post single vtx
		unsigned int nTriplets = selectedvtx.size();
		histos[name+"_nTriplets"]->Fill(nTriplets,wgt);
		
		_DEBUG("");
	}
	
	_DEBUG("");
	
	// get total events processed:
	int nevts = getCounter("nPassing_evt_all");
	nevents[name] = (name==bname) ? nevts : nevents[name]+nevts;
	
	// fill cutflow histo
	fillCutFlowHisto(name,histos);
	fillCountersHisto(name,histos);
	
	// finalize
	string scounters = "";
	scounters = printCounters("cutflow",pname);
	writeCoutners(ftxtname,scounters);
	if(!isdata)
	{
		scounters = printCounters("cutflow",pname,true); // weighted cutflow
		writeCoutners(ftxtname,scounters);
	}
	scounters = printCounters("objects",pname);
	writeCoutners(ftxtname,scounters);
	writeCoutners(ftxtname,"\n\n"); // add 2 blanc lines at the end

	if(name=="Wtaunu_3mu") finalizePileup();
}

// ///////////////////////////////////////////
// ///////////////////////////////////////////
// ///////////////////////////////////////////
// TCanvas* makeEffCnv(TString name, TString title, TString titlex, TH1* htruacc, TH1* htrumatmu, TH1* htrumatcal, TH1* htrumatnmumcal)
// {
// 	TH1D* hTruAcc        = (TH1D*)htruacc->Clone();
// 	TH1D* hTruMatMu      = (TH1D*)htrumatmu->Clone();
// 	TH1D* hTruMatCal     = (TH1D*)htrumatcal->Clone();
// 	TH1D* hTruMatNmuMcal = (TH1D*)htrumatnmumcal->Clone();
// 	
// 	Double_t xmin = hTruAcc->GetXaxis()->GetXmin();
// 	Double_t xmax = hTruAcc->GetXaxis()->GetXmax();
// 
// 	TCanvas* cnv = new TCanvas("cnv_"+name,name,1200,400);
// 	cnv->Divide(2,1);
// 	
// 	cnv->cd(1);
// 	hTruAcc->SetLineColor(kGreen);      hTruAcc->SetMarkerColor(kGreen);      hTruAcc->DrawCopy();
// 	hTruMatMu->SetLineColor(kBlack);    hTruMatMu->SetMarkerColor(kBlack);    hTruMatMu->DrawCopy("same");
// 	hTruMatCal->SetLineColor(kCyan);    hTruMatCal->SetMarkerColor(kCyan);    hTruMatCal->DrawCopy("same");
// 	hTruMatNmuMcal->SetLineColor(kRed); hTruMatNmuMcal->SetMarkerColor(kRed); hTruMatNmuMcal->DrawCopy("same");
// 	// hTruMatNmuMcal->SetLineColor(kOrange); hTruMatNmuMcal->SetMarkerColor(kOrange); hTruMatNmuMcal->DrawCopy("same");
// 	
// 	cnv->cd(2);
// 	TGraphAsymmErrors* gTruMatCal = getDivided((TH1D*)hTruMatCal->Clone(),(TH1D*)hTruAcc->Clone(),"Reconstruction efficiency "+title,titlex,"Efficiency");
// 	gTruMatCal->SetLineColor(kCyan); gTruMatCal->SetMarkerColor(kCyan); gTruMatCal->GetXaxis()->SetLimits(xmin,xmax);
// 	TGraphAsymmErrors* gTruMatMu  = getDivided((TH1D*)hTruMatMu->Clone(),(TH1D*)hTruAcc->Clone(),"Reconstruction efficiency "+title,titlex,"Efficiency");
// 	gTruMatMu->SetLineColor(kBlack); gTruMatMu->SetMarkerColor(kBlack); gTruMatMu->GetXaxis()->SetLimits(xmin,xmax);
// 	TGraphAsymmErrors* gTruMatAll = getDivided(hTruMatNmuMcal,(TH1D*)hTruAcc->Clone(),"Reconstruction efficiency "+title,titlex,"Efficiency");
// 	gTruMatAll->SetLineColor(kRed);  gTruMatAll->SetMarkerColor(kRed);  gTruMatAll->GetXaxis()->SetLimits(xmin,xmax);
// 	
// 	TMultiGraph *mgr = new TMultiGraph("gr","Reconstruction efficiency "+title+";"+titlex+";Efficiency");
// 	mgr->Add(gTruMatMu);
// 	mgr->Add(gTruMatCal);
// 	mgr->Add(gTruMatAll);
// 	mgr->SetMaximum(1.1);
// 	mgr->SetMinimum(0.0);
// 	// mgr->GetXaxis()->SetLimits(xmin,xmax);
// 	TLegend* leg = new TLegend(0.55,0.12,0.83,0.35,NULL,"brNDC");
// 	leg->SetFillStyle(4000); //will be transparent
// 	leg->SetFillColor(0);
// 	leg->SetTextFont(42);
// 	leg->SetBorderSize(0);
// 	leg->AddEntry(gTruMatAll,"All","p");
// 	leg->AddEntry(gTruMatMu,"Muons","p");
// 	leg->AddEntry(gTruMatCal,"Calo muons","p");
// 	mgr->Draw("AP1");
// 	leg->Draw("same");
// 	saveEffCnv(name,mgr,leg);
// 
// 	return cnv;
// }
// TCanvas* makeEffCnv(TString name, TString title, TString titlex, TH1* href, TH1* hpass, TString type="Efficiency")
// {
// 	TH1D* hRef  = (TH1D*)href->Clone();
// 	TH1D* hPass = (TH1D*)hpass->Clone();
// 
// 	Double_t xmin = hRef->GetXaxis()->GetXmin();
// 	Double_t xmax = hRef->GetXaxis()->GetXmax();
// 
// 	TCanvas* cnv = new TCanvas("cnv_"+name,name,1200,400);
// 	cnv->Divide(2,1);
// 	
// 	cnv->cd(1);
// 	hRef->SetLineColor(kGreen);  hRef->SetMarkerColor(kGreen);  hRef->DrawCopy();
// 	hPass->SetLineColor(kBlack); hPass->SetMarkerColor(kBlack); hPass->DrawCopy("same");
// 	
// 	cnv->cd(2);
// 	TGraphAsymmErrors* gEff = getDivided((TH1D*)hPass->Clone(),(TH1D*)hRef->Clone(),type+" "+title,titlex,"Efficiency");
// 	gEff->SetLineColor(kBlack); gEff->SetMarkerColor(kBlack); gEff->GetXaxis()->SetLimits(xmin,xmax);
// 	
// 	TMultiGraph *mgr = new TMultiGraph("gr",type+" "+title+";"+titlex+";Efficiency");
// 	mgr->Add(gEff);
// 	mgr->SetMaximum(1.1);
// 	mgr->SetMinimum(0.0);
// 	// mgr->GetXaxis()->SetLimits(xmin,xmax);
// 	TLegend* leg = new TLegend(0.55,0.12,0.83,0.35,NULL,"brNDC");
// 	leg->SetFillStyle(4000); //will be transparent
// 	leg->SetFillColor(0);
// 	leg->SetTextFont(42);
// 	leg->SetBorderSize(0);
// 	leg->AddEntry(gEff,type,"p");
// 	mgr->Draw("AP1");
// 	leg->Draw("same");
// 	saveEffCnv(name,mgr,leg);
// 	return cnv;
// }

void makeSingleCnv(bool is2d=false)
{
	c1 = (is2d) ? new TCanvas("","",400,400) : new TCanvas("","",600,400);
}
void closeSingleCnv(TString varname, bool is2d=false, TString channel="")
{
	c1->SetTicks(1,1);
	c1->RedrawAxis();
	c1->Update();
	TString epsfilename = pdffilename;
	TString newsuffix = (channel=="") ? "."+varname+".eps"  :  "."+varname+"."+channel+".eps";
	epsfilename.ReplaceAll(".pdf",newsuffix);
	c1->SaveAs(epsfilename);
	if(!is2d)
	{
		// TH1
		epsfilename.ReplaceAll(".eps",".logy.eps");
		c1->SetLogy();
		c1->SaveAs(epsfilename);
	}
	else
	{
		// TH2
		epsfilename.ReplaceAll(".eps",".logz.eps");
		c1->SetLogz();
		c1->SaveAs(epsfilename);
	}
	delete c1;
}


void drawData(TString varname, TMapTSP2TH1& histos, bool isfirst, TString drawopt="p", Bool_t doymin=0, Bool_t doymax=0, Double_t ymin=-1, Double_t ymax=-1)
{	
	TMapTSP2TH1::iterator it=histos.begin();
	for(; it!=histos.end() ; it++)
	{
		TString hname = it->first;
		if(!isData(hname))           continue;
		if(!hname.EndsWith(varname)) continue;
		break; // if got to here, break
	}
	it->second->SetLineColor(kBlack);
	it->second->SetLineWidth(2);
	it->second->SetLineStyle(1);
	it->second->SetMarkerStyle(20);
	it->second->SetMarkerSize(1);
	it->second->SetMarkerColor(kBlack);
	
	double min = getYmin(it->second);
	double max = getYmax(it->second);
	if(min==max && min==0.)
	{
		min = 1e-2;
		max = 1.;
	}
	
	if(doymin) it->second->SetMinimum(min*0.9);
	if(doymax) it->second->SetMaximum(max*2.5);
	if(ymin!=-1) it->second->SetMinimum(ymin);
	if(ymax!=-1) it->second->SetMaximum(ymax);
	
	if(isfirst) it->second->Draw(drawopt);
	else        it->second->Draw(drawopt+" same");
}
void drawData(TString varname, TMapTSP2TH2& histos, bool isfirst, TString drawopt)
{	
	TMapTSP2TH2::iterator it=histos.begin();
	for(; it!=histos.end() ; it++)
	{
		TString hname = it->first;
		if(!isData(hname))           continue;
		if(!hname.EndsWith(varname)) continue;
		break; // if got to here, break
	}
	it->second->SetLineColor(kBlack);
	it->second->SetFillColor(kBlack);
	it->second->SetMarkerColor(kBlack);
	it->second->SetLineStyle(1);
	
	if(drawopt=="")
	{
		it->second->SetContour(50);
		if(isfirst) it->second->Draw("col");
		else        it->second->Draw("col same");
		exeGray->Draw();
		it->second->Draw("col same");
	}
	else
	{
		if(isfirst) it->second->Draw(drawopt);
		else        it->second->Draw(drawopt+" same");
	}
}

void drawPadSingle(TString channel, TString varname, TVirtualPad* pad, TMapuiTS& channels, TMapTSP2TH1& histos, TString drawopt="hist")
{
	makeSingleCnv();
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		if(isData(cit->second))  continue;
		if(cit->second!=channel) continue; // draw only one channel
		
		TString name = cit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		double min = getYmin(histos[hname]);
		double max = getYmax(histos[hname]);
		
		if(min==max && min==0.)
		{
			min = 1e-2;
			max = 1.;
		}
		
		histos[hname]->SetMinimum(min*0.9);
		histos[hname]->SetMaximum(max*2.5);
		
		c1->cd();
		histos[hname]->Draw(drawopt);
		pad->cd();
		histos[hname]->Draw(drawopt);
	}
	// the data
	c1->cd();
	if(isData(channel)) drawData(varname,histos,false,drawopt,true,true);
	pad->cd();
	if(isData(channel)) drawData(varname,histos,false,drawopt,true,true);
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname,false,channel);
}
void drawPadSingle(TString channel, TString varname, TVirtualPad* pad, TMapuiTS& channels, TMapTSP2TH2& histos, TString drawopt="box")
{
	makeSingleCnv();
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	float max = -1.e20;
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		TString name = cit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		max = (histos[hname]->GetMaximum()>max) ? histos[hname]->GetMaximum() : max;
	}
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		if(isData(cit->second))  continue;
		if(cit->second!=channel) continue; // draw only one channel
		
		TString name = cit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		c1->cd();
		histos[hname]->Draw(drawopt);
		pad->cd();
		histos[hname]->Draw(drawopt);
	}
	// the data
	if(isData(channel))
	{
		c1->cd();
		drawData(varname,histos,true,drawopt);
		pad->cd();
		drawData(varname,histos,true,drawopt);
	}
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname,true,channel);
}

void drawPadAll(TString varname, TVirtualPad* pad, TMapuiTS& channels, TMapTSP2TH1& histos, TLegend* leg, Double_t ymin=-1, Double_t ymax=-1)
{
	makeSingleCnv();
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	unsigned int counter = 0;
	double max = -1.e20;
	double min = +1.e20;
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		TString name = cit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if(isBGsum(name))  continue;
		if(isSIGsum(name)) continue;
		
		if(isWsignal(name) || isData(name))
		{
			Double_t hmin = getYmin(histos[hname]);
			Double_t hmax = getYmax(histos[hname]);
			max = (hmax>max) ? hmax : max;
			min = (hmin<min) ? hmin : min;
		}
	}
	
	if(ymin==-1) ymin = min*0.9;
	if(ymax==-1) ymax = max*2.5;
	if(min==max && min==0.)
	{
		min = 1e-2;
		max = 1.;
	}
	
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		if(isData(cit->second)) continue;
		
		TString name = cit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if(isData(name))   continue;
		if(isBGsum(name))  continue;
		if(isSIGsum(name)) continue;
				
		histos[hname]->SetMinimum(ymin);
		histos[hname]->SetMaximum(ymax);
		
		c1->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("hist same");
		
		pad->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("hist same");
		
		counter++;
	}
	// the data
	c1->cd();
	drawData(varname,histos,false,"p",false,false,ymin,ymax);
	leg->Draw("same");
	pad->cd();
	drawData(varname,histos,false,"p",false,false,ymin,ymax);
	leg->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname);
}
void drawPadAll(TString varname, TVirtualPad* pad, TMapuiTS& channels, TMapTSP2TH2& histos, TLegend* leg, TString drawopt)
{
	makeSingleCnv();
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	
	if(!leg) _INFO("leg is NULL");
	
	// the MSs
	TMapuiTS::reverse_iterator rcit=channels.rbegin();
	TString name = rcit->second;
	TString hname = name+"_"+varname;
	TH2* hbkg     = NULL;
	TH2* hsignal  = NULL;
	TH2* hWsignal = NULL;
	TH2* hdata    = NULL;
	for( ; rcit!=channels.rend() ; rcit++)
	{
		name = rcit->second;
		if(drawchannels[name]==0) continue;
		
		hname = name+"_"+varname;
		
		if(!isBGsum(name) && !isSIGsum(name) && !isData(name)) continue;
		
		if(isData(name))    hdata    = histos[hname];
		if(isBGsum(name))   hbkg     = histos[hname];
		if(isSIGsum(name))  hsignal  = histos[hname];
		if(isWsignal(name)) hWsignal = histos[hname];
	}
	
	double max = -1.e20;
	max = (hdata->GetMaximum()>max)    ? hdata->GetMaximum()    : max;
	max = (hWsignal->GetMaximum()>max) ? hWsignal->GetMaximum() : max;
	// max = (hsignal->GetMaximum()>max) ? hsignal->GetMaximum() : max;
	// max = (hbkg->GetMaximum()>max)    ? hbkg->GetMaximum()    : max;
	
	hdata->SetMaximum(max*2.5);
	hsignal->SetMaximum(max*2.5);
	hbkg->SetMaximum(max*2.5);
	
	// the data
	drawData(varname,histos,true,drawopt);
	
	if(drawopt=="")
	{
		c1->cd();
		hbkg->SetContour(50);
		hbkg->Draw("col same");
		exeRed->Draw();
		hbkg->Draw("col same");
		pad->cd();
		hbkg->SetContour(50);
		hbkg->Draw("col same");
		exeRed->Draw();
		hbkg->Draw("col same");
		
		c1->cd();
		hsignal->SetContour(50);
		hsignal->Draw("col same");
		exeGreen->Draw();
		hsignal->Draw("col same");
		pad->cd();
		hsignal->SetContour(50);
		hsignal->Draw("col same");
		exeGreen->Draw();
		hsignal->Draw("col same");
	}
	else
	{
		c1->cd();
		hbkg->Draw(drawopt+" same");
		hsignal->Draw(drawopt+" same");
		pad->cd();
		hbkg->Draw(drawopt+" same");
		hsignal->Draw(drawopt+" same");
	}
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname,true);
}
void drawPadStack(TString varname, TVirtualPad* pad, TMapuiTS& channels, TMapTSP2TH1& histos, TLegend* leg, Double_t ymin=-1, Double_t ymax=-1)
{
	makeSingleCnv();
	pad->cd();
	
	if(varname=="")
	{
		pad->Update();
		pad->cd();
		return;
	}
	
	TMapuiTS::reverse_iterator rit=channels.rbegin();
	TH1* hsumBkg = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TH1* hWsig   = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TH1* hsumSig = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TH1* hData   = (TH1*)histos[rit->second+"_"+varname]->Clone();
	TString t  = histos[rit->second+"_"+varname]->GetTitle();
	TString tx = histos[rit->second+"_"+varname]->GetXaxis()->GetTitle();
	TString ty = histos[rit->second+"_"+varname]->GetYaxis()->GetTitle();
	TString titles = t+";"+tx+";"+ty;
	THStack* hsBkg = new THStack(varname+"_bkg",titles);
	THStack* hsSig = new THStack(varname+"_sig",titles);
	hsumBkg->Reset();
	hsumSig->Reset();
	hData->Reset();
	hWsig->Reset();
	
	for(rit=channels.rbegin() ; rit!=channels.rend() ; rit++)
	{
		TString name = rit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if(isData(rit->second))   { hData->Add(histos[hname]); continue; }
		if(isBGsum(rit->second))  continue;
		if(isSIGsum(rit->second)) continue;
		if(isSignal(rit->second)) continue;
		
		hsBkg->Add(histos[hname],"hist");
		hsumBkg->Add(histos[hname]);
	}
	for(rit=channels.rbegin() ; rit!=channels.rend() ; rit++)
	{
		TString name = rit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname = name+"_"+varname;
		
		if(!isSignal(name)) continue;
		if(isWsignal(name)) hWsig->Add(histos[hname]);
		
		hsSig->Add(histos[hname],"hist");
		hsumSig->Add(histos[hname]);
	}
	
	// Double_t maxBkg = getYmax(hsumBkg);
	// Double_t maxSig = getYmax(hsumSig);
	Double_t maxSig = getYmax(hWsig);
	Double_t maxDat = getYmax(hData);
	// Double_t max = (maxBkg>maxSig) ? maxBkg : maxSig;
	// max = (maxDat>max) ? maxDat : max;
	Double_t max = (maxDat>maxSig) ? maxDat : maxSig;
	
	TH1* htmp;
	TIter nextBhist((TList*)hsBkg->GetHists());
	while( (htmp=(TH1*)nextBhist())!=NULL )
	{
		if(ymin<0.) htmp->SetMinimum(0.);
		else        htmp->SetMinimum(ymin);
		if(ymax<0.) htmp->SetMaximum(max*2.5);
		else        htmp->SetMaximum(ymax);
	}
	TIter nextShist((TList*)hsSig->GetHists());
	while( (htmp=(TH1*)nextShist())!=NULL )
	{
		if(ymin<0.) htmp->SetMinimum(0.);
		else        htmp->SetMinimum(ymin);
		if(ymax<0.) htmp->SetMaximum(max*2.5);
		else        htmp->SetMaximum(ymax);
	}
	
	if(ymin<0.) hsBkg->SetMinimum(0.);
	else        hsBkg->SetMinimum(ymin);
	if(ymax<0.) hsBkg->SetMaximum(max*2.5);
	else        hsBkg->SetMaximum(ymax);
	if(ymin<0.) hsSig->SetMinimum(0.);
	else        hsSig->SetMinimum(ymin);
	if(ymax<0.) hsSig->SetMaximum(max*2.5);
	else        hsSig->SetMaximum(ymax);
	
	// the data
	c1->cd();
	hsBkg->Draw();
	hsSig->Draw("hist same");
	drawData(varname,histos,false);
	leg->Draw("same");
	pad->cd();
	hsBkg->Draw();
	hsSig->Draw("hist same");
	drawData(varname,histos,false);
	leg->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname);
}
void drawPadSignal(TString varname, TVirtualPad* pad, TMapuiTS& channels, TMapTSP2TH1& histos, TLegend* leg)
{
	makeSingleCnv();
	pad->cd();
	
	int counter = 0;
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		TString name = cit->second;
		if(drawchannels[name]==0) continue;
		if(!isSignal(name)) continue;
		TString hname = name+"_"+varname;
		histos[hname]->SetStats(1);
		
		c1->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("same hist");
		
		pad->cd();
		if(counter==0) histos[hname]->Draw("hist");
		else           histos[hname]->Draw("same hist");
		
		counter++;
	}
	c1->cd();
	leg->Draw("same");
	pad->cd();
	leg->Draw("same");
	
	pad->RedrawAxis();
	pad->Update();
	closeSingleCnv(varname);
}
void divideCuts(TString varname_rat, TString varname_num, TString varname_den, TMapuiTS& channels, TMapTSP2TH1& histos)
{
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		TString name         = cit->second;
		if(drawchannels[name]==0) continue;
		
		TString hname_den = name+"_"+varname_den;
		TString hname_num = name+"_"+varname_num;
		TString hname_rat = name+"_"+varname_rat;
		histos[hname_rat]->Divide(histos[hname_num],histos[hname_den]);
	}
}
int getsizex(int ndivx)
{
	if     (ndivx==1) return 600;
	else if(ndivx==2) return 1200;
	else if(ndivx==3) return 1800;
	else _FATAL("ndivx=4 is not supported");
	return 1000;
}
int getsizey(int ndivy)
{
	if     (ndivy==1) return 400;
	else if(ndivy==2) return 800;
	else if(ndivy==3) return 1200;
	else if(ndivy==4) return 1600;
	else _FATAL("ndivy=5 is not supported");
	return 800;
}
void makeCnv(int ndivx, int ndivy, bool logy=false, TString parity="")
{
	cnv = new TCanvas("c","c",getsizex(ndivx),getsizey(ndivy));
	cnv->Divide(divx,divy);	
	pads.clear();
	for(int i=0 ; i<ndivx*ndivy ; i++)
	{
		pads.push_back(cnv->cd(i+1));
		pads[i]->SetTicks(1,1);
		if(logy)
		{
			if(parity=="" || (parity!="even" && parity!="odd")) pads[i]->SetLogy();
			else
			{
				if     (parity=="even" && i%2==0) pads[i]->SetLogy();
				else if(parity=="odd"  && i%2!=0) pads[i]->SetLogy();
			}
		}
	}
	cnv->Draw();
}
void closeCnv(TString fname)
{
	cnv->Update();
	cnv->SaveAs(fname);
	delete cnv;
}


/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

void NTUPmaker(TString runType, TString outDir, TString chnl, TString master, TString method, TString split="")
{
	msglvl[DBG] = SILENT; // SILENT;
	msglvl[INF] = VISUAL;
	msglvl[WRN] = VISUAL;
	msglvl[ERR] = VISUAL;
	msglvl[FAT] = VISUAL;
	
	cout << "-------------------- arguments --------------------" << endl;
	cout << "runType  = " << runType << endl;
	cout << "outDir   = " << outDir  << endl;
	cout << "chnl     = " << chnl    << endl;
	cout << "master   = " << master  << endl;
	cout << "method   = " << method  << endl;
	cout << "split    = " << split   << endl;
	cout << "---------------------------------------------------" << endl;

	// atlstyle();
	// style(false,false);
	// minimalstyle();
	// gStyle->SetErrorX(0);
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// flag for skim ///////////////////////////////////////////////////////////////////////////////////
	if(runType!="skim" && runType!="analysis") _FATAL("unknown runType="+(string)runType); /////////////
	skim = (runType=="skim"); //////////////////////////////////////////////////////////////////////////
	if(skim && !basepath.Contains("root://eosatlas//eos/")) _FATAL("skim can be ran from EOS only"); ///
	if(!outDir.EndsWith("/")) outDir+="/"; /////////////////////////////////////////////////////////////
	if(!skim) basepath = outDir+"data/bphys/"; // in runType=="analysis" mode //////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	////////////////////////////////////////////////
	//// GLOBALS !!! ///////////////////////////////
	mastername = master; ///////////////////////////
	selectionMethod = method; //////////////////////
	doMVA = (method.Contains("MVA:") && !skim); ////
	smethod = (doMVA) ? "mva" : "cuts"; ////////////
	////////////////////////////////////////////////
	
	//////////////////
	buildTPmus(); ////
	//////////////////
	
	
	TPMERegexp re(":");
	splitstr = split;
	TString ssplit = split;
	if(splitstr!="")
	{
		re.Split(splitstr);
		nChunks = re[0].Atoi(); // there are N chunks (nChunks), where each core takes one chunk
		iChunk  = re[1].Atoi(); // now you are processing the i'th chunk (iChunk)
		nChunksMax = 50;
		
		cout << "--------------------- splitter --------------------" << endl;
		re.Print("all");
		cout << "nChunks = " << nChunks << endl;
		cout << "iChunk  = " << iChunk << endl;
		cout << "---------------------------------------------------" << endl;
		
		if(nChunks==0 && iChunk==0) ssplit = "n0.j0";
		else
		{
			Bool_t isValid = ((nChunks>0 && nChunks<=nChunksMax) && (iChunk>0 && iChunk<=nChunks));
			if(!isValid) _FATAL("split argument: "+(string)split+" is invalid");
		
			if(iChunk<10) ssplit.ReplaceAll(":",".j0");
	 		else          ssplit.ReplaceAll(":",".j");
			ssplit = "n"+ssplit;
		}
	}
	else
	{
		nChunks=0;
		iChunk=0;
		nChunksMax = 50;
		
		cout << "--------------------- splitter --------------------" << endl;
		re.Print("all");
		cout << "nChunks -> set to " << nChunks << endl;
		cout << "iChunk  -> set to " << iChunk << endl;
		cout << "---------------------------------------------------" << endl;
		
		ssplit = "n0.j0";
	}
	
	
	////////////////////////////
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
	////////////////////////////
	

	TMapuiTS       channels;
	TMapTSTS       labels;
	TMapTSTS       legoptions;
	TMapTSTC       colors;
	TMapTSI        patterns;
	// TMapTSI        drawchannels;
	TMapTSTS       types;
	TMapTSP2TH1    histos;
	TMapTSP2TH2    histos2;
	TMapTSP2TF     files;
	TMapTSP2TCHAIN chains;
	TMapTSP2TCHAIN chainfriends;
	TMapTSP2TTREE  otrees;
	
	//////////////////////////////////////////
	// To limit data for max entries /////////
	const Long64_t nentriesMax = 0; // 10000;
	//////////////////////////////////////////
	
	
	///////////////////////////////
	// text file for cutflow //////
	///////////////////////////////
	timestr  = getTime();
	ftxtname = outDir+"cutflow."+chnl+"."+master+"."+smethod+"."+runType+"."+ssplit+".txt";
	writeCoutners(ftxtname,timestr);
	
	
	///////////////////////////////
	// root file for out trees ////
	///////////////////////////////
	foutname = outDir+"out."+chnl+"."+master+"."+runType+"."+ssplit+".root"; // the method is irrelevant for the skim
	TFile*  fout = (skim) ? new TFile(foutname,"RECREATE") : NULL;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// if one channel... ///////////////////////////////////////////////////////////////////////
	bool isAllData    = (chnl.Contains("periodall")); //////////////////////////////////////////
	bool isAllMC      = (chnl.Contains("MC")); /////////////////////////////////////////////////
	bool isSignalMC   = (chnl.Contains("_3mu")); ///////////////////////////////////////////////
	bool isWSignalMC  = (chnl=="Wtaunu_3mu"); //////////////////////////////////////////////////
	bool isBkgroundMC = (!isAllMC && (!isSignalMC && !chnl.Contains("period"))); ///////////////
	cout << "-------------------------- 1 channel flags ---------------------------" << endl;
	cout << "is all data:      " << isAllData << endl;
	cout << "is all mc:        " << isAllMC << endl;
	cout << "is signal mc:     " << isSignalMC << endl;
	cout << "is W signal mc:   " << isWSignalMC << endl;
	cout << "is background mc: " << isBkgroundMC << endl;
	cout << "----------------------------------------------------------------------" << endl;
	////////////////////////////////////////////////////////////////////////////////////////////
	

	//////////////////////////////////////////////////////////////////	
	//////////////////////////////////////////////////////////////////	
	// LUMINOSTIY SETUP //////////////////////////////////////////////	
	//////////////////////////////////////////////////////////////////	
	setLumis(); // fixed database ////////////////////////////////////
	if(isAllData || chnl.Contains("periodA")) enablePeriod("A");
	if(isAllData || chnl.Contains("periodB")) enablePeriod("B");
	if(isAllData || chnl.Contains("periodC")) enablePeriod("C");
	if(isAllData || chnl.Contains("periodD")) enablePeriod("D");
	if(isAllData || chnl.Contains("periodE")) enablePeriod("E");
	if(isAllData || chnl.Contains("periodG")) enablePeriod("G");
	if(isAllData || chnl.Contains("periodH")) enablePeriod("H");
	if(isAllData || chnl.Contains("periodI")) enablePeriod("I");
	if(isAllData || chnl.Contains("periodJ")) enablePeriod("J");
	if(isAllData || chnl.Contains("periodL")) enablePeriod("L");
	// if(isAllData || chnl.Contains("periodM")) enablePeriod("M");
	totalLumi = (chnl.Contains("period")) ? getTotalLumi() : 1.; // inverse femtobarns
	TString periods = "";
	for(TMapTSb::iterator it=periodenable.begin() ; it!=periodenable.end() ; ++it) { periods += it->first; }
	_INFO("Total integrated luminosity(periods="+(string)periods+") is: "+_s(totalLumi)+" ifb");
	TString slumi = (TString)(_s(totalLumi,2)+" fb^{-1}");
	//////////////////////////////////////////////////////////////////
	
	
	
	// //////////////////////////////////////////////////////////////////////
	// // Binned MC SETUP ///////////////////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////
	// if(isAllMC || chnl.Contains("WtaunuNp0")) enableBinnedMC("WtaunuNp0");
	// if(isAllMC || chnl.Contains("WtaunuNp1")) enableBinnedMC("WtaunuNp1");
	// if(isAllMC || chnl.Contains("WtaunuNp2")) enableBinnedMC("WtaunuNp2");
	// if(isAllMC || chnl.Contains("WtaunuNp3")) enableBinnedMC("WtaunuNp3");
	// if(isAllMC || chnl.Contains("WtaunuNp4")) enableBinnedMC("WtaunuNp4");
	// if(isAllMC || chnl.Contains("WtaunuNp5")) enableBinnedMC("WtaunuNp5");
	// //////////////////////////////////////////////////////////////////////
	// if(isAllMC || chnl.Contains("WmunuNp0")) enableBinnedMC("WmunuNp0");
	// if(isAllMC || chnl.Contains("WmunuNp1")) enableBinnedMC("WmunuNp1");
	// if(isAllMC || chnl.Contains("WmunuNp2")) enableBinnedMC("WmunuNp2");
	// if(isAllMC || chnl.Contains("WmunuNp3")) enableBinnedMC("WmunuNp3");
	// if(isAllMC || chnl.Contains("WmunuNp4")) enableBinnedMC("WmunuNp4");
	// if(isAllMC || chnl.Contains("WmunuNp5")) enableBinnedMC("WmunuNp5");
	// //////////////////////////////////////////////////////////////////////
	// if(isAllMC || chnl.Contains("ZmumuNp0")) enableBinnedMC("ZmumuNp0");
	// if(isAllMC || chnl.Contains("ZmumuNp1")) enableBinnedMC("ZmumuNp1");
	// if(isAllMC || chnl.Contains("ZmumuNp2")) enableBinnedMC("ZmumuNp2");
	// if(isAllMC || chnl.Contains("ZmumuNp3")) enableBinnedMC("ZmumuNp3");
	// if(isAllMC || chnl.Contains("ZmumuNp4")) enableBinnedMC("ZmumuNp4");
	// if(isAllMC || chnl.Contains("ZmumuNp5")) enableBinnedMC("ZmumuNp5");
	// //////////////////////////////////////////////////////////////////////
	// if(isAllMC || chnl.Contains("JZ0W")) enableBinnedMC("JZ0W");
	// if(isAllMC || chnl.Contains("JZ1W")) enableBinnedMC("JZ1W");
	if(isAllMC || chnl.Contains("JZ2W")) enableBinnedMC("JZ2W");
	if(isAllMC || chnl.Contains("JZ3W")) enableBinnedMC("JZ3W");
	//////////////////////////////////////////////////////////////////////
	
	
	
	//////////////////////////////////////////////////////////
	// this is the most important map                     ////
	// it controls the order of insertion into the stack, ////
	// the order of appearance in the legend and          ////
	// also the inclusion or removal of certain channels  ////
	//////////////////////////////////////////////////////////
	int counter = 0;
	if(isAllMC ||  chnl.Contains("bbTomu15"))        channels.insert(make_pair(increment(counter),"bbTomu15"));  // either that or bb_mu4mu4
	if(isAllMC ||  chnl.Contains("bb_mu4mu4"))       channels.insert(make_pair(increment(counter),"bb_mu4mu4")); // either that or bbTomu15
	if(isAllMC ||  chnl.Contains("ccTomu15"))        channels.insert(make_pair(increment(counter),"ccTomu15"));
	if(isAllMC ||  chnl.Contains("JZ"))              channels.insert(make_pair(increment(counter),"JZxW"));
	if(isAllMC ||  chnl.Contains("bb_Jpsimu4mu4"))   channels.insert(make_pair(increment(counter),"bb_Jpsimu4mu4"));
	// if(isAllMC  ||  chnl.Contains("ZmumuNp"))     channels.insert(make_pair(increment(counter),"ZmumuNpX"));
	// if(isAllMC  ||  chnl.Contains("WmunuNp"))     channels.insert(make_pair(increment(counter),"WmunuNpX"));
	// if(isAllMC  ||  chnl.Contains("WtaunuNp"))    channels.insert(make_pair(increment(counter),"WtaunuNpX"));
	if(isAllMC  ||  chnl.Contains("bbTotau10_3mu"))  channels.insert(make_pair(increment(counter),"bbTotau10_3mu"));
	if(isAllMC  ||  chnl.Contains("ccTotau10_3mu"))  channels.insert(make_pair(increment(counter),"ccTotau10_3mu"));
	if(isAllMC  ||  chnl.Contains("Wtaunu_3mu"))     channels.insert(make_pair(increment(counter),"Wtaunu_3mu"));
	if(isAllMC  ||  chnl.Contains("period"))         channels.insert(make_pair(increment(counter),"Data"));
	if(isAllMC  ||  isSignalMC)                      channels.insert(make_pair(increment(counter),"Signals"));
	if(isAllMC  ||  isBkgroundMC)                    channels.insert(make_pair(increment(counter),"Backgrounds"));
	
	
	
	cout << "------------------------ enabled data periods ------------------------" << endl;
	for(TMapTSb::iterator it=periodenable.begin() ; it!=periodenable.end() ; ++it) { cout << "period: " << it->first << endl; }
	cout << "----------------------------------------------------------------------" << endl;
	cout << "-------------------------- enabled channels --------------------------" << endl;
	for(TMapuiTS::iterator it=channels.begin() ; it!=channels.end() ; ++it)
	{
		if(isData(it->second)) cout << "Data: " << it->second << endl;
		else                   cout << "MC:   " << it->second << endl;
	}
	cout << "----------------------------------------------------------------------" << endl;
	
	
	
	
	nevents.insert(make_pair("Wtaunu_3mu",          0));
	nevents.insert(make_pair("bbTotau10_3mu",       0));
	nevents.insert(make_pair("ccTotau10_3mu",       0));
	nevents.insert(make_pair("bb_Jpsimu4mu4",       0));
	nevents.insert(make_pair("bbTomu15",            0));
	nevents.insert(make_pair("bb_mu4mu4",           0));
	nevents.insert(make_pair("ccTomu15",            0));
	nevents.insert(make_pair("JZxW",                0));
	nevents.insert(make_pair("WmunuNpX",            0));
	nevents.insert(make_pair("ZmumuNpX",            0));
	nevents.insert(make_pair("WtaunuNpX",           0));
	nevents.insert(make_pair("Data",                0));
	nevents.insert(make_pair("Backgrounds",         0));
	nevents.insert(make_pair("Signals",             0));	
	
	
	types.insert(make_pair("Wtaunu_3mu",          "3mu"));
	types.insert(make_pair("bbTotau10_3mu",       "3mu"));
	types.insert(make_pair("ccTotau10_3mu",       "3mu"));
	types.insert(make_pair("bb_Jpsimu4mu4",       "3mu"));
	types.insert(make_pair("bbTomu15",            "3mu"));
	types.insert(make_pair("bb_mu4mu4",           "3mu"));
	types.insert(make_pair("ccTomu15",            "3mu"));
	types.insert(make_pair("JZxW",                "3mu"));
	types.insert(make_pair("WmunuNpX",            "3mu"));
	types.insert(make_pair("ZmumuNpX",            "3mu"));
	types.insert(make_pair("WtaunuNpX",           "3mu"));
	types.insert(make_pair("Data",                "3mu"));
	types.insert(make_pair("Backgrounds",         "3mu"));
	types.insert(make_pair("Signals",             "3mu"));
	
	
	// rest of the maps shoud include all channels,
	// wheather or not they are switched on.
	properties("Wtaunu_3mu",          "#it{W#rightarrow#nu#tau#rightarrow3#mu}#times"+sBRf,  kBlue,     3001,"f",1, labels,colors,patterns,legoptions);
	properties("bbTotau10_3mu",       "#it{bb#rightarrow#tau10#rightarrow3#mu}#times"+sBRf,  kGreen+2,  3001,"f",1, labels,colors,patterns,legoptions);
	properties("ccTotau10_3mu",       "#it{cc#rightarrow#tau10#rightarrow3#mu}#times"+sBRf,  kGreen-6,  3001,"f",1, labels,colors,patterns,legoptions);
	properties("bb_Jpsimu4mu4",       "#it{J/#psi#rightarrow#mu4#mu4}",                      kRed-9,    3015,"f",0, labels,colors,patterns,legoptions);
	properties("bbTomu15",            "#it{bb#rightarrow#mu15} (FONLL)",                     kRed-3,    3005,"f",0, labels,colors,patterns,legoptions);
	properties("bb_mu4mu4",           "#it{bb#rightarrow#mu4#mu4}",                          kRed-3,    3005,"f",1, labels,colors,patterns,legoptions);
	properties("ccTomu15",            "#it{cc#rightarrow#mu15} (FONLL)",                     kRed+3,    3004,"f",1, labels,colors,patterns,legoptions);
	properties("JZxW",                "JZxW (dijet)",                                        kGray+2,   3022,"f",0, labels,colors,patterns,legoptions);
	properties("WmunuNpX",            "#it{W#rightarrow#mu#nu+jets}",                        kViolet+5, 3022,"f",0, labels,colors,patterns,legoptions);
	properties("ZmumuNpX",            "#it{Z#rightarrow#mu#mu+jets}",                        kViolet+7, 3023,"f",0, labels,colors,patterns,legoptions);
	properties("WtaunuNpX",           "#it{W#rightarrow#tau#nu+jets}",                       kAzure,    3011,"f",0, labels,colors,patterns,legoptions);
	properties("Data",                "Data "+slumi,                                         kBlack,    3004,"p",1, labels,colors,patterns,legoptions);
	properties("Backgrounds",         "#SigmaBg",                                            kRed-3,    3005,"f",1, labels,colors,patterns,legoptions);
	properties("Signals",             "#SigmaSig",                                           kGreen+2,  3001,"f",1, labels,colors,patterns,legoptions);
	
	
	///////////////////////////////
	// Initialize trigger bits ////
	buildTriggerbits(); ///////////
	///////////////////////////////
	
	
	///////////////////////////////////////////////////////////////////////////////
	// add the histos - no need to add the binned samples - just one is enough ////
	///////////////////////////////////////////////////////////////////////////////
	for(TMapuiTS::iterator it=channels.begin() ; it!=channels.end() ; it++)
	{
		TString name = it->second;
		_INFO("Booking "+(string)name);
		
		olddir->cd();
		
		clearCounters();
		initCounters();
		makeCutflowHisto(name,histos,labels);
		makeCountersHisto(name,histos,labels);
		makeCategoriesHisto(name,histos,labels);
		clearCounters();
		
		_DEBUG("");
		
		addHist(histos,name,"met_calibration_et_reldiff", ";MET E_{T}^{calib.}/E_{T}^{uncalib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"met_calibration_phi_reldiff", ";MET #phi^{calib.}/#phi^{uncalib.}-1;Normalized", labels,100,-3.,3.);

		addHist(histos,name,"jet_calibration_pt_reldiff", ";Jets p_{T}^{calib.}/p_{T}^{uncalib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_calibration_E_reldiff", ";Jets E^{calib.}/E^{uncalib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_calibration_m_reldiff", ";Jets m^{calib.}/m^{uncalib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_calibration_eta_reldiff", ";Jets #eta^{calib.}/#eta^{uncalib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_calibration_phi_reldiff", ";Jets #phi^{calib.}/#phi^{uncalib.}-1;Normalized", labels,100,-3.,3.);
		
		addHist(histos,name,"jet_resolution_pt_reldiff", ";Jets p_{T}^{smearup.}/p_{T}^{calib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_resolution_E_reldiff", ";Jets E^{smearup.}/E^{calib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_resolution_m_reldiff", ";Jets m^{smearup.}/m^{calib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_resolution_eta_reldiff", ";Jets #eta^{smearup.}/#eta^{calib.}-1;Normalized", labels,100,-3.,3.);
		addHist(histos,name,"jet_resolution_phi_reldiff", ";Jets #phi^{smearup.}/#phi^{calib.}-1;Normalized", labels,100,-3.,3.);
		
		addHist(histos2,name,"jet_calibration_reldiff_vs_pt",  "Jets p_{T}^{calib.}/p_{T}^{uncalib.}-1 vs p_{T}^{uncalib.};p_{T}^{uncalib.} [MeV];Jets p_{T}^{calib.}/p_{T}^{uncalib.}-1;Normalized",labels, 100,0.,200.*GeV2MeV, 100,-3.,+3.);
		addHist(histos2,name,"jet_calibration_reldiff_vs_E",   "Jets E^{calib.}/E^{uncalib.}-1 vs E^{uncalib.};E^{uncalib.} [MeV];Jets E^{calib.}/E^{uncalib.}-1;Normalized",labels, 100,0.,200.*GeV2MeV, 100,-3.,+3.);
		addHist(histos2,name,"jet_calibration_reldiff_vs_m",   "Jets m^{calib.}/m^{uncalib.}-1 vs m^{uncalib.};m^{uncalib.} [MeV];Jets m^{calib.}/m^{uncalib.}-1;Normalized",labels, 100,0.,20.*GeV2MeV, 100,-3.,+3.);
		addHist(histos2,name,"jet_calibration_reldiff_vs_eta", "Jets #eta^{calib.}/#eta^{uncalib.}-1 vs #eta^{uncalib.};#eta^{uncalib.};Jets #eta^{calib.}/#eta^{uncalib.}-1;Normalized",labels, 100,-5.,+5., 100,-3.,+3.);
		addHist(histos2,name,"jet_calibration_reldiff_vs_phi", "Jets #phi^{calib.}/#phi^{uncalib.}-1 vs #phi^{uncalib.};#phi^{uncalib.};Jets #phi^{calib.}/#phi^{uncalib.}-1;Normalized",labels, 100,-TMath::Pi(),+TMath::Pi(), 100,-3.,+3.);
		
		addHist(histos2,name,"jet_resolution_reldiff_vs_pt",  "Jets p_{T}^{smearup.}/p_{T}^{calib.}-1 vs p_{T}^{calib.};p_{T}^{calib.} [MeV];Jets p_{T}^{smearup.}/p_{T}^{calib.}-1;Normalized",labels, 100,0.,200.*GeV2MeV, 100,-3.,+3.);
		addHist(histos2,name,"jet_resolution_reldiff_vs_E",   "Jets E^{smearup.}/E^{calib.}-1 vs E^{calib.};E^{calib.} [MeV];Jets E^{smearup.}/E^{calib.}-1;Normalized",labels, 100,0.,200.*GeV2MeV, 100,-3.,+3.);
		addHist(histos2,name,"jet_resolution_reldiff_vs_m",   "Jets m^{smearup.}/m^{calib.}-1 vs m^{calib.};m^{calib.} [MeV];Jets m^{smearup.}/m^{calib.}-1;Normalized",labels, 100,0.,20.*GeV2MeV, 100,-3.,+3.);
		addHist(histos2,name,"jet_resolution_reldiff_vs_eta", "Jets #eta^{smearup.}/#eta^{calib.}-1 vs #eta^{calib.};#eta^{calib.};Jets #eta^{smearup.}/#eta^{calib.}-1;Normalized",labels, 100,-5.,+5., 100,-3.,+3.);
		addHist(histos2,name,"jet_resolution_reldiff_vs_phi", "Jets #phi^{smearup.}/#phi^{calib.}-1 vs #phi^{calib.};#phi^{calib.};Jets #phi^{smearup.}/#phi^{calib.}-1;Normalized",labels, 100,-TMath::Pi(),+TMath::Pi(), 100,-3.,+3.);

		_DEBUG("");
		
		addHist(histos,name,"MVA_score_all",            ";Inclusive MVA score;Normalized", labels,100,-1.,1.);
		addHist(histos,name,"MVA_score_left_sideband",  ";Left sideband MVA score;Normalized", labels,100,-1.,1.);
		addHist(histos,name,"MVA_score_right_sideband", ";Right sideband MVA score;Normalized", labels,100,-1.,1.);
		addHist(histos2,name,"m3mu_vs_MVA_score",  "m_{3body} vs MVA score;MVA score;m_{3body} [MeV];Normalized",labels,         100,-1.,+1., 50,0.,4000.);
		addHist(histos2,name,"pT3mu_vs_MVA_score", "p_{T}^{3body} vs MVA score;MVA score;p_{T}^{3body} [MeV];Normalized",labels, 100,-1.,+1., 50,0.,4000.);
		addHist(histos2,name,"mOS1_vs_MVA_score",  "m_{OS1} vs MVA score;MVA score;m_{OS1} [MeV];Normalized",labels,             100,-1.,+1., 100,0.,3000.);
		addHist(histos2,name,"mOS2_vs_MVA_score",  "m_{OS2} vs MVA score;MVA score;m_{OS2} [MeV];Normalized",labels,             100,-1.,+1., 100,0.,3000.);
		addHist(histos2,name,"mSS_vs_MVA_score",   "m_{SS}  vs MVA score;MVA score;m_{SS} [MeV];Normalized",labels,              100,-1.,+1., 100,0.,3000.);
		
		_DEBUG("");
		
		addHist(histos,name,"MVA_vars_pT3body_before_vtxmet",        ";p_{T}^{3body} [MeV] (before vtx-met);Normalized",                                     labels,50,20.*GeV2MeV,70.*GeV2MeV);
		addHist(histos,name,"MVA_vars_mOS1_before_vtxmet",           ";m_{OS1} [MeV] (before vtx-met);Normalized",                                           labels,100,0.,2100.);
		addHist(histos,name,"MVA_vars_mOS2_before_vtxmet",           ";m_{OS2} [MeV] (before vtx-met);Normalized",                                           labels,100,0.,2100.);
		addHist(histos,name,"MVA_vars_mSS_before_vtxmet",            ";m_{SS} [MeV] (before vtx-met);Normalized",                                            labels,100,0.,2100.);
		addHist(histos,name,"MVA_vars_isolation_before_vtxmet",      ";#Sigmap_{T}^{trk}/p_{T}^{3body} in #DeltaR_{max}+#delta (before vtx-met);Normalized", labels,100,0.,0.5);
		addHist(histos,name,"MVA_vars_pvalue_before_vtxmet",         ";#it{p}-value (before vtx-met);Normalized",                                            labels,80,0.,1.);
		addHist(histos,name,"MVA_vars_Lxy_before_vtxmet",            ";#it{L}_{xy} (before vtx-met);Normalized",                                             labels,50,-5.,15.);
		addHist(histos,name,"MVA_vars_a0xy_before_vtxmet",           ";#it{a}^{0}_{xy} (before vtx-met);Normalized",                                         labels,100,0.,1.);
		addHist(histos,name,"MVA_vars_cosTxy_before_vtxmet",         ";cos#it{#theta}_{xy} (before vtx-met);Normalized",                                     labels,100,0.,1.);
		addHist(histos,name,"MVA_vars_MET_before_vtxmet",            ";#it{E}_{T}^{miss} [MeV] (before vtx-met);Normalized",                                 labels,50,0.,100.*GeV2MeV);
		addHist(histos,name,"MVA_vars_dPhi3bodyMET_before_vtxmet",   ";#Delta#phi(3body,#it{E}_{T}^{miss}) (before vtx-met);Normalized",                     labels,64,0.,TMath::Pi());
		addHist(histos,name,"MVA_vars_mT3bodyMET_before_vtxmet",     ";m_{T}(3body,#it{E}_{T}^{miss}) [Mev] (before vtx-met);Normalized",                    labels,60,0.,120.*GeV2MeV);
		addHist(histos,name,"MVA_vars_dPhi3bodyJ1_before_vtxmet",    ";#Delta#phi(3body,jet_{1}) (before vtx-met);Normalized",                               labels,64,0.,TMath::Pi());		
		addHist(histos,name,"MVA_vars_dPhiJ1J2_before_vtxmet",       ";#Delta#phi(jet_{1},jet_{2}) (before vtx-met);Normalized",                             labels,64,0.,TMath::Pi());
		addHist(histos,name,"MVA_vars_pTJ1_before_vtxmet",           ";Leading jet #it{p}_{T} [MeV] (before vtx-met);Normalized",                            labels,25,0,100*GeV2MeV);
		addHist(histos,name,"MVA_vars_pTJ2_before_vtxmet",           ";Subleading jet #it{p}_{T} [MeV] (before vtx-met);Normalized",                         labels,25,0,100*GeV2MeV);
		addHist(histos,name,"MVA_vars_mupt12Fraction_before_vtxmet", ";#Deltap_{T}^{#mu1-#mu2}/p_{T}^{3body} (before vtx-met);Normalized",                   labels, 50,0.3,1.00);			
		addHist(histos,name,"MVA_vars_mupt23Fraction_before_vtxmet", ";#Deltap_{T}^{#mu2-#mu3}/p_{T}^{3body} (before vtx-met);Normalized",                   labels, 50,0.3,1.00);			
		addHist(histos,name,"MVA_vars_mupt13Fraction_before_vtxmet", ";#Deltap_{T}^{#mu1-#mu3}/p_{T}^{3body} (before vtx-met);Normalized",                   labels, 50,0.0,0.35);
		addHist(histos,name,"MVA_vars_m3mu_before_vtxmet",           ";m_{3body} [MeV] (before vtx-met);Events",                                              labels, 50,0.,4000.);
		addHist(histos,name,"MVA_vars_m3mu_lin_before_vtxmet",       ";m_{3body} [MeV] (before vtx-met);Events",                                              labels, 50,500.,4000.);
		addHist(histos,name,"MVA_vars_m3mu_lin_zoom_before_vtxmet",  ";m_{3body} [MeV] (before vtx-met);Events",                                              labels, 100,800.,2800.);
		addHist(histos,name,"MVA_vars_m3mu_sigregion_before_vtxmet", ";m_{3body} [MeV] (before vtx-met);Events",                                              labels, 50,1300.,2300.);
		addHist(histos2,name,"MVA_vars_dPhi3muJet1_vs_pTjet1_before_vtxmet",      "#Delta#phi(3#mu,jet_{1}) vs pT jet_{1} (before vtx-met);pT jet_{1};#Delta#phi(3#mu,jet_{1});Normalized",                      labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"MVA_vars_dPhiJet1Jet2_vs_sumpTjet12_before_vtxmet", "#Delta#phi(jet_{1},jet_{2}) vs #SigmapT jet_{1,2} (before vtx-met);#SigmapT jet_{1,2};#Delta#phi(jet_{1},jet_{2});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		
		addHist(histos,name,"MVA_vars_pT3body_after_vtxmet",         ";p_{T}^{3body} [MeV] (after vtx-met);Normalized",                                     labels,50,20.*GeV2MeV,70.*GeV2MeV);
		addHist(histos,name,"MVA_vars_mOS1_after_vtxmet",            ";m_{OS1} [MeV] (after vtx-met);Normalized",                                           labels,100,0.,2500.);
		addHist(histos,name,"MVA_vars_mOS2_after_vtxmet",            ";m_{OS2} [MeV] (after vtx-met);Normalized",                                           labels,100,0.,2500.);
		addHist(histos,name,"MVA_vars_mSS_after_vtxmet",             ";m_{SS} [MeV] (after vtx-met);Normalized",                                            labels,100,0.,2500.);
		addHist(histos,name,"MVA_vars_isolation_after_vtxmet",       ";#Sigmap_{T}^{trk}/p_{T}^{3body} in #DeltaR_{max}+#delta (after vtx-met);Normalized", labels,100,0.,0.5);
		addHist(histos,name,"MVA_vars_pvalue_after_vtxmet",          ";#it{p}-value (after vtx-met);Normalized",                                            labels,80,0.,1.);
		addHist(histos,name,"MVA_vars_Lxy_after_vtxmet",             ";#it{L}_{xy} (after vtx-met);Normalized",                                             labels,50,-5.,15.);
		addHist(histos,name,"MVA_vars_a0xy_after_vtxmet",            ";#it{a}^{0}_{xy} (after vtx-met);Normalized",                                         labels,100,0.,1.);
		addHist(histos,name,"MVA_vars_cosTxy_after_vtxmet",          ";cos#it{#theta}_{xy} (after vtx-met);Normalized",                                     labels,100,0.,1.);
		addHist(histos,name,"MVA_vars_MET_after_vtxmet",             ";#it{E}_{T}^{miss} [MeV] (after vtx-met);Normalized",                                 labels,50,0.,100.*GeV2MeV);
		addHist(histos,name,"MVA_vars_dPhi3bodyMET_after_vtxmet",    ";#Delta#phi(3body,#it{E}_{T}^{miss}) (after vtx-met);Normalized",                     labels,64,0.,TMath::Pi());
		addHist(histos,name,"MVA_vars_mT3bodyMET_after_vtxmet",      ";m_{T}(3body,#it{E}_{T}^{miss}) [Mev] (after vtx-met);Normalized",                    labels,60,0.,120.*GeV2MeV);
		addHist(histos,name,"MVA_vars_dPhi3bodyJ1_after_vtxmet",     ";#Delta#phi(3body,jet_{1}) (after vtx-met);Normalized",                               labels,64,0.,TMath::Pi());		
		addHist(histos,name,"MVA_vars_dPhiJ1J2_after_vtxmet",        ";#Delta#phi(jet_{1},jet_{2}) (after vtx-met);Normalized",                             labels,64,0.,TMath::Pi());
		addHist(histos,name,"MVA_vars_pTJ1_after_vtxmet",            ";Leading jet #it{p}_{T} [MeV] (after vtx-met);Normalized",                            labels,25,0,100*GeV2MeV);
		addHist(histos,name,"MVA_vars_pTJ2_after_vtxmet",            ";Subleading jet #it{p}_{T} [MeV] (after vtx-met);Normalized",                         labels,25,0,100*GeV2MeV);
		addHist(histos,name,"MVA_vars_mupt12Fraction_after_vtxmet",  ";#Deltap_{T}^{#mu1-#mu2}/p_{T}^{3body} (after vtx-met);Normalized",                   labels, 50,0.3,1.00);			
		addHist(histos,name,"MVA_vars_mupt23Fraction_after_vtxmet",  ";#Deltap_{T}^{#mu2-#mu3}/p_{T}^{3body} (after vtx-met);Normalized",                   labels, 50,0.0,0.35);
		addHist(histos,name,"MVA_vars_mupt13Fraction_after_vtxmet",  ";#Deltap_{T}^{#mu1-#mu3}/p_{T}^{3body} (after vtx-met);Normalized",                   labels, 50,0.0,0.35);
		addHist(histos,name,"MVA_vars_m3mu_after_vtxmet",            ";m_{3body} [MeV] (after vtx-met);Events",                                              labels, 50,0.,4000.);
		addHist(histos,name,"MVA_vars_m3mu_lin_after_vtxmet",        ";m_{3body} [MeV] (after vtx-met);Events",                                              labels, 50,500.,4000.);
		addHist(histos,name,"MVA_vars_m3mu_lin_zoom_after_vtxmet",   ";m_{3body} [MeV] (after vtx-met);Events",                                              labels, 100,800.,2800.);
		addHist(histos,name,"MVA_vars_m3mu_sigregion_after_vtxmet",  ";m_{3body} [MeV] (after vtx-met);Events",                                              labels, 50,1300.,2300.);
		addHist(histos2,name,"MVA_vars_dPhi3muJet1_vs_pTjet1_after_vtxmet",      "#Delta#phi(3#mu,jet_{1}) vs pT jet_{1} (after vtx-met);pT jet_{1};#Delta#phi(3#mu,jet_{1});Normalized",                       labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"MVA_vars_dPhiJet1Jet2_vs_sumpTjet12_after_vtxmet", "#Delta#phi(jet_{1},jet_{2}) vs #SigmapT jet_{1,2} (after vtx-met);#SigmapT jet_{1,2};#Delta#phi(jet_{1},jet_{2});Normalized", labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		
		
		addHist(histos2,name,"flwMV1_vs_mjet", "MV1 flavor weight vs #it{m}_{jet};#it{m}_{jet} [MeV];MV1 flavor weight;Normalized",labels,   20,0,20*GeV2MeV, 20,0.,1.);
		addHist(histos2,name,"flwMV1_vs_pTjet","MV1 flavor weight vs #it{pT}_{jet};#it{pT}_{jet} [MeV];MV1 flavor weight;Normalized",labels, 20,0,100*GeV2MeV, 20,0.,1.);

		addHist(histos2,name,"dPhiJet1Jet2_vs_sumpTjet12_before_looseveto", "#Delta#phi(jet_{1},jet_{2}) vs #SigmapT jet_{1,2} before loose veto;#SigmapT jet_{1,2};#Delta#phi(jet_{1},jet_{2});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"dPhiJet1Jet2_vs_sumpTjet12_before_tightveto", "#Delta#phi(jet_{1},jet_{2}) vs #SigmapT jet_{1,2} before tight veto;#SigmapT jet_{1,2};#Delta#phi(jet_{1},jet_{2});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"dPhiJet1Jet2_vs_sumpTjet12_after_looseveto",  "#Delta#phi(jet_{1},jet_{2}) vs #SigmapT jet_{1,2} after loose veto;#SigmapT jet_{1,2};#Delta#phi(jet_{1},jet_{2});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"dPhiJet1Jet2_vs_sumpTjet12_after_tightveto",  "#Delta#phi(jet_{1},jet_{2}) vs #SigmapT jet_{1,2} after tight veto;#SigmapT jet_{1,2};#Delta#phi(jet_{1},jet_{2});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
	
		addHist(histos2,name,"dPhi3muJet1_vs_pTjet1_before_looseveto", "#Delta#phi(3body,jet_{1}) vs pT jet_{1} before loose veto;pT jet_{1};#Delta#phi(3body,jet_{1});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"dPhi3muJet1_vs_pTjet1_before_tightveto", "#Delta#phi(3body,jet_{1}) vs pT jet_{1} before tight veto;pT jet_{1};#Delta#phi(3body,jet_{1});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"dPhi3muJet1_vs_pTjet1_after_looseveto",  "#Delta#phi(3body,jet_{1}) vs pT jet_{1} after loose veto;pT jet_{1};#Delta#phi(3body,jet_{1});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		addHist(histos2,name,"dPhi3muJet1_vs_pTjet1_after_tightveto",  "#Delta#phi(3body,jet_{1}) vs pT jet_{1} after tight veto;pT jet_{1};#Delta#phi(3body,jet_{1});Normalized",labels,20,0,100*GeV2MeV ,20,0.,TMath::Pi());
		
		
		addHist(histos,name,"triplet_isolation_zoom_before_isolation",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+margins) before cutting;Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"triplet_isolation_zoom",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}+margins);Normalized",labels, 80,0.,0.5);
		
		addHist(histos,name,"mu_ptcone10.1_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.1);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone10.2_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.1);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone10.3_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.1);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone20.1_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.2);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone20.2_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.2);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone20.3_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.2);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone30.1_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.3);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone30.2_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.3);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone30.3_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.3);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone40.1_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.4);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone40.2_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.4);Normalized",labels, 80,0.,0.5);
		addHist(histos,name,"mu_ptcone40.3_zoom",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.4);Normalized",labels, 80,0.,0.5);
		
		addHist(histos,name,"triggers",               "Trigger rates (fractional) after 2nd skim;;Normalized rate",labels,triggerorder.size(),0,triggerorder.size(),false);
		addHist(histos,name,"triggers_after_vtxmet",  "Trigger rates (fractional) after vtx-met;;Normalized rate",labels,triggerorder.size(),0,triggerorder.size(),false);
		addHist(histos,name,"triggers_after_hadclean","Trigger rates (fractional) after hadronic cleaning;;Normalized rate",labels,triggerorder.size(),0,triggerorder.size(),false);
		
		addHist(histos,name,"unique_triggers",                 "Unique trigger rates (fractional) after 2nd skim;;Normalized rate",labels,triggerorder.size(),0,triggerorder.size(),false);
		addHist(histos,name,"unique_triggers_after_vtxmet",  "Unique trigger rates (fractional) after vtx-met;;Normalized rate",labels,triggerorder.size(),0,triggerorder.size(),false);
		addHist(histos,name,"unique_triggers_after_hadclean", "Unique trigger rates (fractional) after hadronic cleaning;;Normalized rate",labels,triggerorder.size(),0,triggerorder.size(),false);
		
		///////// btagging 1d
		addHist(histos,name,"mu_m3mu_beforeBtag_norm", "3#mu mass before #it{b}-tag;#it{m}_{3#mu} [MeV];Normalized",labels,30,500,3500,false);
		addHist(histos,name,"mu_m3mu_afterBtag_norm",  "3#mu mass after #it{b}-tag;#it{m}_{3#mu} [MeV];Normalized",labels,30,500,3500,false);
		addHist(histos,name,"mu_m3mu_beforeBtag",      "3#mu mass before #it{b}-tag;#it{m}_{3#mu} [MeV];",labels,30,500,3500,false);
		addHist(histos,name,"mu_m3mu_afterBtag",       "3#mu mass after #it{b}-tag;#it{m}_{3#mu} [MeV];",labels,30,500,3500,false);
		addHist(histos,name,"mu_m3mu_ratioBtag",        "Ratio of 3#mu mass after/before #it{b}-tag;#it{m}_{3#mu} [MeV];Ratio (before/after #it{b}-tag)",labels,30,500,3500,false);

		addHist(histos,name,"bjet_n",  ";#it{n}_{#it{b}-jet};Normalized",labels,10,0,10);
		addHist(histos,name,"bjet_n10",";#it{n}_{#it{b}-jet} for #it{p}_{T}>10 GeV;Normalized",labels,8,0,8);
		addHist(histos,name,"bjet_n20",";#it{n}_{#it{b}-jet} for #it{p}_{T}>20 GeV;Normalized",labels,5,0,5);
		addHist(histos,name,"bjet_E",  ";#it{b}-jet #it{E} [MeV];Normalized",labels,25,0,100*GeV2MeV);
		addHist(histos,name,"bjet_pt", ";#it{b}-jet #it{p}_{T} [MeV];Normalized",labels,25,0,100*GeV2MeV);
		addHist(histos,name,"bjet_m",  ";#it{b}-jet #it{m} [MeV];Normalized",labels,25,0,20*GeV2MeV);
		addHist(histos,name,"bjet_eta",";#it{b}-jet #eta;Normalized",labels,32,-2.7,+2.7);
		addHist(histos,name,"bjet_phi",";#it{b}-jet #phi;Normalized",labels,32,-TMath::Pi(),+TMath::Pi());

		addHist(histos,name,"jets_dPhi3bodyJ1_before_looseveto", ";#Delta#phi(3body,jet_{1}) before loose veto;Normalized",labels,64,0.,TMath::Pi());		
		addHist(histos,name,"jets_dPhiJ1J2_before_looseveto",    ";#Delta#phi(jet_{1},jet_{2}) before loose veto;Normalized",labels,64,0.,TMath::Pi());
		addHist(histos,name,"jets_dPhi3bodyJ1_before_tightveto", ";#Delta#phi(3body,jet_{1}) before tight veto;Normalized",labels,64,0.,TMath::Pi());		
		addHist(histos,name,"jets_dPhiJ1J2_before_tightveto",    ";#Delta#phi(jet_{1},jet_{2}) before tight veto;Normalized",labels,64,0.,TMath::Pi());
		addHist(histos,name,"jets_dPhi3bodyJ1_after_looseveto", ";#Delta#phi(3body,jet_{1}) after loose veto;Normalized",labels,64,0.,TMath::Pi());		
		addHist(histos,name,"jets_dPhiJ1J2_after_looseveto",    ";#Delta#phi(jet_{1},jet_{2}) after loose veto;Normalized",labels,64,0.,TMath::Pi());
		addHist(histos,name,"jets_dPhi3bodyJ1_after_tightveto", ";#Delta#phi(3body,jet_{1}) after tight veto;Normalized",labels,64,0.,TMath::Pi());		
		addHist(histos,name,"jets_dPhiJ1J2_after_tightveto",    ";#Delta#phi(jet_{1},jet_{2}) after tight veto;Normalized",labels,64,0.,TMath::Pi());
		
		addHist(histos,name,"bjet_dR3muLeadingJet",        ";#DeltaR(3#mu,jet) Leading;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muLeadingJet",      ";#Delta#phi(3#mu,jet) Leading;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muSubleadingJet",     ";#DeltaR(3#mu,jet) Subleading;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muSubleadingJet",   ";#Delta#phi(3#mu,jet) Subleading;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muBjet_btagwgt1",     ";#DeltaR(3#mu,#it{b}-jet) highest flavor weight;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muBjet_btagwgt1",   ";#Delta#phi(3#mu,#it{b}-jet) highest flavor weight;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muBjet_btagwgt2",     ";#DeltaR(3#mu,#it{b}-jet) 2nd highest flavor weight;Normalized",labels,50,0,5);
		addHist(histos,name,"bjet_dPhi3muBjet_btagwgt2",   ";#Delta#phi(3#mu,#it{b}-jet) 2nd highest flavor weight;Normalized",labels,32,0.,TMath::Pi());
		addHist(histos,name,"bjet_dPhi_jet1_jet2",         ";#Delta#phi(jet_{pT1},jet_{pT2});Normalized",labels,32,0.,TMath::Pi());
		addHist(histos,name,"bjet_dPhi_btagwgt1_btagwgt2", ";#Delta#phi(jet_{fw1},jet_{fw2});Normalized",labels,32,0.,TMath::Pi());
		
		addHist(histos,name,"bjet_dR3muLeadingJet_noOvrlp",        ";#DeltaR(3#mu,jet) Leading, no 3#mu overlap;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muLeadingJet_noOvrlp",      ";#Delta#phi(3#mu,jet) Leading, no 3#mu overlap;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muSubleadingJet_noOvrlp",     ";#DeltaR(3#mu,jet) Subleading, no 3#mu overlap;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muSubleadingJet_noOvrlp",   ";#Delta#phi(3#mu,jet) Subleading, no 3#mu overlap;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muBjet_btagwgt1_noOvrlp",     ";#DeltaR(3#mu,#it{b}-jet) highest flavor weight, no 3#mu overlap;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muBjet_btagwgt1_noOvrlp",   ";#Delta#phi(3#mu,#it{b}-jet) highest flavor weight, no 3#mu overlap;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muBjet_btagwgt2_noOvrlp",     ";#DeltaR(3#mu,#it{b}-jet) 2nd highest flavor weight, no 3#mu overlap;Normalized",labels,50,0,5);
		addHist(histos,name,"bjet_dPhi3muBjet_btagwgt2_noOvrlp",   ";#Delta#phi(3#mu,#it{b}-jet) 2nd highest flavor weight, no 3#mu overlap;Normalized",labels,32,0.,TMath::Pi());
		addHist(histos,name,"bjet_dPhi_jet1_jet2_noOvrlp",         ";#Delta#phi(jet_{pT1},jet_{pT2}), no 3#mu overlap;Normalized",labels,32,0.,TMath::Pi());
		addHist(histos,name,"bjet_dPhi_btagwgt1_btagwgt2_noOvrlp", ";#Delta#phi(jet_{fw1},jet_{fw2}), no 3#mu overlap;Normalized",labels,32,0.,TMath::Pi());

		addHist(histos,name,"bjet_dR3muLeadingJet_stdMV1",        ";#DeltaR(3#mu,jet) Leading, Standard MV1 only;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muLeadingJet_stdMV1",      ";#Delta#phi(3#mu,jet) Leading, Standard MV1 only;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muSubleadingJet_stdMV1",     ";#DeltaR(3#mu,jet) Subleading, Standard MV1 only;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muSubleadingJet_stdMV1",   ";#Delta#phi(3#mu,jet) Subleading, Standard MV1 only;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muBjet_btagwgt1_stdMV1",     ";#DeltaR(3#mu,#it{b}-jet) highest flavor weight, Standard MV1 only;Normalized",labels,50,0,5);		
		addHist(histos,name,"bjet_dPhi3muBjet_btagwgt1_stdMV1",   ";#Delta#phi(3#mu,#it{b}-jet) highest flavor weight, Standard MV1 only;Normalized",labels,32,0.,TMath::Pi());		
		addHist(histos,name,"bjet_dR3muBjet_btagwgt2_stdMV1",     ";#DeltaR(3#mu,#it{b}-jet) 2nd highest flavor weight, Standard MV1 only;Normalized",labels,50,0,5);
		addHist(histos,name,"bjet_dPhi3muBjet_btagwgt2_stdMV1",   ";#Delta#phi(3#mu,#it{b}-jet) 2nd highest flavor weight, Standard MV1 only;Normalized",labels,32,0.,TMath::Pi());
		addHist(histos,name,"bjet_dPhi_jet1_jet2_stdMV1",         ";#Delta#phi(jet_{pT1},jet_{pT2}), Standard MV1 only;Normalized",labels,32,0.,TMath::Pi());
		addHist(histos,name,"bjet_dPhi_btagwgt1_btagwgt2_stdMV1", ";#Delta#phi(jet_{fw1},jet_{fw2}), Standard MV1 only;Normalized",labels,32,0.,TMath::Pi());
		
		addHist(histos,name,"bjet_maxFlvWgt_E",  ";#it{b}-jet (max MV1 flavor weight) #it{E} [MeV];Normalized",labels,25,0,100*GeV2MeV);
		addHist(histos,name,"bjet_maxFlvWgt_pt", ";#it{b}-jet (max MV1 flavor weight) #it{p}_{T} [MeV];Normalized",labels,25,0,100*GeV2MeV);
		addHist(histos,name,"bjet_maxFlvWgt_m",  ";#it{b}-jet (max MV1 flavor weight) #it{m} [MeV];Normalized",labels,25,0,20*GeV2MeV);
		addHist(histos,name,"bjet_maxFlvWgt_eta",";#it{b}-jet (max MV1 flavor weight) #eta;Normalized",labels,32,-2.7,+2.7);
		addHist(histos,name,"bjet_maxFlvWgt_phi",";#it{b}-jet (max MV1 flavor weight) #phi;Normalized",labels,32,-TMath::Pi(),+TMath::Pi());
		
		addHist(histos,name,"jet_n",  ";Jet #it{n};Normalized",labels,30,0,30);
		addHist(histos,name,"jet_n10",";Jet #it{n} for #it{p}_{T}>10 GeV;Normalized",labels,20,0,20);
		addHist(histos,name,"jet_n20",";Jet #it{n} for #it{p}_{T}>20 GeV;Normalized",labels,15,0,15);
		addHist(histos,name,"jet_E",  ";Jet #it{E} [MeV];Normalized",labels,25,0,100*GeV2MeV);
		addHist(histos,name,"jet_pt", ";Jet #it{p}_{T} [MeV];Normalized",labels,25,0,100*GeV2MeV);
		addHist(histos,name,"jet_m",  ";Jet #it{m} [MeV];Normalized",labels,25,0,20*GeV2MeV);
		addHist(histos,name,"jet_eta",";Jet #eta;Normalized",labels,32,-2.7,+2.7);
		addHist(histos,name,"jet_phi",";Jet #phi;Normalized",labels,32,-TMath::Pi(),+TMath::Pi());
		
		addHist(histos,name,"jet_maxFlvWgt_MV1",";Jet #it{MAX} flavor weight #it{MV1};Normalized",labels,40,0.,+1.);	
		addHist(histos,name,"jet_leadingJetFlvWgt_MV1",";Leading jet flavor weight #it{MV1};Normalized",labels,40,0.,+1.);
		addHist(histos,name,"jet_subleadingJetFlvWgt_MV1",";Subleading jet flavor weight #it{MV1};Normalized",labels,40,0.,+1.);
		
		addHist(histos,name,"met_sumet",    ";#Sigma#it{E}_{T} [MeV];Normalized",labels,50,100.*GeV2MeV,1000.*GeV2MeV);
		addHist(histos,name,"met_et",       ";#it{E}_{T}^{miss} [MeV];Normalized",labels,50,0.,100.*GeV2MeV);
		addHist(histos,name,"met_mt_et3mu", ";m_{T}(3#mu,#it{E}_{T}^{miss}) [MeV];Normalized",labels,60,0,120.*GeV2MeV);
		addHist(histos,name,"met_dphi3mu",  ";#Delta#phi(3#mu,#it{E}_{T}^{miss});Normalized",labels,32,0.,TMath::Pi());
		

		///////// btagging 2d
		addHist(histos2,name,"nBtag_vs_nJet","#it{n}_{#it{b}-jet} vs #it{n}_{jet};#it{n}_{jet};#it{n}_{#it{b}-jet};Normalized",labels,20,0,20, 20,0,20);

		addHist(histos2,name,"bjet_dR3muBjet_btagwgt1_vs_btagwgt1","#DeltaR(3#mu,#it{b}_{1}) vs Flavor weight #it{b}_{1};Flavor weight #it{b}_{1};#DeltaR(3#mu,#it{b}_{1});Normalized",labels,40,0.,+1., 50,0,5);
		addHist(histos2,name,"bjet_dR3muBjet_btagwgt2_vs_btagwgt2","#DeltaR(3#mu,#it{b}_{2}) vs Flavor weight #it{b}_{2};Flavor weight #it{b}_{2};#DeltaR(3#mu,#it{b}_{2});Normalized",labels,40,0.,+1., 50,0,5);
		addHist(histos2,name,"bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1","#Delta#phi(3#mu,#it{b}_{1}) vs Flavor weight #it{b}_{1};Flavor weight #it{b}_{1};#Delta#phi(3#mu,#it{b}_{1});Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2","#Delta#phi(3#mu,#it{b}_{2}) vs Flavor weight #it{b}_{2};Flavor weight #it{b}_{2};#Delta#phi(3#mu,#it{b}_{2});Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhib1b2_vs_b1","#Delta#phi(#it{b}_{1},#it{b}_{2}) vs Flavor weight #it{b}_{1};Flavor weight #it{b}_{1};#Delta#phi(#it{b}_{1},#it{b}_{2});Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhib1b2_vs_b2","#Delta#phi(#it{b}_{1},#it{b}_{2}) vs Flavor weight #it{b}_{2};Flavor weight #it{b}_{2};#Delta#phi(#it{b}_{1},#it{b}_{2});Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_wb1_vs_wb2","Flavor weight #it{b}_{1} vs Flavor weight #it{b}_{2};Flavor weight #it{b}_{2};Flavor weight #it{b}_{1};Normalized",labels,40,0.,+1., 40,0.,+1.);
		addHist(histos2,name,"jet_wj1_vs_wj2","Flavor weight #it{jet}_{1} vs Flavor weight #it{jet}_{2};Flavor weight #it{jet}_{2};Flavor weight #it{jet}_{1};Normalized",labels,40,0.,+1., 40,0.,+1.);
		
		addHist(histos2,name,"bjet_dR3muBjet_btagwgt1_vs_btagwgt1_noOvrlp","#DeltaR(3#mu,#it{b}_{1}) vs Flavor weight #it{b}_{1}, no 3#mu overlap;Flavor weight #it{b}_{1}, no 3#mu overlap;#DeltaR(3#mu,#it{b}_{1}), no 3#mu overlap;Normalized",labels,40,0.,+1., 50,0,5);
		addHist(histos2,name,"bjet_dR3muBjet_btagwgt2_vs_btagwgt2_noOvrlp","#DeltaR(3#mu,#it{b}_{2}) vs Flavor weight #it{b}_{2}, no 3#mu overlap;Flavor weight #it{b}_{2}, no 3#mu overlap;#DeltaR(3#mu,#it{b}_{2}), no 3#mu overlap;Normalized",labels,40,0.,+1., 50,0,5);
		addHist(histos2,name,"bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_noOvrlp","#Delta#phi(3#mu,#it{b}_{1}) vs Flavor weight #it{b}_{1}, no 3#mu overlap;Flavor weight #it{b}_{1}, no 3#mu overlap;#Delta#phi(3#mu,#it{b}_{1}), no 3#mu overlap;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_noOvrlp","#Delta#phi(3#mu,#it{b}_{2}) vs Flavor weight #it{b}_{2}, no 3#mu overlap;Flavor weight #it{b}_{2}, no 3#mu overlap;#Delta#phi(3#mu,#it{b}_{2}), no 3#mu overlap;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhib1b2_vs_b1_noOvrlp","#Delta#phi(#it{b}_{1},#it{b}_{2}) vs Flavor weight #it{b}_{1}, no 3#mu overlap;Flavor weight #it{b}_{1}, no 3#mu overlap;#Delta#phi(#it{b}_{1},#it{b}_{2}), no 3#mu overlap;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhib1b2_vs_b2_noOvrlp","#Delta#phi(#it{b}_{1},#it{b}_{2}) vs Flavor weight #it{b}_{2}, no 3#mu overlap;Flavor weight #it{b}_{2}, no 3#mu overlap;#Delta#phi(#it{b}_{1},#it{b}_{2}), no 3#mu overlap;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_wb1_vs_wb2_noOvrlp","Flavor weight #it{b}_{1} vs Flavor weight #it{b}_{2}, no 3#mu overlap;Flavor weight #it{b}_{2}, no 3#mu overlap;Flavor weight #it{b}_{1}, no 3#mu overlap;Normalized",labels,40,0.,+1., 40,0.,+1.);
		addHist(histos2,name,"jet_wj1_vs_wj2_noOvrlp","Flavor weight #it{jet}_{1} vs Flavor weight #it{jet}_{2}, no 3#mu overlap;Flavor weight #it{jet}_{2}, no 3#mu overlap;Flavor weight #it{jet}_{1}, no 3#mu overlap;Normalized",labels,40,0.,+1., 40,0.,+1.);
		
		addHist(histos2,name,"bjet_dR3muBjet_btagwgt1_vs_btagwgt1_stdMV1","#DeltaR(3#mu,#it{b}_{1}) vs Flavor weight #it{b}_{1}, Standard MV1 only;Flavor weight #it{b}_{1}, Standard MV1 only;#DeltaR(3#mu,#it{b}_{1}), Standard MV1 only;Normalized",labels,40,0.,+1., 50,0,5);
		addHist(histos2,name,"bjet_dR3muBjet_btagwgt2_vs_btagwgt2_stdMV1","#DeltaR(3#mu,#it{b}_{2}) vs Flavor weight #it{b}_{2}, Standard MV1 only;Flavor weight #it{b}_{2}, Standard MV1 only;#DeltaR(3#mu,#it{b}_{2}), Standard MV1 only;Normalized",labels,40,0.,+1., 50,0,5);
		addHist(histos2,name,"bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_stdMV1","#Delta#phi(3#mu,#it{b}_{1}) vs Flavor weight #it{b}_{1}, Standard MV1 only;Flavor weight #it{b}_{1}, Standard MV1 only;#Delta#phi(3#mu,#it{b}_{1}), Standard MV1 only;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_stdMV1","#Delta#phi(3#mu,#it{b}_{2}) vs Flavor weight #it{b}_{2}, Standard MV1 only;Flavor weight #it{b}_{2}, Standard MV1 only;#Delta#phi(3#mu,#it{b}_{2}), Standard MV1 only;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhib1b2_vs_b1_stdMV1","#Delta#phi(#it{b}_{1},#it{b}_{2}) vs Flavor weight #it{b}_{1}, Standard MV1 only;Flavor weight #it{b}_{1}, Standard MV1 only;#Delta#phi(#it{b}_{1},#it{b}_{2}), Standard MV1 only;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_dPhib1b2_vs_b2_stdMV1","#Delta#phi(#it{b}_{1},#it{b}_{2}) vs Flavor weight #it{b}_{2}, Standard MV1 only;Flavor weight #it{b}_{2}, Standard MV1 only;#Delta#phi(#it{b}_{1},#it{b}_{2}), Standard MV1 only;Normalized",labels,40,0.,+1., 32,0.,TMath::Pi());
		addHist(histos2,name,"bjet_wb1_vs_wb2_stdMV1","Flavor weight #it{b}_{1} vs Flavor weight #it{b}_{2}, Standard MV1 only;Flavor weight #it{b}_{2}, Standard MV1 only;Flavor weight #it{b}_{1}, Standard MV1 only;Normalized",labels,40,0.,+1., 40,0.,+1.);
		addHist(histos2,name,"jet_wj1_vs_wj2_stdMV1","Flavor weight #it{jet}_{1} vs Flavor weight #it{jet}_{2}, Standard MV1 only;Flavor weight #it{jet}_{2}, Standard MV1 only;Flavor weight #it{jet}_{1}, Standard MV1 only;Normalized",labels,40,0.,+1., 40,0.,+1.);


		///////// vertexing 2d
		addHist(histos2,name,"m3mu_vs_dRmin",                "m_{3body} vs min[#DeltaR(3#mu,#mu_{i})];min[#DeltaR(3#mu,#mu_{i})];m_{3body} [MeV];Normalized",labels, 50,0.,0.25, 50,0.,4000.);
		addHist(histos2,name,"m3mu_vs_dRmax",                "m_{3body} vs max[#DeltaR(3#mu,#mu_{i})];max[#DeltaR(3#mu,#mu_{i})];m_{3body} [MeV];Normalized",labels, 50,0.,0.50, 50,0.,4000.);
		addHist(histos2,name,"m3mu_vs_dRmin_after_vtxclean", "m_{3body} vs min[#DeltaR(3#mu,#mu_{i})] after vtx cleaning;min[#DeltaR(3#mu,#mu_{i})] after vtx clean;m_{3body} [MeV] after vtx clean;Normalized",labels, 50,0.,0.25, 50,0.,4000.);
		addHist(histos2,name,"m3mu_vs_dRmax_after_vtxclean", "m_{3body} vs max[#DeltaR(3#mu,#mu_{i})] after vtx cleaning;max[#DeltaR(3#mu,#mu_{i})] after vtx clean;m_{3body} [MeV] after vtx clean;Normalized",labels, 50,0.,0.50, 50,0.,4000.);
		addHist(histos2,name,"m3mu_vs_dRmin_after_hadclean", "m_{3body} vs min[#DeltaR(3#mu,#mu_{i})] after hadronic cleaning;min[#DeltaR(3#mu,#mu_{i})] after had clean;m_{3body} [MeV] after isolation;Normalized",labels, 50,0.,0.25, 50,0.,4000.);
		addHist(histos2,name,"m3mu_vs_dRmax_after_hadclean", "m_{3body} vs max[#DeltaR(3#mu,#mu_{i})] after hadronic cleaning;max[#DeltaR(3#mu,#mu_{i})] after had clean;m_{3body} [MeV] after isolationafter isolation;Normalized",labels, 50,0.,0.50, 50,0.,4000.);

		addHist(histos2,name,"pvalVtx3_vs_pvalVtx2",    ";#it{p}-value (vtx 2#mu);#it{p}-value (vtx 3#mu);Normalized",labels, 30,0.,1. ,30,0.,1.);
		addHist(histos2,name,"pvalVtx3_vs_pvalVtx4",    ";#it{p}-value (vtx 3#mu+1trk);#it{p}-value (vtx 3#mu);Normalized",labels, 30,0.,1. ,30,0.,1.);
		addHist(histos2,name,"pvalVtx3_vs_pvalVtx4max", ";max[#it{p}-value] (vtx 3#mu+1trk);#it{p}-value (vtx 3#mu);Normalized",labels, 30,0.,1. ,30,0.,1.);

		/*The best variable to check if this is really 4-track vertex or 3-track vertex is to look on *chi2*  (not probability) of added 4th track.
		This can be a difference in chi2(4trk) and chi2(3trk) [no probabilities!] or it can be calculated chi2 of 4th track given fitted parameters of 3trk vertex.
		He said both values would give us the same number.
		chi2(4trk) = chi2(3trk)+chi2(added 4th track)
		so the value of merit here is difference between 2 chi2.
		It is important that it is not normalized on degrees of freedom [well he says one can normalize difference on n.d.o.f of 4th track which is 2, but after differences in chi2s is taken]
		So the best would be to plot chi2(3track) vs chi2(added 4th track) on 2d plot and see if there is any clear separation between background and signal.
		*/

		addHist(histos2,name,"type1_tau_vs_pval", ";#it{p}-value;#tau [??] (PV type 1);Normalized",labels,        20,0.,1. ,20,-5.,+10.);
		addHist(histos2,name,"type3_tau_vs_pval", ";#it{p}-value;#tau [??] (PV type 3);Normalized",labels,        20,0.,1. ,20,-5.,+10.);
		addHist(histos2,name,"type1_lxy_vs_pval", ";#it{p}-value;#it{L}_{xy} [??] (PV type 1);Normalized",labels, 20,0.,1. ,20,-10.,+40.);
		addHist(histos2,name,"type3_lxy_vs_pval", ";#it{p}-value;#it{L}_{xy} [??] (PV type 3);Normalized",labels, 20,0.,1. ,20,-10.,+40.);

		addHist(histos2,name,"mOS_vs_vtxSigDiffT", ";(T_{3#mu}-T_{2#mu})/#Delta_{T};m_{OS};Normalized",labels, 30,  0.,+10. ,30,0.,3500.);
		addHist(histos2,name,"mOS_vs_vtxSigDiffx", ";(x_{3#mu}-x_{2#mu})/#Delta_{x};m_{OS};Normalized",labels, 30,-10.,+10. ,30,0.,3500.);
		addHist(histos2,name,"mOS_vs_vtxSigDiffy", ";(y_{3#mu}-y_{2#mu})/#Delta_{y};m_{OS};Normalized",labels, 30,-10.,+10. ,30,0.,3500.);
		addHist(histos2,name,"mOS_vs_vtxSigDiffz", ";(z_{3#mu}-z_{2#mu})/#Delta_{z};m_{OS};Normalized",labels, 30,-10.,+10. ,30,0.,3500.);		

		addHist(histos2,name,"mupt1Fraction_vs_mupt1", "p_{T}(#mu1) fraction vs p_{T}(#mu1);p_{T}(#mu1) [MeV];p_{T}(#mu1) fraction;Normalized",labels, 50,2000.,52000. ,50,0.3,1.00);		
		addHist(histos2,name,"mupt2Fraction_vs_mupt2", "p_{T}(#mu2) fraction vs p_{T}(#mu2);p_{T}(#mu2) [MeV];p_{T}(#mu2) fraction;Normalized",labels, 50,2000.,52000. ,50,0.0,0.50);		
		addHist(histos2,name,"mupt3Fraction_vs_mupt3", "p_{T}(#mu3) fraction vs p_{T}(#mu3);p_{T}(#mu3) [MeV];p_{T}(#mu3) fraction;Normalized",labels, 50,2000.,52000. ,50,0.0,0.35);
		
		addHist(histos2,name,"mupt1Fraction_vs_mupt1_afterCut", "p_{T}(#mu1) fraction vs p_{T}(#mu1) after cutting;p_{T}(#mu1) [MeV];p_{T}(#mu1) fraction;Normalized",labels, 50,2000.,52000. ,50,0.3,1.00);		
		addHist(histos2,name,"mupt2Fraction_vs_mupt2_afterCut", "p_{T}(#mu2) fraction vs p_{T}(#mu2) after cutting;p_{T}(#mu2) [MeV];p_{T}(#mu2) fraction;Normalized",labels, 50,2000.,52000. ,50,0.0,0.50);		
		addHist(histos2,name,"mupt3Fraction_vs_mupt3_afterCut", "p_{T}(#mu3) fraction vs p_{T}(#mu3) after cutting;p_{T}(#mu3) [MeV];p_{T}(#mu3) fraction;Normalized",labels, 50,2000.,52000. ,50,0.0,0.35);

		addHist(histos2,name,"mupt12Fraction_vs_mupt12", "#Deltap_{T}(#mu12) fraction vs #Deltap_{T}(#mu12);#Deltap_{T}(#mu12) [MeV];#Deltap_{T}(#mu12) fraction;Normalized",labels, 50,2000.,52000. ,50,0.0,1.0);		
		addHist(histos2,name,"mupt23Fraction_vs_mupt23", "#Deltap_{T}(#mu23) fraction vs #Deltap_{T}(#mu23);#Deltap_{T}(#mu23) [MeV];#Deltap_{T}(#mu23) fraction;Normalized",labels, 50,2000.,52000. ,50,0.0,1.0);		
		addHist(histos2,name,"mupt13Fraction_vs_mupt13", "#Deltap_{T}(#mu13) fraction vs #Deltap_{T}(#mu13);#Deltap_{T}(#mu13) [MeV];#Deltap_{T}(#mu13) fraction;Normalized",labels, 50,2000.,52000. ,50,0.0,1.0);
		
		addHist(histos2,name,"mupt23Fraction_vs_mupt12Fraction", "#Deltap_{T}(#mu23) fraction vs #Deltap_{T}(#mu12) fraction;#Deltap_{T}(#mu12) fraction;#Deltap_{T}(#mu23) fraction;Normalized",labels, 50,0.0,1.0 ,50,0.0,1.0);		
		addHist(histos2,name,"mupt13Fraction_vs_mupt12Fraction", "#Deltap_{T}(#mu13) fraction vs #Deltap_{T}(#mu12) fraction;#Deltap_{T}(#mu12) fraction;#Deltap_{T}(#mu13) fraction;Normalized",labels, 50,0.0,1.0 ,50,0.0,1.0);		
		
		addHist(histos2,name,"mupt1_type_vs_mupt1", "#mu1 type vs p_{T}(#mu1);p_{T}(#mu 1) [MeV];;Normalized",labels, 40,0.,40000., 4,0,4);
		addHist(histos2,name,"mupt2_type_vs_mupt2", "#mu2 type vs p_{T}(#mu2);p_{T}(#mu 2) [MeV];;Normalized",labels, 40,0.,40000., 4,0,4);
		addHist(histos2,name,"mupt3_type_vs_mupt3", "#mu3 type vs p_{T}(#mu3);p_{T}(#mu 3) [MeV];;Normalized",labels, 40,0.,40000., 4,0,4);
		
		addHist(histos,name,"mupt1Fraction", ";p_{T}(#mu1) fraction;Normalized",labels, 50,0.3,1.00);		
		addHist(histos,name,"mupt2Fraction", ";p_{T}(#mu2) fraction;Normalized",labels, 50,0.0,0.50);		
		addHist(histos,name,"mupt3Fraction", ";p_{T}(#mu3) fraction;Normalized",labels, 50,0.0,0.35);
		
		addHist(histos,name,"mupt12Fraction", ";#Deltap_{T}(#mu12) fraction;Normalized",labels, 50,0.3,1.00);		
		addHist(histos,name,"mupt23Fraction", ";#Deltap_{T}(#mu23) fraction;Normalized",labels, 50,0.0,0.50);		
		addHist(histos,name,"mupt13Fraction", ";#Deltap_{T}(#mu13) fraction;Normalized",labels, 50,0.0,0.35);
		
		addHist(histos,name,"mupt1Fraction_afterCut", ";p_{T}(#mu1) fraction;Normalized",labels, 50,0.3,1.00);		
		addHist(histos,name,"mupt2Fraction_afterCut", ";p_{T}(#mu2) fraction;Normalized",labels, 50,0.0,0.50);		
		addHist(histos,name,"mupt3Fraction_afterCut", ";p_{T}(#mu3) fraction;Normalized",labels, 50,0.0,0.35);
		
	
		///////// vertexing 1d
		addHist(histos,name,"vtxSigDiffT", ";(T_{3#mu}-T_{2#mu})/#Delta_{T};Normalized",labels, 40,  0.,+10.);
		addHist(histos,name,"vtxSigDiffx", ";(x_{3#mu}-x_{2#mu})/#Delta_{x};Normalized",labels, 40,-10.,+10.);
		addHist(histos,name,"vtxSigDiffy", ";(y_{3#mu}-y_{2#mu})/#Delta_{y};Normalized",labels, 40,-10.,+10.);
		addHist(histos,name,"vtxSigDiffz", ";(z_{3#mu}-z_{2#mu})/#Delta_{z};Normalized",labels, 40,-10.,+10.);
		
		addHist(histos,name,"pT3mu_before_triplet","p_{T}^{3body} (before triplet cuts);p_{T}^{3body} [MeV] (before triplet cuts);Normalized",labels,80,0.,100000.,false);


		addHist(histos2,name,"m3mu_vs_mOS_rhoomegaphi",   "On #rho/#omega, #phi m_{3body} vs m_{OS};m_{OS} [MeV];m_{3body} [MeV];Events",labels, 100,680.,1100., 100,600.,3500.);
		addHist(histos2,name,"pT3mu_vs_mOS_rhoomegaphi",  "On #rho/#omega, #phi p_{T}^{3body} vs m_{OS};m_{OS} [MeV];p_{T}^{3body} [MeV];Events",labels, 100,680.,1100., 50,10000.,60000.);
		addHist(histos2,name,"mSS_vs_mOS_rhoomegaphi",    "On #rho/#omega, #phi m_{SS} vs m_{OS};m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,680.,1100., 100,200.,3000.);		
		addHist(histos,name,"m3mu_rhoomegaphi",           "On #rho/#omega, #phi;m_{3body} [MeV];Events",labels, 50,600.,3500.,false);
		addHist(histos,name,"pT3mu_rhoomegaphi",          "On #rho/#omega, #phi;p_{T}^{3body} [MeV];Normalized",labels,50,10000.,60000.,false);
		addHist(histos,name,"mOS_rhoomegaphi",            "On #rho/#omega, #phi;m_{OS} [MeV];Events",labels, 100,680.,1100.,false);
		addHist(histos,name,"mSS_rhoomegaphi",            "On #rho/#omega, #phi;m_{SS} [MeV];Normalized",labels, 100,200.,3000.,false);
		addHist(histos,name,"pvalue_rhoomegaphi",         "On #rho/#omega, #phi;#it{p}-value;Normalized",labels, 200,0.,1.,false);
		addHist(histos,name,"mu_chi2ndf_failMedium_rhoomegaphi", "On #rho/#omega, #phi (fail isMedium);#mu track-fit #chi^{2}/N_{DOF} (fail isMedium);Normalized",labels, 100,0.,10.,false);
		addHist(histos,name,"mu_chi2ndf_passMedium_rhoomegaphi", "On #rho/#omega, #phi (pass isMedium);#mu track-fit #chi^{2}/N_{DOF} (pass isMedium);Normalized",labels, 100,0.,10.,false);
		addHist(histos,name,"mu_pvalue_failMedium_rhoomegaphi",  "On #rho/#omega, #phi (fail isMedium);#mu track-fit p-value (fail isMedium);Normalized",labels, 100,0.,1.,false);
		addHist(histos,name,"mu_pvalue_zoom_failMedium_rhoomegaphi",  "On #rho/#omega, #phi (fail isMedium);#mu track-fit p-value (fail isMedium);Normalized",labels, 100,0.,2.e-5,false);
		addHist(histos,name,"mu_pvalue_passMedium_rhoomegaphi",  "On #rho/#omega, #phi (pass isMedium);#mu track-fit p-value (pass isMedium);Normalized",labels, 100,0.,1.,false);
		addHist(histos,name,"mu_pvalue_zoom_passMedium_rhoomegaphi",  "On #rho/#omega, #phi (pass isMedium);#mu track-fit p-value (pass isMedium);Normalized",labels, 100,0.,2.e-5,false);
		addHist(histos,name,"tpmu_chi2ndf_rhoomegaphi",   "On #rho/#omega, #phi;TP#mu track-fit #chi^{2}/N_{DOF};Normalized",labels, 100,0.,10.,false);
		addHist(histos,name,"tpmu_pvalue_rhoomegaphi",    "On #rho/#omega, #phi;TP#mu track-fit p-value;Normalized",labels, 100,0.,1.,false);
		addHist(histos,name,"tpmu_pvalue_zoom_rhoomegaphi",    "On #rho/#omega, #phi;TP#mu track-fit p-value;Normalized",labels, 100,0.,2.e-5,false);
		
		addHist(histos2,name,"m3mu_vs_mOS_Jpsi",   "On J/#psi m_{3body} vs m_{OS};m_{OS} [MeV];m_{3body} [MeV];Events",labels, 100,2900.,3300., 100,3000.,4000.);
		addHist(histos2,name,"pT3mu_vs_mOS_Jpsi",  "On J/#psi p_{T}^{3body} vs m_{OS};m_{OS} [MeV];p_{T}^{3body} [MeV];Events",labels, 100,2900.,3300., 50,10000.,60000.);
		addHist(histos2,name,"mSS_vs_mOS_Jpsi",    "On J/#psi m_{SS} vs m_{OS};m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,2900.,3300., 100,200.,3000.);		
		addHist(histos,name,"m3mu_Jpsi",           "On J/#psi;m_{3body} [MeV];Events",labels, 100,3000.,4000.,false);
		addHist(histos,name,"pT3mu_Jpsi",          "On J/#psi;p_{T}^{3body} [MeV];Normalized",labels,50,10000.,60000.,false);
		addHist(histos,name,"mOS_Jpsi",            "On J/#psi;m_{OS} [MeV];Events",labels, 100,2900.,3300.,false);
		addHist(histos,name,"mSS_Jpsi",            "On J/#psi;m_{SS} [MeV];Normalized",labels, 100,200.,3000.,false);
		addHist(histos,name,"pvalue_Jpsi",         "On J/#psi;#it{p}-value;Normalized",labels, 200,0.,1.,false);
		addHist(histos,name,"mu_chi2ndf_failMedium_Jpsi","On J/#psi (fail isMedium);#mu track-fit #chi^{2}/N_{DOF} (fail isMedium);Normalized",labels, 100,0.,10.,false);
		addHist(histos,name,"mu_chi2ndf_passMedium_Jpsi","On J/#psi (pass isMedium);#mu track-fit #chi^{2}/N_{DOF} (pass isMedium);Normalized",labels, 100,0.,10.,false);
		addHist(histos,name,"mu_pvalue_failMedium_Jpsi", "On J/#psi (fail isMedium);#mu track-fit p-value (fail isMedium);Normalized",labels, 100,0.,1.,false);
		addHist(histos,name,"mu_pvalue_zoom_failMedium_Jpsi", "On J/#psi (fail isMedium);#mu track-fit p-value (fail isMedium);Normalized",labels, 100,0.,2.e-5,false);
		addHist(histos,name,"mu_pvalue_passMedium_Jpsi", "On J/#psi (pass isMedium);#mu track-fit p-value (pass isMedium);Normalized",labels, 100,0.,1.,false);
		addHist(histos,name,"mu_pvalue_zoom_passMedium_Jpsi", "On J/#psi (pass isMedium);#mu track-fit p-value (pass isMedium);Normalized",labels, 100,0.,2.e-5,false);
		addHist(histos,name,"tpmu_chi2ndf_Jpsi",   "On J/#psi;TP#mu track-fit #chi^{2}/N_{DOF};Normalized",labels, 100,0.,10.,false);
		addHist(histos,name,"tpmu_pvalue_Jpsi",    "On J/#psi;TP#mu track-fit p-value;Normalized",labels, 100,0.,1.,false);
		addHist(histos,name,"tpmu_pvalue_zoom_Jpsi",    "On J/#psi;TP#mu track-fit p-value;Normalized",labels, 100,0.,2.e-5,false);
		
		

		addHist(histos,name,"m3mu",           ";m_{3body} [MeV] (skim levle);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin",       ";m_{3body} [MeV] (skim level);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom",  ";m_{3body} [MeV] (skim level);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion", ";m_{3body} [MeV] (skim level);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu",          ";p_{T}^{3body} [MeV] (skim level);Normalized",labels,80,0.,100000.);
		addHist(histos,name,"mOS",            ";m_{OS} [MeV] (skim level);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS",            ";m_{SS} [MeV] (skim level);Normalized",labels, 100,0.,4000.);
		
		addHist(histos,name,"m3mu_before_MVA",           ";m_{3body} [MeV] (before MVA-score cut);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_before_MVA",       ";m_{3body} [MeV] (before MVA-score cut);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_before_MVA",  ";m_{3body} [MeV] (before MVA-score cut);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_before_MVA", ";m_{3body} [MeV] (before MVA-score cut);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_before_MVA",          ";p_{T}^{3body} [MeV] (before MVA-score cut);Normalized",labels,80,0.,100000.);
		addHist(histos,name,"mOS_before_MVA",            ";m_{OS} [MeV] (before MVA-score cut);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_before_MVA",            ";m_{SS} [MeV] (before MVA-score cut);Normalized",labels, 100,0.,4000.);
		
		addHist(histos,name,"m3mu_after_MVA",           ";m_{3body} [MeV] (after MVA-score cut);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_MVA",       ";m_{3body} [MeV] (after MVA-score cut);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_MVA",  ";m_{3body} [MeV] (after MVA-score cut);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_MVA", ";m_{3body} [MeV] (after MVA-score cut);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_MVA",          ";p_{T}^{3body} [MeV] (after MVA-score cut);Normalized",labels,80,0.,100000.);
		addHist(histos,name,"mOS_after_MVA",            ";m_{OS} [MeV] (after MVA-score cut);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_after_MVA",            ";m_{SS} [MeV] (after MVA-score cut);Normalized",labels, 100,0.,4000.);

		addHist(histos,name,"m3mu_after_muons",           ";m_{3body} [MeV] (after muons cuts);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_muons",       ";m_{3body} [MeV] (after muons cuts);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_muons",  ";m_{3body} [MeV] (after muons cuts);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_muons", ";m_{3body} [MeV] (after muons cuts);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_muons",          ";p_{T}^{3body} [MeV] (after muons cuts);Normalized",labels, 80,0.,100000.);
		addHist(histos,name,"mOS_after_muons",            ";m_{OS} [MeV] (after muons cuts);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_after_muons",            ";m_{SS} [MeV] (after muons cuts);Normalized",labels, 100,0.,4000.);
		
		addHist(histos,name,"m3mu_after_triplet",           ";m_{3body} [MeV] (after triplet cuts);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_triplet",       ";m_{3body} [MeV] (after triplet cuts);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_triplet",  ";m_{3body} [MeV] (after triplet cuts);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_triplet", ";m_{3body} [MeV] (after triplet cuts);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_triplet",          ";p_{T}^{3body} [MeV] (after triplet cuts);Normalized",labels, 80,0.,100000.);
		addHist(histos,name,"mOS_after_triplet",            ";m_{OS} [MeV] (after triplet cuts);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_after_triplet",            ";m_{SS} [MeV] (after triplet cuts);Normalized",labels, 100,0.,4000.);

		addHist(histos,name,"m3mu_after_hadclean",           ";m_{3body} [MeV] (after hadronic cleaning);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_hadclean",       ";m_{3body} [MeV] (after hadronic cleaning);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_hadclean",  ";m_{3body} [MeV] (after hadronic cleaning);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_hadclean", ";m_{3body} [MeV] (after hadronic cleaning);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_hadclean",          ";p_{T}^{3body} [MeV] (after hadronic cleaning);Normalized",labels, 80,0.,100000.);
		addHist(histos,name,"mOS_after_hadclean",            ";m_{OS} [MeV] (after hadronic cleaning);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_after_hadclean",            ";m_{SS} [MeV] (after hadronic cleaning);Normalized",labels, 100,0.,4000.);
		
		addHist(histos,name,"m3mu_after_vtxmet",           ";m_{3body} [MeV] (after vtx+met cuts);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_vtxmet",       ";m_{3body} [MeV] (after vtx+met cuts);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_vtxmet",  ";m_{3body} [MeV] (after vtx+met cuts);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_vtxmet", ";m_{3body} [MeV] (after vtx+met cuts);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_vtxmet",          ";p_{T}^{3body} [MeV] (after vtx+met cuts);Normalized",labels, 80,0.,100000.);
		addHist(histos,name,"mOS_after_vtxmet",            ";m_{OS} [MeV] (after vtx+met cuts);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_after_vtxmet",            ";m_{SS} [MeV] (after vtx+met cuts);Normalized",labels, 100,0.,4000.);
		

		addHist(histos,name,"m3mu_after_pval",           ";m_{3body} [MeV] (after vtx #it{p}-value);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_pval",       ";m_{3body} [MeV] (after vtx #it{p}-value);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_pval",  ";m_{3body} [MeV] (after vtx #it{p}-value);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_pval", ";m_{3body} [MeV] (after vtx #it{p}-value);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_pval",          ";p_{T}^{3body} [MeV] (after vtx #it{p}-value);Normalized",labels, 80,0.,100000.);
		addHist(histos,name,"mOS_after_pval",            ";m_{OS} [MeV] (after vtx #it{p}-value);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_after_pval",            ";m_{SS} [MeV] (after vtx #it{p}-value);Normalized",labels, 100,0.,4000.);
		
		addHist(histos,name,"m3mu_after_vtxclean",           ";m_{3body} (after vtx cleaning);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_vtxclean",       ";m_{3body} (after vtx cleaning);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_vtxclean",  ";m_{3body} (after vtx cleaning);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_vtxclean", ";m_{3body} (after vtx cleaning);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_vtxclean",          ";p_{T}^{3body} (after vtx cleaning);Normalized",labels, 80,0.,100000.);
		addHist(histos,name,"mOS_after_vtxclean",            ";m_{OS} (after vtx cleaning);Events",labels, 100,0.,4000.);
		addHist(histos,name,"mSS_after_vtxclean",            ";m_{SS} (after vtx cleaning);Normalized",labels, 100,0.,4000.);
		
		addHist(histos,name,"m3mu_before_W",           ";m_{3body} [MeV] (before #it{W} cuts);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_before_W",       ";m_{3body} [MeV] (before #it{W} cuts);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_before_W",  ";m_{3body} [MeV] (before #it{W} cuts);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_before_W", ";m_{3body} [MeV] (before #it{W} cuts);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_before_W",          ";p_{T}^{3body} [MeV] (before #it{W} cuts);Normalized",labels, 80,0.,100000.);
		
		addHist(histos,name,"m3mu_lin_zoom_ratio_MET",  ";m_{3body} [MeV];R=before/after MET",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_ratio_MET", ";m_{3body} [MeV];R=before/after MET",labels, 50,1300.,2300.);
		addHist(histos,name,"m3mu_after_MET",           ";m_{3body} [MeV] (after MET);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_MET",       ";m_{3body} [MeV] (after MET);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_MET",  ";m_{3body} [MeV] (after MET);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_MET", ";m_{3body} [MeV] (after MET);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_MET",          ";p_{T}^{3body} [MeV] (after MET);Normalized",labels, 80,0.,100000.);
		
		addHist(histos,name,"m3mu_lin_zoom_ratio_dphi3muMET",  ";m_{3body} [MeV];R=before/after #Delta#phi(3#mu,MET)",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_ratio_dphi3muMET", ";m_{3body} [MeV];R=before/after #Delta#phi(3#mu,MET)",labels, 50,1300.,2300.);
		addHist(histos,name,"m3mu_after_dphi3muMET",           ";m_{3body} [MeV] (after #Delta#phi(3#mu,MET));Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_dphi3muMET",       ";m_{3body} [MeV] (after #Delta#phi(3#mu,MET));Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_dphi3muMET",  ";m_{3body} [MeV] (after #Delta#phi(3#mu,MET));Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_dphi3muMET", ";m_{3body} [MeV] (after #Delta#phi(3#mu,MET));Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_dphi3muMET",          ";p_{T}^{3body} [MeV] (after #Delta#phi(3#mu,MET));Normalized",labels, 80,0.,100000.);
		
		addHist(histos,name,"m3mu_lin_zoom_ratio_mT",  ";m_{3body} [MeV];R=before/after m_{T}",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_ratio_mT", ";m_{3body} [MeV];R=before/after m_{T}",labels, 50,1300.,2300.);
		addHist(histos,name,"m3mu_after_mT",           ";m_{3body} [MeV] (after m_{T});Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_mT",       ";m_{3body} [MeV] (after m_{T});Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_mT",  ";m_{3body} [MeV] (after m_{T});Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_mT", ";m_{3body} [MeV] (after m_{T});Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_mT",          ";p_{T}^{3body} [MeV] (after m_{T});Normalized",labels, 80,0.,100000.);
				
		addHist(histos,name,"m3mu_lin_zoom_ratio_isolation",  ";m_{3body} [MeV];R=before/after isolation",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_ratio_isolation", ";m_{3body} [MeV];R=before/after isolation",labels, 50,1300.,2300.);
		addHist(histos,name,"m3mu_after_isolation",           ";m_{3body} [MeV] (after isolation);Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_isolation",       ";m_{3body} [MeV] (after isolation);Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_isolation",  ";m_{3body} [MeV] (after isolation);Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_isolation", ";m_{3body} [MeV] (after isolation);Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_isolation",          ";p_{T}^{3body} [MeV] (after isolation);Normalized",labels, 80,0.,100000.);

		addHist(histos,name,"m3mu_lin_zoom_ratio_dRmax",  ";m_{3body} [MeV];R=before/after #DeltaR_{max}",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_ratio_dRmax", ";m_{3body} [MeV];R=before/after #DeltaR_{max}",labels, 50,1300.,2300.);
		addHist(histos,name,"m3mu_after_dRmax",           ";m_{3body} [MeV] (after #DeltaR_{max});Events",labels, 50,0.,4000.);
		addHist(histos,name,"m3mu_lin_after_dRmax",       ";m_{3body} [MeV] (after #DeltaR_{max});Events",labels, 50,500.,4000.);
		addHist(histos,name,"m3mu_lin_zoom_after_dRmax",  ";m_{3body} [MeV] (after #DeltaR_{max});Events",labels, 100,800.,2800.);
		addHist(histos,name,"m3mu_sigregion_after_dRmax", ";m_{3body} [MeV] (after #DeltaR_{max});Events",labels, 50,1300.,2300.);
		addHist(histos,name,"pT3mu_after_dRmax",          ";p_{T}^{3body} [MeV] (after #DeltaR_{max});Normalized",labels, 80,0.,100000.);
		
		
		
		addHist(histos,name,"chi2",           ";#chi^{2};Normalized",labels, 50,0.,50.);
		addHist(histos,name,"pvalue",         ";#it{p}-value (removed candidates consistent with 0);Normalized",labels, 200,0.,1.);
		addHist(histos,name,"pvalue_all",     ";#it{p}-value (all);Normalized",labels, 200,0.,1.);
		addHist(histos,name,"type1_lxy",      ";#it{L}_{xy} [??] (PV type 1);Normalized",labels, 60,-3.,+23.);
		addHist(histos,name,"type3_lxy",      ";#it{L}_{xy} [??] (PV type 3);Normalized",labels, 60,-3.,+23.);
		addHist(histos,name,"type1_tau",      ";#tau [??] (PV type 1);Normalized",labels, 60,-1.,+5.);
		addHist(histos,name,"type3_tau",      ";#tau [??] (PV type 3);Normalized",labels, 60,-1.,+5.);
		addHist(histos,name,"type1_a0",       ";#it{a}_{0} [??];Normalized",labels, 100,0.,+150.);
		addHist(histos,name,"type1_a0XY",     ";#it{a}_{0}^{xy} [??];Normalized",labels, 100,0.,+5.);
		addHist(histos,name,"type1_cosTh",    ";cos#theta;Normalized",labels, 100,-1.,+1.);
		addHist(histos,name,"type1_cosThXY",  ";cos#theta_{xy};Normalized",labels, 100,.0,+1.);
		addHist(histos,name,"type3_a0",       ";#it{a}_{0} [??];Normalized",labels, 100,0.,+300.);
		addHist(histos,name,"type3_a0XY",     ";#it{a}_{0}^{xy} [??];Normalized",labels, 100,0.,+7.);
		addHist(histos,name,"type3_cosTh",    ";cos#theta;Normalized",labels, 100,-1.,+1.);
		addHist(histos,name,"type3_cosThXY",  ";cos#theta_{xy};Normalized",labels, 100,-1.,+1.);
		
		addHist(histos,name,"vtx_trumatch_x", ";(x^{tru}-x^{fit})/#sigma_{x}^{fit};Normalized",labels, 40,-5.,+5.);
		addHist(histos,name,"vtx_trumatch_y", ";(y^{tru}-y^{fit})/#sigma_{y}^{fit};Normalized",labels, 40,-5.,+5.);
		addHist(histos,name,"vtx_trumatch_z", ";(z^{tru}-z^{fit})/#sigma_{z}^{fit};Normalized",labels, 40,-5.,+5.);

		addHist(histos,name,"actualInteractionsPerXing", ";<#mu>;Normalized", labels,50,0.,45.);
		addHist(histos,name,"vtx_npv",        ";n_{PV} (from vtx);Normalized",labels, 35,0,35);
		addHist(histos,name,"vtx_npv_type1",  ";n_{PV(type1)} (from vtx);Normalized",labels, 5,0,5);
		addHist(histos,name,"vtx_npv_type3",  ";n_{PV(type3)} (from vtx);Normalized",labels, 35,0,35);
		addHist(histos,name,"pv_npv",         ";n_{PV} (from PV);Normalized",labels, 35,0,35);
		addHist(histos,name,"pv_ntrks",       ";PV n_{trk} (from PV);Normalized",labels, 100,0,200);
		addHist(histos,name,"pv_npv_type1",   ";n_{PV(type1)} (from PV);Normalized",labels, 5,0,5);
		addHist(histos,name,"pv_npv_type3",   ";n_{PV(type3)} (from PV);Normalized",labels, 35,0,35);
		addHist(histos,name,"pv_type",        ";type_{PV} (from PV);Normalized",labels, 5,0,5);

		addHist(histos2,name,"mu_chi2_vs_sctangsig",   ";#mu Scattering curvature significance;Track-fit #chi^{2};Normalized",labels,     60,-20.,30., 200,0.,400.);
		addHist(histos2,name,"mu_pvalue_vs_sctangsig", ";#mu Scattering curvature significance;Track-fit #it{p}-value;Normalized",labels, 60,-20.,30., 100,0.,1.);

		addHist(histos,name,"mu_sctangsig",   ";#mu Scattering curvature significance;Normalized",labels,  60,-20.,30.);
		addHist(histos,name,"mu_sctngbsig",   ";#mu Scattering neighbours significance;Normalized",labels, 60,-20.,30.);
		addHist(histos,name,"mu_pbalsig",     ";#mu Momentum balance significance;Normalized",labels,      50,-10.,20.);
		addHist(histos,name,"tpmu_pbalsig",   ";TP#mu Momentum balance significance;Normalized",labels,    50,-10.,20.);
		

		addHist(histos,name,"mu_matchchi2ndf",";#mu Match #chi^{2}_{DOF};Normalized", labels,100,0.,20.);
		addHist(histos,name,"mu_qoverp_me",   ";#mu MS extrapolated Q/p [1/MeV];Normalized", labels,100,-7.e-4,+7.e-4);
		addHist(histos,name,"mu_qoverp_ie",   ";#mu ID extrapolated Q/p [1/MeV];Normalized", labels,100,-7.e-4,+7.e-4);
		addHist(histos,name,"mu_chi2",        ";#mu track-fit #chi^{2};Normalized",labels, 200,0.,400.);
		addHist(histos,name,"mu_chi2_zoom",   ";#mu track-fit #chi^{2};Normalized",labels, 200,0.,150.);
		addHist(histos,name,"mu_ndf",         ";#mu track-fit NDF;Normalized",labels, 120,0.,120.);
		addHist(histos,name,"mu_chi2ndf",     ";#mu track-fit #chi^{2}/N_{DOF};Normalized",labels, 100,0.,10.);
		addHist(histos,name,"mu_pvalue",      ";#mu track-fit p-value;Normalized",labels, 100,0.,1.);
		addHist(histos,name,"mu_pvalue_zoom", ";#mu track-fit p-value;Normalized",labels, 100,0.,0.002);
		addHist(histos,name,"mu_isLoose",     ";#mu isLoose;Normalized",labels, 2,0.,2);
		// addHist(histos,name,"mu_isMedium",    ";#mu isMedium;Normalized",labels, 2,0.,2);
		addHist(histos,name,"mu_isTight",     ";#mu isTight;Normalized",labels, 2,0.,2);
		addHist(histos2,name,"mu_pvalue_vs_chi2ndf", ";#mu track-fit #chi^{2}/N_{DOF};#mu track-fit p-value;Normalized",labels, 100,0.,10., 100,0.,1.);	
		
		addHist(histos,name,"mu_chi2_ismedium",        ";#mu track-fit #chi^{2} (passing isMedium);Normalized",labels, 200,0.,400.);
		addHist(histos,name,"mu_chi2_zoom_ismedium",   ";#mu track-fit #chi^{2} (passing isMedium);Normalized",labels, 200,0.,150.);
		addHist(histos,name,"mu_ndf_ismedium",         ";#mu track-fit NDF (passing isMedium);Normalized",labels, 120,0.,120.);
		addHist(histos,name,"mu_chi2ndf_ismedium",     ";#mu track-fit #chi^{2}/N_{DOF} (passing isMedium);Normalized",labels, 100,0.,10.);
		addHist(histos,name,"mu_pvalue_ismedium",      ";#mu track-fit p-value (passing isMedium);Normalized",labels, 100,0.,1.);
		addHist(histos,name,"mu_pvalue_zoom_ismedium", ";#mu track-fit p-value (passing isMedium);Normalized",labels, 100,0.,2.e-5);
		addHist(histos2,name,"mu_pvalue_vs_chi2ndf_ismedium", ";#mu track-fit #chi^{2}/N_{DOF} (passing isMedium);#mu track-fit p-value (passing isMedium);Normalized",labels, 100,0.,10., 100,0.,1.);	

		addHist(histos,name,"mu_chi2ndf_failMedium",     ";#mu track-fit #chi^{2}/N_{DOF} (failing isMedium);Normalized",labels, 100,0.,10.);
		addHist(histos,name,"mu_pvalue_failMedium",      ";#mu track-fit p-value (failing isMedium);Normalized",labels, 100,0.,1.);
		addHist(histos,name,"mu_pvalue_zoom_failMedium", ";#mu track-fit p-value (failing isMedium);Normalized",labels, 100,0.,2.e-5);
		addHist(histos2,name,"mu_pvalue_vs_chi2ndf_failMedium", ";#mu track-fit #chi^{2}/N_{DOF} (passing isMedium);#mu track-fit p-value (passing isMedium);Normalized",labels, 100,0.,10., 100,0.,1.);
		addHist(histos,name,"mu_chi2ndf_passMedium",     ";#mu track-fit #chi^{2}/N_{DOF} (passing isMedium);Normalized",labels, 100,0.,10.);
		addHist(histos,name,"mu_pvalue_passMedium",      ";#mu track-fit p-value (passing isMedium);Normalized",labels, 100,0.,1.);
		addHist(histos,name,"mu_pvalue_zoom_passMedium", ";#mu track-fit p-value (passing isMedium);Normalized",labels, 100,0.,2.e-5);
		addHist(histos2,name,"mu_pvalue_vs_chi2ndf_passMedium", ";#mu track-fit #chi^{2}/N_{DOF} (passing isMedium);#mu track-fit p-value (passing isMedium);Normalized",labels, 100,0.,10., 100,0.,1.);
		
		addHist(histos,name,"tpmu_qoverp",      ";TP#mu Q/p [1/MeV];Normalized", labels,100,-7.e-4,+7.e-4);
		addHist(histos,name,"tpmu_chi2",        ";TP#mu track-fit #chi^{2};Normalized",labels, 200,0.,400.);
		addHist(histos,name,"tpmu_chi2_zoom",   ";TP#mu track-fit #chi^{2};Normalized",labels, 200,0.,150.);
		addHist(histos,name,"tpmu_ndf",         ";TP#mu track-fit NDF;Normalized",labels, 120,0.,120.);
		addHist(histos,name,"tpmu_chi2ndf",     ";TP#mu track-fit #chi^{2}/N_{DOF};Normalized",labels, 100,0.,10.);
		addHist(histos,name,"tpmu_pvalue",      ";TP#mu track-fit p-value;Normalized",labels, 100,0.,1.);
		addHist(histos,name,"tpmu_pvalue_zoom", ";TP#mu track-fit p-value;Normalized",labels, 100,0.,2.e-5);
		addHist(histos2,name,"tpmu_pvalue_vs_chi2ndf", ";TP#mu track-fit #chi^{2}/N_{DOF};TP#mu track-fit p-value;Normalized",labels, 100,0.,10., 100,0.,1.);	
		
		addHist(histos,name,"tpA_qoverp",      ";TPa Q/p [1/MeV];Normalized", labels,100,-7.e-4,+7.e-4);
		addHist(histos,name,"tpA_chi2",        ";TPa track-fit #chi^{2};Normalized",labels, 100,0.,400.);
		addHist(histos,name,"tpA_chi2_zoom",   ";TPa track-fit #chi^{2};Normalized",labels, 150,0.,150.);
		addHist(histos,name,"tpA_chi2ndf",     ";TPa track-fit #chi^{2}/N_{DOF};Normalized",labels, 100,0.,10.);
		addHist(histos,name,"tpA_pvalue",      ";TPa track-fit p-value;Normalized",labels, 100,0.,1.);
		addHist(histos,name,"tpA_pvalue_zoom", ";TPa track-fit p-value;Normalized",labels, 100,0.,2.e-5);
		addHist(histos2,name,"tpA_pvalue_vs_chi2ndf", ";TPa track-fit #chi^{2}/N_{DOF};TPa track-fit p-value;Normalized",labels, 100,0.,10., 100,0.,1.);
		
		addHist(histos,name,"tpB_qoverp",      ";TPb Q/p [1/MeV];Normalized", labels,100,-7.e-4,+7.e-4);
		addHist(histos,name,"tpB_chi2",        ";TPb track-fit #chi^{2};Normalized",labels, 50,0.,25.);
		addHist(histos,name,"tpB_chi2_zoom",   ";TPb track-fit #chi^{2};Normalized",labels, 80,0.,20.);
		addHist(histos,name,"tpB_chi2ndf",     ";TPb track-fit #chi^{2}/N_{DOF};Normalized",labels, 100,0.,10.);
		addHist(histos,name,"tpB_pvalue",      ";TPb track-fit p-value;Normalized",labels, 100,0.,1.);
		addHist(histos,name,"tpB_pvalue_zoom", ";TPb track-fit p-value;Normalized",labels, 100,0.,2.e-5);
		addHist(histos2,name,"tpB_pvalue_vs_chi2ndf", ";TPb track-fit #chi^{2}/N_{DOF};TPb track-fit p-value;Normalized",labels, 100,0.,10., 100,0.,1.);		
		
		addHist(histos,name,"mu_isMedium",     ";All #mu isMedium;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isLoose1",     ";#mu1 isLoose;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isLoose2",     ";#mu2 isLoose;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isLoose3",     ";#mu3 isLoose;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isMedium1",    ";#mu1 isMedium;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isMedium2",    ";#mu2 isMedium;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isMedium3",    ";#mu3 isMedium;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isTight1",     ";#mu1 isTight;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isTight2",     ";#mu2 isTight;Normalized",labels, 2,0.,2,false);
		addHist(histos,name,"mu_isTight3",     ";#mu3 isTight;Normalized",labels, 2,0.,2,false);
		
		addHist(histos,name,"mu_eta1", ";#eta^{(#mu 1)};Normalized",labels,   64,-3.0,+3.0,false);
		addHist(histos,name,"mu_eta2", ";#eta^{(#mu 2)} ;Normalized",labels, 64,-3.0,+3.0,false);
		addHist(histos,name,"mu_eta3", ";#eta^{(#mu 3)};Normalized",labels,   64,-3.0,+3.0,false);
		addHist(histos,name,"mu_phi1", ";#phi^{(#mu 1)};Normalized",labels,   64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos,name,"mu_phi2", ";#phi^{(#mu 2)} ;Normalized",labels, 64,-TMath::Pi(),+TMath::Pi(),false);
		addHist(histos,name,"mu_phi3", ";#phi^{(#mu 3)};Normalized",labels,   64,-TMath::Pi(),+TMath::Pi(),false);
		
		addHist(histos,name,"mu_pt1",                   ";p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		addHist(histos,name,"mu_pt2",                   ";p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		addHist(histos,name,"mu_pt3",                   ";p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		addHist(histos,name,"mu_pt1_3mu4T",             "EF_3mu4T: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		addHist(histos,name,"mu_pt2_3mu4T",             "EF_3mu4T: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		addHist(histos,name,"mu_pt3_3mu4T",             "EF_3mu4T: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		addHist(histos,name,"mu_pt1_2mu8_EFxe30_tclcw", "EF_2mu8_EFxe30_tclcw: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		addHist(histos,name,"mu_pt2_2mu8_EFxe30_tclcw", "EF_2mu8_EFxe30_tclcw: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		addHist(histos,name,"mu_pt3_2mu8_EFxe30_tclcw", "EF_2mu8_EFxe30_tclcw: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		// addHist(histos,name,"mu_pt1_mu24i_tight",       "EF_mu24i_tight: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		// addHist(histos,name,"mu_pt2_mu24i_tight",       "EF_mu24i_tight: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		// addHist(histos,name,"mu_pt3_mu24i_tight",       "EF_mu24i_tight: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		// addHist(histos,name,"mu_pt1_mu36_tight",        "EF_mu36_tight: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		// addHist(histos,name,"mu_pt2_mu36_tight",        "EF_mu36_tight: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		// addHist(histos,name,"mu_pt3_mu36_tight",        "EF_mu36_tight: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		addHist(histos,name,"mu_pt1_unique_3mu4T",             "EF_3mu4T unique: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		addHist(histos,name,"mu_pt2_unique_3mu4T",             "EF_3mu4T unique: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		addHist(histos,name,"mu_pt3_unique_3mu4T",             "EF_3mu4T unique: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		addHist(histos,name,"mu_pt1_unique_2mu8_EFxe30_tclcw", "EF_2mu8_EFxe30_tclcw unique: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		addHist(histos,name,"mu_pt2_unique_2mu8_EFxe30_tclcw", "EF_2mu8_EFxe30_tclcw unique: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		addHist(histos,name,"mu_pt3_unique_2mu8_EFxe30_tclcw", "EF_2mu8_EFxe30_tclcw unique: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		// addHist(histos,name,"mu_pt1_unique_mu24i_tight",       "EF_mu24i_tight unique: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		// addHist(histos,name,"mu_pt2_unique_mu24i_tight",       "EF_mu24i_tight unique: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		// addHist(histos,name,"mu_pt3_unique_mu24i_tight",       "EF_mu24i_tight unique: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		// addHist(histos,name,"mu_pt1_unique_mu36_tight",        "EF_mu36_tight unique: First #mu;p_{T}^{(#mu 1)} [MeV];Normalized",labels,   50,2000.,102000.,false);
		// addHist(histos,name,"mu_pt2_unique_mu36_tight",        "EF_mu36_tight unique: Second #mu;p_{T}^{(#mu 2)} [MeV];Normalized",labels,   50,2000.,42000.,false);
		// addHist(histos,name,"mu_pt3_unique_mu36_tight",        "EF_mu36_tight unique: Third #mu;p_{T}^{(#mu 3)} [MeV];Normalized",labels,   50,2000.,27000.,false);
		
			
		addHist(histos,name,"muOnly_pt", ";#mu p_{T} [MeV];Normalized",labels, 120,2000.,122000.);
		addHist(histos,name,"tpA_pt",    ";TPa p_{T} [MeV];Normalized",labels, 60,2000.,62000.);
		addHist(histos,name,"tpB_pt",    ";TPb p_{T} [MeV];Normalized",labels, 60,2000.,32000.);
		addHist(histos,name,"muOnly_eta", ";#mu #eta;Normalized",labels,   64,-3.0,+3.0,false);
		addHist(histos,name,"tpA_eta",    ";TPa #eta;Normalized",labels,   64,-3.0,+3.0,false);
		addHist(histos,name,"tpB_eta",    ";TPb #eta;Normalized",labels,   64,-3.0,+3.0,false);
		
		
		addHist(histos,name,"mu_pt",          ";p_{T}^{(#mu all)} [MeV];Normalized",labels, 120,2000.,122000.);
		addHist(histos,name,"mu_eta",         ";#eta^{(#mu all)};Normalized",labels,        40,-2.7,+2.7);
		addHist(histos,name,"mu_MCP",         ";#mu Pass ID hits;Normalized",labels, 2,0,2);
		addHist(histos,name,"mu_nHighThresholdTRTHits",  ";#mu nHighThresholdTRTHits;Normalized",labels, 20,0,20);
		addHist(histos,name,"tpA_nHighThresholdTRTHits", ";TPa nHighThresholdTRTHits;Normalized",labels, 20,0,20);
		addHist(histos,name,"tpB_nHighThresholdTRTHits", ";TPb nHighThresholdTRTHits;Normalized",labels, 20,0,20);
		addHist(histos,name,"mu_type",        ";#mu Type;Normalized",labels,         4,0,4);
		addHist(histos,name,"triplet_isolation",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max});Normalized",labels, 50,0.,15.);
		addHist(histos,name,"triplet_isolation_before_isolation",";#Sigmap_{T}^{trk}/p_{T}^{3body} (cone #DeltaR_{max}) before cutting;Normalized",labels, 50,0.,15.);
		addHist(histos,name,"mu_ptcone10.1",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.1);Normalized",labels, 50,0.,15.);
		addHist(histos,name,"mu_ptcone10.2",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.1);Normalized",labels, 50,0.,25.);
		addHist(histos,name,"mu_ptcone10.3",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.1);Normalized",labels, 50,0.,40.);
		addHist(histos,name,"mu_ptcone20.1",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.2);Normalized",labels, 50,0.,20.);
		addHist(histos,name,"mu_ptcone20.2",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.2);Normalized",labels, 50,0.,30.);
		addHist(histos,name,"mu_ptcone20.3",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.2);Normalized",labels, 50,0.,50.);
		addHist(histos,name,"mu_ptcone30.1",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.3);Normalized",labels, 50,0.,25.);
		addHist(histos,name,"mu_ptcone30.2",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.3);Normalized",labels, 50,0.,40.);
		addHist(histos,name,"mu_ptcone30.3",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.3);Normalized",labels, 50,0.,60.);
		addHist(histos,name,"mu_ptcone40.1",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 1} (cone 0.4);Normalized",labels, 50,0.,25.);
		addHist(histos,name,"mu_ptcone40.2",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 2} (cone 0.4);Normalized",labels, 50,0.,40.);
		addHist(histos,name,"mu_ptcone40.3",  ";Corrected #Sigmap_{T}^{trk}/p_{T}^{#mu 3} (cone 0.4);Normalized",labels, 50,0.,60.);

		addHist(histos,name,"triplet_dRmin_3muons",       ";min[#DeltaR(3#mu,#mu_{i})];Normalized",labels,  100,0.,0.25);
		addHist(histos,name,"triplet_dRmin_2muons_1tpmu", ";min[#DeltaR(3#mu,#mu_{i})];Normalized",labels,  100,0.,0.25);
		addHist(histos,name,"triplet_dRmin_2muons_1calo", ";min[#DeltaR(3#mu,#mu_{i})];Normalized",labels,  100,0.,0.25);
		addHist(histos,name,"triplet_dRmin_1muons_2tpmu", ";min[#DeltaR(3#mu,#mu_{i})];Normalized",labels,  100,0.,0.25);
		addHist(histos,name,"triplet_dRmin_0muons_3tpmu", ";min[#DeltaR(3#mu,#mu_{i})];Normalized",labels,  100,0.,0.25);		
		
		
		addHist(histos,name,"triplet_mQuad_muons", ";m(4body);Events",labels,  100,(91.-25.)*GeV2MeV,(91.+25.)*GeV2MeV);
		addHist(histos,name,"triplet_mQuad_norm_muons", ";m(4body);Normalized",labels,  100,(91.-25.)*GeV2MeV,(91.+25.)*GeV2MeV);
		addHist(histos,name,"triplet_mTrip_muons", ";m(3body);Events",labels,  100,0.,50000.);
		addHist(histos,name,"triplet_mTrip_norm_muons", ";m(3body);Normalized",labels,  100,0.,50000.);
		addHist(histos,name,"triplet_pTQuad_muons",  ";p_{T}(4body);Normalized",labels,  150,0.,150000.);
		addHist(histos,name,"triplet_pTQuad1_muons", ";p_{T}(#mu1);Normalized",labels,  100,0.,100000.);
		addHist(histos,name,"triplet_pTQuad2_muons", ";p_{T}(#mu2);Normalized",labels,  70,0.,70000.);
		addHist(histos,name,"triplet_pTQuad3_muons", ";p_{T}(#mu3);Normalized",labels,  50,0.,50000.);
		addHist(histos,name,"triplet_pTQuad4_muons", ";p_{T}(#mu4);Normalized",labels,  40,0.,40000.);
		
		
		
		addHist(histos,name,"triplet_dRmin",                  ";min[#DeltaR(3#mu,#mu_{i})];Normalized",labels,  100,0.,0.25);
		addHist(histos,name,"triplet_dRmax",                  ";max[#DeltaR(3#mu,#mu_{i})];Normalized",labels,  100,0.,0.50);
		addHist(histos,name,"triplet_dRmin_after_vtxclean",    ";min[#DeltaR(3#mu,#mu_{i})] after vtx clean;Normalized",labels,  100,0.,0.25);
		addHist(histos,name,"triplet_dRmax_after_vtxclean",    ";max[#DeltaR(3#mu,#mu_{i})] after vtx clean;Normalized",labels,  100,0.,0.50);
		addHist(histos,name,"triplet_dRmin_after_hadclean", ";min[#DeltaR(3#mu,#mu_{i})] after had clean;Normalized",labels,  100,0.,0.25);
		addHist(histos,name,"triplet_dRmax_after_hadclean", ";max[#DeltaR(3#mu,#mu_{i})] after had clean;Normalized",labels,  100,0.,0.50);
		
		addHist(histos,name,"cosThOSQ_lightMesons",";cos#theta_{OS}(Q);Normalized",labels, 50,-1.,+1.);
		addHist(histos2,name,"cosThOSQ_vs_mOS","cos#theta_{OS}(Q) vs m_{OS};m_{OS} [MeV];cos#theta_{OS}(Q);Normalized",labels, 100,0.,2000., 50,-1.,+1.);

		addHist(histos,name,"cosThOS_lightMesons",";cos#theta_{OS};Normalized",labels, 50,-1.,+1.);
		addHist(histos2,name,"cosThOS_vs_mOS","cos#theta_{OS} vs m_{OS};m_{OS} [MeV];cos#theta_{OS};Normalized",labels, 100,0.,2000., 50,-1.,+1.);
		
		addHist(histos,name,"cosThOStrk3_lightMesons",";cos#theta_{OS}(trk3);Normalized",labels, 50,-1.,+1.);
		addHist(histos2,name,"cosThOStrk3_vs_mOS","cos#theta_{OS}(trk3) vs m_{OS};m_{OS} [MeV];cos#theta_{OS}(trk3);Normalized",labels, 100,0.,2000., 50,-1.,+1.);

		addHist(histos,name,"dPhiOS_lightMesons",";#Delta#phi_{OS};Normalized",labels, 64,0.,TMath::Pi());
		addHist(histos2,name,"dPhiOS_vs_mOS","#Delta#phi_{OS} vs m_{OS};m_{OS} [MeV];#Delta#phi_{OS};Normalized",labels, 100,0.,2000., 50,0.,TMath::Pi());

		addHist(histos,name,"dmOS_lightMesons",";#Deltam_{OS};Normalized",labels, 100,0.,2000.);
		addHist(histos2,name,"dmOS_vs_mOS","#Deltam_{OS} vs m_{OS};m_{OS} [MeV];#Deltam_{OS};Normalized",labels, 100,0.,2000., 100,0.,2000.);
		
		addHist(histos,name,"triplet_cosTheta",";cos#theta(3#mu,#muQ);Normalized",labels, 50,-1.,+1.);
		
		addHist(histos2,name,"m3mu_vs_mOS",               "m_{3body} vs m_{OS} skim level;m_{OS} [MeV];m_{3body} [MeV];Events",labels, 100,0.,4000., 100,500.,3500.);
		addHist(histos2,name,"m3mu_vs_mOS_after_pval",    "m_{3body} vs m_{OS} after #it{p}-value;m_{OS} [MeV];m_{3body} [MeV] after #it{p}-value;Events",labels, 100,0.,4000., 100,500.,3500.);
		addHist(histos2,name,"m3mu_vs_mOS_after_vtxclean","m_{3body} vs m_{OS} after vtx cleaning;m_{OS} [MeV];m_{3body} [MeV] after vtx cleaning;Events",labels, 100,0.,4000., 100,500.,3500.);
		addHist(histos2,name,"m3mu_vs_mOS_after_triplet", "m_{3body} vs m_{OS} after triplet cuts;m_{OS} [MeV];m_{3body} [MeV] after triplet cuts;Events",labels, 100,0.,4000., 100,500.,3500.);
		addHist(histos2,name,"m3mu_vs_mOS_after_vtxmet",  "m_{3body} vs m_{OS} after vtx+met cuts;m_{OS} [MeV];m_{3body} [MeV] after vtx+met cuts;Events",labels, 100,0.,4000., 100,500.,3500.);
		addHist(histos2,name,"m3mu_vs_mOS_after_muons",   "m_{3body} vs m_{OS} after muons cuts;m_{OS} [MeV];m_{3body} [MeV] after muons cuts;Events",labels, 100,0.,4000., 100,500.,3500.);
		addHist(histos2,name,"m3mu_vs_mOS_after_hadclean","m_{3body} vs m_{OS} after hadronic cleaning;m_{OS} [MeV];m_{3body} [MeV] after hadronic cleaning;Events",labels, 100,0.,4000., 100,500.,3500.);

		addHist(histos2,name,"pT3mu_vs_mOS",               "p_{T}^{3body} vs m_{OS} skim level;m_{OS} [MeV];p_{T}^{3body} [MeV];Events",labels, 100,0.,4000., 80,0.,80000.);
		addHist(histos2,name,"pT3mu_vs_mOS_after_pval",    "p_{T}^{3body} vs m_{OS} after #it{p}-value;m_{OS} [MeV];p_{T}^{3body} [MeV] after #it{p}-value;Events",labels, 100,0.,4000., 80,0.,80000.);
		addHist(histos2,name,"pT3mu_vs_mOS_after_vtxclean","p_{T}^{3body} vs m_{OS} after vtx cleaning;m_{OS} [MeV];p_{T}^{3body} [MeV] after vtx cleaning;Events",labels, 100,0.,4000., 80,0.,80000.);
		addHist(histos2,name,"pT3mu_vs_mOS_after_triplet", "p_{T}^{3body} vs m_{OS} after triplet cuts;m_{OS} [MeV];p_{T}^{3body} [MeV] after triplet cuts;Events",labels, 100,0.,4000., 80,0.,80000.);
		addHist(histos2,name,"pT3mu_vs_mOS_after_vtxmet",  "p_{T}^{3body} vs m_{OS} after vtx+met cuts;m_{OS} [MeV];p_{T}^{3body} [MeV] after vtx+met cuts;Events",labels, 100,0.,4000., 80,0.,80000.);
		addHist(histos2,name,"pT3mu_vs_mOS_after_muons",   "p_{T}^{3body} vs m_{OS} after muons cuts;m_{OS} [MeV];p_{T}^{3body} [MeV] after muons cuts;Events",labels, 100,0.,4000., 80,0.,80000.);
		addHist(histos2,name,"pT3mu_vs_mOS_after_hadclean","p_{T}^{3body} vs m_{OS} after hadronic cleaning;m_{OS} [MeV];p_{T}^{3body} [MeV] after hadronic cleaning;Events",labels, 100,0.,4000., 80,0.,80000.);
		
		addHist(histos2,name,"mOS2_vs_mOS1",             "m_{OS2} vs m_{OS1} at 2nd skim level;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS",               "m_{SS} vs m_{OS} at 2nd skim level;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_beforeCut",   "m_{OS2} vs m_{OS1} before cutting;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_beforeCut",     "m_{SS} vs m_{OS} before cutting;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_afterCut",    "m_{OS2} vs m_{OS1} after cutting;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_afterCut",      "m_{SS} vs m_{OS} after cutting;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_after_pval",  "m_{OS2} vs m_{OS1} after #it{p}-value;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_after_pval",    "m_{SS} vs m_{OS} after #it{p}-value;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_after_vtxclean","m_{OS2} vs m_{OS1} after vtx cleaning;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_after_vtxclean","m_{SS} vs m_{OS} after vtx cleaning;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_after_vtxmet","m_{OS2} vs m_{OS1} after vtx+met cuts;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_after_vtxmet",  "m_{SS} vs m_{OS} after vtx+met cuts;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_after_triplet","m_{OS2} vs m_{OS1} after triplet cuts;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_after_triplet", "m_{SS} vs m_{OS} after triplet cuts;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_after_muons","m_{OS2} vs m_{OS1} after muons cuts;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_after_muons",   "m_{SS} vs m_{OS} after muons cuts;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mOS2_vs_mOS1_after_hadclean","m_{OS2} vs m_{OS1} after triplet cuts;m_{OS1} [MeV];m_{OS2} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		addHist(histos2,name,"mSS_vs_mOS_after_hadclean","m_{SS} vs m_{OS} after hadronic cleaning;m_{OS} [MeV];m_{SS} [MeV];Events",labels, 100,0.,4000., 100,0.,4000.);
		
		addHist(histos,name,"nTriplets", ";N_{triplets};Normalized",labels, 6,0,6);
	
		///////////////////////
		//// force Sumw2 //////
		///////////////////////
		for(TMapTSP2TH1::iterator ii=histos.begin()  ; ii!=histos.end()  ; ++ii) if(ii->first.Contains(name)) ii->second->Sumw2();
		for(TMapTSP2TH2::iterator ii=histos2.begin() ; ii!=histos2.end() ; ++ii) if(ii->first.Contains(name)) ii->second->Sumw2();
	
	}
	
	
	
	///////////////////////////
	///////// ANALYSIS ////////
	///////////////////////////
	for(TMapuiTS::iterator it=channels.begin() ; it!=channels.end() ; it++)
	{
		TString name = it->second;
		_INFO("Starting "+(string)name);
		
		////////////////////////////////////////////////////
		// ignored the dummy channes (sum bg or sum sig) ///
		if(isBGsum(name) || isSIGsum(name)) continue; //////
		////////////////////////////////////////////////////
		
		// initialization and clearance
		clearCounters();
		initCounters();
		clearCounters();
		
		
		// global flag to note if MC or Data
		glob_isMC = (!isData(name));
		glob_isSig = (isSignal(name));
		
		
		///////////////////////////////////////
		// initialize the flatout tree ////////
		if(!skim) initFlatoutTree(name); //////
		///////////////////////////////////////
		
		
		///////////////////////////
		// do the analysis ////////
		///////////////////////////
		
		TString mastertree = getMasterTreeName();
		if(!isBinned(name)) 
		{
			prepareChains(name,mastertree,chains,chainfriends,fout,otrees);
			if(isData(name)) analysis(name,chains,otrees,histos,histos2,nentriesMax); // limit chain only for data and only if nentriesMax!=0 (see up)
			else             analysis(name,chains,otrees,histos,histos2,0);
		}
		else
		{
			if(name.Contains("NpX") || name.Contains("JZxW"))
			{
				for(TMapTSb::iterator imc=binnedmcenable.begin() ; imc!=binnedmcenable.end() ; ++imc)
				{
					TString bname = imc->first;
					if(name.Contains("WmunuNpX")  && !bname.Contains("WmunuNp"))  continue;
					if(name.Contains("ZmumuNpX")  && !bname.Contains("ZmumuNp"))  continue;
					if(name.Contains("WtaunuNpX") && !bname.Contains("WtaunuNp")) continue;
					if(name.Contains("JZxW")      && !bname.Contains("JZ"))       continue;
					prepareChains(bname,mastertree,chains,chainfriends,fout,otrees);
					analysis(bname,chains,otrees,histos,histos2,0);
				}
			}
			else if(name.Contains("Data"))
			{
				for(TMapTSb::iterator ip=periodenable.begin() ; ip!=periodenable.end() ; ++ip)
				{
					TString period = ip->first;
					TString bname  = "period"+period;
					prepareChains(bname,mastertree,chains,chainfriends,fout,otrees);
					analysis(bname,chains,otrees,histos,histos2,0);
				}
			}
		}
		_INFO("Done "+(string)name+"\n\n");	
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	
	_INFO("");
	
	
	///////////////////////////////////////////////////////////////
	// reblind all mass histos 1d and 2d //////////////////////////
	///////////////////////////////////////////////////////////////
	reBlindAllMassHists(histos,  mBlindMinGlob, mBlindMaxGlob);
	reBlindAllMassHists(histos2, mBlindMinGlob, mBlindMaxGlob);



	_INFO("");


	//////////////////////////////////////////////////////////////////////////
	// sum up the signals and the backgrounds and signalse (separately). /////
	// this has to happen before normalizing !!! /////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	for(TMapTSP2TH1::iterator hit=histos.begin() ; hit!=histos.end() ; hit++)
	{	
		TString hname = hit->first;
		if(!isBGsum(hname) && !isSIGsum(hname)) continue;
		
		for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
		{
			TString name = cit->second;
			if(drawchannels[name]==0) continue;
			if(isData(name)) continue;
			
			TString hbarename = hname;
			if(isBGsum(hname)  && !isSignal(name)) { hbarename.ReplaceAll("Backgrounds",""); hit->second->Add(histos[name+hbarename]); }
			if(isSIGsum(hname) && isSignal(name))  { hbarename.ReplaceAll("Signals","");     hit->second->Add(histos[name+hbarename]); }
		}
	}
	for(TMapTSP2TH2::iterator hit=histos2.begin() ; hit!=histos2.end() ; hit++)
	{	
		TString hname = hit->first;
		if(!isBGsum(hname) && !isSIGsum(hname)) continue;
		
		for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
		{
			TString name = cit->second;
			if(drawchannels[name]==0) continue;
			if(isData(name)) continue;
			
			TString hbarename = hname;
			if(isBGsum(hname) && !isSignal(name)) { hbarename.ReplaceAll("Backgrounds",""); hit->second->Add(histos2[name+hbarename]); }
			if(isSIGsum(hname) && isSignal(name)) { hbarename.ReplaceAll("Signals","");     hit->second->Add(histos2[name+hbarename]); }
		}
	}
	
	
	_INFO("");
	
	// histos decoration
	for(TMapTSP2TH1::iterator hit=histos.begin() ; hit!=histos.end() ; hit++)
	{
		TString hname = hit->first;

		
		if(!skim)
		{
			////////////////////////
			// normalize to unity //
			////////////////////////
			TString ytitle = hit->second->GetYaxis()->GetTitle();
			if(ytitle=="Normalized")            NormToEntries(hit->second);
			if(ytitle=="Normalized rate")       NormToEntries(hit->second);
			if(ytitle=="Normalized to 1st bin") NormTo1stBin(hit->second);
			if(ytitle.Contains("Efficiency"))
			{
				for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
				{
					TString cname = cit->second;
					if(!hname.Contains(cname)) continue;
					hit->second->Scale(1./nevents[cname]);
				}
			}
		}
		
		// set bin labels:
		// if(hname.EndsWith("mu_type"))
		if(hname.Contains("isTight"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"!isTight");
			hit->second->GetXaxis()->SetBinLabel(2,"isTight");
		}
		if(hname.Contains("isMedium"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"!isMedium");
			hit->second->GetXaxis()->SetBinLabel(2,"isMedium");
		}
		if(hname.Contains("isLoose"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"!isLoose");
			hit->second->GetXaxis()->SetBinLabel(2,"isLoose");
		}
		if(hname.Contains("_type_") || hname.EndsWith("mu_type"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"CB muon");
			hit->second->GetXaxis()->SetBinLabel(2,"TPmuonA");
			hit->second->GetXaxis()->SetBinLabel(3,"TPmuonB");
			hit->second->GetXaxis()->SetBinLabel(4,"CaloMuon");
		}
		if(hname.EndsWith("mu_MCP"))
		{
			hit->second->GetXaxis()->SetBinLabel(1,"Fail");
			hit->second->GetXaxis()->SetBinLabel(2,"Pass");
		}
		
		if(hname.Contains("triggers"))
		{
			buildTriggerbits();
			for(TMapuiTS::iterator it=triggerorder.begin() ; it!=triggerorder.end() ; ++it)
			{
				unsigned int bin = it->first;
				TString name      = it->second;
				hit->second->GetXaxis()->SetBinLabel(bin,name);
			}
		}
		
		// styles, colors and min/max
		for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
		{
			TString name  = cit->second;	
			Color_t color = colors[name];
			Int_t pattern = patterns[name];
			if(hname.Contains(name))
			{
				// styles, colors
				hit->second->SetLineColor(color);
				hit->second->SetFillColor(color);
				hit->second->SetMarkerColor(color);
				hit->second->SetFillStyle(pattern);
			}
		}
	}
	
	_INFO("");
	
	for(TMapTSP2TH2::iterator hit=histos2.begin() ; hit!=histos2.end() ; hit++)
	{
		TString hname = hit->first;

		if(splitstr=="")
		{
			////////////////////////
			// normalize to unity //
			////////////////////////
			TString ztitle = hit->second->GetZaxis()->GetTitle();
			if(ztitle.Contains("Normalized")) NormToEntries(hit->second);
		}
		
		// set bin labels:
		if(hname.Contains("_type_") || hname.EndsWith("mu_type"))
		{
			hit->second->GetYaxis()->SetBinLabel(1,"CB muon");
			hit->second->GetYaxis()->SetBinLabel(2,"TPmuonA");
			hit->second->GetYaxis()->SetBinLabel(3,"TPmuonB");
			hit->second->GetYaxis()->SetBinLabel(4,"CaloMuon");
		}

		// styles, colors
		for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
		{
			TString name  = cit->second;	
			Color_t color = colors[name];
			if(hname.Contains(name))
			{
				// styles, colors
				hit->second->SetLineColor(color);
				hit->second->SetMarkerColor(color);
				hit->second->SetFillColor(color);
			}
		}
	}
	
	
	_INFO("");
	
	
	// divide stuff...	
	divideCuts("m3mu_lin_zoom_ratio_MET",            "m3mu_lin_zoom_after_MET",            "m3mu_lin_zoom_before_W",             channels,histos);
	divideCuts("m3mu_lin_zoom_ratio_dphi3muMET",     "m3mu_lin_zoom_after_dphi3muMET",     "m3mu_lin_zoom_after_MET",            channels,histos);
	divideCuts("m3mu_lin_zoom_ratio_mT",             "m3mu_lin_zoom_after_mT",             "m3mu_lin_zoom_after_dphi3muMET",     channels,histos);;
	                                                                                                                            
	divideCuts("m3mu_sigregion_ratio_MET",           "m3mu_sigregion_after_MET",           "m3mu_sigregion_before_W",           channels,histos);
	divideCuts("m3mu_sigregion_ratio_dphi3muMET",    "m3mu_sigregion_after_dphi3muMET",    "m3mu_sigregion_after_MET",          channels,histos);
	divideCuts("m3mu_sigregion_ratio_mT",            "m3mu_sigregion_after_mT",            "m3mu_sigregion_after_dphi3muMET",   channels,histos);
	
	_INFO("");
	 

	TString filename = outDir+"../figs/fig."+chnl+"."+master+"."+smethod+"."+runType+"."+ssplit+".root";
	TFile* rootfile = new TFile(filename, "RECREATE");
	for(TMapTSP2TH1::iterator it=histos.begin()  ; it!=histos.end()  ; ++it) it->second->Write();
	for(TMapTSP2TH2::iterator it=histos2.begin() ; it!=histos2.end() ; ++it) it->second->Write();
	rootfile->Write();
	rootfile->Close();
	delete rootfile;
	
	
	_INFO("");
	
	/////////////////////////////////////////
	if(skim) finitOutTrees(fout,otrees); ////
	if(!skim) flatout_finit(olddir); ////////
	/////////////////////////////////////////
	
	int nProcessed = getCounter("nPassing_evt_all");
	if(nProcessed) cout << nProcessed << " EVENTS PROCESSED" << endl;
	else           cout << "NONE PROCESSED" << endl;
	_INFO("<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>");
	_INFO("<<<<<<<<<<<< ALL DONE CORRECTLY >>>>>>>>>>>");
	_INFO("<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>");

	
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	// plot everything only when all MC is processed ////
	// or when the W signal alone is processed //////////
	// bool goOn = (isAllMC || (!skim && isWSignalMC));
	bool goOn = (isAllMC); 
	if(!goOn) return; //////////////////////////////////
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	
	
	

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	ofstr->close();
	ofstr1->close();



	_INFO("");
	
	// pdf file with all the plots
	pdffilename = outDir+"../figs/fig."+chnl+"."+master+"."+smethod+"."+runType+"."+ssplit+".pdf";
	const Int_t nlines = 20;
	TString lines[nlines];
	resetLines(nlines,lines);
	int pcounter = -1;
	
	_INFO("");
	
	// legend...
	TLegend* leg     = new TLegend(0.55,0.6,0.83,0.83,NULL,"brNDC");
	TLegend* legleft = new TLegend(0.15,0.6,0.43,0.83,NULL,"brNDC");
	TLegend* legsig  = new TLegend(0.15,0.6,0.43,0.83,NULL,"brNDC");
	setLegendDefaults(leg);
	setLegendDefaults(legleft);
	setLegendDefaults(legsig);
	for(TMapuiTS::iterator cit=channels.begin() ; cit!=channels.end() ; cit++)
	{
		TString name = cit->second;
		if(drawchannels[name]==0) continue;
		if(isBGsum(name))  continue;
		if(isSIGsum(name)) continue;
		
		TString hname = name+"_chi2";
		leg->AddEntry(histos[hname],labels[name],legoptions[name]);
		legleft->AddEntry(histos[hname],labels[name],legoptions[name]);
		if(isSignal(name)) legsig->AddEntry(histos[hname],labels[name],legoptions[name]);
	}
	
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("cutflow_normalized", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename+"(");
	/////////////////////////////////////////////////////////
	

	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu",    "vtx_trumatch_x", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",    "vtx_trumatch_y", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",    "vtx_trumatch_z", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("bbTotau10_3mu", "vtx_trumatch_x", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("bbTotau10_3mu", "vtx_trumatch_y", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("bbTotau10_3mu", "vtx_trumatch_z", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("ccTotau10_3mu", "vtx_trumatch_x", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("ccTotau10_3mu", "vtx_trumatch_y", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("ccTotau10_3mu", "vtx_trumatch_z", pads[increment(pcounter)], channels, histos, "hist");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_pt_reldiff",  pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_E_reldiff",   pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_m_reldiff",   pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_eta_reldiff", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_phi_reldiff", pads[increment(pcounter)], channels, histos, "hist");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_reldiff_vs_pt",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_reldiff_vs_E",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_reldiff_vs_m",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_reldiff_vs_eta", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_calibration_reldiff_vs_phi", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_pt_reldiff",  pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_E_reldiff",   pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_m_reldiff",   pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_eta_reldiff", pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_phi_reldiff", pads[increment(pcounter)], channels, histos, "hist");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_reldiff_vs_pt",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_reldiff_vs_E",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_reldiff_vs_m",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_reldiff_vs_eta", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "jet_resolution_reldiff_vs_phi", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Wtaunu_3mu",  "met_calibration_et_reldiff",  pads[increment(pcounter)], channels, histos, "hist");
	drawPadSingle("Wtaunu_3mu",  "met_calibration_phi_reldiff",   pads[increment(pcounter)], channels, histos, "hist");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////


	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy);
	pcounter = -1;
	drawPadAll("tripletCategories_noVertexing",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_afterVertexing", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_after_muons",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_after_triplet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_after_hadclean", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_endOfSelection", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy);
	pcounter = -1;
	drawPadAll("tripletCategories_norm_noVertexing",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_norm_afterVertexing", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_norm_after_muons",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_norm_after_triplet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_norm_after_hadclean", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tripletCategories_norm_endOfSelection", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",          "cutflow_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	drawPadSingle("bb_mu4mu4",     "cutflow_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	drawPadSingle("Wtaunu_3mu",    "cutflow_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	drawPadSingle("bbTotau10_3mu", "cutflow_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",          "cutflow_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	drawPadSingle("bb_mu4mu4",     "cutflow_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	drawPadSingle("Wtaunu_3mu",    "cutflow_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	drawPadSingle("bbTotau10_3mu", "cutflow_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=1;
	// divy=1;
	// makeCnv(divx,divy,true);
	// pcounter = -1;
	// drawPadAll("counters_normalized", pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,true);
	// pcounter = -1;
	// drawPadSingle("Data",          "counters_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// drawPadSingle("bb_mu4mu4",     "counters_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// drawPadSingle("Wtaunu_3mu",    "counters_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// drawPadSingle("bbTotau10_3mu", "counters_absolute", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,true);
	// pcounter = -1;
	// drawPadSingle("Data",          "counters_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// drawPadSingle("bb_mu4mu4",     "counters_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// drawPadSingle("Wtaunu_3mu",    "counters_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// drawPadSingle("bbTotau10_3mu", "counters_weighted", pads[increment(pcounter)], channels, histos, "hist TEXT45");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////

	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("triggers",                       pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("unique_triggers",                pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triggers_after_hadclean",        pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("unique_triggers_after_hadclean", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triggers_after_vtxmet",          pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("unique_triggers_after_vtxmet",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("triplet_mQuad_muons", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_mQuad_norm_muons", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_mTrip_muons", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_mTrip_norm_muons", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("triplet_pTQuad1_muons", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_pTQuad2_muons", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_pTQuad3_muons", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_pTQuad4_muons", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("pvalue_all", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pvalue",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("chi2",       pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("",           pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type1_lxy",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type3_lxy",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////	
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("type1_a0",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type1_a0XY",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type1_cosTh",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type1_cosThXY", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("type3_a0",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type3_a0XY",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type3_cosTh",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("type3_cosThXY", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("vtx_npv_type1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pv_npv_type1",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("vtx_npv_type3", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pv_npv_type3",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pv_type",       pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("pv_ntrks",      pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("actualInteractionsPerXing", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////





	/*
	histos2[name+"_tpB_pvalue_vs_chi2ndf"]->Fill(fitchi2TP3/fitndfTP3,fitpvalueTP3,weight);
	*/


	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",       "mu_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadAll("mu_chi2ndf",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue",    pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pvalue_zoom", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",       "mu_pvalue_vs_chi2ndf_ismedium",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_pvalue_vs_chi2ndf_ismedium",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadAll("mu_chi2ndf_ismedium",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_ismedium",    pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pvalue_zoom_ismedium", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_chi2ndf_failMedium",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_failMedium",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_zoom_failMedium", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_chi2ndf_passMedium",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_passMedium",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_zoom_passMedium", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",       "mu_pvalue_vs_chi2ndf_failMedium",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_pvalue_vs_chi2ndf_failMedium",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mu_pvalue_vs_chi2ndf_passMedium",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mu_pvalue_vs_chi2ndf_passMedium",    pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",       "tpmu_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "tpmu_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadAll("tpmu_chi2ndf",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpmu_pvalue",    pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("tpmu_pvalue_zoom", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",       "tpA_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "tpA_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadAll("tpA_chi2ndf",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpA_pvalue",    pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("tpA_pvalue_zoom", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadSingle("Data",       "tpB_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "tpB_pvalue_vs_chi2ndf",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadAll("tpB_chi2ndf",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpB_pvalue",    pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("tpB_pvalue_zoom", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_matchchi2ndf", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mu_qoverp_me", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_qoverp_ie", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpA_qoverp",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpB_qoverp",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////



	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_sctangsig", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_sctngbsig", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pbalsig",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt",        pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_eta",       pads[increment(pcounter)], channels, histos, leg,-1,1);
	drawPadAll("mu_type",      pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_MCP", pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_nHighThresholdTRTHits",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("",                          pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpA_nHighThresholdTRTHits", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("",                          pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpB_nHighThresholdTRTHits", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("",                          pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",       "mu_chi2_vs_sctangsig",   pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu", "mu_chi2_vs_sctangsig",   pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",       "mu_pvalue_vs_sctangsig", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu", "mu_pvalue_vs_sctangsig", pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mu_isLoose1",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isMedium1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isTight1",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isLoose2",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isMedium2", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isTight2",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isLoose3",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isMedium3", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_isTight3",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mu_isMedium", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_pt1",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_eta1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_phi1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt2",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_eta2", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_phi2", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt3",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_eta3", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_phi3", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_pt1_3mu4T", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt2_3mu4T", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt3_3mu4T", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt1_2mu8_EFxe30_tclcw", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt2_2mu8_EFxe30_tclcw", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt3_2mu8_EFxe30_tclcw", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pt1_mu36_tight", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pt2_mu36_tight", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pt3_mu36_tight", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_pt1_unique_3mu4T", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt2_unique_3mu4T", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt3_unique_3mu4T", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt1_unique_2mu8_EFxe30_tclcw", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt2_unique_2mu8_EFxe30_tclcw", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pt3_unique_2mu8_EFxe30_tclcw", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pt1_unique_mu36_tight", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pt2_unique_mu36_tight", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("mu_pt3_unique_mu36_tight", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",       "mupt1_type_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mupt1_type_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mupt2_type_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mupt2_type_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mupt3_type_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mupt3_type_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true,"even");
	pcounter = -1;
	drawPadAll("muOnly_pt",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("muOnly_eta", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpA_pt",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpA_eta",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpB_pt",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpB_eta",    pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////	
	
	
	// _INFO("");
	// divx=1;
	// divy=1;
	// makeCnv(divx,divy);
	// pcounter = -1;
	// drawPadAll("triplet_cosTheta", pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy);
	// pcounter = -1;
	// drawPadAll("cosThOSQ_lightMesons",             pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("",                                 pads[increment(pcounter)], channels, histos, leg);
	// drawPadSingle("Data",       "cosThOSQ_vs_mOS", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu", "cosThOSQ_vs_mOS", pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=2;
	// divy=3;
	// makeCnv(divx,divy);
	// pcounter = -1;
	// drawPadAll("cosThOSQ_lightMesons",     pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("cosThOS_lightMesons",      pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("dPhiOS_lightMesons",       pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("cosThOStrk3_lightMesons",  pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("dmOS_lightMesons",         pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("",                         pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=2;
	// divy=3;
	// makeCnv(divx,divy);
	// pcounter = -1;
	// drawPadSingle("Data",       "cosThOS_vs_mOS",      pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu", "cosThOS_vs_mOS",      pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",       "dPhiOS_vs_mOS",       pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu", "dPhiOS_vs_mOS",       pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",       "cosThOStrk3_vs_mOS",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu", "cosThOStrk3_vs_mOS",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=2;
	// divy=1;
	// makeCnv(divx,divy);
	// pcounter = -1;
	// drawPadSingle("Data",       "dmOS_vs_mOS", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu", "dmOS_vs_mOS", pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("triplet_isolation",                       pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_isolation_zoom",                  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_isolation_before_isolation",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_isolation_zoom_before_isolation", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("triplet_dRmin",                   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_dRmax",                   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_dRmin_after_vtxclean",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_dRmax_after_vtxclean",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_dRmin_after_hadclean", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("triplet_dRmax_after_hadclean", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_muons",            pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "m3mu_vs_mOS_after_muons",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_mOS_after_muons",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_triplet",          pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "m3mu_vs_mOS_after_triplet",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_mOS_after_triplet",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_hadclean",           pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "m3mu_vs_mOS_after_hadclean",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_mOS_after_hadclean",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_vtxmet",           pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "m3mu_vs_mOS_after_vtxmet",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_mOS_after_vtxmet",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_muons",            pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "pT3mu_vs_mOS_after_muons",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3mu_vs_mOS_after_muons",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_triplet",          pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "pT3mu_vs_mOS_after_triplet",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3mu_vs_mOS_after_triplet",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_hadclean",           pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "pT3mu_vs_mOS_after_hadclean",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3mu_vs_mOS_after_hadclean",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_vtxmet",           pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "pT3mu_vs_mOS_after_vtxmet",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3mu_vs_mOS_after_vtxmet",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack("mSS_after_muons",     pads[increment(pcounter)], channels, histos,  leg);
	drawPadStack("mSS_after_triplet",   pads[increment(pcounter)], channels, histos,  leg);
	drawPadStack("mSS_after_hadclean",  pads[increment(pcounter)], channels, histos,  leg);
	drawPadStack("mSS_after_vtxmet",    pads[increment(pcounter)], channels, histos,  leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_muons",           pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mSS_vs_mOS_after_muons",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_mOS_after_muons",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_triplet",         pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mSS_vs_mOS_after_triplet",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_mOS_after_triplet",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_hadclean",          pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mSS_vs_mOS_after_hadclean",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_mOS_after_hadclean",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_vtxmet",          pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mSS_vs_mOS_after_vtxmet",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_mOS_after_vtxmet",   pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_muons",           pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mOS2_vs_mOS1_after_muons",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mOS2_vs_mOS1_after_muons",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_triplet",         pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mOS2_vs_mOS1_after_triplet",pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mOS2_vs_mOS1_after_triplet",pads[increment(pcounter)], channels, histos2, "col");
	drawPadStack(               "mOS_after_hadclean",          pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mOS2_vs_mOS1_after_hadclean", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mOS2_vs_mOS1_after_hadclean", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack(               "mOS_after_vtxmet",          pads[increment(pcounter)], channels, histos,  leg);
	drawPadSingle("Data",       "mOS2_vs_mOS1_after_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mOS2_vs_mOS1_after_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////





	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////


	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVA_score_all",            pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("",                         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_score_left_sideband",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_score_right_sideband", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",       "mOS1_vs_MVA_score", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mOS1_vs_MVA_score", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mOS2_vs_MVA_score", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mOS2_vs_MVA_score", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "mSS_vs_MVA_score",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_MVA_score",  pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",       "m3mu_vs_MVA_score",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_MVA_score",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "pT3mu_vs_MVA_score", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3mu_vs_MVA_score", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	

	
	
	// _INFO("");
	// divx=1;
	// divy=1;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadStack("pT3mu_before_triplet", pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false); //cnv = new TCanvas("c","c",1200,400);
	pcounter = -1;
	drawPadSingle("Data",        "flwMV1_vs_mjet",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Data",        "flwMV1_vs_pTjet", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Wtaunu_3mu",  "flwMV1_vs_mjet",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Wtaunu_3mu",  "flwMV1_vs_pTjet", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "flwMV1_vs_mjet",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "flwMV1_vs_pTjet", pads[increment(pcounter)], channels, histos2, "colz");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("jets_dPhi3bodyJ1_before_looseveto", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("jets_dPhi3bodyJ1_before_tightveto", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("jets_dPhi3bodyJ1_after_looseveto",  pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("jets_dPhi3bodyJ1_after_tightveto",  pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false); //cnv = new TCanvas("c","c",1200,400);
	pcounter = -1;
	drawPadSingle("Data",        "dPhi3muJet1_vs_pTjet1_before_looseveto", pads[increment(pcounter)], channels, histos2, "colz");	
	drawPadSingle("Wtaunu_3mu",  "dPhi3muJet1_vs_pTjet1_before_looseveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhi3muJet1_vs_pTjet1_before_looseveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Data",        "dPhi3muJet1_vs_pTjet1_after_looseveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Wtaunu_3mu",  "dPhi3muJet1_vs_pTjet1_after_looseveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhi3muJet1_vs_pTjet1_after_looseveto",  pads[increment(pcounter)], channels, histos2, "colz");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false); //cnv = new TCanvas("c","c",1200,400);
	pcounter = -1;
	drawPadSingle("Data",        "dPhi3muJet1_vs_pTjet1_before_tightveto", pads[increment(pcounter)], channels, histos2, "colz");	
	drawPadSingle("Wtaunu_3mu",  "dPhi3muJet1_vs_pTjet1_before_tightveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhi3muJet1_vs_pTjet1_before_tightveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Data",        "dPhi3muJet1_vs_pTjet1_after_tightveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Wtaunu_3mu",  "dPhi3muJet1_vs_pTjet1_after_tightveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhi3muJet1_vs_pTjet1_after_tightveto",  pads[increment(pcounter)], channels, histos2, "colz");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
		
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("jets_dPhiJ1J2_before_looseveto", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("jets_dPhiJ1J2_before_tightveto", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("jets_dPhiJ1J2_after_looseveto",  pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("jets_dPhiJ1J2_after_tightveto",  pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false); //cnv = new TCanvas("c","c",1200,400);
	pcounter = -1;
	drawPadSingle("Data",        "dPhiJet1Jet2_vs_sumpTjet12_before_looseveto", pads[increment(pcounter)], channels, histos2, "colz");	
	drawPadSingle("Wtaunu_3mu",  "dPhiJet1Jet2_vs_sumpTjet12_before_looseveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhiJet1Jet2_vs_sumpTjet12_before_looseveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Data",        "dPhiJet1Jet2_vs_sumpTjet12_after_looseveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Wtaunu_3mu",  "dPhiJet1Jet2_vs_sumpTjet12_after_looseveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhiJet1Jet2_vs_sumpTjet12_after_looseveto",  pads[increment(pcounter)], channels, histos2, "colz");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false); //cnv = new TCanvas("c","c",1200,400);
	pcounter = -1;
	drawPadSingle("Data",        "dPhiJet1Jet2_vs_sumpTjet12_before_tightveto", pads[increment(pcounter)], channels, histos2, "colz");	
	drawPadSingle("Wtaunu_3mu",  "dPhiJet1Jet2_vs_sumpTjet12_before_tightveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhiJet1Jet2_vs_sumpTjet12_before_tightveto", pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Data",        "dPhiJet1Jet2_vs_sumpTjet12_after_tightveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Wtaunu_3mu",  "dPhiJet1Jet2_vs_sumpTjet12_after_tightveto",  pads[increment(pcounter)], channels, histos2, "colz");
	drawPadSingle("Backgrounds", "dPhiJet1Jet2_vs_sumpTjet12_after_tightveto",  pads[increment(pcounter)], channels, histos2, "colz");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("met_sumet",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("met_et",       pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("met_mt_et3mu", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("met_dphi3mu",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	///////////////////////////////////////////////////////
	
	
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack("m3mu_lin_zoom_after_muons",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pT3mu_after_muons",               pads[increment(pcounter)], channels, histos, leg);
	drawPadStack("m3mu_lin_zoom_after_triplet",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pT3mu_after_triplet",             pads[increment(pcounter)], channels, histos, leg);
	drawPadStack("m3mu_lin_zoom_after_hadclean",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pT3mu_after_hadclean",            pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack("m3mu_lin_zoom_after_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadStack("pT3mu_after_vtxmet",          pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack("m3mu_sigregion_after_muons",   pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("m3mu_sigregion_after_triplet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("m3mu_sigregion_after_hadclean",pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("m3mu_sigregion_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVA_vars_mSS_before_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mOS1_before_vtxmet", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mOS2_before_vtxmet", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mSS_after_vtxmet",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mOS1_after_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mOS2_after_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=2;
	divy=1;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVA_vars_pT3body_before_vtxmet",        pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_pT3body_after_vtxmet",         pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVA_vars_mupt12Fraction_before_vtxmet", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mupt23Fraction_before_vtxmet", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mupt13Fraction_before_vtxmet", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mupt12Fraction_after_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mupt23Fraction_after_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_mupt13Fraction_after_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("MVA_vars_isolation_before_vtxmet", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_pvalue_before_vtxmet",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_isolation_after_vtxmet",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_pvalue_after_vtxmet",     pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("MVA_vars_Lxy_before_vtxmet",    pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_a0xy_before_vtxmet",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_cosTxy_before_vtxmet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_Lxy_after_vtxmet",     pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_a0xy_after_vtxmet",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("MVA_vars_cosTxy_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVA_vars_MET_before_vtxmet",          pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_dPhi3bodyMET_before_vtxmet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_mT3bodyMET_before_vtxmet",   pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_MET_after_vtxmet",           pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_dPhi3bodyMET_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_mT3bodyMET_after_vtxmet",    pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVA_vars_pTJ1_before_vtxmet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_pTJ2_before_vtxmet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_pTJ1_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_pTJ2_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("MVA_vars_dPhi3bodyJ1_before_vtxmet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_dPhiJ1J2_before_vtxmet",    pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_dPhi3bodyJ1_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	drawPadAll("MVA_vars_dPhiJ1J2_after_vtxmet",     pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",        "MVA_vars_dPhi3muJet1_vs_pTjet1_before_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "MVA_vars_dPhi3muJet1_vs_pTjet1_before_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "MVA_vars_dPhi3muJet1_vs_pTjet1_after_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "MVA_vars_dPhi3muJet1_vs_pTjet1_after_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",        "MVA_vars_dPhiJet1Jet2_vs_sumpTjet12_before_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "MVA_vars_dPhiJet1Jet2_vs_sumpTjet12_before_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "MVA_vars_dPhiJet1Jet2_vs_sumpTjet12_after_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "MVA_vars_dPhiJet1Jet2_vs_sumpTjet12_after_vtxmet", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack("MVA_vars_m3mu_before_vtxmet",     pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("MVA_vars_m3mu_lin_before_vtxmet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("MVA_vars_m3mu_after_vtxmet",      pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("MVA_vars_m3mu_lin_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadStack("MVA_vars_m3mu_lin_zoom_before_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("MVA_vars_m3mu_sigregion_before_vtxmet", pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("MVA_vars_m3mu_lin_zoom_after_vtxmet",   pads[increment(pcounter)], channels, histos, legleft);
	drawPadStack("MVA_vars_m3mu_sigregion_after_vtxmet",  pads[increment(pcounter)], channels, histos, legleft);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////	
	/////////////////////////////////////////////////////////////////////////////	
	/////////////////////////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",        "mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",        "mupt1Fraction_vs_mupt1_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt1Fraction_vs_mupt1_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt1Fraction_vs_mupt1_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "mupt2Fraction_vs_mupt2_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt2Fraction_vs_mupt2_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt2Fraction_vs_mupt2_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "mupt3Fraction_vs_mupt3_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt3Fraction_vs_mupt3_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt3Fraction_vs_mupt3_afterCut", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",         "mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("bbTotau10_3mu","mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds",  "mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",         "mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("bbTotau10_3mu","mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds",  "mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",         "mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("bbTotau10_3mu","mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds",  "mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",         "mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("ccTotau10_3mu","mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds",  "mupt1Fraction_vs_mupt1", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",         "mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("ccTotau10_3mu","mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds",  "mupt2Fraction_vs_mupt2", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",         "mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("ccTotau10_3mu","mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds",  "mupt3Fraction_vs_mupt3", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",        "mupt12Fraction_vs_mupt12", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt12Fraction_vs_mupt12", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt12Fraction_vs_mupt12", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "mupt23Fraction_vs_mupt23", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt23Fraction_vs_mupt23", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt23Fraction_vs_mupt23", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "mupt13Fraction_vs_mupt13", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt13Fraction_vs_mupt13", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt13Fraction_vs_mupt13", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",        "mupt23Fraction_vs_mupt12Fraction", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt23Fraction_vs_mupt12Fraction", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt23Fraction_vs_mupt12Fraction", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",        "mupt13Fraction_vs_mupt12Fraction", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu",  "mupt13Fraction_vs_mupt12Fraction", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Backgrounds", "mupt13Fraction_vs_mupt12Fraction", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=4;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_ptcone10.1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone10.2", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone10.3", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone20.1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone20.2", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone20.3", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone30.1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone30.2", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone30.3", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone40.1", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone40.2", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone40.3", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=4;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_ptcone10.1_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone10.2_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone10.3_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone20.1_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone20.2_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone20.3_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone30.1_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone30.2_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone30.3_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone40.1_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone40.2_zoom", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_ptcone40.3_zoom", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	_INFO("");
	divx=3;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mupt1Fraction",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mupt2Fraction",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mupt3Fraction",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mupt1Fraction_afterCut",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mupt2Fraction_afterCut",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mupt3Fraction_afterCut",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=3;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mupt12Fraction",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mupt23Fraction",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mupt13Fraction",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",       "m3mu_vs_dRmin",                   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_dRmin",                   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "m3mu_vs_dRmin_after_vtxclean",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_dRmin_after_vtxclean",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "m3mu_vs_dRmin_after_hadclean", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_dRmin_after_hadclean", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",       "m3mu_vs_dRmax",                   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_dRmax",                   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "m3mu_vs_dRmax_after_vtxclean",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_dRmax_after_vtxclean",    pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "m3mu_vs_dRmax_after_hadclean", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_dRmax_after_hadclean", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true,"");
	pcounter = -1;
	drawPadAll("bjet_dR3muLeadingJet",            pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muSubleadingJet",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muLeadingJet_noOvrlp",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muSubleadingJet_noOvrlp", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muLeadingJet_stdMV1",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muSubleadingJet_stdMV1",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true,"");
	pcounter = -1;
	drawPadAll("bjet_dR3muBjet_btagwgt1",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muBjet_btagwgt2",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muBjet_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muBjet_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muBjet_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dR3muBjet_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true,"");
	pcounter = -1;
	drawPadAll("bjet_dPhi3muLeadingJet",            pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muSubleadingJet",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muLeadingJet_noOvrlp",    pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muSubleadingJet_noOvrlp", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muLeadingJet_stdMV1",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muSubleadingJet_stdMV1",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true,"");
	pcounter = -1;
	drawPadAll("bjet_dPhi3muBjet_btagwgt1",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muBjet_btagwgt2",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muBjet_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muBjet_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muBjet_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi3muBjet_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true,"");
	pcounter = -1;
	drawPadAll("bjet_dPhi_jet1_jet2",                 pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi_btagwgt1_btagwgt2",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi_jet1_jet2_noOvrlp",         pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi_btagwgt1_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi_jet1_jet2_stdMV1",          pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("bjet_dPhi_btagwgt1_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "bjet_dR3muBjet_btagwgt1_vs_btagwgt1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dR3muBjet_btagwgt1_vs_btagwgt1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dR3muBjet_btagwgt1_vs_btagwgt1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dR3muBjet_btagwgt1_vs_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dR3muBjet_btagwgt1_vs_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dR3muBjet_btagwgt1_vs_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dR3muBjet_btagwgt1_vs_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dR3muBjet_btagwgt1_vs_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dR3muBjet_btagwgt1_vs_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "bjet_dR3muBjet_btagwgt2_vs_btagwgt2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dR3muBjet_btagwgt2_vs_btagwgt2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dR3muBjet_btagwgt2_vs_btagwgt2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dR3muBjet_btagwgt2_vs_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dR3muBjet_btagwgt2_vs_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dR3muBjet_btagwgt2_vs_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dR3muBjet_btagwgt2_vs_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dR3muBjet_btagwgt2_vs_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dR3muBjet_btagwgt2_vs_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhi3muBjet_btagwgt1_vs_btagwgt1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhi3muBjet_btagwgt2_vs_btagwgt2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "bjet_dPhib1b2_vs_b1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhib1b2_vs_b1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhib1b2_vs_b1",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhib1b2_vs_b1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhib1b2_vs_b1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhib1b2_vs_b1_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhib1b2_vs_b1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhib1b2_vs_b1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhib1b2_vs_b1_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "bjet_dPhib1b2_vs_b2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhib1b2_vs_b2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhib1b2_vs_b2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhib1b2_vs_b2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhib1b2_vs_b2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhib1b2_vs_b2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_dPhib1b2_vs_b2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_dPhib1b2_vs_b2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_dPhib1b2_vs_b2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "jet_wj1_vs_wj2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "jet_wj1_vs_wj2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "jet_wj1_vs_wj2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "jet_wj1_vs_wj2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "jet_wj1_vs_wj2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "jet_wj1_vs_wj2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "jet_wj1_vs_wj2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "jet_wj1_vs_wj2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "jet_wj1_vs_wj2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=3;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Data",        "bjet_wb1_vs_wb2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_wb1_vs_wb2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_wb1_vs_wb2",         pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_wb1_vs_wb2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_wb1_vs_wb2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_wb1_vs_wb2_noOvrlp", pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Data",        "bjet_wb1_vs_wb2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Wtaunu_3mu",  "bjet_wb1_vs_wb2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// drawPadSingle("Backgrounds", "bjet_wb1_vs_wb2_stdMV1",  pads[increment(pcounter)], channels, histos2, "col");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	// _INFO("");
	// divx=3;
	// divy=1;
	// makeCnv(divx,divy,true,"");
	// pcounter = -1;
	// drawPadAll("jet_maxFlvWgt_MV1",                pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("jet_leadingJetFlvWgt_MV1",         pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("jet_subleadingJetFlvWgt_MV1",      pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,true,"");
	pcounter = -1;
	drawPadAll("jet_n",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("jet_E",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("jet_pt",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("jet_m",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("jet_eta", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("jet_phi", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	
	
	
	
	// _INFO("");
	// divx=2;
	// divy=3;
	// makeCnv(divx,divy,true,"");
	// pcounter = -1;
	// drawPadAll("bjet_n",   pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_E",   pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_pt",  pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_m",   pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_eta", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_phi", pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=3;
	// divy=2;
	// makeCnv(divx,divy,true,"");
	// pcounter = -1;
	// drawPadAll("bjet_maxFlvWgt_E",   pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_maxFlvWgt_pt",  pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_maxFlvWgt_m",   pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("",                   pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_maxFlvWgt_eta", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_maxFlvWgt_phi", pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=2;
	// divy=2;
	// makeCnv(divx,divy,true,"");
	// pcounter = -1;
	// drawPadAll("jet_n10",  pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_n10", pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("jet_n20",  pads[increment(pcounter)], channels, histos, leg);
	// drawPadAll("bjet_n20", pads[increment(pcounter)], channels, histos, leg);
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	// 
	// _INFO("");
	// divx=2;
	// divy=3;
	// makeCnv(divx,divy,false);
	// pcounter = -1;
	// drawPadSingle("Wtaunu_3mu",    "nBtag_vs_nJet", pads[increment(pcounter)], channels, histos2, "colz");
	// drawPadSingle("ccTotau10_3mu", "nBtag_vs_nJet", pads[increment(pcounter)], channels, histos2, "colz");
	// drawPadSingle("bbTotau10_3mu", "nBtag_vs_nJet", pads[increment(pcounter)], channels, histos2, "colz");
	// drawPadSingle("bbTomu15",      "nBtag_vs_nJet", pads[increment(pcounter)], channels, histos2, "colz");
	// drawPadSingle("ccTomu15",      "nBtag_vs_nJet", pads[increment(pcounter)], channels, histos2, "colz");
	// drawPadSingle("Data",          "nBtag_vs_nJet", pads[increment(pcounter)], channels, histos2, "colz");
	// closeCnv(pdffilename);
	// /////////////////////////////////////////////////////////
	
	
	
	
	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mOS_rhoomegaphi",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mSS_rhoomegaphi",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("m3mu_rhoomegaphi",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pT3mu_rhoomegaphi", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_chi2ndf_failMedium_rhoomegaphi",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_failMedium_rhoomegaphi",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_zoom_failMedium_rhoomegaphi", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_chi2ndf_passMedium_rhoomegaphi",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_passMedium_rhoomegaphi",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_zoom_passMedium_rhoomegaphi", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpmu_chi2ndf_rhoomegaphi",              pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpmu_pvalue_rhoomegaphi",               pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpmu_pvalue_zoom_rhoomegaphi",          pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",       "mSS_vs_mOS_rhoomegaphi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_mOS_rhoomegaphi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "m3mu_vs_mOS_rhoomegaphi",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_mOS_rhoomegaphi",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "pT3mu_vs_mOS_rhoomegaphi", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3mu_vs_mOS_rhoomegaphi", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("pvalue_rhoomegaphi",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	
	
	
	
	
	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=2;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadAll("mOS_Jpsi",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mSS_Jpsi",   pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("m3mu_Jpsi",  pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("pT3mu_Jpsi",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=3;
	divy=3;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("mu_chi2ndf_failMedium_Jpsi",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_failMedium_Jpsi",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_zoom_failMedium_Jpsi", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_chi2ndf_passMedium_Jpsi",     pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_passMedium_Jpsi",      pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("mu_pvalue_zoom_passMedium_Jpsi", pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpmu_chi2ndf_Jpsi",              pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpmu_pvalue_Jpsi",               pads[increment(pcounter)], channels, histos, leg);
	drawPadAll("tpmu_pvalue_zoom_Jpsi",          pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	_INFO("");
	divx=2;
	divy=3;
	makeCnv(divx,divy,false);
	pcounter = -1;
	drawPadSingle("Data",       "mSS_vs_mOS_Jpsi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "mSS_vs_mOS_Jpsi",   pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "m3mu_vs_mOS_Jpsi",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "m3mu_vs_mOS_Jpsi",  pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Data",       "pT3mu_vs_mOS_Jpsi", pads[increment(pcounter)], channels, histos2, "col");
	drawPadSingle("Wtaunu_3mu", "pT3mu_vs_mOS_Jpsi", pads[increment(pcounter)], channels, histos2, "col");
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////

	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,true);
	pcounter = -1;
	drawPadAll("pvalue_Jpsi",   pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename);
	/////////////////////////////////////////////////////////
	
	

	
	
	// blank page
	_INFO("");
	divx=1;
	divy=1;
	makeCnv(divx,divy,false,"");
	pcounter = -1;
	drawPadAll("", pads[increment(pcounter)], channels, histos, leg);
	closeCnv(pdffilename+")");
	/////////////////////////////////////////////////////////
	
	
	
	_INFO("");
	
	
	// Normalize on the J/psi peak (m_Jpsi = 3.1 GeV)
	// where the J/psi + HF components are dominant.
	// This is valid only for the 3mu analysis.
	float xminHFnorm = 2.95*GeV2MeV;
	float xmaxHFnorm = 3.20*GeV2MeV;
	getHFnormalization(channels,histos,xminHFnorm,xmaxHFnorm,"mOS_after_vtxclean"); // this is for the histo filled after the 3mu p-value and the Lxy+pvalueOS+mindchi2 cuts), with the J/psi->mu4mu4 component in and just before cutting on mOS (removing J/psi's)
	// getHFnormalization(channels,histos,xminHFnorm,xmaxHFnorm,"mOS");
}
