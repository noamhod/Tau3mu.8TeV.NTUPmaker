#include "rawROOT.h"
//#include "enums.h"

static const unsigned int nMaxTracks = 6; // including the triplet, i.e. 3+3 here

/*
enum vtx_systshift
{
	VTX_SHIFTDWN=-1,
	VTX_NOSHIFT,
	VTX_SHIFTUP
};
enum vtx_uncertainties
{
	VTX_NOJES,
	VTX_NOMINAL,
	VTX_JESUP, VTX_JESDWN,
	VTX_JERUP, VTX_JERDWN,
	VTX_JETTRKNOM, VTX_JETTRKUP, VTX_JETTRKDWN
	VTX_SOFTTRKNOM, VTX_SOFTTRKUP, VTX_SOFTTRKDWN
	VTX_SOFTTRKRESPARA, VTX_SOFTTRKRESPREP, VTX_SOFTTRKRESCORR
};
*/

enum met_types
{
	METSTACO,METMUONS,METTRACK
};


class vertex
{
public:
	vertex() {};
	~vertex() {};
	
public:
	//// setter
	void           set(unsigned int vtx);
	
	//// classificators
	int            vtxIndex() { return m_index; };
	int            vtxCode()  { return m_code;  };
	TString        vtxType()  { return m_type;  };
	
	//// vertex
	TLorentzVector vtxP()              { return m_p4;           };
	TLorentzVector vtxPOS1()           { return m_pOS1;         };
	TLorentzVector vtxPOS2()           { return m_pOS2;         };
	TLorentzVector vtxPSS()            { return m_pSS;          };
	TLorentzVector vtxPquad(int i)     { return m_pQuad[i];     };
	double         vtxPt()             { return m_p4.Pt();      };
	double         vtxM()              { return m_p4.M();       };
	double         vtxMOS1()           { return m_pOS1.M();     };
	double         vtxMOS2()           { return m_pOS2.M();     };
	double         vtxMSS()            { return m_pSS.M();      };
	double         vtxMquad(int i)     { return m_pQuad[i].M(); };
	int            vtxNdf()            { return m_ndf;          };
	double         vtxChi2()           { return m_chi2;         };
	double         vtxChi2Ndf()        { return m_chi2ndf;      };
	double         vtxPvalue()         { return m_pvalue;       };
	double         vtxLxy()            { return m_lxy;          };
	double         vtxLxyErr()         { return m_lxyErr;       };
	double         vtxTau()            { return m_tau;          };
	double         vtxTauErr()         { return m_tauErr;       };
	double         vtxA0()             { return m_a0;           };
	double         vtxA0Err()          { return m_a0Err;        };
	double         vtxA0xy()           { return m_a0xy;         };
	double         vtxA0xyErr()        { return m_a0xyErr;      };
	double         vtxCosT()           { return m_cosT;         };
	double         vtxCosTxy()         { return m_cosTxy;       };
	double         vtxDPt12()          { return m_dpt12;        };
	double         vtxDPt23()          { return m_dpt23;        };
	double         vtxDPt13()          { return m_dpt13;        };
	double         vtxPtFrac12()       { return m_ptFrac12;     };
	double         vtxPtFrac23()       { return m_ptFrac23;     };
	double         vtxPtFrac13()       { return m_ptFrac13;     };
	double         vtxQ()              { return m_charge;       };
	double         vtxIsolation(int i) { return m_isolation[i]; }; // margins=i*0.01
	double         vtxDRmax()          { return m_drmax;        };
	double         vtxDRmin()          { return m_drmin;        };
	int            vtxPVntrk()         { return m_pvNtrk;       };
	int            vtxNelectrons()     { return m_nElectrons;   };
	
	//// tracks
 	string         trkSrc(int i)            { return m_src[iSorted(i)];               };
	int            trkSrcIndex(int i)       { return m_isrc[iSorted(i)];              };
 	int            trkTrkIndex(int i)       { return m_itrk[iSorted(i)];              };
	int            trkType(int i)           { return m_trktype[iSorted(i)];           };
	int            trkOrder(int i)          { return m_order[iSorted(i)];             };
	bool           trkIsMuon(int i)         { return m_ismuon[iSorted(i)];            };
	bool           trkIsTP(int i)           { return m_istp[iSorted(i)];              };
	bool           trkIsTP1(int i)          { return m_istpa[iSorted(i)];             };
	bool           trkIsTP2(int i)          { return m_istpb[iSorted(i)];             };
	bool           trkIsCalo(int i)         { return m_iscalo[iSorted(i)];            };
	bool           trkIsMediumMu(int i)     { return m_ismedium[iSorted(i)];          };
	bool           trkIsTightMu(int i)      { return m_istight[iSorted(i)];           };
	bool           trkIsLooseMu(int i)      { return m_isloose[iSorted(i)];           };
	bool           trkIsCBMu(int i)         { return m_iscb[iSorted(i)];              };
	double         trkMatchChi2NdfMu(int i) { return m_trkMuMatchChi2Ndf[iSorted(i)]; };
	TLorentzVector trkP(int i)              { return m_trkP[iSorted(i)];              };
	double         trkChi2(int i)           { return m_trkChi2[iSorted(i)];           };
	int            trkNdf(int i)            { return m_trkNdf[iSorted(i)];            };
	double         trkChi2Ndf(int i)        { return m_trkChi2Ndf[iSorted(i)];        };
	double         trkPval(int i)           { return m_trkPval[iSorted(i)];           };
	double         trkPtfrac(int i)         { return m_trkptfrac[iSorted(i)];         };
	double         trkPbalMu(int i)         { return m_trkpbal[iSorted(i)];           };
	double         trkSctAngMu(int i)       { return m_trksctang[iSorted(i)];         };
	double         trkSctNgbMu(int i)       { return m_trksctngb[iSorted(i)];         };
	double         trkQoverPtrk(int i)      { return m_trkQoverP[iSorted(i)];         };
	double         trkQoverPsrc(int i)      { return m_srcQoverP[iSorted(i)];         };
	double         trkPixeldEdx(int i)      { return m_trkPixeldEdx[iSorted(i)];      };
	int            trkPIXhits(int i)        { return m_trkPIXhits[iSorted(i)];        };
	int            trkDeadPIX(int i)        { return m_trkDeadPIX[iSorted(i)];        };
	int            trkPIXholes(int i)       { return m_trkPIXholes[iSorted(i)];       };	
	int            trkSCThits(int i)        { return m_trkSCThits[iSorted(i)];        };
	int            trkDeadSCT(int i)        { return m_trkDeadSCT[iSorted(i)];        };
	int            trkSCTholes(int i)       { return m_trkSCTholes[iSorted(i)];       };
	int            trkTRThits(int i)        { return m_trkTRThits[iSorted(i)];        };
	int            trkTRToutliers(int i)    { return m_trkTRToutliers[iSorted(i)];    };
	int            trkHtTRThits(int i)      { return m_trkHtTRThits[iSorted(i)];      };
	int            trkUsedHitsdEdx(int i)   { return m_trkUsedHitsdEdx[iSorted(i)];   };
	
	int            trkMDThits(int i)         { return m_trkMDThits[iSorted(i)];         };
	int            trkTGCPhiHits(int i)      { return m_trkTGCPhiHits[iSorted(i)];      };
	int            trkTGCEtaHits(int i)      { return m_trkTGCEtaHits[iSorted(i)];      };
	int            trkCSCPhiHits(int i)      { return m_trkCSCPhiHits[iSorted(i)];      };
	int            trkCSCEtaHits(int i)      { return m_trkCSCEtaHits[iSorted(i)];      };
	int            trkRPCPhiHits(int i)      { return m_trkRPCPhiHits[iSorted(i)];      };
	int            trkRPCEtaHits(int i)      { return m_trkRPCEtaHits[iSorted(i)];      };
	int            trkCSCEtaHoles(int i)     { return m_trkCSCEtaHoles[iSorted(i)];     };
	int            trkCSCPhiHoles(int i)     { return m_trkCSCPhiHoles[iSorted(i)];     };
	int            trkRPCEtaHoles(int i)     { return m_trkRPCEtaHoles[iSorted(i)];     };
	int            trkRPCPhiHoles(int i)     { return m_trkRPCPhiHoles[iSorted(i)];     };
	int            trkMDTholes(int i)        { return m_trkMDTholes[iSorted(i)];        };
	int            trkTGCEtaHoles(int i)     { return m_trkTGCEtaHoles[iSorted(i)];     };
	int            trkTGCPhiHoles(int i)     { return m_trkTGCPhiHoles[iSorted(i)];     };
	int            trkOutliersOnTrack(int i) { return m_trkOutliersOnTrack[iSorted(i)]; };
	int            trkStdDevOfChi2OS(int i)  { return m_trkStdDevOfChi2OS[iSorted(i)];  };

	int            trkPrecisionHits(int i)        { return m_trkPrecisionHits[iSorted(i)];        };
	int            trkPhiLayers(int i)            { return m_trkPhiLayers[iSorted(i)];            };
	int            trkEtaPhiLayers(int i)         { return m_trkEtaPhiLayers[iSorted(i)];         };
	int            trkPrecisionHoles(int i)       { return m_trkPrecisionHoles[iSorted(i)];       };
	int            trkEtaTriggerHoleLayers(int i) { return m_trkEtaTriggerHoleLayers[iSorted(i)]; };
	int            trkPhiHoleLayers(int i)        { return m_trkPhiHoleLayers[iSorted(i)];        };
	int            trkPrecisionOutliers(int i)    { return m_trkPrecisionOutliers[iSorted(i)];    };
	
	
	//// jets
	int            jetN(int mode)              { return m_njets[mode]; };
	TLorentzVector jetPEall(int mode, int i)   { return m_jetPEall[mode][i];   };
	TLorentzVector jetPMall(int mode, int i)   { return m_jetPMall[mode][i];   };
	double         jetMV1all(int mode, int i)  { return m_jetMV1all[mode][i];  };
	double         jetVtxFall(int mode, int i) { return m_jetVtxFall[mode][i]; };
	int            jetNtrkall(int mode, int i) { return m_jetNtrkall[mode][i]; };
	double         jetSumPtAll(int mode)       { return m_jetSumpt12[mode];    };
	double         jetDphi12All(int mode)      { return m_jetDphi12[mode];     };
	double         jetDR12All(int mode)        { return m_jetDR12[mode];       };
	double         jetDphi3bodyAll(int mode)   { return m_jetDphi3body[mode];  };
	double         jetDR3bodyAll(int mode)     { return m_jetDR3body[mode];    };
	double         jetShiftJES(int i) { return m_jet_shiftJES[i]; };
	double         jetShiftJER(int i) { return m_jet_shiftJER[i]; };
	
	//// MET
	double         met(int type, int mode)          { return m_met[type][mode];          };
	double         metPhi(int type, int mode)       { return m_metPhi[type][mode];       };
	double         metDphi3body(int type, int mode) { return m_metDphi3body[type][mode]; };
	double         metMt(int type, int mode)        { return m_metMt[type][mode];        };
	


private:
	int iSorted(int i);
	
public:
	
private:
	int m_index, m_code;
	TString m_type;
	TLorentzVector m_p4, m_pOS1, m_pOS2, m_pSS, m_pQuad[nMaxTracks-3];
	TLorentzVector m_trkP[nMaxTracks];
	int m_njets[6];
	int m_pvNtrk;
	int m_nElectrons;
	TLorentzVector m_jetPEall[6][4], m_jetPMall[6][4];
	double m_jet_shiftJES[4], m_jet_shiftJER[4];
	double m_jetMV1all[6][4], m_jetVtxFall[6][4];
	int m_jetNtrkall[6][4];
	double m_jetSumpt12[6], m_jetDphi12[6], m_jetDphi3body[6], m_jetDR12[6], m_jetDR3body[6];
	bool m_ismuon[nMaxTracks], m_istp[nMaxTracks], m_istpa[nMaxTracks], m_istpb[nMaxTracks], m_iscalo[nMaxTracks];
	bool m_ismedium[nMaxTracks], m_istight[nMaxTracks], m_isloose[nMaxTracks], m_iscb[nMaxTracks];
	int m_isrc[nMaxTracks], m_itrk[nMaxTracks], m_trktype[nMaxTracks], m_order[nMaxTracks];
	string m_src[nMaxTracks];
	double m_trkChi2[nMaxTracks], m_trkChi2Ndf[nMaxTracks], m_trkPval[nMaxTracks], m_trkMuMatchChi2Ndf[nMaxTracks];
	double m_trkptfrac[nMaxTracks], m_ptFrac12, m_ptFrac13, m_ptFrac23, m_dpt12, m_dpt13, m_dpt23;
	double m_trkpbal[nMaxTracks], m_trksctang[nMaxTracks], m_trksctngb[nMaxTracks];
	double m_trkQoverP[nMaxTracks], m_srcQoverP[nMaxTracks], m_trkPixeldEdx[nMaxTracks];
	int m_trkNdf[nMaxTracks];
	int m_trkPIXhits[nMaxTracks], m_trkDeadPIX[nMaxTracks], m_trkPIXholes[nMaxTracks];
	int m_trkSCThits[nMaxTracks], m_trkDeadSCT[nMaxTracks], m_trkSCTholes[nMaxTracks];
	int m_trkTRThits[nMaxTracks], m_trkTRToutliers[nMaxTracks], m_trkHtTRThits[nMaxTracks];
	int m_trkUsedHitsdEdx[nMaxTracks];
	int m_trkMDThits[nMaxTracks], m_trkTGCPhiHits[nMaxTracks], m_trkTGCEtaHits[nMaxTracks], m_trkCSCPhiHits[nMaxTracks], m_trkCSCEtaHits[nMaxTracks], m_trkRPCPhiHits[nMaxTracks], m_trkRPCEtaHits[nMaxTracks];
	int m_trkCSCEtaHoles[nMaxTracks], m_trkCSCPhiHoles[nMaxTracks], m_trkRPCEtaHoles[nMaxTracks], m_trkRPCPhiHoles[nMaxTracks], m_trkMDTholes[nMaxTracks], m_trkTGCEtaHoles[nMaxTracks], m_trkTGCPhiHoles[nMaxTracks];
	int m_trkOutliersOnTrack[nMaxTracks], m_trkStdDevOfChi2OS[nMaxTracks];
	int m_trkPrecisionHits[nMaxTracks], m_trkPhiLayers[nMaxTracks], m_trkEtaPhiLayers[nMaxTracks], m_trkPrecisionHoles[nMaxTracks], m_trkEtaTriggerHoleLayers[nMaxTracks], m_trkPhiHoleLayers[nMaxTracks], m_trkPrecisionOutliers[nMaxTracks];
	
	double m_met[3][15], m_metPhi[3][15], m_metDphi3body[3][15], m_metMt[3][15];
	
	double m_chi2, m_ndf, m_chi2ndf, m_pvalue, m_lxy, m_lxyErr, m_tau, m_tauErr;
	double m_a0, m_a0Err, m_a0xy, m_a0xyErr, m_cosT, m_cosTxy, m_charge, m_isolation[21], m_drmax, m_drmin;
};

int vertex::iSorted(int i)
{
	unsigned int j = 999;
	if(i<0) _FATAL("i<0: "+_s(i));
	j = (unsigned int)i;
	if(j>2 && j<nMaxTracks) return j;
	if((m_order[0]-1)==i)   return 0;
	if((m_order[1]-1)==i)   return 1;
	if((m_order[2]-1)==i)   return 2;
	return -1;
}
