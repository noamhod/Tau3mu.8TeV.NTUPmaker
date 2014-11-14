#include "rawROOT.h"
//#include "enums.h"

static const unsigned int nMaxTracks = 6; // including the triplet, i.e. 3+3 here

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
	VTX_JERUP, VTX_JERDWN
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
	double         vtxA0()             { return m_a0;           };
	double         vtxA0xy()           { return m_a0xy;         };
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
	int            jetN()             { return m_njets;           };
	TLorentzVector jetPE(int i)       { return m_jetPE[i];        };
	TLorentzVector jetPM(int i)       { return m_jetPM[i];        };
	double         jetShiftJES(int i) { return m_jet_shiftJES[i]; };
	double         jetShiftJER(int i) { return m_jet_shiftJER[i]; };
	TLorentzVector jetPE_JES(int i, int shift=VTX_NOSHIFT) { return m_jetPE[i]*(1+shift*jetShiftJES(i)); };
	TLorentzVector jetPE_JER(int i, int shift=VTX_NOSHIFT) { return m_jetPE[i]*(1+shift*jetShiftJER(i)); };
	TLorentzVector jetPM_JES(int i, int shift=VTX_NOSHIFT) { return m_jetPM[i]*(1+shift*jetShiftJES(i)); };
	TLorentzVector jetPM_JER(int i, int shift=VTX_NOSHIFT) { return m_jetPM[i]*(1+shift*jetShiftJER(i)); };
	double         jetMV1(int i)      { return m_jetMV1[i];      };
	double         jetVtxFrac(int i)  { return m_jetVtxFrac[i];  };
	
	//// jets+X
	// double         jetSumPt()         { return m_jetSumpt12;     };
	// double         jetDphi12()        { return m_jetDphi12;      };
	// double         jetDR12()          { return m_jetDR12;        };
	// double         jetDphi3body()     { return m_jetDphi3body;   };
	// double         jetDR3body()       { return m_jetDR3body;     };
	double jetSumPt(int shift=VTX_NOMINAL)         
	{
		switch(shift)
		{
			case VTX_JESUP:   return m_jetSumpt12_jes_up;
			case VTX_JESDWN:  return m_jetSumpt12_jes_dwn;
			case VTX_JERUP:   return m_jetSumpt12_jer_up;
			case VTX_JERDWN:  return m_jetSumpt12_jer_dwn;
			case VTX_NOMINAL: return m_jetSumpt12;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	double jetDphi12(int shift=VTX_NOMINAL)         
	{
		switch(shift)
		{
			case VTX_JESUP:   return m_jetDphi12_jes_up;
			case VTX_JESDWN:  return m_jetDphi12_jes_dwn;
			case VTX_JERUP:   return m_jetDphi12_jer_up;
			case VTX_JERDWN:  return m_jetDphi12_jer_dwn;
			case VTX_NOMINAL: return m_jetDphi12;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	double jetDR12(int shift=VTX_NOMINAL)         
	{
		switch(shift)
		{
			case VTX_JESUP:   return m_jetDR12_jes_up;
			case VTX_JESDWN:  return m_jetDR12_jes_dwn;
			case VTX_JERUP:   return m_jetDR12_jer_up;
			case VTX_JERDWN:  return m_jetDR12_jer_dwn;
			case VTX_NOMINAL: return m_jetDR12;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	double jetDphi3body(int shift=VTX_NOMINAL)         
	{
		switch(shift)
		{
			case VTX_JESUP:   return m_jetDphi3body_jes_up;
			case VTX_JESDWN:  return m_jetDphi3body_jes_dwn;
			case VTX_JERUP:   return m_jetDphi3body_jer_up;
			case VTX_JERDWN:  return m_jetDphi3body_jer_dwn;
			case VTX_NOMINAL: return m_jetDphi3body;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	double jetDR3body(int shift=VTX_NOMINAL)         
	{
		switch(shift)
		{
			case VTX_JESUP:   return m_jetDR3body_jes_up;
			case VTX_JESDWN:  return m_jetDR3body_jes_dwn;
			case VTX_JERUP:   return m_jetDR3body_jer_up;
			case VTX_JERDWN:  return m_jetDR3body_jer_dwn;
			case VTX_NOMINAL: return m_jetDR3body;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	
	
	//// MET
	// double         met()          { return m_met;          };
	// double         metPhi()       { return m_metPhi;       };
	// double         metDphi3body() { return m_metDphi3body; };
	// double         metMt()        { return m_metMt;        };
	double met(int shift=VTX_NOJES)
	{
		switch(shift)
		{
			case VTX_NOJES:   return m_met;
			case VTX_NOMINAL: return m_met_jes_nominal;
			case VTX_JESUP:   return m_met_jes_up;
			case VTX_JESDWN:  return m_met_jes_dwn;
			case VTX_JERUP:   return m_met_jer_up;
			case VTX_JERDWN:  return m_met_jer_dwn;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	double metPhi(int shift=VTX_NOJES)
	{
		switch(shift)
		{
			case VTX_NOJES:   return m_metPhi;
			case VTX_NOMINAL: return m_metPhi_jes_nominal;
			case VTX_JESUP:   return m_metPhi_jes_up;
			case VTX_JESDWN:  return m_metPhi_jes_dwn;
			case VTX_JERUP:   return m_metPhi_jer_up;
			case VTX_JERDWN:  return m_metPhi_jer_dwn;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	double metDphi3body(int shift=VTX_NOJES)
	{
		switch(shift)
		{
			case VTX_NOJES:   return m_metDphi3body;
			case VTX_NOMINAL: return m_metDphi3body_jes_nominal;
			case VTX_JESUP:   return m_metDphi3body_jes_up;
			case VTX_JESDWN:  return m_metDphi3body_jes_dwn;
			case VTX_JERUP:   return m_metDphi3body_jer_up;
			case VTX_JERDWN:  return m_metDphi3body_jer_dwn;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};
	double metMt(int shift=VTX_NOJES)
	{
		switch(shift)
		{
			case VTX_NOJES:   return m_metMt;
			case VTX_NOMINAL: return m_metMt_jes_nominal;
			case VTX_JESUP:   return m_metMt_jes_up;
			case VTX_JESDWN:  return m_metMt_jes_dwn;
			case VTX_JERUP:   return m_metMt_jer_up;
			case VTX_JERDWN:  return m_metMt_jer_dwn;
			default: break; // _FATAL("Unsupported shit code: "+_s(shit));
		}
		return -999.;
	};

private:
	int iSorted(int i);
	
public:
	
private:
	int m_index, m_code;
	TString m_type;
	TLorentzVector m_p4, m_pOS1, m_pOS2, m_pSS, m_pQuad[nMaxTracks-3];
	TLorentzVector m_trkP[nMaxTracks];
	int m_njets;
	TLorentzVector m_jetPE[4], m_jetPM[4];
	double m_jet_shiftJES[4], m_jet_shiftJER[4];
	double m_jetMV1[4], m_jetVtxFrac[4];
	double m_jetSumpt12, m_jetDphi12, m_jetDphi3body, m_jetDR12, m_jetDR3body;
	double m_jetSumpt12_jes_up,  m_jetDphi12_jes_up,  m_jetDphi3body_jes_up,  m_jetDR12_jes_up,  m_jetDR3body_jes_up;
	double m_jetSumpt12_jes_dwn, m_jetDphi12_jes_dwn, m_jetDphi3body_jes_dwn, m_jetDR12_jes_dwn, m_jetDR3body_jes_dwn;
	double m_jetSumpt12_jer_up,  m_jetDphi12_jer_up,  m_jetDphi3body_jer_up,  m_jetDR12_jer_up,  m_jetDR3body_jer_up;
	double m_jetSumpt12_jer_dwn, m_jetDphi12_jer_dwn, m_jetDphi3body_jer_dwn, m_jetDR12_jer_dwn, m_jetDR3body_jer_dwn;
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
	
	double m_met, m_metPhi, m_metDphi3body, m_metMt;
	double m_met_jes_nominal,  m_metPhi_jes_nominal,  m_metDphi3body_jes_nominal,  m_metMt_jes_nominal;
	double m_met_jes_up,  m_metPhi_jes_up,  m_metDphi3body_jes_up,  m_metMt_jes_up;
	double m_met_jes_dwn, m_metPhi_jes_dwn, m_metDphi3body_jes_dwn, m_metMt_jes_dwn;
	double m_met_jer_up,  m_metPhi_jer_up,  m_metDphi3body_jer_up,  m_metMt_jer_up;
	double m_met_jer_dwn, m_metPhi_jer_dwn, m_metDphi3body_jer_dwn, m_metMt_jer_dwn;
	double m_chi2, m_ndf, m_chi2ndf, m_pvalue, m_lxy, m_lxyErr, m_tau;
	double m_a0, m_a0xy, m_cosT, m_cosTxy, m_charge, m_isolation[21], m_drmax, m_drmin;
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
