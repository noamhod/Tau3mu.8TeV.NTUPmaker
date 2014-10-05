#include "rawROOT.h"
//#include "enums.h"


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
	double         vtxPt()             { return m_p4.Pt();      };
	double         vtxM()              { return m_p4.M();       };
	double         vtxMOS1()           { return m_pOS1.M();     };
	double         vtxMOS2()           { return m_pOS2.M();     };
	double         vtxMSS()            { return m_pSS.M();      };
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
	int            jetN()         { return m_njets;        };
	TLorentzVector jetPE(int i)   { return m_jetPE[i];     };
	TLorentzVector jetPM(int i)   { return m_jetPM[i];     };
	double         jetMV1(int i)  { return m_jetMV1[i];    };
	double         jetSumPt()     { return m_jetSumpt12;   };
	double         jetDphi12()    { return m_jetDphi12;    };
	double         jetDR12()      { return m_jetDR12;      };
	double         jetDphi3body() { return m_jetDphi3body; };
	double         jetDR3body()   { return m_jetDR3body;   };
	
	//// MET
	double         met()          { return m_met;          };
	double         metPhi()       { return m_metPhi;       };
	double         metDphi3body() { return m_metDphi3body; };
	double         metMt()        { return m_metMt;        };

private:
	int iSorted(int i);
	
public:
	
private:
	int m_index, m_code;
	TString m_type;
	TLorentzVector m_p4, m_pOS1, m_pOS2, m_pSS;
	TLorentzVector m_trkP[6];
	int m_njets;
	TLorentzVector m_jetPE[4], m_jetPM[4];
	double m_jetMV1[4];
	double m_jetSumpt12, m_jetDphi12, m_jetDphi3body, m_jetDR12, m_jetDR3body;
	bool m_ismuon[6], m_istp[6], m_istpa[6], m_istpb[6], m_iscalo[6];
	bool m_ismedium[6], m_istight[6], m_isloose[6], m_iscb[6];
	int m_isrc[6], m_itrk[6], m_trktype[6], m_order[6];
	string m_src[6];
	double m_trkChi2[6], m_trkChi2Ndf[6], m_trkPval[6], m_trkMuMatchChi2Ndf[6];
	double m_trkptfrac[6], m_ptFrac12, m_ptFrac13, m_ptFrac23, m_dpt12, m_dpt13, m_dpt23;
	double m_trkpbal[6], m_trksctang[6], m_trksctngb[6];
	double m_trkQoverP[6], m_srcQoverP[6], m_trkPixeldEdx[6];
	int m_trkNdf[6];
	int m_trkPIXhits[6], m_trkDeadPIX[6], m_trkPIXholes[6];
	int m_trkSCThits[6], m_trkDeadSCT[6], m_trkSCTholes[6];
	int m_trkTRThits[6], m_trkTRToutliers[6], m_trkHtTRThits[6];
	int m_trkUsedHitsdEdx[6];
	int m_trkMDThits[6], m_trkTGCPhiHits[6], m_trkTGCEtaHits[6], m_trkCSCPhiHits[6], m_trkCSCEtaHits[6], m_trkRPCPhiHits[6], m_trkRPCEtaHits[6];
	int m_trkCSCEtaHoles[6], m_trkCSCPhiHoles[6], m_trkRPCEtaHoles[6], m_trkRPCPhiHoles[6], m_trkMDTholes[6], m_trkTGCEtaHoles[6], m_trkTGCPhiHoles[6];
	int m_trkOutliersOnTrack[6], m_trkStdDevOfChi2OS[6];
	int m_trkPrecisionHits[6], m_trkPhiLayers[6], m_trkEtaPhiLayers[6], m_trkPrecisionHoles[6], m_trkEtaTriggerHoleLayers[6], m_trkPhiHoleLayers[6], m_trkPrecisionOutliers[6];
	
	double m_met, m_metPhi, m_metDphi3body, m_metMt;
	double m_chi2, m_ndf, m_chi2ndf, m_pvalue, m_lxy, m_lxyErr, m_tau;
	double m_a0, m_a0xy, m_cosT, m_cosTxy, m_charge, m_isolation[21], m_drmax, m_drmin;
};

int vertex::iSorted(int i)
{
	if(i>3)               return i;
	if((m_order[0]-1)==i) return 0;
	if((m_order[1]-1)==i) return 1;
	if((m_order[2]-1)==i) return 2;
	return -1;
}
