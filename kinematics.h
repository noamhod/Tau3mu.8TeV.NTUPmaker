#ifndef kinematics_H
#define kinematics_H

namespace
{
#include "rawROOT.h"

TLorentzVector m_pTmp;

inline float pT(   TLorentzVector* p )     { return p->Perp(); }
inline float pT(   float qOp, float theta) { return (qOp!=0) ? fabs(1./qOp)*sin(theta) : 0.; }
inline float eta(  TLorentzVector* p )     { return p->PseudoRapidity(); }
inline float pAbs( TLorentzVector* p )     { return p->Mag(); }   
inline float y(    TLorentzVector* p )     { return p->Rapidity(); }

inline float chi(   TLorentzVector* pa, TLorentzVector* pb ) { return TMath::Exp(fabs(pa->Rapidity()-pb->Rapidity())); }
inline float yStar( TLorentzVector* pa, TLorentzVector* pb ) { return (pa->Rapidity()-pb->Rapidity())/2.; }
inline float yB(    TLorentzVector* pa, TLorentzVector* pb ) { return (pa->Rapidity()+pb->Rapidity())/2.; }
inline float dR(    TLorentzVector* pa, TLorentzVector* pb ) { return pa->DeltaR(*pb); }

inline TVector3 systemBoostVector(TLorentzVector* pa, TLorentzVector* pb)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	return pTmp.BoostVector();
}
inline TLorentzVector* Boost(TLorentzVector* pBoost, TLorentzVector* p)
{
	TLorentzVector* pBoosted = (TLorentzVector*)p->Clone("");
	pBoosted->Boost(-1.*pBoost->BoostVector());
	return pBoosted;
}
inline TLorentzVector* Boost(TLorentzVector* pa, TLorentzVector* pb, TLorentzVector* p)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	TLorentzVector* pBoosted = (TLorentzVector*)p->Clone("");
	pBoosted->Boost(-1.*pTmp.BoostVector());
	return pBoosted;
}
inline float cosThetaDimu(TLorentzVector* pa, TLorentzVector* pb)
{
	float paMag = sqrt( pa->Px()*pa->Px() + pa->Py()*pa->Py() + pa->Pz()*pa->Pz() );
	float pbMag = sqrt( pb->Px()*pb->Px() + pb->Py()*pb->Py() + pb->Pz()*pb->Pz() );	
	float costh = (pa->Px()*pb->Px() + pa->Py()*pb->Py() + pa->Pz()*pb->Pz()) / (paMag*pbMag);
	return costh;
}
inline float imass(TLorentzVector* pa, TLorentzVector* pb)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	return pTmp.M();
}
inline float QT(TLorentzVector* pa, TLorentzVector* pb)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	return pTmp.Perp();
} 
inline float ySystem(TLorentzVector* pa, TLorentzVector* pb)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	return pTmp.Rapidity();
}
inline float etaSystem(TLorentzVector* pa, TLorentzVector* pb)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	return pTmp.PseudoRapidity();
}
inline float magBetaSystem(TLorentzVector* pa, TLorentzVector* pb)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	return pTmp.Beta();
}
inline float betaiSystem(TLorentzVector* pa, TLorentzVector* pb, int i)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	float E = pTmp.E();
	if(E<=0.) return -99999999.;

	if     (i==1) return pTmp.Px()/E;
	else if(i==2) return pTmp.Py()/E;
	else if(i==3) return pTmp.Pz()/E;
	else _WARNING("i needs to be 1,2,3 (x,y,z), returning -99999999.");
	return -99999999.;
}
inline float betaTSystem(TLorentzVector* pa, TLorentzVector* pb)
{
	TLorentzVector pTmp = (*pa)+(*pb);
	return (pTmp.E()!=0.) ? pTmp.Pt()/pTmp.E() : -999.;
}
inline float magBetaSystem(float pa, float pb, float ea, float eb)
{
	return ((ea+eb)>0.) ? (pa+pb)/(ea+eb) : -999999.;
}
inline float cosThetaCollinsSoper(TLorentzVector* pa, float ca, TLorentzVector* pb, float cb)
{
	// this will work only for leptons e, mu and tau
	// by default it is assumed that pa is the lepton
	// if instead pb is the lepton, then the result is
	// reclassified by a (-) sign - see line 4.
	float mass2 = imass(pa,pb)*imass(pa,pb);
	float QT2 = QT(pa,pb)*QT(pa,pb);
	//float cosThetaCS = 2.*( pa->Plus()*pb->Minus() - pa->Minus()*pb->Plus() ) / sqrt( mass2 * (mass2 + QT2) );
	float cosThetaCS = 2.*( pa->Pz()*pb->E() - pa->E()*pb->Pz() ) / sqrt( mass2 * (mass2 + QT2) );
	if (ca>0. && cb<0.)     cosThetaCS *= -1.; // if pb is the lepton
	if (ySystem(pa,pb) < 0) cosThetaCS *= -1.; // reclassify
	return cosThetaCS;
}
inline float cosThetaBoost(TLorentzVector* pa, float ca, TLorentzVector* pb, float cb)
{
	// http://xrootd.slac.stanford.edu/BFROOT/www/doc/workbook_backup_010108/analysis/analysis.html
	// A useful quantity in many analyses is the helicity angle.
	// In the reaction Y -> X -> a + b, the helicity angle of 
	// particle a is the angle measured in the rest frame of the
	//decaying parent particle, X, between the direction of the
	// decay daughter a and the direction of the grandparent particle Y.

	TLorentzVector pTmp = (*pa)+(*pb); // this is the mumu system (Z) 4vector
	TVector3 ZboostVector = pTmp.BoostVector(); // this is the 3vector of the Z
	TLorentzVector p; // this is the muon 4vector
	
	if(ca<0)      p.SetPxPyPzE(pa->Px(),pa->Py(),pa->Pz(),pa->E());
	else if(cb<0) p.SetPxPyPzE(pb->Px(),pb->Py(),pb->Pz(),pb->E());
	p.Boost( -ZboostVector ); // boost p to the dimuon CM (rest) frame
	float cosThetaB = p.Vect()*pTmp.Vect()/(p.P()*pTmp.P());
	//if (ySystem(pa,pb) < 0) cosThetaB *= -1.; // reclassify ???
	return cosThetaB;
}
inline float cosThetaBoost(vector<TLorentzVector>& p, unsigned int index)
{
	TLorentzVector pTmp = p[0];
	for(unsigned int n=1 ; n<p.size() ; n++) pTmp = pTmp + p[n];
	TVector3 boost = pTmp.BoostVector();
	
	TLorentzVector pStar; 
	pStar.SetPtEtaPhiM(p[index].Pt(),p[index].Eta(),p[index].Phi(),p[index].M()); // not yet in CMF
	pStar.Boost( -boost ); // boosted to CMF
	
	float cosTheta = pStar.Vect()*pTmp.Vect()/(pStar.P()*pTmp.P());
	return cosTheta;
}

}

#endif

