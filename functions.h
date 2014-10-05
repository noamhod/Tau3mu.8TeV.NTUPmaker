#ifndef functions_H
#define functions_H

namespace
{

Double_t fdijet(Double_t* v, Double_t* par)
{
	Double_t x = v[0];
	
	// DIJET: a0 * x^a1 * x^[a2*log(x)]
	return par[0]*TMath::Power(x, par[1]) * TMath::Power(x, par[2]*TMath::Log(x));
}

}

#endif
