#ifndef constants_H
#define constants_H

static const double luminosity   = 3.2531; // periodD fb-1
static const double Ltot         = 20.3;   // full 2012 fb-1
static const double tau3muBR     = 2.1e-8; // Belle limit

static const float MeV2TeV = 1.e-6;
static const float MeV2GeV = 1.e-3;
static const float GeV2TeV = 1.e-3;
static const float GeV2MeV = 1.e+3;
static const float TeV2MeV = 1.e+6;
static const float TeV2GeV = 1.e+3;

static const double mb2fb = 1.e12;
static const double nb2fb = 1.e+6;
static const double pb2fb = 1.e+3;
static const double nb2mb = 1.e-6;

static const float muonMass    = 0.105658367; // GeV
static const float muonMassMeV = 105.658367; // MeV
static const float tauMassGeV  = 1.77682; // GeV
static const float tauMassMeV  = 1776.82; // MeV
static const float elecMassGeV = 0.00051099891; // GeV
static const float elecMassMeV = 0.51099891; // MeV

static const double pi   = 3.14159265;
static const double sq2 = sqrt(2.);
static const double f12 = 1./2.;
static const double f13 = 1./3.;
static const double f23 = 2./3.;

// static const double mZ0 = 91.187; // GeV ->PYTHIA6
static const double mZ0    = 91.1876;   // GeV ->PDG
// static const double wZ0 = 2.489;  // GeV ->PYTHIA6
static const double wZ0    = 2.4952;    // GeV ->PDG

static const double mZ2 = mZ0*mZ0;   // GeV^2

// static const double sw2  = 0.23;      // ->PYTHIA6
static const double sw2     = 0.2312015; // ->PDG
// static const double cw2  = 0.74;      // ->PYTHIA6
static const double cw2     = 0.7687985; // ->PDG
// static const double alphaEM = 0.00729735; // ->PYTHIA6
static const double alphaEM    = 0.00729735; // ->PDG

static const double Gmu     = 0.0000116633980690699; // GeV^-2

static const double GeV2mb  = 0.38937930419;      // 1/GeV^2 to mb (mili-barn) conversion constant
static const double GeV2nb  = 1e6*0.38937930419;  // 1/GeV^2 to nb (nano-barn) conversion constant
static const double GeV2pb  = 1e9*0.38937930419;  // 1/GeV^2 to pb (pico-barn) conversion constant
static const double GeV2fb  = 1e12*0.38937930419;  // 1/GeV^2 to fb (femto-barn) conversion constant
/*-----------------------------------------------------------------------------------------*/
/*                                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* [1] == (hbar*c)^2 = 0.38937930419 [GeV^2 * mb] -> !!! [1/GeV^2] = 0.38937930419[mb] !!! */
/*                                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* sigma ~ [mb] ~ [1/GeV^2]                                                                */ 
/* d(sigms)/ds ~ [mb/GeV^2] ~ [1/GeV^4]                                                    */
/* d(sigma)/d(sqrt(s) == d(sigms)/ds * 2sqrt(s) ~ [mb/GeV] ~ [1/GeV^3]                     */
/* 1[mb] = 1e6[nb]                                                                         */
/*-----------------------------------------------------------------------------------------*/

#endif

