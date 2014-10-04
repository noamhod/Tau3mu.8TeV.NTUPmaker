#ifndef TauLFVCommonTools_enums_H
#define TauLFVCommonTools_enums_H

enum pdtEnum { // *** note that this is not a complete list of pdt's particles ***
	PDTDWN=1, PDTUP=2, PDTSTR=3, PDTCHM=4, PDTBOT=5, PDTTOP=6,
	PDTBOTPRIME=7, PDTTOPPRIME=8,
	PDTE=11, PDTNUE=12, PDTMU=13, PDTNUMU=14, PDTTAU=15, PDTNUTAU=16,
	PDTTAUPRIME=17, PDTNUTAUPRIME=18,
	PDTGLU=21, PDTGAMMA=22, PDTZ=23, PDTWPLUS=24, PDTZPRIME0=32,PDTGRV=5000039,PDTKK=5000023,
	PDTPI0=111, PDTPIPLUS=211, PDTETA=221, PDTRHO0=113, PDTRHOPLUS=213,
	PDTK0L=130, PDTK0S=310, PDTK0=311, PDTKPLUS=321,
	PDTOMEGA=223, PDTPHI=333,
	PDTJPSI1S=443,
	PDTFIRSTMESON=111, PDTLASTMESON=555,
	PDTFIRSTBARYON=1103, PDTLASTBARYON=5554
};

enum skimEnum {
	SKIMOBJECT, SKIMSYSTEM, SKIMTRUTH
};

enum objecttype
{
	INT,FLT,DBL,STR,
	VINT,VFLT,VDBL,VSTR,
	VVINT,VVFLT,VVDBL,VVSTR
};

enum genealsourcecode
{
	MUON, TPA, TPB, CALOMUON, TPC // general TPs
};

enum sourcecode
{
	CALO, MUONS, MUID, // muon...
	MUONSTPA, MUONSTPB, MUONSTPC, // muons chain TPs
	MUIDTPA, MUIDTPB, MUIDTPC // muid chain TPs
};

enum categorycodes
{
	MUONS3,
	MUONS2TPA1, MUONS2TPB1, MUONS2CALO1,
	MUONS1TPA2, MUONS1TPB2,
	MUONS1TPA1TPB1,
	MUONS0TPA3, MUONS0TPB3,
	MUONS0TPA2TPB1, MUONS0TPA2CALO1,
	MUONS0TPA1TPB2,
	
	MUID3,
	MUID2TPA1, MUID2TPB1, MUID2CALO1,
	MUID1TPA2, MUID1TPB2,
	MUID1TPA1TPB1,
	MUID0TPA3, MUID0TPB3,
	MUID0TPA2TPB1, MUID0TPA2CALO1,
	MUID0TPA1TPB2
};

inline bool isBhadron(unsigned int pdgId)
{	
	bool isB = false;
	switch(pdgId)
	{
		// b quark
		//case 5:   isB = true; break;
		
		// B mesons
		case 511:   isB = true; break;
		case 521:   isB = true; break;
		case 10511: isB = true; break;
		case 10521: isB = true; break;
		case 513:   isB = true; break;
		case 523:   isB = true; break;
		case 10513: isB = true; break;
		case 10523: isB = true; break;
		case 20513: isB = true; break;
		case 20523: isB = true; break;
		case 515:   isB = true; break;
		case 525:   isB = true; break;
		case 531:   isB = true; break;
		case 10531: isB = true; break;
		case 533:   isB = true; break;
		case 10533: isB = true; break;
		case 20533: isB = true; break;
		case 535:   isB = true; break;
		case 541:   isB = true; break;
		case 10541: isB = true; break;
		case 543:   isB = true; break;
		case 10543: isB = true; break;
		case 20543: isB = true; break;
		case 545:   isB = true; break;
		
		// bbar mesons
		case 551:     isB = true; break;
		case 10551:   isB = true; break;
		case 100551:  isB = true; break;
		case 110551:  isB = true; break;
		case 200551:  isB = true; break;
		case 210551:  isB = true; break;
		case 553:     isB = true; break;
		case 10553:   isB = true; break;
		case 20553:   isB = true; break;
		case 30553:   isB = true; break;
		case 100553:  isB = true; break;
		case 110553:  isB = true; break;
		case 120553:  isB = true; break;
		case 130553:  isB = true; break;
		case 200553:  isB = true; break;
		case 210553:  isB = true; break;
		case 220553:  isB = true; break;
		case 300553:  isB = true; break;
		case 9000553: isB = true; break;
		case 9010553: isB = true; break;
		case 555:     isB = true; break;
		case 10555:   isB = true; break;
		case 20555:   isB = true; break;
		case 100555:  isB = true; break;
		case 110555:  isB = true; break;
		case 120555:  isB = true; break;
		case 200555:  isB = true; break;
		case 557:     isB = true; break;
		case 100557:  isB = true; break;
		
		// B baryons
		case 5122: isB = true; break;
		case 5112: isB = true; break;
		case 5212: isB = true; break;
		case 5222: isB = true; break;
		case 5114: isB = true; break;
		case 5214: isB = true; break;
		case 5224: isB = true; break;
		case 5132: isB = true; break;
		case 5232: isB = true; break;
		case 5312: isB = true; break;
		case 5322: isB = true; break;
		case 5314: isB = true; break;
		case 5324: isB = true; break;
		case 5332: isB = true; break;
		case 5334: isB = true; break;
		case 5142: isB = true; break;
		case 5242: isB = true; break;
		case 5412: isB = true; break;
		case 5422: isB = true; break;
		case 5414: isB = true; break;
		case 5424: isB = true; break;
		case 5342: isB = true; break;
		case 5432: isB = true; break;
		case 5434: isB = true; break;
		case 5442: isB = true; break;
		case 5444: isB = true; break;
		case 5512: isB = true; break;
		case 5522: isB = true; break;
		case 5514: isB = true; break;
		case 5524: isB = true; break;
		case 5532: isB = true; break;
		case 5534: isB = true; break;
		case 5542: isB = true; break;
		case 5544: isB = true; break;
		case 5554: isB = true; break;
		
		default:   isB = false; break;
	}
	return isB;
}


inline bool isDhadron(unsigned int pdgId)
{	
	bool isD = false;
	switch(pdgId)
	{
		// c quark
		//case 4:   isD = true; break;
		
		// D mesons
		case 511:   isD = true; break;
		case 411:   isD = true; break;
		case 421:   isD = true; break;
		case 10411: isD = true; break;
		case 10421: isD = true; break;
		case 413:   isD = true; break;
		case 423:   isD = true; break;
		case 10413: isD = true; break;
		case 10423: isD = true; break;
		case 20413: isD = true; break;
		case 20423: isD = true; break;
		case 415:   isD = true; break;
		case 425:   isD = true; break;
		case 431:   isD = true; break;
		case 10431: isD = true; break;
		case 433:   isD = true; break;
		case 10433: isD = true; break;
		case 20433: isD = true; break;
		case 435:   isD = true; break;
		
		// ccbar mesons
		case 441:     isD = true; break;
		case 10441:   isD = true; break;
		case 100441:  isD = true; break;
		case 443:     isD = true; break;
		case 10443:   isD = true; break;
		case 20443:   isD = true; break;
		case 100443:  isD = true; break;
		case 30443:   isD = true; break;
		case 9000443: isD = true; break;
		case 9010443: isD = true; break;
		case 9020443: isD = true; break;
		case 445:     isD = true; break;
		case 9000445: isD = true; break;
		
		// D baryons
		case 4122: isD = true; break;
		case 4222: isD = true; break;
		case 4212: isD = true; break;
		case 4112: isD = true; break;
		case 4224: isD = true; break;
		case 4214: isD = true; break;
		case 4114: isD = true; break;
		case 4232: isD = true; break;
		case 4132: isD = true; break;
		case 4322: isD = true; break;
		case 4312: isD = true; break;
		case 4324: isD = true; break;
		case 4314: isD = true; break;
		case 4332: isD = true; break;
		case 4334: isD = true; break;
		case 4412: isD = true; break;
		case 4422: isD = true; break;
		case 4414: isD = true; break;
		case 4424: isD = true; break;
		case 4432: isD = true; break;
		case 4434: isD = true; break;
		case 4444: isD = true; break;
		
		default:   isD = false; break;
	}
	return isD;
}




#endif

