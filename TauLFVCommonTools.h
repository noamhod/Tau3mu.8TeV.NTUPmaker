#ifndef TauLFVCommonTools_TauLFVCommonTools_H
#define TauLFVCommonTools_TauLFVCommonTools_H

#include "rawStd.h"
#include "rawROOT.h"
#include "types.h"
#include "logs.h"
#include "enums.h"
#include "constants.h"
#include "colors.h"
#include "style.h"
#include "histos.h"
#include "kinematics.h"
#include "functions.h"

class TauLFVCommonTools
{
	public:
		int msgDBG;
		int msgINF;
		int msgWRN;
		int msgERR;
		int msgFAT;


		// this is a standard constructor
		TauLFVCommonTools ();
		virtual ~TauLFVCommonTools() { }
		static TauLFVCommonTools* getInstance(); // singleton
		
		// initializers
		virtual void initialize();
		virtual void initializeLOG();
		virtual void initializeBIN();
		
		// helpers
		virtual int  increment(int& counter);
		virtual bool findInVector(vector<unsigned int>* v, unsigned int x);
		virtual bool findInVector(vector<unsigned int>& v, unsigned int x);
	
	private:
		static TauLFVCommonTools* theInstance; // singleton
		
	// this is needed to distribute the algorithm to the workers
	ClassDef(TauLFVCommonTools, 1);
};

#endif
