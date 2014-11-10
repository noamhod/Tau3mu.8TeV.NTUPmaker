//////////////////////////////////////////////
//// run ./execute here and follow ///////////
//// the instructions ////////////////////////
//////////////////////////////////////////////
{
	TString p = gSystem->pwd();
	
	TString mode = ""; // flag to build/build+run/load+run the process
	
	TString type   = "";
	TString run    = "";
	TString outDir = "";
	TString chnl   = "";
	TString master = "";
	TString method = "";
	TString split  = "";
	Int_t nArgs    = 13;
	
	for(int i=0 ; i<gApplication->Argc() ; i++) printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
	if(gApplication->Argc()<nArgs)
	{
		cout << "gApplication->Argc() = " << gApplication->Argc() << endl;
		cout << "should be: " << nArgs  << " -> exit." << endl;
		exit(-1);
	}
	else
	{
		for(int i=0 ; i<gApplication->Argc() ; i++)
		{
			TString arg = gApplication->Argv(i);
			if(arg.Contains("--mode="))   mode   = arg.ReplaceAll("--mode=","");
			if(arg.Contains("--outDir=")) outDir = arg.ReplaceAll("--outDir=","");
			if(arg.Contains("--run="))    run    = arg.ReplaceAll("--run=","");
			if(arg.Contains("--type="))   type   = arg.ReplaceAll("--type=","");
			if(arg.Contains("--chnl="))   chnl   = arg.ReplaceAll("--chnl=","");
			if(arg.Contains("--master=")) master = arg.ReplaceAll("--master=","");
			if(arg.Contains("--method=")) method = arg.ReplaceAll("--method=","");
			if(arg.Contains("--split="))  split  = arg.ReplaceAll("--split=","");
		}
	}
	
	Int_t modeflag = -1;
	if     (mode=="build")    modeflag = 0;
	else if(mode=="load")     modeflag = 1;
	else if(mode=="buildrun") modeflag = 2;
	else if(mode=="loadrun")  modeflag = 3;
	else
	{
		cout << "build mode=" << mode << " is unknown" << endl;
		cout << "should be: build/buildrun/loadrun -> exit." << endl;
		exit(-1);
	}
	
	gROOT->ProcessLine(".include "+p+"/../PileupReweighting/");
	gROOT->ProcessLine(".include "+p+"/../PileupReweighting/PileupReweighting");
	gROOT->ProcessLine(".include "+p+"/../ApplyJetCalibration/");
	gROOT->ProcessLine(".include "+p+"/../ApplyJetCalibration/ApplyJetCalibration");
	gROOT->ProcessLine(".include "+p+"/../JetUncertainties/");
	gROOT->ProcessLine(".include "+p+"/../JetUncertainties/JetUncertainties");
	gROOT->ProcessLine(".include "+p+"/../JetResolution/");
	gROOT->ProcessLine(".include "+p+"/../JetResolution/JetResolution");
	gROOT->ProcessLine(".include "+p+"/../ApplyJetResolutionSmearing/");
	gROOT->ProcessLine(".include "+p+"/../ApplyJetResolutionSmearing/ApplyJetResolutionSmearing");
	gROOT->ProcessLine(".include "+p+"/../METAnalysisCommon/");
	gROOT->ProcessLine(".include "+p+"/../METAnalysisCommon/METAnalysisCommon");
	gROOT->ProcessLine(".include "+p+"/../METSystematics/");
	gROOT->ProcessLine(".include "+p+"/../METSystematics/METSystematics");
	gROOT->ProcessLine(".include "+p+"/../MissingETUtility/");
	gROOT->ProcessLine(".include "+p+"/../MissingETUtility/MissingETUtility");
	gROOT->ProcessLine(".include "+p+"/../PATCore/");
	gROOT->ProcessLine(".include "+p+"/../PATCore/PATCore/");
	gROOT->ProcessLine(".include "+p+"/../TileTripReader/");
	gROOT->ProcessLine(".include "+p+"/../TileTripReader/TileTripReader");
	gROOT->ProcessLine(".include "+p+"/../BCHCleaningTool/");
	gROOT->ProcessLine(".include "+p+"/../BCHCleaningTool/BCHCleaningTool");
	gROOT->ProcessLine(".include "+p);
	gSystem->Load( "libCintex.so" );
	Cintex::Cintex::Enable();
	
	if(modeflag%2==0)
	{
		cout << "===> Building" << endl;
		gROOT->ProcessLine(".L Loader.C+");
		gROOT->ProcessLine(".L ../PileupReweighting/StandAlone/PileupReweightingLib.so");
		gROOT->ProcessLine(".L ../ApplyJetCalibration/StandAlone/libApplyJetCalibration.so");
		gROOT->ProcessLine(".L ../JetUncertainties/StandAlone/libJetUncertainties.so");
		gROOT->ProcessLine(".L ../ApplyJetResolutionSmearing/StandAlone/libApplyJetResolutionSmearing.so");
		gROOT->ProcessLine(".L ../JetResolution/StandAlone/libJERProvider.so");
		gROOT->ProcessLine(".L ../METAnalysisCommon/StandAlone/libMETAnalysisCommon.so");
		gROOT->ProcessLine(".L ../METSystematics/StandAlone/libMETSystematics.so");
		gROOT->ProcessLine(".L ../MissingETUtility/StandAlone/libMETUtility.so");
		gROOT->ProcessLine(".L ../TileTripReader/StandAlone/libTTileTripReader.so");
		gROOT->ProcessLine(".L ../BCHCleaningTool/StandAlone/libBCHCleaningTool.so");
		gROOT->ProcessLine(".L "+type+".C++");
	}
	if(modeflag%2!=0)
	{
		cout << "===> Loading" << endl;
		gROOT->ProcessLine(".L Loader_C.so");
		gROOT->ProcessLine(".L ../PileupReweighting/StandAlone/PileupReweightingLib.so");
		gROOT->ProcessLine(".L ../ApplyJetCalibration/StandAlone/libApplyJetCalibration.so");
		gROOT->ProcessLine(".L ../JetUncertainties/StandAlone/libJetUncertainties.so");
		gROOT->ProcessLine(".L ../ApplyJetResolutionSmearing/StandAlone/libApplyJetResolutionSmearing.so");
		gROOT->ProcessLine(".L ../JetResolution/StandAlone/libJERProvider.so");
		gROOT->ProcessLine(".L ../METAnalysisCommon/StandAlone/libMETAnalysisCommon.so");
		gROOT->ProcessLine(".L ../METSystematics/StandAlone/libMETSystematics.so");
		gROOT->ProcessLine(".L ../MissingETUtility/StandAlone/libMETUtility.so");
		gROOT->ProcessLine(".L ../TileTripReader/StandAlone/libTTileTripReader.so");
		gROOT->ProcessLine(".L ../BCHCleaningTool/StandAlone/libBCHCleaningTool.so");
		gROOT->ProcessLine(".L "+type+"_C.so");
	}
	if(modeflag>1)
	{
		cout << "===> Running" << endl;
		TString args = "";
		args = run+"\",\""+outDir+"\",\""+chnl+"\",\""+master+"\",\""+method;
		args += (split!="") ? "\",\""+split : "";
		TString proc = type+"(\""+args+"\")";
		cout << proc << endl;
		gROOT->ProcessLine(proc);
	}
}
