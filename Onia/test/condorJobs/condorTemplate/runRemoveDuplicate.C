 #include "TSystem.h"
 #include "TChain.h"
 #include "../RemoveDuplicate.C"
 void runRemoveDuplicate_NUMBER()
 {
    TChain * chain = new TChain("rootuple/oniaTree","");
	 chain->Add("INPUTPATH/INPUTFILE");
	 RemoveDuplicate a(chain);
    a.Loop("/eos/uscms/store/user/zhenhu/DoubleMuon/Onia2MuMuRootuple-Run2015B-DoubleMuon-v3/f64ae4aafa8965110120322f8d8de3c5/RootupleAll.root","RemoveDuplicate_NUMBER.root");

 }
