// -*- C++ -*-
//
// Package:    MuMuGammaRootupler
// Class:      MuMuGammaRootupler
// 
// Description: Dump mumugamma decays
//
// Author:  Zhen Hu
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

//
// class declaration
//

class MuMuGammaRootupler:public edm::EDAnalyzer {
	public:
		explicit MuMuGammaRootupler(const edm::ParameterSet &);
		~MuMuGammaRootupler();

		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

	private:
		UInt_t getTriggerBits(const edm::Event &);
		bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
		const  reco::Candidate* GetAncestor(const reco::Candidate *);

		virtual void beginJob();
		virtual void analyze(const edm::Event &, const edm::EventSetup &);
		virtual void endJob();

		virtual void beginRun(edm::Run const &, edm::EventSetup const &);
		virtual void endRun(edm::Run const &, edm::EventSetup const &);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
		virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);

		// ----------member data ---------------------------
		std::string file_name;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> conversion_Label;
		edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
		edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
		int  pdgid_;
		std::vector<double> OniaMassCuts_;
		bool isMC_;
		bool OnlyBest_;
		bool OnlyGen_;
		double upsilon_mass_;
		uint32_t triggerCuts_;

		std::vector<TH1F*> myRefitY1sConversionmassVeryLoose;
		std::vector<TH1F*> myRefitY1sConversionmassLoose;
		std::vector<TH1F*> myRefitY1sConversionmassTight;

		TH1F* myConversionPt_all;
		TH1F* myConversionMass_all;
		TH1F* myConversionRadius_all;
		TH2F* myConversionVertex_all;
		TH1F* muGammaDeltaP;
		TH1F* muGammaDeltaR;
		TH1F* myDimuonMass_all;
		TH1F* myY1SFitMass_all;
		TH1F* mydz_photon_refit_all;
		TH1F* mydz_before_refit_all;
		TH1F* myChiBVtxP_all;

		TH1F* myConversionPt_used;
		TH1F* myConversionMass_used;
		TH1F* myConversionRadius_used;
		TH2F* myConversionVertex_used;
		TH1F* myDimuonMass_used;
		TH1F* myY1SFitMass_used;
		TH1F* mydz_photon_refit_used;
		TH1F* mydz_before_refit_used;
		TH1F* myChiBVtxP_used;

		UInt_t run;
		UInt_t lumi;
		UInt_t event;
		Int_t  irank;
		UInt_t trigger; 

		TLorentzVector dimuon_p4;
		TLorentzVector muonP_p4;
		TLorentzVector muonM_p4;
		Float_t MassErr;
		Float_t vProb;
		Float_t DCA;
		Float_t ppdlPV;
		Float_t ppdlErrPV;
		Float_t cosAlpha;

		Int_t numPrimaryVertices;

		Float_t JmumuMass;
		Float_t JmumuMassErr;
		Float_t JmumuPx;
		Float_t JmumuPy;
		Float_t JmumuPz;
		Float_t JmumuVtxCL;
		Float_t JmumuVtxCL2;
		Float_t JmumuDecayVtxX;
		Float_t JmumuDecayVtxY;
		Float_t JmumuDecayVtxZ;
		Float_t JmumuDecayVtxXE;
		Float_t JmumuDecayVtxYE;
		Float_t JmumuDecayVtxZE;
		TLorentzVector Jmumu_p4;

		Float_t deltaR;
		Float_t deltaP;
		Float_t ChiBM_fit;
		Float_t ChiBM_fitErr;
		Float_t ChiBPx_fit;
		Float_t ChiBPy_fit;
		Float_t ChiBPz_fit;
		Float_t ChiBVtxX_fit;
		Float_t ChiBVtxY_fit;
		Float_t ChiBVtxZ_fit;
		Float_t ChiBVtxP_fit;
		Float_t dz_photon_refit;
		Float_t dz_before_refit;
		Float_t mu1M_fit;
		Float_t mu1Px_fit;
		Float_t mu1Py_fit;
		Float_t mu1Pz_fit;
		Float_t mu2M_fit;
		Float_t mu2Px_fit;
		Float_t mu2Py_fit;
		Float_t mu2Pz_fit;
		Float_t gammaM_fit;
		Float_t gammaPx_fit;
		Float_t gammaPy_fit;
		Float_t gammaPz_fit;
		TLorentzVector mu1_p4_fit;
		TLorentzVector mu2_p4_fit;
		TLorentzVector gamma_p4_fit;
		TLorentzVector gamma_scaledp4_fit;
		Float_t mymumugammamassY1s;

		Float_t conversionMass;
		Float_t conversionPt;
		Float_t conversionEta;
		Float_t conversionVertexRho;
		Float_t conversionVertexX;
		Float_t conversionVertexY;
		Float_t conversionVertexZ;
		TLorentzVector conversion_p4;

		TTree *onia_tree;

		Int_t mother_pdgId;
		Int_t dimuon_pdgId;
		TLorentzVector gen_dimuon_p4;
		TLorentzVector gen_muonP_p4;
		TLorentzVector gen_muonM_p4;

		edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
		edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constructors and destructor
//

MuMuGammaRootupler::MuMuGammaRootupler(const edm::ParameterSet & iConfig):
	dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
	conversion_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("conversions"))),
	primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
	triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
	OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
	isMC_(iConfig.getParameter<bool>("isMC")),
	OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
	OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
	upsilon_mass_(iConfig.getParameter<double>("upsilon_mass")),
	triggerCuts_(iConfig.getParameter<uint32_t>("triggerCuts"))
{
	edm::Service < TFileService > fs;
	onia_tree = fs->make < TTree > ("oniaTree", "Tree of MuMuGamma");

	if (!OnlyGen_) {
		onia_tree->Branch("run",     &run,     "run/I");
		onia_tree->Branch("lumi",     &lumi,     "lumi/I");
		onia_tree->Branch("event",   &event,   "event/I");
		onia_tree->Branch("irank",   &irank,   "irank/I");
		onia_tree->Branch("trigger", &trigger, "trigger/I");

		onia_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
		onia_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
		onia_tree->Branch("muonM_p4",  "TLorentzVector", &muonM_p4);
		onia_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
		onia_tree->Branch("vProb",     &vProb,      "vProb/F");
		onia_tree->Branch("DCA",       &DCA,        "DCA/F");
		onia_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
		onia_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
		onia_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");

		onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");

		onia_tree->Branch("JmumuMass",&JmumuMass,"JmumuMass/F");
		onia_tree->Branch("JmumuMassErr",&JmumuMassErr,"JmumuMassErr/F");
		onia_tree->Branch("JmumuPx",&JmumuPx,"JmumuPx/F");
		onia_tree->Branch("JmumuPy",&JmumuPy,"JmumuPy/F");
		onia_tree->Branch("JmumuPz",&JmumuPz,"JmumuPz/F");
		onia_tree->Branch("JmumuVtxCL",&JmumuVtxCL,"JmumuVtxCL/F");
		onia_tree->Branch("JmumuVtxCL2",&JmumuVtxCL2,"JmumuVtxCL2/F");
		onia_tree->Branch("JmumuDecayVtxX",&JmumuDecayVtxX,"JmumuDecayVtxX/F");
		onia_tree->Branch("JmumuDecayVtxY",&JmumuDecayVtxY,"JmumuDecayVtxY/F");
		onia_tree->Branch("JmumuDecayVtxZ",&JmumuDecayVtxZ,"JmumuDecayVtxZ/F");
		onia_tree->Branch("JmumuDecayVtxXE",&JmumuDecayVtxXE,"JmumuDecayVtxXE/F");
		onia_tree->Branch("JmumuDecayVtxYE",&JmumuDecayVtxYE,"JmumuDecayVtxYE/F");
		onia_tree->Branch("JmumuDecayVtxZE",&JmumuDecayVtxZE,"JmumuDecayVtxZE/F");
		onia_tree->Branch("Jmumu_p4",  "TLorentzVector", &Jmumu_p4);

		onia_tree->Branch("ChiBM_fit",&ChiBM_fit,"ChiBM_fit/F");
		onia_tree->Branch("ChiBM_fitErr",&ChiBM_fitErr,"ChiBM_fitErr/F");
		onia_tree->Branch("ChiBPx_fit",&ChiBPx_fit,"ChiBPx_fit/F");
		onia_tree->Branch("ChiBPy_fit",&ChiBPy_fit,"ChiBPy_fit/F");
		onia_tree->Branch("ChiBPz_fit",&ChiBPz_fit,"ChiBPz_fit/F");
		onia_tree->Branch("ChiBVtxX_fit",&ChiBVtxX_fit,"ChiBVtxX_fit/F");
		onia_tree->Branch("ChiBVtxY_fit",&ChiBVtxY_fit,"ChiBVtxY_fit/F");
		onia_tree->Branch("ChiBVtxZ_fit",&ChiBVtxZ_fit,"ChiBVtxZ_fit/F");
		onia_tree->Branch("ChiBVtxP_fit",&ChiBVtxP_fit,"ChiBVtxP_fit/F");
		onia_tree->Branch("dz_photon_refit",&dz_photon_refit,"dz_photon_refit/F");
		onia_tree->Branch("dz_before_refit",&dz_before_refit,"dz_before_refit/F");
		onia_tree->Branch("mu1M_fit",&mu1M_fit,"mu1M_fit/F");
		onia_tree->Branch("mu1Px_fit",&mu1Px_fit,"mu1Px_fit/F");
		onia_tree->Branch("mu1Py_fit",&mu1Py_fit,"mu1Py_fit/F");
		onia_tree->Branch("mu1Pz_fit",&mu1Pz_fit,"mu1Pz_fit/F");
		onia_tree->Branch("mu2M_fit",&mu2M_fit,"mu2M_fit/F");
		onia_tree->Branch("mu2Px_fit",&mu2Px_fit,"mu2Px_fit/F");
		onia_tree->Branch("mu2Py_fit",&mu2Py_fit,"mu2Py_fit/F");
		onia_tree->Branch("mu2Pz_fit",&mu2Pz_fit,"mu2Pz_fit/F");
		onia_tree->Branch("gammaM_fit",&gammaM_fit,"gammaM_fit/F");
		onia_tree->Branch("gammaPx_fit",&gammaPx_fit,"gammaPx_fit/F");
		onia_tree->Branch("gammaPy_fit",&gammaPy_fit,"gammaPy_fit/F");
		onia_tree->Branch("gammaPz_fit",&gammaPz_fit,"gammaPz_fit/F");
		onia_tree->Branch("mu1_p4_fit",  "TLorentzVector", &mu1_p4_fit);
		onia_tree->Branch("mu2_p4_fit",  "TLorentzVector", &mu2_p4_fit);
		onia_tree->Branch("gamma_p4_fit",  "TLorentzVector", &gamma_p4_fit);
		onia_tree->Branch("gamma_scaledp4_fit",  "TLorentzVector", &gamma_scaledp4_fit);

		onia_tree->Branch("conversionMass",&conversionMass,"conversionMass/F");
		onia_tree->Branch("conversionPt",&conversionPt,"conversionPt/F");
		onia_tree->Branch("conversionEta",&conversionEta,"conversionEta/F");
		onia_tree->Branch("conversionVertexRho",&conversionVertexRho,"conversionVertexRho/F");
		onia_tree->Branch("conversionVertexX",&conversionVertexX,"conversionVertexX/F");
		onia_tree->Branch("conversionVertexY",&conversionVertexY,"conversionVertexY/F");
		onia_tree->Branch("conversionVertexZ",&conversionVertexZ,"conversionVertexZ/F");
		onia_tree->Branch("conversion_p4",  "TLorentzVector", &conversion_p4);


	}

	if (isMC_ || OnlyGen_) {
		onia_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
		onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
		onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
		onia_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
		onia_tree->Branch("gen_muonM_p4",  "TLorentzVector",  &gen_muonM_p4);
	}
	genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
	packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");


	TFileDirectory hists = fs->mkdir( "hists" );
	for (int m=0; m<15; m++) { 
		char str[100];
		sprintf(str,"refitY1sGammaMass_VeryLoose_Pt%02d",m);
		myRefitY1sConversionmassVeryLoose.push_back(hists.make<TH1F>(str,str,50000,0,500.0));
		myRefitY1sConversionmassVeryLoose[m]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}#gamma} [GeV/c^{2}]");
		sprintf(str,"refitY1sGammaMass_Loose_Pt%02d",m);
		myRefitY1sConversionmassLoose.push_back(hists.make<TH1F>(str,str,50000,0,500.0));
		myRefitY1sConversionmassLoose[m]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}#gamma} [GeV/c^{2}]");
		sprintf(str,"refitY1sGammaMass_Tight_Pt%02d",m);
		myRefitY1sConversionmassTight.push_back(hists.make<TH1F>(str,str,50000,0,500.0));
		myRefitY1sConversionmassTight[m]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}#gamma} [GeV/c^{2}]");
	}
	myConversionPt_all   = hists.make<TH1F>("myConversionPt_all","myConversionPt_all",5000,0,50.0);
	myConversionPt_all->GetXaxis()->SetTitle("Conversion p_{T} [GeV/c]");
	myConversionMass_all = hists.make<TH1F>("myConversionMass_all","myConversionMass_all",10000,0,1.0);;
	myConversionMass_all->GetXaxis()->SetTitle("Conversion mass [GeV/c^{2}]");
	myConversionRadius_all = hists.make<TH1F>("myConversionRadius_all","myConversionRadius_all",100,0,50.0);
	myConversionRadius_all->GetXaxis()->SetTitle("Conversion radius [cm]");
	myConversionVertex_all = hists.make<TH2F>("myConversionVertex_all","myConversionVertex_all",150,-15,15,150,-15,15);
	myConversionVertex_all->GetXaxis()->SetTitle("x [cm]");
	myConversionVertex_all->GetYaxis()->SetTitle("y [cm]");
	muGammaDeltaP = hists.make<TH1F>("muGammaDeltaP","muGammaDeltaP",2000,0,20);
	muGammaDeltaR = hists.make<TH1F>("muGammaDeltaR","muGammaDeltaR",2000,0,10);
	myDimuonMass_all = hists.make<TH1F>("myDimuonMass_all","myDimuonMass_all",2000,0,20.0);
	myDimuonMass_all->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
	myY1SFitMass_all = hists.make<TH1F>("myY1SFitMass_all","myY1SFitMass_all",2000,0,20.0);
	myY1SFitMass_all->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
	mydz_photon_refit_all = hists.make<TH1F>("mydz_photon_refit_all","mydz_photon_refit_all",4000,-10,10);
	mydz_photon_refit_all->GetXaxis()->SetTitle("dz[cm] after photon refit");
	mydz_before_refit_all = hists.make<TH1F>("mydz_before_refit_all","mydz_before_refit_all",4000,-10,10);
	mydz_before_refit_all->GetXaxis()->SetTitle("dz[cm] before photon refit");
	myChiBVtxP_all = hists.make<TH1F>("myChiBVtxP_all","myChiBVtxP_all",1000,0,1);
	myChiBVtxP_all->GetXaxis()->SetTitle("#mu^{+}#mu^{-}#gamma vertex fit #chi^{2} probability");

	myConversionPt_used   = hists.make<TH1F>("myConversionPt_used","myConversionPt_used",5000,0,50.0);
	myConversionPt_used->GetXaxis()->SetTitle("Conversion p_{T} [GeV/c]");
	myConversionMass_used = hists.make<TH1F>("myConversionMass_used","myConversionMass_used",10000,0,1.0);;
	myConversionMass_used->GetXaxis()->SetTitle("Conversion mass [GeV/c^{2}]");
	myConversionRadius_used = hists.make<TH1F>("myConversionRadius_used","myConversionRadius_used",100,0,50.0);
	myConversionRadius_used->GetXaxis()->SetTitle("Conversion radius [cm]");
	myConversionVertex_used = hists.make<TH2F>("myConversionVertex_used","myConversionVertex_used",150,-15,15,150,-15,15);
	myConversionVertex_used->GetXaxis()->SetTitle("x [cm]");
	myConversionVertex_used->GetYaxis()->SetTitle("y [cm]");
	myDimuonMass_used = hists.make<TH1F>("myDimuonMass_used","myDimuonMass_used",2000,0,20.0);
	myDimuonMass_used->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
	myY1SFitMass_used = hists.make<TH1F>("myY1SFitMass_used","myY1SFitMass_used",2000,0,20.0);
	myY1SFitMass_used->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
	mydz_photon_refit_used = hists.make<TH1F>("mydz_photon_refit_used","mydz_photon_refit_used",4000,-10,10);
	mydz_photon_refit_used->GetXaxis()->SetTitle("dz[cm] after photon refit");
	mydz_before_refit_used = hists.make<TH1F>("mydz_before_refit_used","mydz_before_refit_used",4000,-10,10);
	mydz_before_refit_used->GetXaxis()->SetTitle("dz[cm] before photon refit");
	myChiBVtxP_used = hists.make<TH1F>("myChiBVtxP_used","myChiBVtxP_used",1000,0,1);
	myChiBVtxP_used->GetXaxis()->SetTitle("#mu^{+}#mu^{-}#gamma vertex fit #chi^{2} probability");
}

MuMuGammaRootupler::~MuMuGammaRootupler() {}

//
// member functions
//

const reco::Candidate* MuMuGammaRootupler::GetAncestor(const reco::Candidate* p) {
	if (p->numberOfMothers()) {
		if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
		else return p->mother(0);
	}
	//std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
	return p;
}

//Check recursively if any ancestor of particle is the given one
bool MuMuGammaRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
	if (ancestor == particle ) return true;
	for (size_t i=0; i< particle->numberOfMothers(); i++) {
		if (isAncestor(ancestor, particle->mother(i))) return true;
	}
	return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an int between 0 and 127, in binary it is:
	(pass 2)(pass 1)(pass 0)
	ex. 7 = pass 0, 1 and 2
	ex. 6 = pass 2, 3
	ex. 1 = pass 0
	*/

UInt_t MuMuGammaRootupler::getTriggerBits(const edm::Event& iEvent ) {
	UInt_t itrigger = 0;
	edm::Handle<edm::TriggerResults> triggerResults_handle;
	iEvent.getByToken(triggerResults_Label, triggerResults_handle);
	if ( triggerResults_handle.isValid() ) { 
		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
		std::vector <unsigned int> bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7, bits_8, bits_9;
		for ( int version = 1; version<3; version ++ ) {
			std::stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9;
			ss0<<"HLT_Dimuon16_Jpsi_v"<<version;
			bits_0.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss0.str()).label().c_str()));
			ss1<<"HLT_Dimuon13_PsiPrime_v"<<version;
			bits_1.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss1.str()).label().c_str()));
			ss2<<"HLT_Dimuon13_Upsilon_v"<<version;
			bits_2.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss2.str()).label().c_str()));
			ss3<<"HLT_Dimuon10_Jpsi_Barrel_v"<<version;
			bits_3.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss3.str()).label().c_str()));
			ss4<<"HLT_Dimuon8_PsiPrime_Barrel_v"<<version;
			bits_4.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss4.str()).label().c_str()));
			ss5<<"HLT_Dimuon8_Upsilon_Barrel_v"<<version;
			bits_5.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss5.str()).label().c_str()));
			ss6<<"HLT_Dimuon20_Jpsi_v"<<version;
			bits_6.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss6.str()).label().c_str()));
			ss7<<"HLT_Dimuon0_Phi_Barrel_v"<<version;
			bits_7.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss7.str()).label().c_str()));
			ss8<<"HLT_DoubleMu4_3_Bs_v"<<version;
			bits_8.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss8.str()).label().c_str()));
			ss9<<"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"<<version;
			bits_9.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss9.str()).label().c_str()));
		}
		for (unsigned int i=0; i<bits_0.size(); i++) {
			unsigned int bit = bits_0[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 1;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_1.size(); i++) {
			unsigned int bit = bits_1[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 2;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_2.size(); i++) {
			unsigned int bit = bits_2[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 4;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_3.size(); i++) {
			unsigned int bit = bits_3[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 8;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_4.size(); i++) {
			unsigned int bit = bits_4[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 16;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_5.size(); i++) {
			unsigned int bit = bits_5[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 32;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_6.size(); i++) {
			unsigned int bit = bits_6[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 64;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_7.size(); i++) {
			unsigned int bit = bits_7[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 128;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_8.size(); i++) {
			unsigned int bit = bits_8[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 256;
					break;
				}   
			}   
		}   
		for (unsigned int i=0; i<bits_9.size(); i++) {
			unsigned int bit = bits_9[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 512;
					break;
				}   
			}   
		}   
	}
	return itrigger;
}

// ------------ method called for each event  ------------
void MuMuGammaRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

	edm::Handle<pat::CompositeCandidateCollection> dimuons;
	iEvent.getByToken(dimuon_Label,dimuons);

	edm::Handle<pat::CompositeCandidateCollection> conversions;
	iEvent.getByToken(conversion_Label,conversions);


	edm::Handle<reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

	/*edm::Handle<reco::ConversionCollection> conversionHandle;
	  iEvent.getByLabel(conversion_Label,conversionHandle);

	  edm::Handle<reco::PFCandidateCollection> pfcandidates;
	  iEvent.getByLabel("particleFlow",pfcandidates);
	  const reco::PFCandidateCollection pfphotons = selectPFPhotons(*pfcandidates);
	  */

	if (!OnlyGen_) {
		numPrimaryVertices = -1;
		if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();

		trigger = getTriggerBits(iEvent);
		run     = iEvent.id().run();
		lumi    = iEvent.id().luminosityBlock();
		event   = iEvent.id().event();
		//std::cout<<"trigger:"<<trigger<<std::endl;
	}

	dimuon_pdgId = 0;
	mother_pdgId = 0;
	irank = 0;

	// Pruned particles are the one containing "important" stuff
	edm::Handle<reco::GenParticleCollection> pruned;
	iEvent.getByToken(genCands_, pruned);

	// Packed particles are all the status 1, so usable to remake jets
	// The navigation from status 1 to pruned is possible (the other direction should be made by hand)
	edm::Handle<pat::PackedGenParticleCollection> packed;
	iEvent.getByToken(packCands_,  packed);

	if ((isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid()) {
		dimuon_pdgId  = 0;
		gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		int foundit   = 0;

		for (size_t i=0; i<pruned->size(); i++) {
			int p_id = abs((*pruned)[i].pdgId());
			const reco::Candidate *aonia = &(*pruned)[i];
			if (( p_id == pdgid_ ) && (aonia->status() == 2)) {
				dimuon_pdgId = p_id;
				foundit++;
				for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
					const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
					const reco::Candidate * d = &(*packed)[j];
					if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ){
						gen_muonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
						foundit++;
					} 
					if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
						gen_muonP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
						foundit++;
					}
					if ( foundit == 3 ) break;               
				}
				if ( foundit == 3 ) {
					gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;   // this should take into account FSR
					mother_pdgId  = GetAncestor(aonia)->pdgId();
					break;
				} else {
					foundit = 0;
					dimuon_pdgId = 0;
					mother_pdgId = 0;
					gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
				}            
			}  // if ( p_id
		} // for (size

		// sanity check
		//if ( ! dimuon_pdgId ) std::cout << "MuMuGammaRootupler: does not found the given decay " << run << "," << event << std::endl;
	}  // end if isMC

	float OniaMassMax_ = OniaMassCuts_[1];
	float OniaMassMin_ = OniaMassCuts_[0];

	// Kinematic fit

	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
	//edm::ESHandle<MagneticField> bFieldHandle;
	//iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	if ( ! OnlyGen_ && dimuons.isValid() && dimuons->size() > 0) {

		for(pat::CompositeCandidateCollection::const_iterator dimuonCand=dimuons->begin();dimuonCand!= dimuons->end(); ++dimuonCand)
		{
			if (dimuonCand->mass() < OniaMassMin_ || dimuonCand->mass() > OniaMassMax_) continue;
			if (dimuonCand->daughter("muon1")->charge() == dimuonCand->daughter("muon2")->charge() ) continue;

			//raw dimuon and muon
			dimuon_p4.SetPtEtaPhiM(dimuonCand->pt(), dimuonCand->eta(), dimuonCand->phi(), dimuonCand->mass());
			reco::Candidate::LorentzVector vP = dimuonCand->daughter("muon1")->p4();
			reco::Candidate::LorentzVector vM = dimuonCand->daughter("muon2")->p4();
			//std::cout<<"muon charge"<<dimuonCand->daughter("muon1")->charge()<<" "<<dimuonCand->daughter("muon2")->charge()<<std::endl;
			if ( dimuonCand->daughter("muon1")->charge() < 0) {
				vP = dimuonCand->daughter("muon2")->p4();
				vM = dimuonCand->daughter("muon1")->p4();
			}   
			muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
			muonM_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
			MassErr = dimuonCand->userFloat("MassErr");
			vProb = dimuonCand->userFloat("vProb");
			DCA = dimuonCand->userFloat("DCA");
			ppdlPV = dimuonCand->userFloat("ppdlPV");
			ppdlErrPV = dimuonCand->userFloat("ppdlErrPV");
			cosAlpha = dimuonCand->userFloat("cosAlpha");
			myDimuonMass_all->Fill(dimuon_p4.M());

			//dimuon refit. 
			//Here we use the KinematicParticleVertexFitter with muon mass. But in the Onia2MuMu skim, it was just KalmanVertexFitter. 
			//Fitted Vertex from both methods are the same (dimuonCand->userData<reco::Vertex>("commonVertex")->z() == mumu_vFit_vertex_noMC->position().z()), 
			//but fitted p4 and mass are different. 
			reco::TrackRef JpsiTk[2]={  //this is from Chib code
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1") ) )->innerTrack(),
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2") ) )->innerTrack()
			};  
			std::vector<reco::TransientTrack> MuMuTT;
			MuMuTT.push_back((*theB).build(&JpsiTk[0]));
			MuMuTT.push_back((*theB).build(&JpsiTk[1]));	
			//reco::TransientTrack muonPTT(JpsiTk[0], &(*bFieldHandle));
			//reco::TransientTrack muonMTT(JpsiTk[1], &(*bFieldHandle));

			KinematicParticleFactoryFromTransientTrack pmumuFactory;
			std::vector<RefCountedKinematicParticle> mumuParticles;
			const ParticleMass muonMass(0.1056583);
			float muonSigma = muonMass*1E-6;
			mumuParticles.push_back(pmumuFactory.particle(MuMuTT[0],muonMass,float(0),float(0),muonSigma));
			mumuParticles.push_back(pmumuFactory.particle(MuMuTT[1],muonMass,float(0),float(0),muonSigma));

			KinematicParticleVertexFitter mumufitter;
			RefCountedKinematicTree mumuVertexFitTree;
			try {mumuVertexFitTree = mumufitter.fit(mumuParticles);}
			catch (...) {
				std::cout<<"mumu fit: PerigeeKinematicState::kinematic state passed is not valid!"<<std::endl;
				continue;
			}   
			//mumuVertexFitTree = mumufitter.fit(mumuParticles);
			RefCountedKinematicParticle mumu_vFit_noMC;
			RefCountedKinematicVertex mumu_vFit_vertex_noMC;

			if (mumuVertexFitTree->isValid()) {
				mumuVertexFitTree->movePointerToTheTop();     
				mumu_vFit_noMC = mumuVertexFitTree->currentParticle();    
				mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
				//KinematicParameters mymumupara=  mumu_vFit_noMC->currentState().kinematicParameters();
				//cout<<"mymumupara px="<<mymumupara.momentum().x()<<",py="<<mymumupara.momentum().y()<<", m="<<mumu_vFit_noMC->currentState().mass()<<endl;      
				//float mymumuonlyctau=GetcTau(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);
				//float mymumuonlyctauerr=GetcTauErr(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);     
				//cout<<"mymumuonlyctau="<<mymumuonlyctau<<endl;
				//JmumuCtau->push_back( mymumuonlyctau );
				//JmumuCtauerr->push_back( mymumuonlyctauerr );
				JmumuMass = mumu_vFit_noMC->currentState().mass();
				JmumuMassErr = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) )  ;
				JmumuPx = mumu_vFit_noMC->currentState().globalMomentum().x();
				JmumuPy = mumu_vFit_noMC->currentState().globalMomentum().y();
				JmumuPz = mumu_vFit_noMC->currentState().globalMomentum().z();
				JmumuVtxCL = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) ;
				JmumuVtxCL2 = mumu_vFit_vertex_noMC->chiSquared() ;
				JmumuDecayVtxX = mumu_vFit_vertex_noMC->position().x() ;
				JmumuDecayVtxY = mumu_vFit_vertex_noMC->position().y() ;
				JmumuDecayVtxZ = mumu_vFit_vertex_noMC->position().z() ;
				JmumuDecayVtxXE = mumu_vFit_vertex_noMC->error().cxx() ;
				JmumuDecayVtxYE = mumu_vFit_vertex_noMC->error().cyy() ;
				JmumuDecayVtxZE = mumu_vFit_vertex_noMC->error().czz() ;
				Jmumu_p4.SetXYZM( JmumuPx, JmumuPy, JmumuPz, JmumuMass ); 
				//JmumumupIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon1));
				//JmumumumIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon2));
				//std::cout<<"mass0="<<dimuonCand->mass()<<", mass1="<<JmumuMass<<std::endl;
				//std::cout<<"vProb0="<<vProb<<", vProb1="<<JmumuVtxCL<<std::endl;
				myY1SFitMass_all->Fill(Jmumu_p4.M());


				//mumugamma refit
				if ( ! OnlyGen_ && conversions.isValid() && conversions->size() > 0) {
					for(pat::CompositeCandidateCollection::const_iterator conv=conversions->begin(); conv!= conversions->end(); ++conv){

						conversionMass = conv->mass();
						conversionPt = conv->pt() ;
						conversionEta = conv->eta() ;
						conversionVertexRho = conv->vertex().rho() ;
						conversionVertexX = conv->vertex().x() ;
						conversionVertexY = conv->vertex().y() ;
						conversionVertexZ = conv->vertex().z() ;
						conversion_p4.SetPtEtaPhiM(conv->pt(), conv->eta(), conv->phi(), conv->mass());

						float deltaP_P; float deltaR_P; float deltaP_M; float deltaR_M;
						if ( conv->userData<reco::Track>("track0")->charge() > 0) {
							//std::cout<<"ele charge"<<conv->userData<reco::Track>("track0")->charge()<<" "<<conv->userData<reco::Track>("track1")->charge()<<std::endl; 
							deltaP_P = sqrt(pow(conv->userData<reco::Track>("track0")->px()-muonP_p4.Px(),2) + pow(conv->userData<reco::Track>("track0")->py()-muonP_p4.Py(),2) + pow(conv->userData<reco::Track>("track0")->pz()-muonP_p4.Pz(),2));
							deltaR_P = sqrt(pow(conv->userData<reco::Track>("track0")->eta()-muonP_p4.Eta(),2) + pow(conv->userData<reco::Track>("track0")->phi()-muonP_p4.Phi(),2));
							deltaP_M = sqrt(pow(conv->userData<reco::Track>("track1")->px()-muonM_p4.Px(),2) + pow(conv->userData<reco::Track>("track1")->py()-muonM_p4.Py(),2) + pow(conv->userData<reco::Track>("track1")->pz()-muonM_p4.Pz(),2));
							deltaR_M = sqrt(pow(conv->userData<reco::Track>("track1")->eta()-muonM_p4.Eta(),2) + pow(conv->userData<reco::Track>("track1")->phi()-muonM_p4.Phi(),2));
						}
						else {
							deltaP_P = sqrt(pow(conv->userData<reco::Track>("track1")->px()-muonP_p4.Px(),2) + pow(conv->userData<reco::Track>("track1")->py()-muonP_p4.Py(),2) + pow(conv->userData<reco::Track>("track1")->pz()-muonP_p4.Pz(),2));
							deltaR_P = sqrt(pow(conv->userData<reco::Track>("track1")->eta()-muonP_p4.Eta(),2) + pow(conv->userData<reco::Track>("track1")->phi()-muonP_p4.Phi(),2));
							deltaP_M = sqrt(pow(conv->userData<reco::Track>("track0")->px()-muonM_p4.Px(),2) + pow(conv->userData<reco::Track>("track0")->py()-muonM_p4.Py(),2) + pow(conv->userData<reco::Track>("track0")->pz()-muonM_p4.Pz(),2));
							deltaR_M = sqrt(pow(conv->userData<reco::Track>("track0")->eta()-muonM_p4.Eta(),2) + pow(conv->userData<reco::Track>("track0")->phi()-muonM_p4.Phi(),2));
						}
						deltaP = ((deltaP_P<deltaP_M)?deltaP_P:deltaP_M);
						deltaR = ((deltaR_P<deltaR_M)?deltaR_P:deltaR_M);

						if (mumuVertexFitTree->isValid()) 
							dz_before_refit = (conv->vertex().z()-mumu_vFit_vertex_noMC->position().z()) -
								((conv->vertex().x()-mumu_vFit_vertex_noMC->position().x())*conv->p4().x()+
								 (conv->vertex().y()-mumu_vFit_vertex_noMC->position().y())*conv->p4().y())/conv->p4().rho() *
								conv->p4().z()/conv->p4().rho();
						else dz_before_refit = -999;

						myConversionPt_all->Fill(conversionPt);
						myConversionMass_all->Fill(conversionMass);
						myConversionRadius_all->Fill(conversionVertexRho);
						myConversionVertex_all->Fill(conversionVertexX,conversionVertexY);
						muGammaDeltaP->Fill(deltaP);
						muGammaDeltaR->Fill(deltaR);
						mydz_before_refit_all->Fill(dz_before_refit);

						//dielectron fit
						std::vector<reco::TransientTrack> EETT;
						EETT.push_back((*theB).build(conv->userData<reco::Track>("track0")));
						EETT.push_back((*theB).build(conv->userData<reco::Track>("track1")));

						const ParticleMass zero_mass(0);
						float zero_sigma = 1E-6;

						const ParticleMass eleMass(0.000511);
						float eleSigma = 1E-6;
						KinematicParticleFactoryFromTransientTrack pFactory;
						std::vector<RefCountedKinematicParticle> PhotonParticles;
						PhotonParticles.push_back(pFactory.particle(EETT[0],eleMass,float(0),float(0),eleSigma));
						PhotonParticles.push_back(pFactory.particle(EETT[1],eleMass,float(0),float(0),eleSigma));

						//KinematicParticleVertexFitter fitter(ps);
						KinematicParticleVertexFitter fitter;
						RefCountedKinematicTree photonVertexFitTree;
						try {photonVertexFitTree = fitter.fit(PhotonParticles);}
						catch (...) {
							std::cout<<"Photon fit: PerigeeKinematicState::kinematic state passed is not valid!"<<std::endl;
							continue;
						}
						//photonVertexFitTree = fitter.fit(PhotonParticles);

						if (!photonVertexFitTree->isValid())
						{
							edm::ParameterSet pSet;

							pSet.addParameter<double>("maxDistance", 3);
							pSet.addParameter<int>("maxNbrOfIterations", 10000); //10

							KinematicParticleVertexFitter fitter2(pSet);

							//RefCountedKinematicTree photonVertexFitTree;
							photonVertexFitTree = fitter2.fit(PhotonParticles);
						}
						if (photonVertexFitTree->isValid())
						{
							// now apply Photon mass constraint                                                                                                            
							KinematicParticleFitter csFitterPhoton;
							KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass,zero_sigma);
							// add mass constraint to the photon fit to do a constrained fit:                                                                              

							photonVertexFitTree->movePointerToTheTop();
							photonVertexFitTree = csFitterPhoton.fit(pho_c,photonVertexFitTree);

							if (photonVertexFitTree->isValid())
							{
								photonVertexFitTree->movePointerToTheTop();
								RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
								// compute dz. Note: mumu_vFit_vertex_noMC->position().z() == dimuonVtxZ
								RefCountedKinematicVertex fitPhotonDecayVertex = photonVertexFitTree->currentDecayVertex();
								double MomRho_fit = sqrt(fitPhoton->currentState().globalMomentum().x()*fitPhoton->currentState().globalMomentum().x()+
										fitPhoton->currentState().globalMomentum().y()*fitPhoton->currentState().globalMomentum().y());
								if (mumuVertexFitTree->isValid()) {
									dz_photon_refit = (fitPhotonDecayVertex->position().z()-mumu_vFit_vertex_noMC->position().z()) -
										((fitPhotonDecayVertex->position().x()-mumu_vFit_vertex_noMC->position().x())*fitPhoton->currentState().globalMomentum().x()+
										 (fitPhotonDecayVertex->position().y()-mumu_vFit_vertex_noMC->position().y())*fitPhoton->currentState().globalMomentum().y())/MomRho_fit *
										fitPhoton->currentState().globalMomentum().z()/MomRho_fit;
									//std::cout<<"photonZfit="<<fitPhotonDecayVertex->position().z()<<", dimuonZfit="<<mumu_vFit_vertex_noMC->position().z()<<", dzFit="<<dz_fit<<", photonMomRhoFit="<<MomRho_fit<<std::endl;
									//std::cout<<"photonZ="<<conv->vertex().z()<<", dimuonZfit="<<mumu_vFit_vertex_noMC->position().z()<<", dz="<<dz<<", photonMomRho="<<conv->p4().rho()<<std::endl;
								}
								else dz_photon_refit = -999;
								mydz_photon_refit_all->Fill(dz_photon_refit);

								std::vector<RefCountedKinematicParticle> allChiBDaughters;
								allChiBDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
								allChiBDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
								allChiBDaughters.push_back(fitPhoton);

								KinematicConstrainedVertexFitter constVertexFitter;
								MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
								RefCountedKinematicTree ChiBTree = constVertexFitter.fit(allChiBDaughters,upsilon_mtc);

								//fit w/o any mass constraint
								//RefCountedKinematicTree ChiBTree = constVertexFitter.fit(allChiBDaughters);


								if(!ChiBTree->isEmpty())
								{
									ChiBTree->movePointerToTheTop();
									RefCountedKinematicParticle fitChiB = ChiBTree->currentParticle();
									RefCountedKinematicVertex ChiBDecayVertex = ChiBTree->currentDecayVertex();

									if (fitChiB->currentState().isValid())
									{ //Get chib         
										ChiBM_fit = fitChiB->currentState().mass();
										ChiBM_fitErr = sqrt(fitChiB->currentState().kinematicParametersError().matrix()(6,6));
										ChiBPx_fit = fitChiB->currentState().kinematicParameters().momentum().x();
										ChiBPy_fit = fitChiB->currentState().kinematicParameters().momentum().y();
										ChiBPz_fit = fitChiB->currentState().kinematicParameters().momentum().z();
										ChiBVtxX_fit = ChiBDecayVertex->position().x();
										ChiBVtxY_fit = ChiBDecayVertex->position().y();
										ChiBVtxZ_fit = ChiBDecayVertex->position().z();
										ChiBVtxP_fit=ChiSquaredProbability((double)(ChiBDecayVertex->chiSquared()),(double)(ChiBDecayVertex->degreesOfFreedom()));
										myChiBVtxP_all->Fill(ChiBVtxP_fit);
										//reco::CompositeCandidate recoChib(0,math::XYZTLorentzVector(ChiBPx_fit,ChiBPy_fit,ChiBPz_fit,sqrt(ChiBM_fit*ChiBM_fit+ChiBPx_fit*ChiBPx_fit+ChiBPy_fit*ChiBPy_fit+ChiBPz_fit*ChiBPz_fit)),math::XYZPoint(ChiBVtxX_fit,ChiBVtxY_fit,ChiBVtxZ_fit),50551);

										//pat::CompositeCandidate patChib(recoChib);
										//patChib.addUserFloat("vProb",ChiBVtxP_fit);

										//get first muon
										bool child = ChiBTree->movePointerToTheFirstChild();
										RefCountedKinematicParticle fitMu1 = ChiBTree->currentParticle();
										if(!child) break;

										mu1M_fit = fitMu1->currentState().mass();
										mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
										mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
										mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
										mu1_p4_fit.SetXYZM( mu1Px_fit, mu1Py_fit, mu1Pz_fit, mu1M_fit );
										//reco::CompositeCandidate recoMu1(0,math::XYZTLorentzVector(mu1Px_fit,mu1Py_fit,mu1Pz_fit,sqrt(mu1M_fit*mu1M_fit+mu1Px_fit*mu1Px_fit+mu1Py_fit*mu1Py_fit+mu1Pz_fit*mu1Pz_fit)),math::XYZPoint(ChiBVtxX_fit,ChiBVtxY_fit,ChiBVtxZ_fit),13);
										//pat::CompositeCandidate patMu1(recoMu1);


										//get second muon
										child = ChiBTree->movePointerToTheNextChild();
										RefCountedKinematicParticle fitMu2 = ChiBTree->currentParticle();
										if(!child) break;

										mu2M_fit = fitMu2->currentState().mass();
										mu2Px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
										mu2Py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
										mu2Pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
										mu2_p4_fit.SetXYZM( mu2Px_fit, mu2Py_fit, mu2Pz_fit, mu2M_fit );
										//reco::CompositeCandidate recoMu2(0,math::XYZTLorentzVector(mu2Px_fit,mu2Py_fit,mu2Pz_fit,sqrt(mu2M_fit*mu2M_fit+mu2Px_fit*mu2Px_fit+mu2Py_fit*mu2Py_fit+mu2Pz_fit*mu2Pz_fit)),math::XYZPoint(ChiBVtxX_fit,ChiBVtxY_fit,ChiBVtxZ_fit),13);
										//pat::CompositeCandidate patMu2(recoMu2);

										//Define Upsilon from two muons
										//pat::CompositeCandidate ups;
										//ups.addDaughter(patMu1, "muon1");
										//ups.addDaughter(patMu2,"muon2");
										//ups.setP4(patMu1.p4()+patMu2.p4());


										//get photon
										child = ChiBTree->movePointerToTheNextChild();
										RefCountedKinematicParticle fitGamma = ChiBTree->currentParticle();
										if(!child) break;

										gammaM_fit = fitGamma->currentState().mass();
										gammaPx_fit = fitGamma->currentState().kinematicParameters().momentum().x();
										gammaPy_fit = fitGamma->currentState().kinematicParameters().momentum().y();
										gammaPz_fit = fitGamma->currentState().kinematicParameters().momentum().z();
										gamma_p4_fit.SetXYZM( gammaPx_fit, gammaPy_fit, gammaPz_fit, gammaM_fit );
										gamma_scaledp4_fit.SetPtEtaPhiM( 1.04*gamma_p4_fit.Pt(), gamma_p4_fit.Eta(), gamma_p4_fit.Phi(),  gamma_p4_fit.M() );
										//reco::CompositeCandidate recoGamma(0,math::XYZTLorentzVector(gammaPx_fit,gammaPy_fit,gammaPz_fit,sqrt(gammaM_fit*gammaM_fit+gammaPx_fit*gammaPx_fit+gammaPy_fit*gammaPy_fit+gammaPz_fit*gammaPz_fit)),math::XYZPoint(ChiBVtxX_fit,ChiBVtxY_fit,ChiBVtxZ_fit),13);
										//pat::CompositeCandidate patGamma(recoGamma);

										//patChib.addDaughter(ups,"dimuon");
										//patChib.addDaughter(patGamma,"photon");

										//chicCompCandRefitColl->push_back(patChib);

										conversionMass = conv->mass();
										conversionPt = conv->pt() ;
										conversionEta = conv->eta() ;
										conversionVertexRho = conv->vertex().rho() ;
										conversionVertexX = conv->vertex().x() ;
										conversionVertexY = conv->vertex().y() ;
										conversionVertexZ = conv->vertex().z() ;
										conversion_p4.SetPtEtaPhiM(conv->pt(), conv->eta(), conv->phi(), conv->mass());

										//mass with photon energy scaled by X3872
										mymumugammamassY1s = ((mu1_p4_fit + mu2_p4_fit + gamma_scaledp4_fit).M() - (mu1_p4_fit + mu2_p4_fit).M() + upsilon_mass_); // (mu1_p4_fit + mu2_p4_fit).M() == upsilon_mass_
										//without scale
										//mymumugammamassY1s= ChiBM_fit - (mu1_p4_fit + mu2_p4_fit).M() + upsilon_mass_;  // (mu1_p4_fit + mu2_p4_fit).M() == upsilon_mass_

										//std::cout<<"trigger:"<<trigger<<"->"<<(trigger&36)<<", deltaP="<<deltaP<<", deltaR="<<deltaR<<", flags:"<<(conv->userInt("flags")&16)<<"->"<<(conv->userInt("flags")&16)<<std::endl;
										//apply final cuts and fill histogram
										if (1
												//apply trigger cut: 36 for upsilon, 73 for Jpsi
												&& (trigger&triggerCuts_)!=0

												//apply single muon cuts
												//&& myNumGoodLooseMuon>=2 //not needed anymore. applied in PAT
												&& fabs(muonP_p4.Eta())<=2.4 && fabs(muonM_p4.Eta())<=2.4  //
												&& muonP_p4.Pt()>=3.5 && muonM_p4.Pt()>=3.5 

												//apply dimuon mass cuts
												&& JmumuVtxCL>0.005
												&& fabs(JmumuMass-upsilon_mass_)<3*1.105*(JmumuMassErr) //another method: using dimuon instead of Jmumu
												&& dimuon_p4.Pt() >= 7 //another method: using Jmumu instead of dimuon

												//apply photon cuts //not needed anymore. applied in PAT
												//&& !(*IsPi0Bkg_photon_refit)[iRefit] 
												//&& (*conversionVertexChi2Prob)[rawPhotonIdx]>=0.00001   //default 0.0005, or 0.00001
												//&& (*conversionDaughter1qd0)[rawPhotonIdx]!=0
												//&& (*conversionDaughter2qd0)[rawPhotonIdx]!=0
												//&& (*conversionDaughter1NormalizedChi2)[rawPhotonIdx]<=10
												//&& (*conversionDaughter2NormalizedChi2)[rawPhotonIdx]<=10
												//&& ((*conversionDaughter1NHits)[rawPhotonIdx]+(*conversionDaughter2NHits)[rawPhotonIdx])>=7  //track hits cut at 4,3
												//&& fabs((*conversionDeltaPhi)[rawPhotonIdx])<0.2
												//apply pi_zero veto: 
											&& (conv->userInt("flags")&16) == 0 // for large window
												//&& (conv->userInt("flags")&8) == 0 // for small window

												//apply muon-ele iso cuts
												&& deltaP>0.01
												//&& deltaR>0.4

												//apply Y(1S)gamma cuts
												&& ChiBM_fitErr<10
												&& ChiBVtxP_fit>0.001
												&& fabs(dz_photon_refit)<0.2
												)
												{
													for (int pt=0; pt<15; pt++) {
														if (conversion_p4.Pt()>pt) myRefitY1sConversionmassVeryLoose[pt]->Fill(mymumugammamassY1s);
													}

													if (ChiBVtxP_fit>0.001 && fabs(dz_photon_refit)<0.1 ) {  
														for (int pt=0; pt<15; pt++) {
															if (conversion_p4.Pt()>pt) myRefitY1sConversionmassLoose[pt]->Fill(mymumugammamassY1s);
														}
														myConversionPt_used->Fill(conversionPt);
														myConversionMass_used->Fill(conversionMass);
														myConversionRadius_used->Fill(conversionVertexRho);
														myConversionVertex_used->Fill(conversionVertexX,conversionVertexY);;
														myDimuonMass_used->Fill(dimuon_p4.M());
														myY1SFitMass_used->Fill(Jmumu_p4.M());
														mydz_photon_refit_used->Fill(dz_photon_refit);
														mydz_before_refit_used->Fill(dz_before_refit);
														myChiBVtxP_used->Fill(ChiBVtxP_fit);


														irank++;
														onia_tree->Fill();
														//myoutfilesignal<<mymumugammamassY1s<<" "<<conversion_p4.Pt()<<" "<<ChiBM_fitErr<<" "<<muonP_p4.Pt()<<" "<<muonM_p4.Pt()<<" "<<deltaR<<" "<<dimuon_p4.Pt()<<endl;
													}

													if (ChiBVtxP_fit>0.02 && fabs(dz_photon_refit)<0.1 ) {  
														for (int pt=0; pt<15; pt++) {
															if (conversion_p4.Pt()>pt) myRefitY1sConversionmassTight[pt]->Fill(mymumugammamassY1s);
														}  
													}
												}
									}
								}
							}
						}
					}
				}
			}
			if (OnlyBest_) break; //only use the best dimuon pair
		}
	}  // End kinematic fit
	else {
		if (dimuon_pdgId && OnlyGen_) onia_tree->Fill();
		//else std::cout << "MuMuGammaRootupler: does not find a valid dimuon combination " << run << "," << event << std::endl;
	}
	/*
		dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		dimuon_pdgId = 0;
		mother_pdgId = 0;
		vProb=-1;
		*/
}

// ------------ method called once each job just before starting event loop  ------------
void MuMuGammaRootupler::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void MuMuGammaRootupler::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void MuMuGammaRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void MuMuGammaRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void MuMuGammaRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void MuMuGammaRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuMuGammaRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuGammaRootupler);
