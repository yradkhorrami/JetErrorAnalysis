#include "JetErrorAnalysis.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TPad.h"
#include "TPaveStats.h"


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;
using namespace std ;

JetErrorAnalysis aJetErrorAnalysis ;

JetErrorAnalysis::JetErrorAnalysis() : Processor("JetErrorAnalysis"),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0),
m_nTrueJets(0),
m_nTrueIsoLeps(0),
m_nRecoJets(0),
m_nIsoLeptons(0),
m_pTFile(NULL),
m_pTTree(NULL)
{

	// modify processor description
	_description = "JetErrorAnalysis does whatever it does ..." ;


	// register steering parameters: name, description, class-variable, default value


	// Inputs: MC-particles, Reco-particles, the link between the two

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"isolatedleptonCollection" ,
					"Name of the Isolated Lepton collection"  ,
					m_inputIsolatedleptonCollection ,
					std::string("ISOLeptons")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"RecoJetCollection" ,
					"Name of the input Reconstructed Jet collection"  ,
					m_recoJetCollectionName ,
					std::string("Durham_nJets")
				);

	registerProcessorParameter(	"jetMatchingMethod",
					"method for matching TRUE-RECO jets: [1] finding leading particle of reco jet in true jets, [2] direction of jets",
					m_jetMatchingMethod,
					int(1)
				);

	registerProcessorParameter(	"outputFilename",
					"name of output file",
					m_outputFile,
					std::string("")
				);

	registerProcessorParameter(	"HistogramsName",
					"name of histograms",
					m_histName,
					std::string("Std. Reco")
				);

	registerProcessorParameter(	"HistogramsColour",
					"color of histograms",
					m_histColour,
					int(1)
				);

	registerProcessorParameter(	"fillRootTree",
					"Whether create root file for new pfo collection or not true:create / false: do not create",
					m_fillRootTree,
					bool(true)
				);

	registerProcessorParameter(	"normalizeHistograms",
					"Whether normalize histograms to 1 or not true:normalize / false: do not normalize",
					m_normalizeHistograms,
					bool(true)
				);

	// Inputs: True jets (as a recoparticle, will be the sum of the _reconstructed particles_
	// created by the true particles in each true jet, in the RecoMCTruthLink sense.
	// link jet-to-reco particles, link jet-to-MC-particles.

	registerInputCollection( 	LCIO::LCRELATION,
					"MCTruthRecoLink",
					"Name of the MCTruthRecoLink input collection"  ,
					_MCTruthRecoLink,
					std::string("MCTruthRecoLink")
				);

	registerInputCollection( 	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					_MCParticleColllectionName ,
					std::string("MCParticlesSkimmed")
					);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"RecoParticleCollection" ,
					"Name of the ReconstructedParticles input collection"  ,
					_recoParticleCollectionName ,
					std::string("PandoraPFOs")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"RecoMCTruthLink",
					"Name of the RecoMCTruthLink input collection"  ,
					_recoMCTruthLink,
					std::string("RecoMCTruthLink")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"TrueJets" ,
					"Name of the TrueJetCollection input collection",
					_trueJetCollectionName ,
					std::string("TrueJets")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"FinalColourNeutrals" ,
					"Name of the FinalColourNeutralCollection input collection"  ,
					_finalColourNeutralCollectionName ,
					std::string("FinalColourNeutrals")
				);

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"InitialColourNeutrals" ,
					"Name of the InitialColourNeutralCollection input collection"  ,
					_initialColourNeutralCollectionName ,
					std::string("InitialColourNeutrals")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"TrueJetPFOLink" ,
					"Name of the TrueJetPFOLink input collection"  ,
					_trueJetPFOLink,
					std::string("TrueJetPFOLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"TrueJetMCParticleLink" ,
					"Name of the TrueJetMCParticleLink input collection"  ,
					_trueJetMCParticleLink,
					std::string("TrueJetMCParticleLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"FinalElementonLink rueJetMCParticleLink" ,
					"Name of the  FinalElementonLink input collection"	,
					_finalElementonLink,
					std::string("FinalElementonLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"InitialElementonLink" ,
					"Name of the  InitialElementonLink input collection"  ,
					_initialElementonLink,
					std::string("InitialElementonLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"FinalColourNeutralLink" ,
					"Name of the  FinalColourNeutralLink input collection"  ,
					_finalColourNeutralLink,
					std::string("FinalColourNeutralLink")
				);

	registerInputCollection( 	LCIO::LCRELATION,
					"InitialColourNeutralLink" ,
					"Name of the  InitialColourNeutralLink input collection"  ,
					_initialColourNeutralLink,
					std::string("InitialColourNeutralLink")
				);


}


void JetErrorAnalysis::init()
{

	streamlog_out(DEBUG6) << "   init called  " << std::endl ;

	printParameters() ;

	m_nRun = 0 ;
	m_nEvt = 0 ;
	if ( m_fillRootTree )
	{
		m_pTFile = new TFile(m_outputFile.c_str(),"recreate");
		m_pTTree = new TTree("eventTree","eventTree");
		m_pTTree->SetDirectory(m_pTFile);
		m_pTTree->Branch("run", &m_nRun, "run/I");
		m_pTTree->Branch("event", &m_nEvt, "event/I");
		m_pTTree->Branch("nTrueJets",&m_nTrueJets,"nTrueJets/I") ;
		m_pTTree->Branch("nTrueIsoLeps",&m_nTrueIsoLeps,"nTrueIsoLeps/I") ;
		m_pTTree->Branch("nRecoJets",&m_nRecoJets,"nRecoJets/I") ;
		m_pTTree->Branch("nIsoLeptons",&m_nIsoLeptons,"nIsoLeptons/I") ;
		m_pTTree->Branch("nSeenISRs",&m_nSeenISRs) ;
		m_pTTree->Branch( "nSLDecayBHadron" , &m_nSLDecayBHadron , "nSLDecayBHadron/I" );
		m_pTTree->Branch( "nSLDecayCHadron" , &m_nSLDecayCHadron , "nSLDecayCHadron/I" );
		m_pTTree->Branch( "nSLDecayTauLepton" , &m_nSLDecayTauLepton , "nSLDecayTauLepton/I" );
		m_pTTree->Branch( "nSLDecayTotal" , &m_nSLDecayTotal , "nSLDecayTotal/I" );
		m_pTTree->Branch( "nCorrectedSLD" , &m_nCorrectedSLD , "nCorrectedSLD/I" );
		m_pTTree->Branch("quarkPx",&m_quarkPx) ;
		m_pTTree->Branch("quarkPy",&m_quarkPy) ;
		m_pTTree->Branch("quarkPz",&m_quarkPz) ;
		m_pTTree->Branch("quarkE",&m_quarkE) ;
		m_pTTree->Branch("quarkTheta",&m_quarkTheta) ;
		m_pTTree->Branch("quarkPhi",&m_quarkPhi) ;
		m_pTTree->Branch("trueJetPx",&m_trueJetPx) ;
		m_pTTree->Branch("trueJetPy",&m_trueJetPy) ;
		m_pTTree->Branch("trueJetPz",&m_trueJetPz) ;
		m_pTTree->Branch("trueJetE",&m_trueJetE) ;
		m_pTTree->Branch("trueJetTheta",&m_trueJetTheta) ;
		m_pTTree->Branch("trueJetPhi",&m_trueJetPhi) ;
		m_pTTree->Branch("trueSeenJetPx",&m_trueSeenJetPx) ;
		m_pTTree->Branch("trueSeenJetPy",&m_trueSeenJetPy) ;
		m_pTTree->Branch("trueSeenJetPz",&m_trueSeenJetPz) ;
		m_pTTree->Branch("trueSeenJetE",&m_trueSeenJetE) ;
		m_pTTree->Branch("trueSeenJetTheta",&m_trueSeenJetTheta) ;
		m_pTTree->Branch("trueSeenJetPhi",&m_trueSeenJetPhi) ;
		m_pTTree->Branch("seenJetPx",&m_seenJetPx) ;
		m_pTTree->Branch("seenJetPy",&m_seenJetPy) ;
		m_pTTree->Branch("seenJetPz",&m_seenJetPz) ;
		m_pTTree->Branch("seenJetE",&m_seenJetE) ;
		m_pTTree->Branch("seenJetTheta",&m_seenJetTheta) ;
		m_pTTree->Branch("seenJetPhi",&m_seenJetPhi) ;
		m_pTTree->Branch("recoJetPx",&m_recoJetPx) ;
		m_pTTree->Branch("recoJetPy",&m_recoJetPy) ;
		m_pTTree->Branch("recoJetPz",&m_recoJetPz) ;
		m_pTTree->Branch("recoJetE",&m_recoJetE) ;
		m_pTTree->Branch("recoJetTheta",&m_recoJetTheta) ;
		m_pTTree->Branch("recoJetPhi",&m_recoJetPhi) ;
		m_pTTree->Branch("jetSigmaPx2",&m_jetSigmaPx2) ;
		m_pTTree->Branch("jetSigmaPxPy",&m_jetSigmaPxPy) ;
		m_pTTree->Branch("jetSigmaPy2",&m_jetSigmaPy2) ;
		m_pTTree->Branch("jetSigmaPxPz",&m_jetSigmaPxPz) ;
		m_pTTree->Branch("jetSigmaPyPz",&m_jetSigmaPyPz) ;
		m_pTTree->Branch("jetSigmaPz2",&m_jetSigmaPz2) ;
		m_pTTree->Branch("jetSigmaPxE",&m_jetSigmaPxE) ;
		m_pTTree->Branch("jetSigmaPyE",&m_jetSigmaPyE) ;
		m_pTTree->Branch("jetSigmaPzE",&m_jetSigmaPzE) ;
		m_pTTree->Branch("jetSigmaE2",&m_jetSigmaE2) ;
		m_pTTree->Branch("jetSigmaTheta2",&m_jetSigmaTheta2) ;
		m_pTTree->Branch("jetSigmaPhi2",&m_jetSigmaPhi2) ;
		m_pTTree->Branch("ResidualJetPx",&m_ResidualJetPx) ;
		m_pTTree->Branch("ResidualJetPy",&m_ResidualJetPy) ;
		m_pTTree->Branch("ResidualJetPz",&m_ResidualJetPz) ;
		m_pTTree->Branch("ResidualJetE",&m_ResidualJetE) ;
		m_pTTree->Branch("ResidualJetTheta",&m_ResidualJetTheta) ;
		m_pTTree->Branch("ResidualJetPhi",&m_ResidualJetPhi) ;
		m_pTTree->Branch("NormalizedResidualJetPx",&m_NormalizedResidualJetPx) ;
		m_pTTree->Branch("NormalizedResidualJetPy",&m_NormalizedResidualJetPy) ;
		m_pTTree->Branch("NormalizedResidualJetPz",&m_NormalizedResidualJetPz) ;
		m_pTTree->Branch("NormalizedResidualJetE",&m_NormalizedResidualJetE) ;
		m_pTTree->Branch("NormalizedResidualJetTheta",&m_NormalizedResidualJetTheta) ;
		m_pTTree->Branch("NormalizedResidualJetPhi",&m_NormalizedResidualJetPhi) ;
		m_pTTree->Branch("trueJetType",&m_trueJetType) ;
		m_pTTree->Branch("trueJetFlavour",&m_trueJetFlavour) ;
		m_pTTree->Branch("recoJetFlavour",&m_recoJetFlavour) ;
		h_ResidualJetPx = new TH1F( "h_ResidualJetPx" , "; ^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetPy = new TH1F( "h_ResidualJetPy" , "; ^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetPz = new TH1F( "h_ResidualJetPz" , "; ^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetE = new TH1F( "h_ResidualJetE" , "; ^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetTheta = new TH1F( "h_ResidualJetTheta" , "; ^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE} [rad]; Normalized #jet / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_ResidualJetPhi = new TH1F( "h_ResidualJetPhi" , "; ^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE} [rad]; Normalized #jet / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_NormalizedResidualJetPx = new TH1F( "h_NormalizedResidualJetPx" , "; (^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE}) / #sigma_{p_{x,jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetPy = new TH1F( "h_NormalizedResidualJetPy" , "; (^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE}) / #sigma_{p_{y,jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetPz = new TH1F( "h_NormalizedResidualJetPz" , "; (^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE}) / #sigma_{p_{z,jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetE = new TH1F( "h_NormalizedResidualJetE" , "; (^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE}) / #sigma_{E_{jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetTheta = new TH1F( "h_NormalizedResidualJetTheta" , "; (^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE}) / #sigma_{#theta_{jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetPhi = new TH1F( "h_NormalizedResidualJetPhi" , "; (^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE}) / #sigma_{#phi_{jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_ResidualJetPxSeen = new TH1F( "h_ResidualJetPxSeen" , "; ^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetPySeen = new TH1F( "h_ResidualJetPySeen" , "; ^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetPzSeen = new TH1F( "h_ResidualJetPzSeen" , "; ^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetESeen = new TH1F( "h_ResidualJetESeen" , "; ^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualJetThetaSeen = new TH1F( "h_ResidualJetThetaSeen" , "; ^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE} [rad]; Normalized #jet / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_ResidualJetPhiSeen = new TH1F( "h_ResidualJetPhiSeen" , "; ^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE} [rad]; Normalized #jet / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_NormalizedResidualJetPxSeen = new TH1F( "h_NormalizedResidualJetPxSeen" , "; (^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE}) / #sigma_{p_{x,jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetPySeen = new TH1F( "h_NormalizedResidualJetPySeen" , "; (^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE}) / #sigma_{p_{y,jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetPzSeen = new TH1F( "h_NormalizedResidualJetPzSeen" , "; (^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE}) / #sigma_{p_{z,jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetESeen = new TH1F( "h_NormalizedResidualJetESeen" , "; (^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE}) / #sigma_{E_{jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetThetaSeen = new TH1F( "h_NormalizedResidualJetThetaSeen" , "; (^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE}) / #sigma_{#theta_{jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualJetPhiSeen = new TH1F( "h_NormalizedResidualJetPhiSeen" , "; (^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE}) / #sigma_{#phi_{jet}}; Normalized #jet / 0.1" , 220 , -11.0 , 11.0 );

		m_pTTree->Branch("trueIsoLepPx",&m_trueIsoLepPx) ;
		m_pTTree->Branch("trueIsoLepPy",&m_trueIsoLepPy) ;
		m_pTTree->Branch("trueIsoLepPz",&m_trueIsoLepPz) ;
		m_pTTree->Branch("trueIsoLepE",&m_trueIsoLepE) ;
		m_pTTree->Branch("trueIsoLepTheta",&m_trueIsoLepTheta) ;
		m_pTTree->Branch("trueIsoLepPhi",&m_trueIsoLepPhi) ;
		m_pTTree->Branch("seenIsoLepPx",&m_seenIsoLepPx) ;
		m_pTTree->Branch("seenIsoLepPy",&m_seenIsoLepPy) ;
		m_pTTree->Branch("seenIsoLepPz",&m_seenIsoLepPz) ;
		m_pTTree->Branch("seenIsoLepE",&m_seenIsoLepE) ;
		m_pTTree->Branch("seenIsoLepTheta",&m_seenIsoLepTheta) ;
		m_pTTree->Branch("seenIsoLepPhi",&m_seenIsoLepPhi) ;
		m_pTTree->Branch("recoIsoLepPx",&m_recoIsoLepPx) ;
		m_pTTree->Branch("recoIsoLepPy",&m_recoIsoLepPy) ;
		m_pTTree->Branch("recoIsoLepPz",&m_recoIsoLepPz) ;
		m_pTTree->Branch("recoIsoLepE",&m_recoIsoLepE) ;
		m_pTTree->Branch("recoIsoLepTheta",&m_recoIsoLepTheta) ;
		m_pTTree->Branch("recoIsoLepPhi",&m_recoIsoLepPhi) ;
		m_pTTree->Branch("isoLepSigmaPx2",&m_isoLepSigmaPx2) ;
		m_pTTree->Branch("isoLepSigmaPxPy",&m_isoLepSigmaPxPy) ;
		m_pTTree->Branch("isoLepSigmaPy2",&m_isoLepSigmaPy2) ;
		m_pTTree->Branch("isoLepSigmaPxPz",&m_isoLepSigmaPxPz) ;
		m_pTTree->Branch("isoLepSigmaPyPz",&m_isoLepSigmaPyPz) ;
		m_pTTree->Branch("isoLepSigmaPz2",&m_isoLepSigmaPz2) ;
		m_pTTree->Branch("isoLepSigmaPxE",&m_isoLepSigmaPxE) ;
		m_pTTree->Branch("isoLepSigmaPyE",&m_isoLepSigmaPyE) ;
		m_pTTree->Branch("isoLepSigmaPzE",&m_isoLepSigmaPzE) ;
		m_pTTree->Branch("isoLepSigmaE2",&m_isoLepSigmaE2) ;
		m_pTTree->Branch("isoLepSigmaTheta2",&m_isoLepSigmaTheta2) ;
		m_pTTree->Branch("isoLepSigmaPhi2",&m_isoLepSigmaPhi2) ;
		m_pTTree->Branch("ResidualIsoLepPx",&m_ResidualIsoLepPx) ;
		m_pTTree->Branch("ResidualIsoLepPy",&m_ResidualIsoLepPy) ;
		m_pTTree->Branch("ResidualIsoLepPz",&m_ResidualIsoLepPz) ;
		m_pTTree->Branch("ResidualIsoLepE",&m_ResidualIsoLepE) ;
		m_pTTree->Branch("ResidualIsoLepTheta",&m_ResidualIsoLepTheta) ;
		m_pTTree->Branch("ResidualIsoLepPhi",&m_ResidualIsoLepPhi) ;
		m_pTTree->Branch("NormalizedResidualIsoLepPx",&m_NormalizedResidualIsoLepPx) ;
		m_pTTree->Branch("NormalizedResidualIsoLepPy",&m_NormalizedResidualIsoLepPy) ;
		m_pTTree->Branch("NormalizedResidualIsoLepPz",&m_NormalizedResidualIsoLepPz) ;
		m_pTTree->Branch("NormalizedResidualIsoLepE",&m_NormalizedResidualIsoLepE) ;
		m_pTTree->Branch("NormalizedResidualIsoLepTheta",&m_NormalizedResidualIsoLepTheta) ;
		m_pTTree->Branch("NormalizedResidualIsoLepPhi",&m_NormalizedResidualIsoLepPhi) ;
		m_pTTree->Branch("trueIsoLepType",&m_trueIsoLepType) ;
		h_ResidualIsoLepPx = new TH1F( "h_ResidualIsoLepPx" , "; ^{}_{}p_{x,#font[12]{l}^{#pm}}^{REC} - p_{x,#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepPy = new TH1F( "h_ResidualIsoLepPy" , "; ^{}_{}p_{y,#font[12]{l}^{#pm}}^{REC} - p_{y,#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepPz = new TH1F( "h_ResidualIsoLepPz" , "; ^{}_{}p_{z,#font[12]{l}^{#pm}}^{REC} - p_{z,#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepE = new TH1F( "h_ResidualIsoLepE" , "; ^{}_{}E_{#font[12]{l}^{#pm}}^{REC} - E_{#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepTheta = new TH1F( "h_ResidualIsoLepTheta" , "; ^{}_{}#theta_{#font[12]{l}^{#pm}}^{REC} - #theta_{#font[12]{l}^{#pm}}^{TRUE} [rad]; Normalized ##font[12]{l}^{#pm} / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_ResidualIsoLepPhi = new TH1F( "h_ResidualIsoLepPhi" , "; ^{}_{}#phi_{#font[12]{l}^{#pm}}^{REC} - #phi_{#font[12]{l}^{#pm}}^{TRUE} [rad]; Normalized ##font[12]{l}^{#pm} / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_NormalizedResidualIsoLepPx = new TH1F( "h_NormalizedResidualIsoLepPx" , "; (^{}_{}p_{x,#font[12]{l}^{#pm}}^{REC} - p_{x,#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{p_{x,#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepPy = new TH1F( "h_NormalizedResidualIsoLepPy" , "; (^{}_{}p_{y,#font[12]{l}^{#pm}}^{REC} - p_{y,#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{p_{y,#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepPz = new TH1F( "h_NormalizedResidualIsoLepPz" , "; (^{}_{}p_{z,#font[12]{l}^{#pm}}^{REC} - p_{z,#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{p_{z,#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepE = new TH1F( "h_NormalizedResidualIsoLepE" , "; (^{}_{}E_{#font[12]{l}^{#pm}}^{REC} - E_{#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{E_{#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepTheta = new TH1F( "h_NormalizedResidualIsoLepTheta" , "; (^{}_{}#theta_{#font[12]{l}^{#pm}}^{REC} - #theta_{#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{#theta_{#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepPhi = new TH1F( "h_NormalizedResidualIsoLepPhi" , "; (^{}_{}#phi_{#font[12]{l}^{#pm}}^{REC} - #phi_{#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{#phi_{#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_ResidualIsoLepPxSeen = new TH1F( "h_ResidualIsoLepPxSeen" , "; ^{}_{}p_{x,#font[12]{l}^{#pm}}^{REC} - p_{x,#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepPySeen = new TH1F( "h_ResidualIsoLepPySeen" , "; ^{}_{}p_{y,#font[12]{l}^{#pm}}^{REC} - p_{y,#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepPzSeen = new TH1F( "h_ResidualIsoLepPzSeen" , "; ^{}_{}p_{z,#font[12]{l}^{#pm}}^{REC} - p_{z,#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepESeen = new TH1F( "h_ResidualIsoLepESeen" , "; ^{}_{}E_{#font[12]{l}^{#pm}}^{REC} - E_{#font[12]{l}^{#pm}}^{TRUE} [GeV]; Normalized ##font[12]{l}^{#pm} / 0.1 GeV" , 440 , -22.0 , 22.0 );
		h_ResidualIsoLepThetaSeen = new TH1F( "h_ResidualIsoLepThetaSeen" , "; ^{}_{}#theta_{#font[12]{l}^{#pm}}^{REC} - #theta_{#font[12]{l}^{#pm}}^{TRUE} [rad]; Normalized ##font[12]{l}^{#pm} / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_ResidualIsoLepPhiSeen = new TH1F( "h_ResidualIsoLepPhiSeen" , "; ^{}_{}#phi_{#font[12]{l}^{#pm}}^{REC} - #phi_{#font[12]{l}^{#pm}}^{TRUE} [rad]; Normalized ##font[12]{l}^{#pm} / 0.001" , 2000 * 3.15 , -3.15 , 3.15 );
		h_NormalizedResidualIsoLepPxSeen = new TH1F( "h_NormalizedResidualIsoLepPxSeen" , "; (^{}_{}p_{x,#font[12]{l}^{#pm}}^{REC} - p_{x,#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{p_{x,#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepPySeen = new TH1F( "h_NormalizedResidualIsoLepPySeen" , "; (^{}_{}p_{y,#font[12]{l}^{#pm}}^{REC} - p_{y,#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{p_{y,#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepPzSeen = new TH1F( "h_NormalizedResidualIsoLepPzSeen" , "; (^{}_{}p_{z,#font[12]{l}^{#pm}}^{REC} - p_{z,#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{p_{z,#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepESeen = new TH1F( "h_NormalizedResidualIsoLepESeen" , "; (^{}_{}E_{#font[12]{l}^{#pm}}^{REC} - E_{#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{E_{#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepThetaSeen = new TH1F( "h_NormalizedResidualIsoLepThetaSeen" , "; (^{}_{}#theta_{#font[12]{l}^{#pm}}^{REC} - #theta_{#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{#theta_{#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
		h_NormalizedResidualIsoLepPhiSeen = new TH1F( "h_NormalizedResidualIsoLepPhiSeen" , "; (^{}_{}#phi_{#font[12]{l}^{#pm}}^{REC} - #phi_{#font[12]{l}^{#pm}}^{TRUE}) / #sigma_{#phi_{#font[12]{l}^{#pm}}}; Normalized ##font[12]{l}^{#pm} / 0.1" , 220 , -11.0 , 11.0 );
	}
}

void JetErrorAnalysis::Clear()
{
	m_nTrueJets = 0;
	m_nTrueIsoLeps = 0;
	m_nRecoJets = 0;
	m_nIsoLeptons = 0;
	m_nSeenISRs.clear();
	m_nSLDecayBHadron = 0;
	m_nSLDecayCHadron = 0;
	m_nSLDecayTauLepton = 0;
	m_nSLDecayTotal = 0;
	m_nCorrectedSLD = 0;
	m_quarkPx.clear();
	m_quarkPy.clear();
	m_quarkPz.clear();
	m_quarkE.clear();
	m_quarkTheta.clear();
	m_quarkPhi.clear();
	m_trueJetPx.clear();
	m_trueJetPy.clear();
	m_trueJetPz.clear();
	m_trueJetE.clear();
	m_trueJetTheta.clear();
	m_trueJetPhi.clear();
	m_trueSeenJetPx.clear();
	m_trueSeenJetPy.clear();
	m_trueSeenJetPz.clear();
	m_trueSeenJetE.clear();
	m_trueSeenJetTheta.clear();
	m_trueSeenJetPhi.clear();
	m_seenJetPx.clear();
	m_seenJetPy.clear();
	m_seenJetPz.clear();
	m_seenJetE.clear();
	m_seenJetTheta.clear();
	m_seenJetPhi.clear();
	m_recoJetPx.clear();
	m_recoJetPy.clear();
	m_recoJetPz.clear();
	m_recoJetE.clear();
	m_recoJetTheta.clear();
	m_recoJetPhi.clear();
	m_jetSigmaPx2.clear();
	m_jetSigmaPxPy.clear();
	m_jetSigmaPy2.clear();
	m_jetSigmaPxPz.clear();
	m_jetSigmaPyPz.clear();
	m_jetSigmaPz2.clear();
	m_jetSigmaPxE.clear();
	m_jetSigmaPyE.clear();
	m_jetSigmaPzE.clear();
	m_jetSigmaE2.clear();
	m_jetSigmaTheta2.clear();
	m_jetSigmaPhi2.clear();
	m_ResidualJetPx.clear();
	m_ResidualJetPy.clear();
	m_ResidualJetPz.clear();
	m_ResidualJetE.clear();
	m_ResidualJetTheta.clear();
	m_ResidualJetPhi.clear();
	m_NormalizedResidualJetPx.clear();
	m_NormalizedResidualJetPy.clear();
	m_NormalizedResidualJetPz.clear();
	m_NormalizedResidualJetE.clear();
	m_NormalizedResidualJetTheta.clear();
	m_NormalizedResidualJetPhi.clear();
	m_trueJetType.clear();
	m_trueJetFlavour.clear();
	m_recoJetFlavour.clear();

	m_trueIsoLepPx.clear();
	m_trueIsoLepPy.clear();
	m_trueIsoLepPz.clear();
	m_trueIsoLepE.clear();
	m_trueIsoLepTheta.clear();
	m_trueIsoLepPhi.clear();
	m_seenIsoLepPx.clear();
	m_seenIsoLepPy.clear();
	m_seenIsoLepPz.clear();
	m_seenIsoLepE.clear();
	m_seenIsoLepTheta.clear();
	m_seenIsoLepPhi.clear();
	m_recoIsoLepPx.clear();
	m_recoIsoLepPy.clear();
	m_recoIsoLepPz.clear();
	m_recoIsoLepE.clear();
	m_recoIsoLepTheta.clear();
	m_recoIsoLepPhi.clear();
	m_isoLepSigmaPx2.clear();
	m_isoLepSigmaPxPy.clear();
	m_isoLepSigmaPy2.clear();
	m_isoLepSigmaPxPz.clear();
	m_isoLepSigmaPyPz.clear();
	m_isoLepSigmaPz2.clear();
	m_isoLepSigmaPxE.clear();
	m_isoLepSigmaPyE.clear();
	m_isoLepSigmaPzE.clear();
	m_isoLepSigmaE2.clear();
	m_isoLepSigmaTheta2.clear();
	m_isoLepSigmaPhi2.clear();
	m_ResidualIsoLepPx.clear();
	m_ResidualIsoLepPy.clear();
	m_ResidualIsoLepPz.clear();
	m_ResidualIsoLepE.clear();
	m_ResidualIsoLepTheta.clear();
	m_ResidualIsoLepPhi.clear();
	m_NormalizedResidualIsoLepPx.clear();
	m_NormalizedResidualIsoLepPy.clear();
	m_NormalizedResidualIsoLepPz.clear();
	m_NormalizedResidualIsoLepE.clear();
	m_NormalizedResidualIsoLepTheta.clear();
	m_NormalizedResidualIsoLepPhi.clear();
	m_trueIsoLepType.clear();
}

void JetErrorAnalysis::processRunHeader()
{
    m_nRun++ ;
}


void JetErrorAnalysis::processEvent( LCEvent* pLCEvent)
{
	LCCollection *recoJetCol{};
	LCCollection *isoLepCol{};
	//LCCollection *trueJetCol{};
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	this->Clear();
	std::string trueJetType[6]{ "hadronic(string)" , "leptonic" , "hadronic(cluster)" , "ISR" , "overlay" , "M.E. photon" };
	std::string icnType[6]{ "quark pair" , "lepton pair" , "quark pair" , "ISR" , "???" , "M.E. photon" };
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event: 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	try
	{
		LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( _recoMCTruthLink ) );
		LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( _MCTruthRecoLink ) );
		streamlog_out(DEBUG6) << "	Reading collection " << m_recoJetCollectionName << std::endl;
		recoJetCol		= pLCEvent->getCollection( m_recoJetCollectionName );
		m_nRecoJets = recoJetCol->getNumberOfElements();
		streamlog_out(DEBUG6) << "	Collection " << m_recoJetCollectionName << " was red with " << m_nRecoJets << " elements" << std::endl;
		m_nSLDecayBHadron = recoJetCol->getParameters().getIntVal( "nBHadronSLD_found" );
		m_nSLDecayCHadron = recoJetCol->getParameters().getIntVal( "nCHadronSLD_found" );
		m_nSLDecayTauLepton = recoJetCol->getParameters().getIntVal( "nTauLeptonSLD_found" );
		m_nSLDecayTotal = recoJetCol->getParameters().getIntVal( "nTotalSLD_found" );
		m_nCorrectedSLD = recoJetCol->getParameters().getIntVal( "nSolvedSLD" );

		streamlog_out(DEBUG6) << "	Reading collection " << m_inputIsolatedleptonCollection << std::endl;
		isoLepCol		= pLCEvent->getCollection( m_inputIsolatedleptonCollection );
		m_nIsoLeptons = isoLepCol->getNumberOfElements();
		streamlog_out(DEBUG6) << "	Collection " << m_inputIsolatedleptonCollection << " was red with " << m_nIsoLeptons << " elements" << std::endl;
		//trueJetCol		= pLCEvent->getCollection( _trueJetCollectionName );
		TrueJet_Parser* trueJet	= this;
		trueJet->getall( pLCEvent );

		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;

		int njets = trueJet->njets();
		
		//_COUNT_FSR = 1;
		std::vector<int> seenISRinRecoJetIndex{};for ( int i = 0; i < m_nRecoJets; ++i ) seenISRinRecoJetIndex.push_back( -1 );
		streamlog_out(DEBUG8) << "	Number of True Jets: " << njets << std::endl;
		for ( int i_trueJet = 0 ; i_trueJet < njets ; ++i_trueJet )
		{
			streamlog_out(DEBUG7) << "Type of trueJet[ " << i_trueJet << " ] = " << type_jet( i_trueJet ) << std::endl;
			if ( type_jet( i_trueJet ) != 5 )
			{
				streamlog_out(DEBUG4) << "		initial_element (Px,Py,Pz,E):	" << initial_elementon( i_trueJet )->getMomentum()[ 0 ] << " , " << initial_elementon( i_trueJet )->getMomentum()[ 1 ] << " , " << initial_elementon( i_trueJet )->getMomentum()[ 2 ] << " , " << initial_elementon( i_trueJet )->getEnergy() << std::endl;
				streamlog_out(DEBUG4) << "		final_element (Px,Py,Pz,E):	" << final_elementon( i_trueJet )->getMomentum()[ 0 ] << " , " << final_elementon( i_trueJet )->getMomentum()[ 1 ] << " , " << final_elementon( i_trueJet )->getMomentum()[ 2 ] << " , " << final_elementon( i_trueJet )->getEnergy() << std::endl;
				streamlog_out(DEBUG4) << "		Quark (Px,Py,Pz,E):		" << pquark( i_trueJet )[ 0 ] << " , " << pquark( i_trueJet )[ 1 ] << " , " << pquark( i_trueJet )[ 2 ] << " , " << Equark( i_trueJet ) << std::endl;
			}
			streamlog_out(DEBUG4) << "		trueJet (Px,Py,Pz,E):		" << ptrue( i_trueJet )[ 0 ] << " , " << ptrue( i_trueJet )[ 1 ] << " , " << ptrue( i_trueJet )[ 2 ] << " , " << Etrue( i_trueJet ) << std::endl;
			streamlog_out(DEBUG4) << "		trueOfSeenJet (Px,Py,Pz,E):	" << ptrueseen( i_trueJet )[ 0 ] << " , " << ptrueseen( i_trueJet )[ 1 ] << " , " << ptrueseen( i_trueJet )[ 2 ] << " , " << Etrueseen( i_trueJet ) << std::endl;
			streamlog_out(DEBUG4) << "		seenJet (Px,Py,Pz,E):		" << pseen( i_trueJet )[ 0 ] << " , " << pseen( i_trueJet )[ 1 ] << " , " << pseen( i_trueJet )[ 2 ] << " , " << Eseen( i_trueJet ) << std::endl;
			if ( type_jet( i_trueJet ) == 4 )
			{
				m_nSeenISRs.push_back( seen_partics( i_trueJet ).size() );
				streamlog_out(DEBUG4) << "		This jet has " << seen_partics( i_trueJet ).size() << " reconstructed ISR photon(s)" << std::endl;
				TVector3 seenISRDirecion( pseen( i_trueJet ) ); seenISRDirecion.SetMag( 1.0 );
				double cosMaxAngle = -1.0;
				int i_matchedRecoJet = -1;
				for ( int i_recoJet = 0 ; i_recoJet < m_nRecoJets ; ++i_recoJet )
				{
					ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( i_recoJet ) );
					TVector3 recoJetDirection( recoJet->getMomentum() ); recoJetDirection.SetMag( 1.0 );
					if ( recoJetDirection.Dot( seenISRDirecion ) > cosMaxAngle )
					{
						cosMaxAngle = recoJetDirection.Dot( seenISRDirecion );
						i_matchedRecoJet = i_recoJet;
					}
				}
				seenISRinRecoJetIndex[ i_matchedRecoJet ] = i_trueJet;
			}
			
			streamlog_out(DEBUG4) << "" << std::endl;
		}
		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(DEBUG8) << "	Number of Reconstructed Jets: " << m_nRecoJets << std::endl;
		for ( int i_recoJet = 0 ; i_recoJet < m_nRecoJets ; ++i_recoJet )
		{
			ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( i_recoJet ) );
			streamlog_out(DEBUG4) << "		recoJet[ " << i_recoJet << " ] (Px,Py,Pz,E):		" << recoJet->getMomentum()[ 0 ] << " , " << recoJet->getMomentum()[ 1 ] << " , " << recoJet->getMomentum()[ 2 ] << " , " << recoJet->getEnergy() << std::endl;
			if ( seenISRinRecoJetIndex[ i_recoJet ] != -1 ) streamlog_out(DEBUG4) << "		recoJet[ " << i_recoJet << " ] contains one reconstructed ISR photon at index " << seenISRinRecoJetIndex[ i_recoJet ] << std::endl; 
			streamlog_out(DEBUG4) << "" << std::endl;
		}
		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(DEBUG8) << "	Number of Isolated Leptons: " << m_nIsoLeptons << std::endl;
		for ( int i_isoLep = 0 ; i_isoLep < m_nIsoLeptons ; ++i_isoLep )
		{
			ReconstructedParticle *isoLep = dynamic_cast<ReconstructedParticle*>( isoLepCol->getElementAt( i_isoLep ) );
			streamlog_out(DEBUG4) << "		isoLep[ " << i_isoLep << " ] (Px,Py,Pz,E):		" << isoLep->getMomentum()[ 0 ] << " , " << isoLep->getMomentum()[ 1 ] << " , " << isoLep->getMomentum()[ 2 ] << " , " << isoLep->getEnergy() << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
		}
		streamlog_out(DEBUG1) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(DEBUG1) << "	Number of Initial Colour Neutrals: " << nicn() << std::endl;
		for ( int i_icn = 0 ; i_icn < nicn() ; ++i_icn )
		{
			streamlog_out(DEBUG1) << "Initial Colour Neutral[ " << i_icn << " ], Type = " << type_icn_parent( i_icn ) << " , PDG = " << pdg_icn_parent( i_icn ) << " with " << elementons_initial_cn( i_icn ).size() << " MCParticles" << std::endl;
			for ( unsigned int i_mcp = 0 ; i_mcp < elementons_initial_cn( i_icn ).size() ; ++i_mcp )
			{
				streamlog_out(DEBUG1) << "MCParticle[ " << i_mcp << " ]: " << std::endl;
				streamlog_out(DEBUG1) << *elementons_initial_cn( i_icn )[ i_mcp ] << std::endl;
			}
		}
		streamlog_out(DEBUG1) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(DEBUG1) << "	Number of Final Colour Neutrals: " << nfcn() << std::endl;
		for ( int i_fcn = 0 ; i_fcn < nfcn() ; ++i_fcn )
		{
			streamlog_out(DEBUG1) << "Final Colour Neutral[ " << i_fcn << " ], Type = " << type_fcn_parent( i_fcn ) << " , PDG = " << pdg_fcn_parent( i_fcn ) << " with " << elementons_final_cn( i_fcn ).size() << " MCParticles" << std::endl;
			for ( unsigned int i_mcp = 0 ; i_mcp < elementons_final_cn( i_fcn ).size() ; ++i_mcp )
			{
				streamlog_out(DEBUG1) << "MCParticle[ " << i_mcp << " ]: " << std::endl;
				streamlog_out(DEBUG1) << *elementons_final_cn( i_fcn )[ i_mcp ] << std::endl;
			}
		}
		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		if ( m_nRecoJets == 0 || njets == 0 ) return;

		std::vector<ReconstructedParticle*> recoJetVector{};
		std::vector<ReconstructedParticle*> leadingParticles{};
		std::vector<int> trueJetVectorIndex{};

		for ( int i_recoJet = 0 ; i_recoJet < m_nRecoJets ; ++i_recoJet )
		{
			ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( i_recoJet ) );
			if ( m_jetMatchingMethod == 1 )
			{
				ReconstructedParticle* leadingParticle = NULL;
				float leadingEnergy = 0.0;
				streamlog_out(DEBUG2) << "Looking for leading particle in recoJet[ " << i_recoJet << " ] with " << ( recoJet->getParticles() ).size() << " particles" << std::endl;
				for ( unsigned int i_par = 0 ; i_par < ( recoJet->getParticles() ).size() ; ++i_par )
				{
					ReconstructedParticle* pfo = ( ReconstructedParticle* )recoJet->getParticles()[ i_par ];
					streamlog_out(DEBUG0) << "Checking particle " << i_par << " with Energy = " << pfo->getEnergy() << std::endl;
					if ( abs( pfo->getType() ) == 12 || abs( pfo->getType() ) == 14 || abs( pfo->getType() ) == 16 ) continue;
					if ( pfo->getEnergy() > leadingEnergy )
					{
						leadingParticle = pfo;
						leadingEnergy = pfo->getEnergy();
						streamlog_out(DEBUG0) << "Leading particle so far: " << std::endl;
						streamlog_out(DEBUG0) << *leadingParticle << std::endl;
					}
				}
				streamlog_out(DEBUG2) << "Leading particle in recoJet [ " << i_recoJet << " ]:" << std::endl;
				streamlog_out(DEBUG2) << *leadingParticle << std::endl;

				if ( leadingParticle != NULL )
				{
					int trueJetIndex = recojet( leadingParticle );
					//LCObjectVec jetvec = reltjreco->getRelatedFromObjects( leadingParticle );
					//streamlog_out(DEBUG2) << jetvec.size() << " true Jet found for leading particle of jet " << i_recoJet << std::endl;
					//if ( jetvec.size() != 0 && type_jet( jetindex( dynamic_cast<ReconstructedParticle*>( jetvec[ 0 ] ) ) ) == 1 )
					if ( trueJetIndex >= 0 && type_jet( trueJetIndex ) == 1 )
					{
						//trueJetIndex = jetindex( dynamic_cast<ReconstructedParticle*>( jetvec[ 0 ] ) );
						recoJetVector.push_back( recoJet );
						leadingParticles.push_back( leadingParticle );
						trueJetVectorIndex.push_back( trueJetIndex );
						++m_nTrueJets;
						streamlog_out(DEBUG2) << " The index of trueJet containing leading particle of recoJet [ " << i_recoJet << " ]: " << trueJetIndex << " , trueJet Type = " << type_jet( trueJetIndex ) << std::endl;
					}
				}
			}
			else if ( m_jetMatchingMethod == 2 )
			{
				streamlog_out(DEBUG4) << "	Finding closest trueJet to recoJet[ " << i_recoJet << " ]" << std::endl;
				TVector3 recoJetDirection( recoJet->getMomentum() ); recoJetDirection.SetMag( 1.0 );
				streamlog_out(DEBUG4) << "		Direction of recoJet (Ux,Uy,Uz):	" << recoJetDirection.X() << " , " << recoJetDirection.Y() << " , " << recoJetDirection.Z() << std::endl;
				double cosMaxAngle = -1.0;
				int i_matchedtrueJet = -1;
				for ( int i_trueJet = 0 ; i_trueJet < njets ; ++i_trueJet )
				{
					if ( type_jet( i_trueJet ) == 1 )
					{
						TVector3 trueJetDirection( initial_elementon( i_trueJet )->getMomentum() ); trueJetDirection.SetMag( 1.0 );
						streamlog_out(DEBUG4) << "		Checking trueJet[ " << i_trueJet << " ] with Type: " << type_jet( i_trueJet ) << std::endl;
						streamlog_out(DEBUG4) << "		Direction of trueJet (Ux,Uy,Uz):	" << trueJetDirection.X() << " , " << trueJetDirection.Y() << " , " << trueJetDirection.Z() << std::endl;
						if ( recoJetDirection.Dot( trueJetDirection ) > cosMaxAngle )
						{
							bool trueJetExist = false;
							for ( unsigned int i_jet = 0 ; i_jet < trueJetVectorIndex.size() ; ++i_jet )
							{
								if ( initial_elementon( trueJetVectorIndex[ i_jet ] ) == initial_elementon( i_trueJet ) ) trueJetExist = true;
							}
							if ( !trueJetExist )
							{
								cosMaxAngle = recoJetDirection.Dot( trueJetDirection );
								i_matchedtrueJet = i_trueJet;
							}
						}
					}
				}
				if ( i_matchedtrueJet != -1 )
				{
					recoJetVector.push_back( recoJet );
					trueJetVectorIndex.push_back( i_matchedtrueJet );
					++m_nTrueJets;
					streamlog_out(DEBUG2) << " The index of trueJet with lowest angle to recoJet [ " << i_recoJet << " ]: " << i_matchedtrueJet << " , trueJet Type = " << type_jet( i_matchedtrueJet ) << std::endl;
				}
			}
		}
		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(DEBUG4) << "	" << recoJetVector.size() << " recoJets and " << trueJetVectorIndex.size() << " true Jets are found and matched" << std::endl;
		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		for ( unsigned int i_jet = 0 ; i_jet < recoJetVector.size() ; ++i_jet )
		{
			TLorentzVector initial_elementFourMomentum( initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum() , initial_elementon( trueJetVectorIndex[ i_jet ] )->getEnergy() );
			TLorentzVector quarkFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector trueJetFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector trueOfSeenJetFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector seenJetFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector recoJetFourMomentum( recoJetVector[ i_jet ]->getMomentum() , recoJetVector[ i_jet ]->getEnergy() );
			quarkFourMomentum = TLorentzVector( initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum() , initial_elementon( trueJetVectorIndex[ i_jet ] )->getEnergy() );
			for ( int i_trueJet = 0 ; i_trueJet < trueJet->njets() ; ++i_trueJet )
			{
				if ( initial_elementon( trueJetVectorIndex[ i_jet ] ) == initial_elementon( i_trueJet ) )
				{
					//quarkFourMomentum += TLorentzVector( pquark( i_trueJet ) , Equark( i_trueJet ) );
					trueJetFourMomentum += TLorentzVector( ptrue( i_trueJet ) , Etrue( i_trueJet ) );
					trueOfSeenJetFourMomentum += TLorentzVector( ptrueseen( i_trueJet ) , Etrueseen( i_trueJet ) );
					seenJetFourMomentum += TLorentzVector( pseen( i_trueJet ) , Eseen( i_trueJet ) );
				}
			}
			for ( unsigned int i_daughter = 0 ; i_daughter < initial_elementon( trueJetVectorIndex[ i_jet ] )->getDaughters().size() ; ++i_daughter )
			{
				MCParticle *daughter = initial_elementon( trueJetVectorIndex[ i_jet ] )->getDaughters()[ i_daughter ];
				mcpVector trueFSRPhotons{}; trueFSRPhotons.clear();
				pfoVector recoFSRPhotons{}; recoFSRPhotons.clear();
				ReconstructedParticle *FSRPhoton = NULL;
				if ( abs( daughter->getPDG() ) == 22 && daughter->getGeneratorStatus() == 1 )
				{
					trueJetFourMomentum += TLorentzVector( daughter->getMomentum() , daughter->getEnergy() );
					float weightPFOtoMCP = 0.0;
					float weightMCPtoPFO = 0.0;
					FSRPhoton = getLinkedPFO( daughter , RecoMCParticleNav , MCParticleRecoNav , false , true , weightPFOtoMCP , weightMCPtoPFO );
					if ( FSRPhoton != NULL )
					{
						recoFSRPhotons.push_back( FSRPhoton );
					}
					else if ( daughter->getDaughters().size() > 0 )
					{
						ReconstructedParticle *recoDaughter = NULL;
						for ( unsigned int i_daughtersOfDaughters = 0 ; i_daughtersOfDaughters < daughter->getDaughters().size() ; ++i_daughtersOfDaughters )
						{
							streamlog_out(DEBUG1) << "Looking for reconstructed particle linked to daughter of true FSR photon:" << std::endl;
							streamlog_out(DEBUG1) << *( daughter->getDaughters()[ i_daughtersOfDaughters ] ) << std::endl;
							recoDaughter = getLinkedPFO( daughter->getDaughters()[ i_daughtersOfDaughters ] , RecoMCParticleNav , MCParticleRecoNav , false , false , weightPFOtoMCP , weightMCPtoPFO );
							bool linkedRecoParExist = false;
							for ( unsigned int i_recoPar = 0 ; i_recoPar < recoFSRPhotons.size() ; ++i_recoPar )
							{
								if ( recoFSRPhotons[ i_recoPar ] == recoDaughter ) linkedRecoParExist = true;
							}
							if ( recoDaughter != NULL && !linkedRecoParExist )
							{
								streamlog_out(DEBUG1) << "Linked reconstructed particle to daughter of true FSR photon:" << std::endl;
								streamlog_out(DEBUG1) << *recoDaughter << std::endl;
								recoFSRPhotons.push_back( recoDaughter );
							}
						}
					}
					for ( unsigned int i_d = 0 ; i_d < recoFSRPhotons.size() ; ++i_d )
					{
						if ( i_d == 0 )
						{
							trueOfSeenJetFourMomentum += TLorentzVector( daughter->getMomentum() , daughter->getEnergy() );
							streamlog_out(DEBUG1) << "A FSR photon is added to true/seen jets:" << std::endl;
							streamlog_out(DEBUG1) << *daughter << std::endl;
						}
						bool linkedRecoParExist = false;
						for ( int i_trueJet = 0 ; i_trueJet < trueJet->njets() ; ++i_trueJet )
						{
							for ( unsigned int i_recoPar = 0 ; i_recoPar < seen_partics( i_trueJet ).size() ; ++i_recoPar )
							{
								if ( seen_partics( i_trueJet )[ i_recoPar ] == recoFSRPhotons[ i_d ] ) linkedRecoParExist = true;
							}
						}
						if ( !linkedRecoParExist )
						{
							seenJetFourMomentum += TLorentzVector( recoFSRPhotons[ i_d ]->getMomentum() , recoFSRPhotons[ i_d ]->getEnergy() );
							streamlog_out(DEBUG1) << *recoFSRPhotons[ i_d ] << std::endl;
						}
					}
				}
			}
			//if ( seenISRinRecoJetIndex[ i_jet ] != -1 ) seenJetFourMomentum += TLorentzVector( pseen( seenISRinRecoJetIndex[ i_jet ] ) , Eseen( seenISRinRecoJetIndex[ i_jet ] ) );
			/*
			trueJetFourMomentum = TLorentzVector( ptrue( trueJetVectorIndex[ i_jet ] ) , Etrue( trueJetVectorIndex[ i_jet ] ) );
			trueOfSeenJetFourMomentum = TLorentzVector( ptrueseen( trueJetVectorIndex[ i_jet ] ) , Etrueseen( trueJetVectorIndex[ i_jet ] ) );
			seenJetFourMomentum = TLorentzVector( pseen( trueJetVectorIndex[ i_jet ] ) , Eseen( trueJetVectorIndex[ i_jet ] ) );
			*/

			streamlog_out(DEBUG4) << "----------------------------------------------------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG4) << "Reconstructed Jet[ " << i_jet << " ] matches true jet [ " << trueJetVectorIndex[ i_jet ] << " ]" << std::endl;
			streamlog_out(DEBUG4) << "Four-Momenta using siblings" << std::endl;
			streamlog_out(DEBUG4) << "	trueJet TYPE:	" << type_jet( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG4) << "	initial_element (Px,Py,Pz,E):	" << initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 0 ] << " , " << initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 1 ] << " , " << initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 2 ] << " , " << initial_elementon( trueJetVectorIndex[ i_jet ] )->getEnergy() << std::endl;
			streamlog_out(DEBUG4) << "	Quark (Px,Py,Pz,E):		" << quarkFourMomentum.Px() << " , " << quarkFourMomentum.Py() << " , " << quarkFourMomentum.Pz() << " , " << quarkFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "	trueJet (Px,Py,Pz,E):		" << trueJetFourMomentum.Px() << " , " << trueJetFourMomentum.Py() << " , " << trueJetFourMomentum.Pz() << " , " << trueJetFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "	trueOfSeenJet (Px,Py,Pz,E):	" << trueOfSeenJetFourMomentum.Px() << " , " << trueOfSeenJetFourMomentum.Py() << " , " << trueOfSeenJetFourMomentum.Pz() << " , " << trueOfSeenJetFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "	seenJet (Px,Py,Pz,E):		" << seenJetFourMomentum.Px() << " , " << seenJetFourMomentum.Py() << " , " << seenJetFourMomentum.Pz() << " , " << seenJetFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "	recoJet (Px,Py,Pz,E):		" << recoJetFourMomentum.Px() << " , " << recoJetFourMomentum.Py() << " , " << recoJetFourMomentum.Pz() << " , " << recoJetFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "----------------------------------------------------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG4) << "Four-Momenta using Indices" << std::endl;
			streamlog_out(DEBUG4) << "	trueJet TYPE:	" << type_jet( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG4) << "	initial_element (Px,Py,Pz,E):	" << initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 0 ] << " , " << initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 1 ] << " , " << initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 2 ] << " , " << initial_elementon( trueJetVectorIndex[ i_jet ] )->getEnergy() << std::endl;
			streamlog_out(DEBUG4) << "	Quark (Px,Py,Pz,E):		" << pquark( trueJetVectorIndex[ i_jet ] )[ 0 ] << " , " << pquark( trueJetVectorIndex[ i_jet ] )[ 1 ] << " , " << pquark( trueJetVectorIndex[ i_jet ] )[ 2 ] << " , " << Equark( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG4) << "	trueJet (Px,Py,Pz,E):		" << ptrue( trueJetVectorIndex[ i_jet ] )[ 0 ] << " , " << ptrue( trueJetVectorIndex[ i_jet ] )[ 1 ] << " , " << ptrue( trueJetVectorIndex[ i_jet ] )[ 2 ] << " , " << Etrue( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG4) << "	trueOfSeenJet (Px,Py,Pz,E):	" << ptrueseen( trueJetVectorIndex[ i_jet ] )[ 0 ] << " , " << ptrueseen( trueJetVectorIndex[ i_jet ] )[ 1 ] << " , " << ptrueseen( trueJetVectorIndex[ i_jet ] )[ 2 ] << " , " << Etrueseen( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG4) << "	seenJet (Px,Py,Pz,E):		" << pseen( trueJetVectorIndex[ i_jet ] )[ 0 ] << " , " << pseen( trueJetVectorIndex[ i_jet ] )[ 1 ] << " , " << pseen( trueJetVectorIndex[ i_jet ] )[ 2 ] << " , " << Eseen( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG4) << "	recoJet (Px,Py,Pz,E):		" << recoJetVector[ i_jet ]->getMomentum()[ 0 ] << " , " << recoJetVector[ i_jet ]->getMomentum()[ 1 ] << " , " << recoJetVector[ i_jet ]->getMomentum()[ 2 ] << " , " << recoJetVector[ i_jet ]->getEnergy() << std::endl;
			streamlog_out(DEBUG3) << "	recoJet[ " << i_jet << " ]:" << std::endl;
			streamlog_out(DEBUG3) << *recoJetVector[ i_jet ] << std::endl;
			if ( m_jetMatchingMethod == 1 ) streamlog_out(DEBUG3) << "	leadingParticle:" << std::endl;
			if ( m_jetMatchingMethod == 1 ) streamlog_out(DEBUG3) << *leadingParticles[ i_jet ] << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			m_trueJetType.push_back( type_jet( trueJetVectorIndex[ i_jet ] ) );
			/*
			TVector3 quarkMomentum( initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 0 ] , initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 1 ] , initial_elementon( trueJetVectorIndex[ i_jet ] )->getMomentum()[ 2 ] );
			double quarkEnergy = initial_elementon( trueJetVectorIndex[ i_jet ] )->getEnergy();
			//TVector3 quarkMomentum( pquark( trueJetVectorIndex[ i_jet ] )[ 0 ] , pquark( trueJetVectorIndex[ i_jet ] )[ 1 ] , pquark( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			//double quarkEnergy = Equark( trueJetVectorIndex[ i_jet ] );
			TVector3 jetTrueMomentum( ptrue( trueJetVectorIndex[ i_jet ] )[ 0 ] , ptrue( trueJetVectorIndex[ i_jet ] )[ 1 ] , ptrue( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			double jetTrueEnergy = Etrue( trueJetVectorIndex[ i_jet ] );
			TVector3 jetTrueSeenMomentum( ptrueseen( trueJetVectorIndex[ i_jet ] )[ 0 ] , ptrueseen( trueJetVectorIndex[ i_jet ] )[ 1 ] , ptrueseen( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			double jetTrueSeenEnergy = Etrueseen( trueJetVectorIndex[ i_jet ] );
			TVector3 jetSeenMomentum( pseen( trueJetVectorIndex[ i_jet ] )[ 0 ] , pseen( trueJetVectorIndex[ i_jet ] )[ 1 ] , pseen( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			double jetSeenEnergy = Eseen( trueJetVectorIndex[ i_jet ] );
			TVector3 jetRecoMomentum( recoJetVector[ i_jet ]->getMomentum() );
			double jetRecoEnergy = recoJetVector[ i_jet ]->getEnergy();
			*/

			TVector3 quarkMomentum = quarkFourMomentum.Vect(); double quarkEnergy = quarkFourMomentum.E();
			TVector3 jetTrueMomentum = trueJetFourMomentum.Vect(); double jetTrueEnergy = trueJetFourMomentum.E();
			TVector3 jetTrueSeenMomentum = trueOfSeenJetFourMomentum.Vect(); double jetTrueSeenEnergy = trueOfSeenJetFourMomentum.E();
			TVector3 jetSeenMomentum = seenJetFourMomentum.Vect(); double jetSeenEnergy = seenJetFourMomentum.E();
			TVector3 jetRecoMomentum = recoJetFourMomentum.Vect(); double jetRecoEnergy = recoJetFourMomentum.E();
			double jetPxResidual , jetPyResidual , jetPzResidual , jetEnergyResidual , jetThetaResidual , jetPhiResidual;
			double jetPxResidualSeen , jetPyResidualSeen , jetPzResidualSeen , jetEnergyResidualSeen , jetThetaResidualSeen , jetPhiResidualSeen;
			double jetSigmaE , jetSigmaTheta , jetSigmaPhi;
			getResiduals( jetTrueMomentum , jetTrueEnergy , jetRecoMomentum , jetRecoEnergy , jetPxResidual , jetPyResidual , jetPzResidual , jetEnergyResidual , jetThetaResidual , jetPhiResidual );
			getResiduals( jetTrueSeenMomentum , jetTrueSeenEnergy , jetRecoMomentum , jetRecoEnergy , jetPxResidualSeen , jetPyResidualSeen , jetPzResidualSeen , jetEnergyResidualSeen , jetThetaResidualSeen , jetPhiResidualSeen );
			getResolutions( TLorentzVector( jetRecoMomentum , jetRecoEnergy ) , recoJetVector[ i_jet ]->getCovMatrix() , jetSigmaE , jetSigmaTheta , jetSigmaPhi );
			if ( m_fillRootTree )
			{
				m_quarkPx.push_back( quarkMomentum.Px() );		m_quarkPy.push_back( quarkMomentum.Py() );			m_quarkPz.push_back( quarkMomentum.Pz() );
				m_quarkE.push_back( quarkEnergy );			m_quarkTheta.push_back( quarkMomentum.Theta() );		m_quarkPhi.push_back( quarkMomentum.Phi() );
				m_trueJetPx.push_back( jetTrueMomentum.Px() );		m_trueJetPy.push_back( jetTrueMomentum.Py() );			m_trueJetPz.push_back( jetTrueMomentum.Pz() );
				m_trueJetE.push_back( jetTrueEnergy );			m_trueJetTheta.push_back( jetTrueMomentum.Theta() );		m_trueJetPhi.push_back( jetTrueMomentum.Phi() );
				m_trueSeenJetPx.push_back( jetTrueSeenMomentum.Px() );	m_trueSeenJetPy.push_back( jetTrueSeenMomentum.Py() );		m_trueSeenJetPz.push_back( jetTrueSeenMomentum.Pz() );
				m_trueSeenJetE.push_back( jetTrueSeenEnergy );		m_trueSeenJetTheta.push_back( jetTrueSeenMomentum.Theta() );	m_trueSeenJetPhi.push_back( jetTrueSeenMomentum.Phi() );
				m_seenJetPx.push_back( jetSeenMomentum.Px() );		m_seenJetPy.push_back( jetSeenMomentum.Py() );			m_seenJetPz.push_back( jetSeenMomentum.Pz() );
				m_seenJetE.push_back( jetSeenEnergy );			m_seenJetTheta.push_back( jetSeenMomentum.Theta() );		m_seenJetPhi.push_back( jetSeenMomentum.Phi() );
				m_recoJetPx.push_back( jetRecoMomentum.Px() );		m_recoJetPy.push_back( jetRecoMomentum.Py() );			m_recoJetPz.push_back( jetRecoMomentum.Pz() );
				m_recoJetE.push_back( jetRecoEnergy );			m_recoJetTheta.push_back( jetRecoMomentum.Theta() );		m_recoJetPhi.push_back( jetRecoMomentum.Phi() );
				m_jetSigmaPx2.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] );
				m_jetSigmaPxPy.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 1 ] );
				m_jetSigmaPy2.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] );
				m_jetSigmaPxPz.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 3 ] );
				m_jetSigmaPyPz.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 4 ] );
				m_jetSigmaPz2.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] );
				m_jetSigmaPxE.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 6 ] );
				m_jetSigmaPyE.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 7 ] );
				m_jetSigmaPzE.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 8 ] );
				m_jetSigmaE2.push_back( recoJetVector[ i_jet ]->getCovMatrix()[ 9 ] );
				m_jetSigmaTheta2.push_back( pow( jetSigmaTheta , 2 ) );
				m_jetSigmaPhi2.push_back( pow( jetSigmaPhi , 2 ) );
				m_ResidualJetPx.push_back( jetPxResidual ); h_ResidualJetPx->Fill( jetPxResidual );
				m_ResidualJetPy.push_back( jetPyResidual ); h_ResidualJetPy->Fill( jetPyResidual );
				m_ResidualJetPz.push_back( jetPzResidual ); h_ResidualJetPz->Fill( jetPzResidual );
				m_ResidualJetE.push_back( jetEnergyResidual ); h_ResidualJetE->Fill( jetEnergyResidual );
				m_ResidualJetTheta.push_back( jetThetaResidual ); h_ResidualJetTheta->Fill( jetThetaResidual );
				m_ResidualJetPhi.push_back( jetPhiResidual ); h_ResidualJetPhi->Fill( jetPhiResidual );
				m_NormalizedResidualJetPx.push_back( jetPxResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
				m_NormalizedResidualJetPy.push_back( jetPyResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
				m_NormalizedResidualJetPz.push_back( jetPzResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
				m_NormalizedResidualJetE.push_back( jetEnergyResidual / jetSigmaE );
				m_NormalizedResidualJetTheta.push_back( jetThetaResidual / jetSigmaTheta );
				m_NormalizedResidualJetPhi.push_back( jetPhiResidual / jetSigmaPhi );
				h_NormalizedResidualJetPx->Fill( jetPxResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
				h_NormalizedResidualJetPy->Fill( jetPyResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
				h_NormalizedResidualJetPz->Fill( jetPzResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
				h_NormalizedResidualJetE->Fill( jetEnergyResidual / jetSigmaE );
				h_NormalizedResidualJetTheta->Fill( jetThetaResidual / jetSigmaTheta );
				h_NormalizedResidualJetPhi->Fill( jetPhiResidual / jetSigmaPhi );
				m_ResidualJetPxSeen.push_back( jetPxResidualSeen ); h_ResidualJetPxSeen->Fill( jetPxResidualSeen );
				m_ResidualJetPySeen.push_back( jetPyResidualSeen ); h_ResidualJetPySeen->Fill( jetPyResidualSeen );
				m_ResidualJetPzSeen.push_back( jetPzResidualSeen ); h_ResidualJetPzSeen->Fill( jetPzResidualSeen );
				m_ResidualJetESeen.push_back( jetEnergyResidualSeen ); h_ResidualJetESeen->Fill( jetEnergyResidualSeen );
				m_ResidualJetThetaSeen.push_back( jetThetaResidualSeen ); h_ResidualJetThetaSeen->Fill( jetThetaResidualSeen );
				m_ResidualJetPhiSeen.push_back( jetPhiResidualSeen ); h_ResidualJetPhiSeen->Fill( jetPhiResidualSeen );
				m_NormalizedResidualJetPxSeen.push_back( jetPxResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
				m_NormalizedResidualJetPySeen.push_back( jetPyResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
				m_NormalizedResidualJetPzSeen.push_back( jetPzResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
				m_NormalizedResidualJetESeen.push_back( jetEnergyResidualSeen / jetSigmaE );
				m_NormalizedResidualJetThetaSeen.push_back( jetThetaResidualSeen / jetSigmaTheta );
				m_NormalizedResidualJetPhiSeen.push_back( jetPhiResidualSeen / jetSigmaPhi );
				h_NormalizedResidualJetPxSeen->Fill( jetPxResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
				h_NormalizedResidualJetPySeen->Fill( jetPyResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
				h_NormalizedResidualJetPzSeen->Fill( jetPzResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
				h_NormalizedResidualJetESeen->Fill( jetEnergyResidualSeen / jetSigmaE );
				h_NormalizedResidualJetThetaSeen->Fill( jetThetaResidualSeen / jetSigmaTheta );
				h_NormalizedResidualJetPhiSeen->Fill( jetPhiResidualSeen / jetSigmaPhi );
			}
		}

		std::vector<ReconstructedParticle*> recoIsoLepVector{};
		std::vector<int> trueIsoLepVectorIndex{};
		for ( int i_isoLep = 0 ; i_isoLep < m_nIsoLeptons ; ++i_isoLep )
		{
			ReconstructedParticle *recoIsoLep = dynamic_cast<ReconstructedParticle*>( isoLepCol->getElementAt( i_isoLep ) );
			streamlog_out(DEBUG4) << "	Finding closest trueIsoLep to recoIsoLep[ " << i_isoLep << " ]" << std::endl;
			TVector3 recoIsoLepDirection( recoIsoLep->getMomentum() ); recoIsoLepDirection.SetMag( 1.0 );
			streamlog_out(DEBUG4) << "		Direction of recoIsoLep (Ux,Uy,Uz):	" << recoIsoLepDirection.X() << " , " << recoIsoLepDirection.Y() << " , " << recoIsoLepDirection.Z() << std::endl;
			double cosMaxAngle = -1.0;
			int i_matchedtrueIsoLep = -1;
			for ( int i_trueIsoLep = 0 ; i_trueIsoLep < njets ; ++i_trueIsoLep )
			{
				if ( type_jet( i_trueIsoLep ) == 2 )
				{
					TVector3 trueIsoLepDirection( initial_elementon( i_trueIsoLep )->getMomentum() ); trueIsoLepDirection.SetMag( 1.0 );
					streamlog_out(DEBUG4) << "		Checking trueIsoLep[ " << i_trueIsoLep << " ] with Type: " << type_jet( i_trueIsoLep ) << std::endl;
					streamlog_out(DEBUG4) << "		Direction of trueIsoLep (Ux,Uy,Uz):	" << trueIsoLepDirection.X() << " , " << trueIsoLepDirection.Y() << " , " << trueIsoLepDirection.Z() << std::endl;
					if ( recoIsoLepDirection.Dot( trueIsoLepDirection ) > cosMaxAngle )
					{
						bool trueIsoLepExist = false;
						for ( unsigned int i_lep = 0 ; i_lep < trueIsoLepVectorIndex.size() ; ++i_lep )
						{
							if ( initial_elementon( trueIsoLepVectorIndex[ i_lep ] ) == initial_elementon( i_trueIsoLep ) ) trueIsoLepExist = true;
						}
						if ( !trueIsoLepExist )
						{
							cosMaxAngle = recoIsoLepDirection.Dot( trueIsoLepDirection );
							i_matchedtrueIsoLep = i_trueIsoLep;
						}
					}
				}
			}
			if ( i_matchedtrueIsoLep != -1 )
			{
				recoIsoLepVector.push_back( recoIsoLep );
				trueIsoLepVectorIndex.push_back( i_matchedtrueIsoLep );
				++m_nTrueIsoLeps;
				streamlog_out(DEBUG2) << " The index of trueIsoLep with lowest angle to recoIsoLep [ " << i_isoLep << " ]: " << i_matchedtrueIsoLep << " , trueIsoLep Type = " << type_jet( i_matchedtrueIsoLep ) << std::endl;
			}
		}
		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(DEBUG4) << "	" << recoIsoLepVector.size() << " reco IsoLeptons and " << trueIsoLepVectorIndex.size() << " true IsoLeptons are found and matched" << std::endl;
		streamlog_out(DEBUG8) << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		for ( unsigned int i_lep = 0 ; i_lep < recoIsoLepVector.size() ; ++i_lep )
		{
			TLorentzVector initial_elementFourMomentum( initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getMomentum() , initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getEnergy() );
			TLorentzVector trueIsoLepFourMomentum( initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getMomentum() , initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getEnergy() );
			TLorentzVector seenIsoLepFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector recoIsoLepFourMomentum( recoIsoLepVector[ i_lep ]->getMomentum() , recoIsoLepVector[ i_lep ]->getEnergy() );
			for ( int i_trueIsoLep = 0 ; i_trueIsoLep < trueJet->njets() ; ++i_trueIsoLep )
			{
				if ( initial_elementon( trueIsoLepVectorIndex[ i_lep ] ) == initial_elementon( i_trueIsoLep ) )
				{
					for ( unsigned int i_par = 0; i_par < seen_partics( i_trueIsoLep ).size() ; ++i_par )
					{
						seenIsoLepFourMomentum += TLorentzVector( ( seen_partics( i_trueIsoLep )[ i_par ] )->getMomentum() , ( seen_partics( i_trueIsoLep )[ i_par ] )->getEnergy() );
					}
				}
			}

			streamlog_out(DEBUG4) << "----------------------------------------------------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG4) << "Reconstructed IsoLep[ " << i_lep << " ] matches true IsoLep [ " << trueIsoLepVectorIndex[ i_lep ] << " ]" << std::endl;
			streamlog_out(DEBUG4) << "Four-Momenta using siblings" << std::endl;
			streamlog_out(DEBUG4) << "	trueIsoLep TYPE:	" << type_jet( trueIsoLepVectorIndex[ i_lep ] ) << std::endl;
			streamlog_out(DEBUG4) << "	initial_element (Px,Py,Pz,E):	" << initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getMomentum()[ 0 ] << " , " << initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getMomentum()[ 1 ] << " , " << initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getMomentum()[ 2 ] << " , " << initial_elementon( trueIsoLepVectorIndex[ i_lep ] )->getEnergy() << std::endl;
			streamlog_out(DEBUG4) << "	trueIsoLep (Px,Py,Pz,E):		" << trueIsoLepFourMomentum.Px() << " , " << trueIsoLepFourMomentum.Py() << " , " << trueIsoLepFourMomentum.Pz() << " , " << trueIsoLepFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "	seenIsoLep (Px,Py,Pz,E):		" << seenIsoLepFourMomentum.Px() << " , " << seenIsoLepFourMomentum.Py() << " , " << seenIsoLepFourMomentum.Pz() << " , " << seenIsoLepFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "	recoIsoLep (Px,Py,Pz,E):		" << recoIsoLepFourMomentum.Px() << " , " << recoIsoLepFourMomentum.Py() << " , " << recoIsoLepFourMomentum.Pz() << " , " << recoIsoLepFourMomentum.E() << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG3) << "	recoIsoLep[ " << i_lep << " ]:" << std::endl;
			streamlog_out(DEBUG3) << *recoIsoLepVector[ i_lep ] << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			m_trueIsoLepType.push_back( type_jet( trueIsoLepVectorIndex[ i_lep ] ) );

			TVector3 isoLepTrueMomentum = trueIsoLepFourMomentum.Vect(); double isoLepTrueEnergy = trueIsoLepFourMomentum.E();
			TVector3 isoLepSeenMomentum = seenIsoLepFourMomentum.Vect(); double isoLepSeenEnergy = seenIsoLepFourMomentum.E();
			TVector3 isoLepRecoMomentum = recoIsoLepFourMomentum.Vect(); double isoLepRecoEnergy = recoIsoLepFourMomentum.E();
			double isoLepPxResidual , isoLepPyResidual , isoLepPzResidual , isoLepEnergyResidual , isoLepThetaResidual , isoLepPhiResidual;
			double isoLepPxResidualSeen , isoLepPyResidualSeen , isoLepPzResidualSeen , isoLepEnergyResidualSeen , isoLepThetaResidualSeen , isoLepPhiResidualSeen;
			double isoLepSigmaE , isoLepSigmaTheta , isoLepSigmaPhi;
			getResiduals( isoLepTrueMomentum , isoLepTrueEnergy , isoLepRecoMomentum , isoLepRecoEnergy , isoLepPxResidual , isoLepPyResidual , isoLepPzResidual , isoLepEnergyResidual , isoLepThetaResidual , isoLepPhiResidual );
			getResiduals( isoLepSeenMomentum , isoLepSeenEnergy , isoLepRecoMomentum , isoLepRecoEnergy , isoLepPxResidualSeen , isoLepPyResidualSeen , isoLepPzResidualSeen , isoLepEnergyResidualSeen , isoLepThetaResidualSeen , isoLepPhiResidualSeen );
			getResolutions( TLorentzVector( isoLepRecoMomentum , isoLepRecoEnergy ) , recoIsoLepVector[ i_lep ]->getCovMatrix() , isoLepSigmaE , isoLepSigmaTheta , isoLepSigmaPhi );
			if ( m_fillRootTree )
			{
				m_trueIsoLepPx.push_back( isoLepTrueMomentum.Px() );		m_trueIsoLepPy.push_back( isoLepTrueMomentum.Py() );			m_trueIsoLepPz.push_back( isoLepTrueMomentum.Pz() );
				m_trueIsoLepE.push_back( isoLepTrueEnergy );			m_trueIsoLepTheta.push_back( isoLepTrueMomentum.Theta() );		m_trueIsoLepPhi.push_back( isoLepTrueMomentum.Phi() );
				m_seenIsoLepPx.push_back( isoLepSeenMomentum.Px() );		m_seenIsoLepPy.push_back( isoLepSeenMomentum.Py() );			m_seenIsoLepPz.push_back( isoLepSeenMomentum.Pz() );
				m_seenIsoLepE.push_back( isoLepSeenEnergy );			m_seenIsoLepTheta.push_back( isoLepSeenMomentum.Theta() );		m_seenIsoLepPhi.push_back( isoLepSeenMomentum.Phi() );
				m_recoIsoLepPx.push_back( isoLepRecoMomentum.Px() );		m_recoIsoLepPy.push_back( isoLepRecoMomentum.Py() );			m_recoIsoLepPz.push_back( isoLepRecoMomentum.Pz() );
				m_recoIsoLepE.push_back( isoLepRecoEnergy );			m_recoIsoLepTheta.push_back( isoLepRecoMomentum.Theta() );		m_recoIsoLepPhi.push_back( isoLepRecoMomentum.Phi() );
				m_isoLepSigmaPx2.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 0 ] );
				m_isoLepSigmaPxPy.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 1 ] );
				m_isoLepSigmaPy2.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 2 ] );
				m_isoLepSigmaPxPz.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 3 ] );
				m_isoLepSigmaPyPz.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 4 ] );
				m_isoLepSigmaPz2.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 5 ] );
				m_isoLepSigmaPxE.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 6 ] );
				m_isoLepSigmaPyE.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 7 ] );
				m_isoLepSigmaPzE.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 8 ] );
				m_isoLepSigmaE2.push_back( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 9 ] );
				m_isoLepSigmaTheta2.push_back( pow( isoLepSigmaTheta , 2 ) );
				m_isoLepSigmaPhi2.push_back( pow( isoLepSigmaPhi , 2 ) );
				m_ResidualIsoLepPx.push_back( isoLepPxResidual ); h_ResidualIsoLepPx->Fill( isoLepPxResidual );
				m_ResidualIsoLepPy.push_back( isoLepPyResidual ); h_ResidualIsoLepPy->Fill( isoLepPyResidual );
				m_ResidualIsoLepPz.push_back( isoLepPzResidual ); h_ResidualIsoLepPz->Fill( isoLepPzResidual );
				m_ResidualIsoLepE.push_back( isoLepEnergyResidual ); h_ResidualIsoLepE->Fill( isoLepEnergyResidual );
				m_ResidualIsoLepTheta.push_back( isoLepThetaResidual ); h_ResidualIsoLepTheta->Fill( isoLepThetaResidual );
				m_ResidualIsoLepPhi.push_back( isoLepPhiResidual ); h_ResidualIsoLepPhi->Fill( isoLepPhiResidual );
				m_NormalizedResidualIsoLepPx.push_back( isoLepPxResidual / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 0 ] ) );
				m_NormalizedResidualIsoLepPy.push_back( isoLepPyResidual / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 2 ] ) );
				m_NormalizedResidualIsoLepPz.push_back( isoLepPzResidual / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 5 ] ) );
				m_NormalizedResidualIsoLepE.push_back( isoLepEnergyResidual / isoLepSigmaE );
				m_NormalizedResidualIsoLepTheta.push_back( isoLepThetaResidual / isoLepSigmaTheta );
				m_NormalizedResidualIsoLepPhi.push_back( isoLepPhiResidual / isoLepSigmaPhi );
				h_NormalizedResidualIsoLepPx->Fill( isoLepPxResidual / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 0 ] ) );
				h_NormalizedResidualIsoLepPy->Fill( isoLepPyResidual / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 2 ] ) );
				h_NormalizedResidualIsoLepPz->Fill( isoLepPzResidual / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 5 ] ) );
				h_NormalizedResidualIsoLepE->Fill( isoLepEnergyResidual / isoLepSigmaE );
				h_NormalizedResidualIsoLepTheta->Fill( isoLepThetaResidual / isoLepSigmaTheta );
				h_NormalizedResidualIsoLepPhi->Fill( isoLepPhiResidual / isoLepSigmaPhi );
				m_ResidualIsoLepPxSeen.push_back( isoLepPxResidualSeen ); h_ResidualIsoLepPxSeen->Fill( isoLepPxResidualSeen );
				m_ResidualIsoLepPySeen.push_back( isoLepPyResidualSeen ); h_ResidualIsoLepPySeen->Fill( isoLepPyResidualSeen );
				m_ResidualIsoLepPzSeen.push_back( isoLepPzResidualSeen ); h_ResidualIsoLepPzSeen->Fill( isoLepPzResidualSeen );
				m_ResidualIsoLepESeen.push_back( isoLepEnergyResidualSeen ); h_ResidualIsoLepESeen->Fill( isoLepEnergyResidualSeen );
				m_ResidualIsoLepThetaSeen.push_back( isoLepThetaResidualSeen ); h_ResidualIsoLepThetaSeen->Fill( isoLepThetaResidualSeen );
				m_ResidualIsoLepPhiSeen.push_back( isoLepPhiResidualSeen ); h_ResidualIsoLepPhiSeen->Fill( isoLepPhiResidualSeen );
				m_NormalizedResidualIsoLepPxSeen.push_back( isoLepPxResidualSeen / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 0 ] ) );
				m_NormalizedResidualIsoLepPySeen.push_back( isoLepPyResidualSeen / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 2 ] ) );
				m_NormalizedResidualIsoLepPzSeen.push_back( isoLepPzResidualSeen / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 5 ] ) );
				m_NormalizedResidualIsoLepESeen.push_back( isoLepEnergyResidualSeen / isoLepSigmaE );
				m_NormalizedResidualIsoLepThetaSeen.push_back( isoLepThetaResidualSeen / isoLepSigmaTheta );
				m_NormalizedResidualIsoLepPhiSeen.push_back( isoLepPhiResidualSeen / isoLepSigmaPhi );
				h_NormalizedResidualIsoLepPxSeen->Fill( isoLepPxResidualSeen / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 0 ] ) );
				h_NormalizedResidualIsoLepPySeen->Fill( isoLepPyResidualSeen / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 2 ] ) );
				h_NormalizedResidualIsoLepPzSeen->Fill( isoLepPzResidualSeen / std::sqrt( recoIsoLepVector[ i_lep ]->getCovMatrix()[ 5 ] ) );
				h_NormalizedResidualIsoLepESeen->Fill( isoLepEnergyResidualSeen / isoLepSigmaE );
				h_NormalizedResidualIsoLepThetaSeen->Fill( isoLepThetaResidualSeen / isoLepSigmaTheta );
				h_NormalizedResidualIsoLepPhiSeen->Fill( isoLepPhiResidualSeen / isoLepSigmaPhi );
			}
		}

		m_nEvtSum++;
		m_nEvt++ ;
		if( m_fillRootTree ) m_pTTree->Fill();
		delall();
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "	Check : Input collections not found in event " << m_nEvt << std::endl;
	}


}

void JetErrorAnalysis::getResiduals( TVector3 jetTrueMomentum , double jetTrueEnergy , TVector3 jetRecoMomentum , double jetRecoEnergy , double &jetPxResidual , double &jetPyResidual , double &jetPzResidual , double &jetEnergyResidual , double &jetThetaResidual , double &jetPhiResidual )
{
	jetPxResidual = jetRecoMomentum.Px() - jetTrueMomentum.Px();
	jetPyResidual = jetRecoMomentum.Py() - jetTrueMomentum.Py();
	jetPzResidual = jetRecoMomentum.Pz() - jetTrueMomentum.Pz();
	jetEnergyResidual = jetRecoEnergy - jetTrueEnergy;
	TVector3 jetTrueMomentumUnit = jetTrueMomentum; jetTrueMomentumUnit.SetMag( 1.0 );
	TVector3 jetTruePtUnit( jetTrueMomentum.Px() , jetTrueMomentum.Py() , 0.0 ); jetTruePtUnit.SetMag( 1.0 );
	TVector3 jetRecoMomentumUnitRotated = jetRecoMomentum; jetRecoMomentumUnitRotated.SetPhi( jetTrueMomentum.Phi() ); jetRecoMomentumUnitRotated.SetMag( 1.0 );
	TVector3 jetRecoPtUnit( jetRecoMomentum.Px() , jetRecoMomentum.Py() , 0.0 ); jetRecoPtUnit.SetMag( 1.0 );
	jetThetaResidual = ( jetRecoMomentum.Theta() >= jetTrueMomentum.Theta() ? acos( jetTrueMomentumUnit.Dot( jetRecoMomentumUnitRotated ) ) : -1.0 * acos( jetTrueMomentumUnit.Dot( jetRecoMomentumUnitRotated ) ) );
	jetPhiResidual = ( jetRecoMomentum.Phi() >= jetTrueMomentum.Phi() ? acos( jetTruePtUnit.Dot( jetRecoPtUnit ) ) : -1.0 * acos( jetTruePtUnit.Dot( jetRecoPtUnit ) ) );
}

void JetErrorAnalysis::getResolutions(	TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , double &sigmaE , double &sigmaTheta , double &sigmaPhi )
{
	float Px , Py , Pz , P2 , Pt , Pt2;
	float dTheta_dPx , dTheta_dPy , dTheta_dPz , dPhi_dPx , dPhi_dPy;
	float sigmaPx2 , sigmaPy2 , sigmaPz2 , sigmaPxPy , sigmaPxPz , sigmaPyPz;

	Px		= jetFourMomentum.Px();
	Py		= jetFourMomentum.Py();
	Pz		= jetFourMomentum.Pz();
	P2		= ( jetFourMomentum.Vect() ).Mag2();
	Pt2		= pow( Px , 2 ) + pow( Py , 2 );
	Pt		= sqrt( Pt2 );
	sigmaPx2	= jetCovMat[ 0 ];
	sigmaPxPy	= jetCovMat[ 1 ];
	sigmaPy2	= jetCovMat[ 2 ];
	sigmaPxPz	= jetCovMat[ 3 ];
	sigmaPyPz	= jetCovMat[ 4 ];
	sigmaPz2	= jetCovMat[ 5 ];

	dTheta_dPx	= Px * Pz / ( P2 * Pt );
	dTheta_dPy	= Py * Pz / ( P2 * Pt );
	dTheta_dPz	= -Pt / P2;
	dPhi_dPx	= -Py / Pt2;
	dPhi_dPy	= Px / Pt2;

	sigmaE		= std::sqrt( jetCovMat[ 9 ] );
	sigmaTheta	= std::sqrt( std::fabs( std::pow( dTheta_dPx , 2 ) * sigmaPx2 + std::pow( dTheta_dPy , 2 ) * sigmaPy2 + std::pow( dTheta_dPz , 2 ) * sigmaPz2 +
	 					2.0 * dTheta_dPx * dTheta_dPy * sigmaPxPy + 2.0 * dTheta_dPx * dTheta_dPz * sigmaPxPz + 2.0 * dTheta_dPy * dTheta_dPz * sigmaPyPz ) );
	sigmaPhi	= std::sqrt( std::fabs( std::pow( dPhi_dPx , 2 ) * sigmaPx2 + std::pow( dPhi_dPy , 2 ) * sigmaPy2 + 2.0 * dPhi_dPx * dPhi_dPy * sigmaPxPy ) );
	streamlog_out(DEBUG1) << "			E		= " << jetFourMomentum.E() << std::endl ;
	streamlog_out(DEBUG1) << "			Theta		= " << jetFourMomentum.Theta() << std::endl ;
	streamlog_out(DEBUG1) << "			Phi		= " << jetFourMomentum.Phi() << std::endl ;
	streamlog_out(DEBUG1) << "			sigmaE		= " << sigmaE << std::endl ;
	streamlog_out(DEBUG1) << "			sigmaTheta	= " << sigmaTheta << std::endl ;
	streamlog_out(DEBUG1) << "			sigmaPhi	= " << sigmaPhi << std::endl ;

}

EVENT::ReconstructedParticle* JetErrorAnalysis::getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO )
{
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "Look for PFO linked to visible MCParticle:" << std::endl;

	ReconstructedParticle* linkedPFO{};
	bool foundlinkedPFO = false;
	const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( mcParticle );
	const EVENT::FloatVec&	PFOweightvec = MCParticleRecoNav.getRelatedToWeights( mcParticle );
	streamlog_out(DEBUG0) << "Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
	weightPFOtoMCP = 0.0;
	weightMCPtoPFO = 0.0;
	double maxweightPFOtoMCP = 0.;
	double maxweightMCPtoPFO = 0.;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
	{
		double pfo_weight = 0.0;
		double trackWeight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
		double clusterWeight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
		if ( getChargedPFO && !getNeutralPFO )
		{
			pfo_weight = trackWeight;
		}
		else if ( getNeutralPFO && !getChargedPFO )
		{
			pfo_weight = clusterWeight;
		}
		else
		{
			pfo_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
		}
		streamlog_out(DEBUG0) << "Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << trackWeight << " , Cluster: " << clusterWeight << ")" << std::endl;
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
		if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightMCPtoPFO = pfo_weight;
			iMCPtoPFOmax = i_pfo;
			streamlog_out(DEBUG0) << "PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
		}
	}
	if ( getChargedPFO && maxweightMCPtoPFO < 0.8 )
	{
		streamlog_out(DEBUG1) << "MCParticle has link weight lower than 0.8 ( " << maxweightMCPtoPFO << " ), looking for linked PFO in clusters" << std::endl;
		for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
		{
			double pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
			streamlog_out(DEBUG0) << "Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0 << " , Cluster: " << ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0 << ")" << std::endl;
			ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
			if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoPFO = pfo_weight;
				iMCPtoPFOmax = i_pfo;
				streamlog_out(DEBUG0) << "PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
			}
		}
	}
	if ( iMCPtoPFOmax != -1 )
	{
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( iMCPtoPFOmax );
		const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( testPFO );
		const EVENT::FloatVec&	MCPweightvec = RecoMCParticleNav.getRelatedToWeights( testPFO );
		for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
		{
			double mcp_weight = 0.0;
			double trackWeight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
			double clusterWeight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
			if ( getChargedPFO && !getNeutralPFO )
			{
				mcp_weight = trackWeight;
			}
			else if ( getNeutralPFO && !getChargedPFO )
			{
				mcp_weight = clusterWeight;
			}
			else
			{
				mcp_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
			}
			MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
			if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightPFOtoMCP = mcp_weight;
				iPFOtoMCPmax = i_mcp;
				streamlog_out(DEBUG0) << "MCParticle at index: " << testMCP->id() << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
			}
		}
		if ( iPFOtoMCPmax != -1 )
		{
			if ( MCPvec.at( iPFOtoMCPmax ) == mcParticle )
			{
				linkedPFO = testPFO;
				foundlinkedPFO = true;
			}
		}
	}
	if( foundlinkedPFO )
	{
		streamlog_out(DEBUG1) << "Linked PFO to MCParticle found successfully " << std::endl;
		weightPFOtoMCP = maxweightPFOtoMCP;
		weightMCPtoPFO = maxweightMCPtoPFO;
		return linkedPFO;
	}
	else
	{
		streamlog_out(DEBUG1) << "Couldn't Find a PFO linked to MCParticle" << std::endl;
		return NULL;
	}
}


void JetErrorAnalysis::check()
{

}

void JetErrorAnalysis::InitializeHistogram( TH1F *histogram , bool scale , int color , int lineWidth , int markerSize , int markerStyle , bool fitGaussian )
{
	/*
	TPad *pad = new TPad("pad1", "pad1", 0.0 , 0.0 , 1.0 , 1.0 , color , 0 , 0 );
	pad->SetBottomMargin(0.0175);
	pad->SetRightMargin(0.02);
	pad->SetLeftMargin(0.15);
	pad->SetTopMargin(0.015);
	pad->Draw();
	pad->cd();
	*/
	gStyle->SetStatColor( color );
	if ( scale ) histogram->Scale( 1.0 / ( histogram->GetEntries() ) );
	histogram->SetLineColor( color );
	histogram->SetLineWidth( lineWidth );
	histogram->SetMarkerSize( markerSize );
	histogram->SetMarkerStyle( markerStyle );
	histogram->SetMarkerColor( color );
	float fit_range = 2.0;
	float fit_min = -2.0;
	float fit_max = 2.0;
	if ( fitGaussian )
	{
		doProperGaussianFit( histogram , fit_min , fit_max , fit_range );
		histogram->GetFunction("gaus")->SetLineColor( color );
	}
	float y_max = 1.2 * histogram->GetMaximum();
	histogram->GetXaxis()->SetLabelSize(0.06);
	histogram->GetXaxis()->SetLabelOffset(0.005);
	histogram->GetXaxis()->SetTitleSize(0.07);
	histogram->GetXaxis()->SetTitleOffset(1.10);
	histogram->GetYaxis()->SetRangeUser(0.0, y_max);
	histogram->GetYaxis()->SetLabelSize(0.05);
	histogram->GetYaxis()->SetLabelOffset(0.005);
	histogram->GetYaxis()->SetTitleSize(0.07);
	histogram->GetYaxis()->SetTitleOffset(1.24);
	histogram->Write();

	TPaveText *ildFullsim = new TPaveText( 0.2 , 0.5 , 0.2 , 0.5 , "tbNDC" );
	ildFullsim->SetTextSize(0.05);
	ildFullsim->SetBorderSize(1);
	ildFullsim->SetTextAlign(13);
	ildFullsim->SetFillColor(0);
	ildFullsim->SetTextColor(1);
	ildFullsim->AddText("ILD full simulation");
	ildFullsim->Draw("same");

	gPad->Update();
	//TPaveStats *tps = (TPaveStats *)histogram->GetListOfFunctions()->FindObject("stats");
	//TPaveStats *tps = (TPaveStats *)histogram->FindObject("stats");
	//tps->SetTextColor( color );
	//tps->SetLineColor( color );
	//tps->SetX1NDC( 0.65 );
	//tps->SetX2NDC( 0.95 );
	//tps->SetY1NDC( 0.5 );
	//tps->SetY2NDC( 0.95 );

	gPad->SetRightMargin(0.027);
	gPad->SetLeftMargin(0.17);
	gPad->SetTopMargin(0.027);
	gPad->SetBottomMargin(0.15);
	gPad->Update();

}

void JetErrorAnalysis::doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange )
{
	float Chi2 = 0.0;
	float NDF = 0.0;
	for ( int i_fit = 0 ; i_fit < 3 ; ++i_fit )
	{
		histogram->Fit( "gaus" , "" , "" , fitMin , fitMax );
		TF1 *fitFunction = (TF1 *)histogram->GetFunction("gaus");
		double fitMean = fitFunction->GetParameter( 1 );
		double fitSigma = fitFunction->GetParameter( 2 );
		fitMin = fitMean - fitRange * fitSigma;
		fitMax = fitMean + fitRange * fitSigma;
		Chi2 = fitFunction->GetChisquare();
		NDF = fitFunction->GetNDF();
	}
	streamlog_out(DEBUG4) << "	FIT : CHI2(" << Chi2 << ") / NDF(" << NDF << ") = " << Chi2 / NDF << " 	, fitrange = " << fitRange << std::endl;
	streamlog_out(DEBUG4) << "" << std::endl;
	if ( Chi2 != 0.0 && NDF != 0.0 && Chi2 / NDF > 2.0 && fitRange >= 0.5 )
	{
		doProperGaussianFit( histogram , fitMin , fitMax , fitRange - 0.1 );
	}
}

void JetErrorAnalysis::end()
{
	if( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree->Write();
		h_ResidualJetPx->SetName( m_histName.c_str() );
		h_ResidualJetPy->SetName( m_histName.c_str() );
		h_ResidualJetPz->SetName( m_histName.c_str() );
		h_ResidualJetE->SetName( m_histName.c_str() );
		h_ResidualJetTheta->SetName( m_histName.c_str() );
		h_ResidualJetPhi->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPx->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPy->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPz->SetName( m_histName.c_str() );
		h_NormalizedResidualJetE->SetName( m_histName.c_str() );
		h_NormalizedResidualJetTheta->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPhi->SetName( m_histName.c_str() );
		m_histName = ( m_histName + "(without invisibles)" ).c_str();
		h_ResidualJetPxSeen->SetName( m_histName.c_str() );
		h_ResidualJetPySeen->SetName( m_histName.c_str() );
		h_ResidualJetPzSeen->SetName( m_histName.c_str() );
		h_ResidualJetESeen->SetName( m_histName.c_str() );
		h_ResidualJetThetaSeen->SetName( m_histName.c_str() );
		h_ResidualJetPhiSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPxSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPySeen->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPzSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualJetESeen->SetName( m_histName.c_str() );
		h_NormalizedResidualJetThetaSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualJetPhiSeen->SetName( m_histName.c_str() );
		InitializeHistogram( h_ResidualJetPx , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetPy , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetPz , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetE , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetTheta , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetPhi , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualJetPx , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetPy , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetPz , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetE , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetTheta , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetPhi , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_ResidualJetPxSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetPySeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetPzSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetESeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetThetaSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualJetPhiSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualJetPxSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetPySeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetPzSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetESeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetThetaSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualJetPhiSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		m_histName = "IsolatedLepton(Reco-True)";
		h_ResidualIsoLepPx->SetName( m_histName.c_str() );
		h_ResidualIsoLepPy->SetName( m_histName.c_str() );
		h_ResidualIsoLepPz->SetName( m_histName.c_str() );
		h_ResidualIsoLepE->SetName( m_histName.c_str() );
		h_ResidualIsoLepTheta->SetName( m_histName.c_str() );
		h_ResidualIsoLepPhi->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPx->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPy->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPz->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepE->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepTheta->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPhi->SetName( m_histName.c_str() );
		m_histName = "IsolatedLepton(Reco-Seen)";
		h_ResidualIsoLepPxSeen->SetName( m_histName.c_str() );
		h_ResidualIsoLepPySeen->SetName( m_histName.c_str() );
		h_ResidualIsoLepPzSeen->SetName( m_histName.c_str() );
		h_ResidualIsoLepESeen->SetName( m_histName.c_str() );
		h_ResidualIsoLepThetaSeen->SetName( m_histName.c_str() );
		h_ResidualIsoLepPhiSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPxSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPySeen->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPzSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepESeen->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepThetaSeen->SetName( m_histName.c_str() );
		h_NormalizedResidualIsoLepPhiSeen->SetName( m_histName.c_str() );
		InitializeHistogram( h_ResidualIsoLepPx , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepPy , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepPz , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepE , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepTheta , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepPhi , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualIsoLepPx , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualIsoLepPy , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualIsoLepPz , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualIsoLepE , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualIsoLepTheta , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_NormalizedResidualIsoLepPhi , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , true );
		InitializeHistogram( h_ResidualIsoLepPxSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepPySeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepPzSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepESeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepThetaSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_ResidualIsoLepPhiSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualIsoLepPxSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualIsoLepPySeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualIsoLepPzSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualIsoLepESeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualIsoLepThetaSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		InitializeHistogram( h_NormalizedResidualIsoLepPhiSeen , m_normalizeHistograms , m_histColour , 1 , 1.0 , 1 , false );
		m_pTFile->Close();
		delete m_pTFile;
	}


}
