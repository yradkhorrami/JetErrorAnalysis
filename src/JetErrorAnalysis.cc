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
m_nRecoJets(0),
m_pTFile(NULL),
m_pTTree(NULL)
{

	// modify processor description
	_description = "JetErrorAnalysis does whatever it does ..." ;


	// register steering parameters: name, description, class-variable, default value


	// Inputs: MC-particles, Reco-particles, the link between the two

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"RecoJetCollection" ,
					"Name of the input Reconstructed Jet collection"  ,
					m_recoJetCollectionName ,
					std::string("Durham_nJets")
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

	// Inputs: True jets (as a recoparticle, will be the sum of the _reconstructed particles_
	// created by the true particles in each true jet, in the RecoMCTruthLink sense.
	// link jet-to-reco particles, link jet-to-MC-particles.

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

	m_pTFile = new TFile(m_outputFile.c_str(),"recreate");
	m_pTTree = new TTree("eventTree","eventTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("nTrueJets",&m_nTrueJets,"nTrueJets/I") ;
	m_pTTree->Branch("nRecoJets",&m_nRecoJets,"nRecoJets/I") ;
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
	m_pTTree->Branch("ResidualPx",&m_ResidualPx) ;
	m_pTTree->Branch("ResidualPy",&m_ResidualPy) ;
	m_pTTree->Branch("ResidualPz",&m_ResidualPz) ;
	m_pTTree->Branch("ResidualE",&m_ResidualE) ;
	m_pTTree->Branch("ResidualTheta",&m_ResidualTheta) ;
	m_pTTree->Branch("ResidualPhi",&m_ResidualPhi) ;
	m_pTTree->Branch("NormalizedResidualPx",&m_NormalizedResidualPx) ;
	m_pTTree->Branch("NormalizedResidualPy",&m_NormalizedResidualPy) ;
	m_pTTree->Branch("NormalizedResidualPz",&m_NormalizedResidualPz) ;
	m_pTTree->Branch("NormalizedResidualE",&m_NormalizedResidualE) ;
	m_pTTree->Branch("NormalizedResidualTheta",&m_NormalizedResidualTheta) ;
	m_pTTree->Branch("NormalizedResidualPhi",&m_NormalizedResidualPhi) ;
	m_pTTree->Branch("trueJetType",&m_trueJetType) ;
	m_pTTree->Branch("trueJetFlavour",&m_trueJetFlavour) ;
	m_pTTree->Branch("recoJetFlavour",&m_recoJetFlavour) ;
	h_ResidualPx = new TH1F( m_histName.c_str() , "; ^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualPy = new TH1F( m_histName.c_str() , "; ^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualPz = new TH1F( m_histName.c_str() , "; ^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualE = new TH1F( m_histName.c_str() , "; ^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualTheta = new TH1F( m_histName.c_str() , "; ^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE} [rad]; Normalized #jet / 0.01" , 200 * 3.14159265 , -3.14159265 , 3.14159265 );
	h_ResidualPhi = new TH1F( m_histName.c_str() , "; ^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE} [rad]; Normalized #jet / 0.01" , 200 * 3.14159265 , -3.14159265 , 3.14159265 );
	h_NormalizedResidualPx = new TH1F( m_histName.c_str() , "; (^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE}) / #sigma_{p_{x,jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualPy = new TH1F( m_histName.c_str() , "; (^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE}) / #sigma_{p_{y,jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualPz = new TH1F( m_histName.c_str() , "; (^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE}) / #sigma_{p_{z,jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualE = new TH1F( m_histName.c_str() , "; (^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE}) / #sigma_{E_{jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualTheta = new TH1F( m_histName.c_str() , "; (^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE}) / #sigma_{#theta_{jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualPhi = new TH1F( m_histName.c_str() , "; (^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE}) / #sigma_{#phi_{jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	m_histName = ( m_histName + "(without invisibles)" ).c_str();
	h_ResidualPxSeen = new TH1F( m_histName.c_str() , "; ^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualPySeen = new TH1F( m_histName.c_str() , "; ^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualPzSeen = new TH1F( m_histName.c_str() , "; ^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualESeen = new TH1F( m_histName.c_str() , "; ^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE} [GeV]; Normalized #jet / 0.1 GeV" , 400 , -20.0 , 20.0 );
	h_ResidualThetaSeen = new TH1F( m_histName.c_str() , "; ^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE} [rad]; Normalized #jet / 0.01" , 200 * 3.14159265 , -3.14159265 , 3.14159265 );
	h_ResidualPhiSeen = new TH1F( m_histName.c_str() , "; ^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE} [rad]; Normalized #jet / 0.01" , 200 * 3.14159265 , -3.14159265 , 3.14159265 );
	h_NormalizedResidualPxSeen = new TH1F( m_histName.c_str() , "; (^{}_{}p_{x,jet}^{REC} - p_{x,jet}^{TRUE}) / #sigma_{p_{x,jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualPySeen = new TH1F( m_histName.c_str() , "; (^{}_{}p_{y,jet}^{REC} - p_{y,jet}^{TRUE}) / #sigma_{p_{y,jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualPzSeen = new TH1F( m_histName.c_str() , "; (^{}_{}p_{z,jet}^{REC} - p_{z,jet}^{TRUE}) / #sigma_{p_{z,jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualESeen = new TH1F( m_histName.c_str() , "; (^{}_{}E_{jet}^{REC} - E_{jet}^{TRUE}) / #sigma_{E_{jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualThetaSeen = new TH1F( m_histName.c_str() , "; (^{}_{}#theta_{jet}^{REC} - #theta_{jet}^{TRUE}) / #sigma_{#theta_{jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );
	h_NormalizedResidualPhiSeen = new TH1F( m_histName.c_str() , "; (^{}_{}#phi_{jet}^{REC} - #phi_{jet}^{TRUE}) / #sigma_{#phi_{jet}}; Normalized #jet / 0.1" , 200 , -10.0 , 10.0 );

}

void JetErrorAnalysis::Clear()
{
	m_nTrueJets = 0;
	m_nRecoJets = 0;
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
	m_ResidualPx.clear();
	m_ResidualPy.clear();
	m_ResidualPz.clear();
	m_ResidualE.clear();
	m_ResidualTheta.clear();
	m_ResidualPhi.clear();
	m_NormalizedResidualPx.clear();
	m_NormalizedResidualPy.clear();
	m_NormalizedResidualPz.clear();
	m_NormalizedResidualE.clear();
	m_NormalizedResidualTheta.clear();
	m_NormalizedResidualPhi.clear();
	m_trueJetType.clear();
	m_trueJetFlavour.clear();
	m_recoJetFlavour.clear();
}

void JetErrorAnalysis::processRunHeader()
{
    m_nRun++ ;
}


void JetErrorAnalysis::processEvent( LCEvent* pLCEvent)
{
	LCCollection *recoJetCol{};
//	LCCollection *trueJetCol{};
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
		recoJetCol		= pLCEvent->getCollection( m_recoJetCollectionName );
//		trueJetCol		= pLCEvent->getCollection( _trueJetCollectionName );
		TrueJet_Parser* trueJet	= this;
		trueJet->getall( pLCEvent );

		m_nRecoJets = recoJetCol->getNumberOfElements();
		streamlog_out(DEBUG3) << "	Number of Reconstructed Jets: " << m_nRecoJets << std::endl;

		int njets = trueJet->njets();
		streamlog_out(DEBUG3) << "	Number of True Jets: " << njets << std::endl;

		if ( m_nRecoJets == 0 || njets == 0 ) return;

		std::vector<ReconstructedParticle*> recoJetVector{};
		std::vector<ReconstructedParticle*> leadingParticles{};
		std::vector<int> trueJetVectorIndex{};

		for ( int i_recoJet = 0 ; i_recoJet < m_nRecoJets ; ++i_recoJet )
		{
			ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( i_recoJet ) );
			ReconstructedParticle* leadingParticle = NULL;
			float leadingEnergy = 0.0;
			streamlog_out(DEBUG2) << " looking for leading particle in recoJet " << i_recoJet << " with " << ( recoJet->getParticles() ).size() << " particles" << std::endl;
			for ( unsigned int i_par = 0 ; i_par < ( recoJet->getParticles() ).size() ; ++i_par )
			{
				ReconstructedParticle* pfo = ( ReconstructedParticle* )recoJet->getParticles()[ i_par ];
				streamlog_out(DEBUG0) << " checking particle " << i_par << " with Energy = " << pfo->getEnergy() << std::endl;
				if ( abs( pfo->getType() ) == 12 || abs( pfo->getType() ) == 14 || abs( pfo->getType() ) == 16 ) continue;
				if ( pfo->getEnergy() > leadingEnergy )
				{
					leadingParticle = pfo;
					leadingEnergy = pfo->getEnergy();
					streamlog_out(DEBUG0) << "leading particle so far: " << std::endl;
					streamlog_out(DEBUG0) << *leadingParticle << std::endl;
				}
			}
			if ( leadingParticle != NULL )
			{
				int trueJetIndex;
				LCObjectVec jetvec = reltjreco->getRelatedFromObjects( leadingParticle );
				streamlog_out(DEBUG2) << jetvec.size() << " true Jet found for leading particle of jet " << i_recoJet << std::endl;
				if ( jetvec.size() != 0 )
				{
					trueJetIndex = jetindex( dynamic_cast<ReconstructedParticle*>( jetvec[ 0 ] ) );
					recoJetVector.push_back( recoJet );
					leadingParticles.push_back( leadingParticle );
					trueJetVectorIndex.push_back( trueJetIndex );
				}
			}
		}
		streamlog_out(DEBUG2) << "	" << recoJetVector.size() << " recoJets and " << trueJetVectorIndex.size() << " true Jets are find and matched" << std::endl;
		for ( unsigned int i_jet = 0 ; i_jet < recoJetVector.size() ; ++i_jet )
		{
			streamlog_out(DEBUG2) << "----------------------------------------------------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG2) << "Reconstructed Jet[ " << i_jet << " ] matches true jet [ " << trueJetVectorIndex[ i_jet ] << " ] by finding leading particle" << std::endl;
			streamlog_out(DEBUG2) << "	trueJet TYPE:	" << type_jet( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG2) << "	trueJet (Px,Py,Pz,E):	" << ptrue( trueJetVectorIndex[ i_jet ] )[ 0 ] << " , " << ptrue( trueJetVectorIndex[ i_jet ] )[ 1 ] << " , " << ptrue( trueJetVectorIndex[ i_jet ] )[ 2 ] << " , " << Etrue( trueJetVectorIndex[ i_jet ] ) << std::endl;
			streamlog_out(DEBUG2) << "	recoJet:" << std::endl;
			streamlog_out(DEBUG2) << *recoJetVector[ i_jet ] << std::endl;
			streamlog_out(DEBUG2) << "	leadingParticle:" << std::endl;
			streamlog_out(DEBUG2) << *leadingParticles[ i_jet ] << std::endl;
			streamlog_out(DEBUG2) << "" << std::endl;
			TVector3 quarkMomentum( pquark( trueJetVectorIndex[ i_jet ] )[ 0 ] , pquark( trueJetVectorIndex[ i_jet ] )[ 1 ] , pquark( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			double quarkEnergy = Equark( trueJetVectorIndex[ i_jet ] );
			TVector3 jetTrueMomentum( ptrue( trueJetVectorIndex[ i_jet ] )[ 0 ] , ptrue( trueJetVectorIndex[ i_jet ] )[ 1 ] , ptrue( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			double jetTrueEnergy = Etrue( trueJetVectorIndex[ i_jet ] );
			TVector3 jetTrueSeenMomentum( ptrueseen( trueJetVectorIndex[ i_jet ] )[ 0 ] , ptrueseen( trueJetVectorIndex[ i_jet ] )[ 1 ] , ptrueseen( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			double jetTrueSeenEnergy = Etrueseen( trueJetVectorIndex[ i_jet ] );
			TVector3 jetSeenMomentum( pseen( trueJetVectorIndex[ i_jet ] )[ 0 ] , pseen( trueJetVectorIndex[ i_jet ] )[ 1 ] , pseen( trueJetVectorIndex[ i_jet ] )[ 2 ] );
			double jetSeenEnergy = Eseen( trueJetVectorIndex[ i_jet ] );
			TVector3 jetRecoMomentum( recoJetVector[ i_jet ]->getMomentum() );
			double jetRecoEnergy = recoJetVector[ i_jet ]->getEnergy();
			double jetPxResidual , jetPyResidual , jetPzResidual , jetEnergyResidual , jetThetaResidual , jetPhiResidual;
			double jetPxResidualSeen , jetPyResidualSeen , jetPzResidualSeen , jetEnergyResidualSeen , jetThetaResidualSeen , jetPhiResidualSeen;
			double jetSigmaE , jetSigmaTheta , jetSigmaPhi;
			getJetResiduals( jetTrueMomentum , jetTrueEnergy , jetRecoMomentum , jetRecoEnergy , jetPxResidual , jetPyResidual , jetPzResidual , jetEnergyResidual , jetThetaResidual , jetPhiResidual );
			getJetResiduals( jetTrueSeenMomentum , jetTrueSeenEnergy , jetRecoMomentum , jetRecoEnergy , jetPxResidualSeen , jetPyResidualSeen , jetPzResidualSeen , jetEnergyResidualSeen , jetThetaResidualSeen , jetPhiResidualSeen );
			getJetResolutions( TLorentzVector( jetRecoMomentum , jetRecoEnergy ) , recoJetVector[ i_jet ]->getCovMatrix() , jetSigmaE , jetSigmaTheta , jetSigmaPhi );
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
			m_ResidualPx.push_back( jetPxResidual ); h_ResidualPx->Fill( jetPxResidual );
			m_ResidualPy.push_back( jetPyResidual ); h_ResidualPy->Fill( jetPyResidual );
			m_ResidualPz.push_back( jetPzResidual ); h_ResidualPz->Fill( jetPzResidual );
			m_ResidualE.push_back( jetEnergyResidual ); h_ResidualE->Fill( jetEnergyResidual );
			m_ResidualTheta.push_back( jetThetaResidual ); h_ResidualTheta->Fill( jetThetaResidual );
			m_ResidualPhi.push_back( jetPhiResidual ); h_ResidualPhi->Fill( jetPhiResidual );
			m_NormalizedResidualPx.push_back( jetPxResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
			m_NormalizedResidualPy.push_back( jetPyResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
			m_NormalizedResidualPz.push_back( jetPzResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
			m_NormalizedResidualE.push_back( jetEnergyResidual / jetSigmaE );
			m_NormalizedResidualTheta.push_back( jetThetaResidual / jetSigmaTheta );
			m_NormalizedResidualPhi.push_back( jetPhiResidual / jetSigmaPhi );
			h_NormalizedResidualPx->Fill( jetPxResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
			h_NormalizedResidualPy->Fill( jetPyResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
			h_NormalizedResidualPz->Fill( jetPzResidual / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
			h_NormalizedResidualE->Fill( jetEnergyResidual / jetSigmaE );
			h_NormalizedResidualTheta->Fill( jetThetaResidual / jetSigmaTheta );
			h_NormalizedResidualPhi->Fill( jetPhiResidual / jetSigmaPhi );
			m_ResidualPxSeen.push_back( jetPxResidualSeen ); h_ResidualPxSeen->Fill( jetPxResidualSeen );
			m_ResidualPySeen.push_back( jetPyResidualSeen ); h_ResidualPySeen->Fill( jetPyResidualSeen );
			m_ResidualPzSeen.push_back( jetPzResidualSeen ); h_ResidualPzSeen->Fill( jetPzResidualSeen );
			m_ResidualESeen.push_back( jetEnergyResidualSeen ); h_ResidualESeen->Fill( jetEnergyResidualSeen );
			m_ResidualThetaSeen.push_back( jetThetaResidualSeen ); h_ResidualThetaSeen->Fill( jetThetaResidualSeen );
			m_ResidualPhiSeen.push_back( jetPhiResidualSeen ); h_ResidualPhiSeen->Fill( jetPhiResidualSeen );
			m_NormalizedResidualPxSeen.push_back( jetPxResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
			m_NormalizedResidualPySeen.push_back( jetPyResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
			m_NormalizedResidualPzSeen.push_back( jetPzResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
			m_NormalizedResidualESeen.push_back( jetEnergyResidualSeen / jetSigmaE );
			m_NormalizedResidualThetaSeen.push_back( jetThetaResidualSeen / jetSigmaTheta );
			m_NormalizedResidualPhiSeen.push_back( jetPhiResidualSeen / jetSigmaPhi );
			h_NormalizedResidualPxSeen->Fill( jetPxResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 0 ] ) );
			h_NormalizedResidualPySeen->Fill( jetPyResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 2 ] ) );
			h_NormalizedResidualPzSeen->Fill( jetPzResidualSeen / std::sqrt( recoJetVector[ i_jet ]->getCovMatrix()[ 5 ] ) );
			h_NormalizedResidualESeen->Fill( jetEnergyResidualSeen / jetSigmaE );
			h_NormalizedResidualThetaSeen->Fill( jetThetaResidualSeen / jetSigmaTheta );
			h_NormalizedResidualPhiSeen->Fill( jetPhiResidualSeen / jetSigmaPhi );
		}
		m_nEvtSum++;
		m_nEvt++ ;
		m_pTTree->Fill();
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "	Check : Input collections not found in event " << m_nEvt << std::endl;
	}


}

void JetErrorAnalysis::getJetResiduals( TVector3 jetTrueMomentum , double jetTrueEnergy , TVector3 jetRecoMomentum , double jetRecoEnergy , double &jetPxResidual , double &jetPyResidual , double &jetPzResidual , double &jetEnergyResidual , double &jetThetaResidual , double &jetPhiResidual )
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

void JetErrorAnalysis::getJetResolutions(	TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , double &sigmaE , double &sigmaTheta , double &sigmaPhi )
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

void JetErrorAnalysis::check()
{

}

void JetErrorAnalysis::InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle , bool fitGaussian )
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
	histogram->Scale( 1.0 / scale );
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
//	TPaveStats *tps = (TPaveStats *)histogram->GetListOfFunctions()->FindObject("stats");
//	TPaveStats *tps = (TPaveStats *)histogram->FindObject("stats");
//	tps->SetTextColor( color );
//	tps->SetLineColor( color );
//	tps->SetX1NDC( 0.65 );
//	tps->SetX2NDC( 0.95 );
//	tps->SetY1NDC( 0.5 );
//	tps->SetY2NDC( 0.95 );

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
	m_pTFile->cd();
	m_pTTree->Write();
	InitializeHistogram( h_ResidualPx , h_ResidualPx->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualPy , h_ResidualPy->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualPz , h_ResidualPz->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualE , h_ResidualE->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualTheta , h_ResidualTheta->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualPhi , h_ResidualPhi->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_NormalizedResidualPx , h_NormalizedResidualPx->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualPy , h_NormalizedResidualPy->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualPz , h_NormalizedResidualPz->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualE , h_NormalizedResidualE->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualTheta , h_NormalizedResidualTheta->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualPhi , h_NormalizedResidualPhi->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_ResidualPxSeen , h_ResidualPxSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualPySeen , h_ResidualPySeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualPzSeen , h_ResidualPzSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualESeen , h_ResidualESeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualThetaSeen , h_ResidualThetaSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_ResidualPhiSeen , h_ResidualPhiSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , false );
	InitializeHistogram( h_NormalizedResidualPxSeen , h_NormalizedResidualPxSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualPySeen , h_NormalizedResidualPySeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualPzSeen , h_NormalizedResidualPzSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualESeen , h_NormalizedResidualESeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualThetaSeen , h_NormalizedResidualThetaSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	InitializeHistogram( h_NormalizedResidualPhiSeen , h_NormalizedResidualPhiSeen->GetEntries() , m_histColour , 1 , 1.0 , 1 , true );
	m_pTFile->Close();
	delete m_pTFile;


}
