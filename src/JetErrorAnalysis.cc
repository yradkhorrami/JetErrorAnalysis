#include "JetErrorAnalysis.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>


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
m_nTrueLeptons(0),
m_nRecoJets(0),
m_nRecoLeptons(0),
m_HDecayMode(0),
m_nSLDecayBHadron(0),
m_nSLDecayCHadron(0),
m_nSLDecayTotal(0),
m_trueKaonEnergyTotal(0.0),
m_trueProtonEnergyTotal(0.0),
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

	registerProcessorParameter(	"outputFilename",
					"name of output file",
					m_outputFile,
					std::string("")
				);


	// Inputs: True jets (as a recoparticle, will be the sum of the _reconstructed particles_
	// created by the true particles in each true jet, in the RecoMCTruthLink sense.
	// link jet-to-reco particles, link jet-to-MC-particles.

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

	// usually a good idea to
	printParameters() ;

	m_nRun = 0 ;
	m_nEvt = 0 ;

	m_pTFile = new TFile(m_outputFile.c_str(),"recreate");

	m_pTTree = new TTree("eventTree","eventTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("nTrueJets",&m_nTrueJets,"nTrueJets/I") ;
	m_pTTree->Branch("nTrueLeptons",&m_nTrueLeptons,"nTrueLeptons/I") ;
	m_pTTree->Branch("nRecoJets",&m_nRecoJets,"nRecoJets/I") ;
	m_pTTree->Branch("nRecoLeptons",&m_nRecoLeptons,"nRecoLeptons/I") ;
	m_pTTree->Branch("HDecayMode",&m_HDecayMode,"HDecayMode/I") ;
	m_pTTree->Branch("nSLDecayBHadron",&m_nSLDecayBHadron,"nSLDecayBHadron/I") ;
	m_pTTree->Branch("nSLDecayCHadron",&m_nSLDecayCHadron,"nSLDecayCHadron/I") ;
	m_pTTree->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree->Branch("trueKaonEnergy",&m_trueKaonEnergy) ;
	m_pTTree->Branch("trueKaonEnergyTotal",&m_trueKaonEnergyTotal,"trueKaonEnergyTotal/F") ;
	m_pTTree->Branch("trueProtonEnergy",&m_trueProtonEnergy) ;
	m_pTTree->Branch("trueProtonEnergyTotal",&m_trueProtonEnergyTotal,"trueProtonEnergyTotal/F") ;
	m_pTTree->Branch("trueJetType",&m_trueJetType) ;

}

void JetErrorAnalysis::Clear()
{
	m_nTrueJets = 0;
	m_nTrueLeptons = 0;
	m_nTrueLeptons = 0;
	m_nRecoJets = 0;
	m_nRecoLeptons = 0;
	m_HDecayMode = 0;
	m_nSLDecayBHadron = 0;
	m_nSLDecayCHadron = 0;
	m_nSLDecayTotal = 0;
	m_trueKaonEnergyTotal = 0.0;
	m_trueKaonEnergy.clear();
	m_trueProtonEnergyTotal = 0.0;
	m_trueProtonEnergy.clear();
	m_trueJetType.clear();
}

void JetErrorAnalysis::processRunHeader()
{
    m_nRun++ ;
}


void JetErrorAnalysis::processEvent( LCEvent* pLCEvent)
{
	LCCollection *recoJetCol{};
	LCCollection *trueJetCol{};
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event: 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	try
	{
		recoJetCol		= pLCEvent->getCollection( m_recoJetCollectionName );
		trueJetCol		= pLCEvent->getCollection( _trueJetCollectionName );
		TrueJet_Parser* trueJet	= this;
		trueJet->getall(pLCEvent);

		m_nRecoJets = recoJetCol->getNumberOfElements();
		streamlog_out(DEBUG0) << "	Number of Reconstructed Jets: " << m_nRecoJets << std::endl;

		int njets = trueJet->njets();
		std::vector<int> trueHadronicJetIndices; trueHadronicJetIndices.clear();
		std::vector<int> recoJetIndices; recoJetIndices.clear();
		for (int i_jet = 0 ; i_jet < njets ; i_jet++ )
		{
			m_trueJetType.push_back( type_jet( i_jet ) );
			if ( type_jet( i_jet ) == 1 )
			{
				++m_nTrueJets;
				trueHadronicJetIndices.push_back( i_jet );
			}
		}
		streamlog_out(DEBUG0) << "	Number of True Hadronic Jets(type = 1): " << m_nTrueJets << std::endl;
		if ( m_nRecoJets == m_nTrueJets )
		{
			for ( int i_trueJet = 0 ; i_trueJet < m_nTrueJets ; ++i_trueJet )
			{
				TVector3 trueJetMomentum( ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 0 ] , ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 1 ] , ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 2 ] );
				streamlog_out(DEBUG0) << "	True(seen) Jet Momentum[ " << trueHadronicJetIndices[ i_trueJet ] << " ]: (	" << ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 0 ] << " 	, " << ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 1 ] << " 	, " << ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 2 ] << "	)" << std::endl;
				TVector3 trueJetMomentumUnit = trueJetMomentum; trueJetMomentumUnit.SetMag(1.0);
				float widestAngleCosTheta = -1.0;
				int matchedRecoJetIndex = -1;
				for ( int i_recoJet = 0 ; i_recoJet < m_nRecoJets ; ++i_recoJet )
				{
					ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( i_recoJet ) );
					TVector3 recoJetMometnum( recoJet->getMomentum() );
					TVector3 recoJetMometnumUnit = recoJetMometnum; recoJetMometnumUnit.SetMag(1.0);
					streamlog_out(DEBUG0) << "	Reco Jet Momentum[ " << i_recoJet << " ]: (	" << recoJet->getMomentum()[ 0 ] << " 	, " << recoJet->getMomentum()[ 1 ] << " 	, " << recoJet->getMomentum()[ 2 ] << "	)" << std::endl;
					if ( trueJetMomentumUnit.Dot( recoJetMometnumUnit ) >= widestAngleCosTheta )
					{
						widestAngleCosTheta = trueJetMomentumUnit.Dot( recoJetMometnumUnit );
						matchedRecoJetIndex = i_recoJet;
					}
				}
				streamlog_out(DEBUG0) << "	True(seen) Jet [ " << trueHadronicJetIndices[ i_trueJet ] << " ] is matched with RecoJet [ " << matchedRecoJetIndex << " ]" << std::endl;
				recoJetIndices.push_back( matchedRecoJetIndex );
			}
			for ( int i_jet = 0 ; i_jet < m_nTrueJets ; ++i_jet )
			{
				const EVENT::MCParticleVec& mcpVec =  true_partics( trueHadronicJetIndices[ i_jet ] );
				streamlog_out(DEBUG0) << "	Number of all MCParticles in Jet [ " << trueHadronicJetIndices[ i_jet ] << " ] : " << mcpVec.size() << std::endl;
				for ( unsigned int i_mcp = 0 ; i_mcp < mcpVec.size() ; ++i_mcp )
				{
					EVENT::MCParticle *testMCP = mcpVec.at( i_mcp );
					streamlog_out(DEBUG0) << "	MCParticle [ " << i_mcp << " ] : 	GeneratorStatus = " << testMCP->getGeneratorStatus() << " ; 	PDGCode = " << testMCP->getPDG() << std::endl;
					if ( testMCP->getGeneratorStatus() == 1 && abs( testMCP->getPDG() ) == 321 )
					{
						m_trueKaonEnergy.push_back( testMCP->getEnergy() );
						m_trueKaonEnergyTotal += testMCP->getEnergy();
					}
					else if ( testMCP->getGeneratorStatus() == 1 && abs( testMCP->getPDG() ) == 2212 )
					{
						m_trueProtonEnergy.push_back( testMCP->getEnergy() );
						m_trueProtonEnergyTotal += testMCP->getEnergy();
					}
				}
			}

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


void JetErrorAnalysis::check()
{

}

void JetErrorAnalysis::end()
{
	m_pTFile->cd();
	m_pTTree->Write();
	m_pTFile->Close();
	delete m_pTFile;


}
