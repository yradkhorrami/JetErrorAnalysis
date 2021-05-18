#include "JetErrorAnalysis.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
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
m_Bfield(0.f),
c(0.),
mm2m(0.),
eV2GeV(0.),
eB(0.),
m_pion_mass(0.0),
m_proton_mass(0.0),
m_kaon_mass(0.0),
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
m_pionTrackEnergyTotal(0.0),
m_protonTrackEnergyTotal(0.0),
m_kaonTrackEnergyTotal(0.0),
n_NormalizedResidualPx(0.0),
n_NormalizedResidualPy(0.0),
n_NormalizedResidualPz(0.0),
n_NormalizedResidualE(0.0),
n_NormalizedResidualTheta(0.0),
n_NormalizedResidualPhi(0.0),
m_pTFile(NULL),
m_pTTree(NULL)
{

	// modify processor description
	_description = "JetErrorAnalysis does whatever it does ..." ;


	// register steering parameters: name, description, class-variable, default value


	// Inputs: MC-particles, Reco-particles, the link between the two

	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
					"referenceJetCollection" ,
					"Name of the Reference Jet collection"  ,
					m_referenceJetCollection ,
					std::string("Durham_nJets")
				);

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

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollection" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracks,
					std::string("MarlinTrkTracks")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionKaon" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksKAON,
					std::string("MarlinTrkTracksKaon")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionProton" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksPROTON,
					std::string("MarlinTrkTracksProton")
				);

	registerProcessorParameter(	"outputFilename",
					"name of output file",
					m_outputFile,
					std::string("")
				);

	registerProcessorParameter(	"HistogramsName",
					"name of histograms",
					m_histName,
					std::string("Std. Trk")
				);

	registerProcessorParameter(	"HistogramsColour",
					"color of histograms",
					m_histColour,
					int(1)
				);

	registerProcessorParameter(	"minKaonTrackEnergy",
					"min Energy of Kaons Tracks for histograming",
					m_minKaonTrackEnergy,
					float(0.0)
				);

	registerProcessorParameter(	"minProtonTrackEnergy",
					"min Energy of Protons Tracks for histograming",
					m_minProtonTrackEnergy,
					float(0.0)
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
//	m_Bfield = MarlinUtil::getBzAtOrigin();
	m_Bfield = 3.5;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;
	printParameters() ;

	m_pion_mass = 0.13957018;
	m_proton_mass = 0.938272088;
	m_kaon_mass = 0.493677;

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
	m_pTTree->Branch("pionTrackEnergy",&m_pionTrackEnergy) ;
	m_pTTree->Branch("pionTrackEnergyTotal",&m_pionTrackEnergyTotal,"pionTrackEnergyTotal/F") ;
	m_pTTree->Branch("protonTrackEnergy",&m_protonTrackEnergy) ;
	m_pTTree->Branch("protonTrackEnergyinJet",&m_ProtonTrackEnergyinJet) ;
	m_pTTree->Branch("protonTrackEnergyTotal",&m_protonTrackEnergyTotal,"protonTrackEnergyTotal/F") ;
	m_pTTree->Branch("kaonTrackEnergy",&m_kaonTrackEnergy) ;
	m_pTTree->Branch("kaonTrackEnergyinJet",&m_KaonTrackEnergyinJet) ;
	m_pTTree->Branch("kaonTrackEnergyTotal",&m_kaonTrackEnergyTotal,"kaonTrackEnergyTotal/F") ;
	m_pTTree->Branch("NormalizedResidualPx",&m_NormalizedResidualPx) ;
	m_pTTree->Branch("NormalizedResidualPy",&m_NormalizedResidualPy) ;
	m_pTTree->Branch("NormalizedResidualPz",&m_NormalizedResidualPz) ;
	m_pTTree->Branch("NormalizedResidualE",&m_NormalizedResidualE) ;
	m_pTTree->Branch("NormalizedResidualTheta",&m_NormalizedResidualTheta) ;
	m_pTTree->Branch("NormalizedResidualPhi",&m_NormalizedResidualPhi) ;
	m_pTTree->Branch("trueJetType",&m_trueJetType) ;
	m_pTTree->Branch("trueJetFlavour",&m_trueJetFlavour) ;
	h_NormalizedResidualPx = new TH1F( m_histName.c_str() , "; (_{}p_{x,jet}^{REC} - p_{x,jet}^{MC}) / #sigma_{p_{x,jet}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NormalizedResidualPx = 0;
	h_NormalizedResidualPy = new TH1F( m_histName.c_str() , "; (_{}p_{y,jet}^{REC} - p_{y,jet}^{MC}) / #sigma_{p_{y,jet}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NormalizedResidualPy = 0;
	h_NormalizedResidualPz = new TH1F( m_histName.c_str() , "; (_{}p_{z,jet}^{REC} - p_{z,jet}^{MC}) / #sigma_{p_{z,jet}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NormalizedResidualPz = 0;
	h_NormalizedResidualE = new TH1F( m_histName.c_str() , "; (_{}E_{jet}^{REC} - E_{jet}^{MC}) / #sigma_{E_{jet}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NormalizedResidualE = 0;
	h_NormalizedResidualTheta = new TH1F( m_histName.c_str() , "; (_{}#theta_{jet}^{REC} - #theta_{jet}^{MC}) / #sigma_{#theta_{jet}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NormalizedResidualTheta = 0;
	h_NormalizedResidualPhi = new TH1F( m_histName.c_str() , "; (_{}#phi_{jet}^{REC} - #phi_{jet}^{MC}) / #sigma_{#phi_{jet}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NormalizedResidualPhi = 0;

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
	m_pionTrackEnergy.clear();
	m_pionTrackEnergyTotal = 0.0;
	m_protonTrackEnergy.clear();
	m_ProtonTrackEnergyinJet.clear();
	m_protonTrackEnergyTotal = 0.0;
	m_kaonTrackEnergy.clear();
	m_KaonTrackEnergyinJet.clear();
	m_kaonTrackEnergyTotal = 0.0;
	m_NormalizedResidualPx.clear();
	m_NormalizedResidualPy.clear();
	m_NormalizedResidualPz.clear();
	m_NormalizedResidualE.clear();
	m_NormalizedResidualTheta.clear();
	m_NormalizedResidualPhi.clear();
	m_trueJetType.clear();
	m_trueJetFlavour.clear();
}

void JetErrorAnalysis::processRunHeader()
{
    m_nRun++ ;
}


void JetErrorAnalysis::processEvent( LCEvent* pLCEvent)
{
	LCCollection *recoJetCol{};
	LCCollection *trueJetCol{};
	LCCollection *refJetCol{};
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	this->Clear();
	std::string trueJetType[6]{ "hadronic (string)" , "leptonic" , "hadronic(cluster)" , "ISR" , "overlay" , "M.E. photon" };
	std::string icnType[6]{ "quark pair" , "lepton pair" , "quark pair" , "ISR" , "???" , "M.E. photon" };
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event: 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	try
	{
		recoJetCol		= pLCEvent->getCollection( m_recoJetCollectionName );
		trueJetCol		= pLCEvent->getCollection( _trueJetCollectionName );
		refJetCol		= pLCEvent->getCollection( m_referenceJetCollection );
		TrueJet_Parser* trueJet	= this;
		trueJet->getall(pLCEvent);

		m_nRecoJets = recoJetCol->getNumberOfElements();
		streamlog_out(DEBUG3) << "	Number of Reconstructed Jets: " << m_nRecoJets << std::endl;

		int njets = trueJet->njets();
		streamlog_out(DEBUG3) << "	Number of True Jets: " << njets << std::endl;
		std::vector<int> trueHadronicJetIndices; trueHadronicJetIndices.clear();
		std::vector<int> recoJetIndices; recoJetIndices.clear();
		for (int i_jet = 0 ; i_jet < njets ; i_jet++ )
		{
			m_trueJetType.push_back( type_jet( i_jet ) );
			streamlog_out(DEBUG0) << "	Type of True Jet[ " << i_jet << " ]: " << type_jet( i_jet ) << " ( " << trueJetType[ abs( type_jet( i_jet ) ) ] << " ) ; 	PDG of Initial Colour Neutral = " << pdg_icn_parent( initial_cn( i_jet ) ) << " ; 	Type of Final Colour Neutral : " << type_icn_parent( initial_cn( i_jet ) ) << "( " << icnType[ type_icn_parent( initial_cn( i_jet ) ) ] << " )" << std::endl;
			if ( type_jet( i_jet ) == 1 )
			{
				++m_nTrueJets;
				trueHadronicJetIndices.push_back( i_jet );
			}
		}
		streamlog_out(DEBUG3) << "	Number of True Hadronic Jets(type = 1): " << m_nTrueJets << std::endl;
		double KaonTrackEnergyinJet;
		double ProtonTrackEnergyinJet;
		if ( m_nRecoJets == m_nTrueJets )
		{
			for ( int i_trueJet = 0 ; i_trueJet < m_nTrueJets ; ++i_trueJet )
			{
				TVector3 trueJetMomentum( ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 0 ] , ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 1 ] , ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 2 ] );
				streamlog_out(DEBUG2) << "	True(seen) Jet Momentum[ " << trueHadronicJetIndices[ i_trueJet ] << " ]: (	" << ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 0 ] << " 	, " << ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 1 ] << " 	, " << ptrueseen( trueHadronicJetIndices[ i_trueJet ] )[ 2 ] << "	)" << std::endl;
				TVector3 trueJetMomentumUnit = trueJetMomentum; trueJetMomentumUnit.SetMag(1.0);
				float CosWidestAngle = -1.0;
				int matchedRecoJetIndex = -1;
				for ( int i_recoJet = 0 ; i_recoJet < m_nRecoJets ; ++i_recoJet )
				{
					ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( i_recoJet ) );
					TVector3 recoJetMometnum( recoJet->getMomentum() );
					TVector3 recoJetMometnumUnit = recoJetMometnum; recoJetMometnumUnit.SetMag(1.0);
					streamlog_out(DEBUG2) << "	Reco Jet Momentum[ " << i_recoJet << " ]: (	" << recoJet->getMomentum()[ 0 ] << " 	, " << recoJet->getMomentum()[ 1 ] << " 	, " << recoJet->getMomentum()[ 2 ] << "	)" << std::endl;
					if ( trueJetMomentumUnit.Dot( recoJetMometnumUnit ) >= CosWidestAngle )
					{
						CosWidestAngle = trueJetMomentumUnit.Dot( recoJetMometnumUnit );
						matchedRecoJetIndex = i_recoJet;
					}
				}
				streamlog_out(DEBUG2) << "	True(seen) Jet [ " << trueHadronicJetIndices[ i_trueJet ] << " ] is matched with RecoJet [ " << matchedRecoJetIndex << " ]" << std::endl;
				recoJetIndices.push_back( matchedRecoJetIndex );
			}
			for ( int i_jet = 0 ; i_jet < m_nTrueJets ; ++i_jet )
			{
				bool trueJetCriteria = false;
				bool recoJetCriteria = false;
				KaonTrackEnergyinJet = 0.0;
				ProtonTrackEnergyinJet = 0.0;
				const EVENT::MCParticleVec& mcpVec =  true_partics( trueHadronicJetIndices[ i_jet ] );
				streamlog_out(DEBUG3) << "	Number of all MCParticles in trueJet [ " << trueHadronicJetIndices[ i_jet ] << " ] : " << mcpVec.size() << std::endl;
				for ( unsigned int i_mcp = 0 ; i_mcp < mcpVec.size() ; ++i_mcp )
				{
					EVENT::MCParticle *testMCP = mcpVec.at( i_mcp );
					streamlog_out(DEBUG1) << "	MCParticle [ " << i_mcp << " ] : 	GeneratorStatus = " << testMCP->getGeneratorStatus() << " ; 	PDGCode = " << testMCP->getPDG() << std::endl;
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
				ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( recoJetIndices[ i_jet ] ) );
				ReconstructedParticle *refJet = dynamic_cast<ReconstructedParticle*>( refJetCol->getElementAt( recoJetIndices[ i_jet ] ) );
				ReconstructedParticleVec jetRecoPFOs  = recoJet->getParticles();
				ReconstructedParticleVec refjetRecoPFOs  = refJet->getParticles();
				streamlog_out(DEBUG3) << "	Number of all Reconstructed Particles in recoJet [ " << recoJetIndices[ i_jet ] << " ] : " << jetRecoPFOs.size() << std::endl;
				streamlog_out(DEBUG3) << "	Number of all Reconstructed Particles in refJet [ " << recoJetIndices[ i_jet ] << " ] : " << refjetRecoPFOs.size() << std::endl;
				if ( jetRecoPFOs.size() != refjetRecoPFOs.size() ) continue;
				for ( unsigned int i_pfo = 0 ; i_pfo < jetRecoPFOs.size() ; ++i_pfo )
				{
					EVENT::ReconstructedParticle *testPFO = jetRecoPFOs.at( i_pfo );
					EVENT::ReconstructedParticle *refPFO = refjetRecoPFOs.at( i_pfo );
					streamlog_out(DEBUG1) << "	PFO [ " << i_pfo << " ] : 	PFO Type = " << testPFO->getType() << std::endl;
					getTrackInformation( pLCEvent , refPFO , KaonTrackEnergyinJet , ProtonTrackEnergyinJet );
				}
				TLorentzVector trueJetFourMomentum( p4trueseen( trueHadronicJetIndices[ i_jet ] )[ 1 ] , p4trueseen( trueHadronicJetIndices[ i_jet ] )[ 2 ] , p4trueseen( trueHadronicJetIndices[ i_jet ] )[ 3 ] , p4trueseen( trueHadronicJetIndices[ i_jet ] )[ 0 ] );
				m_KaonTrackEnergyinJet.push_back( KaonTrackEnergyinJet );
				m_ProtonTrackEnergyinJet.push_back( ProtonTrackEnergyinJet );
				if ( KaonTrackEnergyinJet >= m_minKaonTrackEnergy && ProtonTrackEnergyinJet >= m_minProtonTrackEnergy ) getJetResiduals( trueJetFourMomentum , recoJet );
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

void JetErrorAnalysis::getJetResiduals( TLorentzVector trueJetFourMomentum , EVENT::ReconstructedParticle *recoJet )
{
	double trueJetPx = trueJetFourMomentum.Px();
	double trueJetPy = trueJetFourMomentum.Py();
	double trueJetPz = trueJetFourMomentum.Pz();
	double trueJetE = trueJetFourMomentum.E();
	double trueJetTheta = trueJetFourMomentum.Theta();
	double trueJetPhi = trueJetFourMomentum.Phi();
	TVector3 trueP( trueJetPx , trueJetPy , trueJetPz );
	TVector3 truePunit = trueP; truePunit.SetMag(1.0);
	TVector3 truePt( trueJetPx , trueJetPy , 0.0 );
	TVector3 truePtunit = truePt; truePtunit.SetMag(1.0);
	TLorentzVector recoJetFourMomentum( recoJet->getMomentum()[ 0 ] , recoJet->getMomentum()[ 1 ] , recoJet->getMomentum()[ 2 ] , recoJet->getEnergy() );
	double recoJetPx = recoJetFourMomentum.Px();
	double recoJetPy = recoJetFourMomentum.Py();
	double recoJetPz = recoJetFourMomentum.Pz();
	double recoJetPt = std::sqrt( pow( recoJetPx , 2 ) + pow( recoJetPy , 2 ) );
	double recoJetPt2 = pow( recoJetPx , 2 ) + pow( recoJetPy , 2 );
	double recoJetP2 = pow( recoJetPx , 2 ) + pow( recoJetPy , 2 ) + pow( recoJetPz , 2 );
	double recoJetE = recoJetFourMomentum.E();
	double recoJetTheta = recoJetFourMomentum.Theta();
	double recoJetPhi = recoJetFourMomentum.Phi();
	TVector3 recoP( recoJetPx , recoJetPy , recoJetPz );
	TVector3 recoPunit = recoP; recoPunit.SetMag(1.0);
	TVector3 recoProtated = recoP; recoProtated.SetMag(1.0); recoProtated.SetPhi( trueJetPhi );
	TVector3 recoPt( recoJetPx , recoJetPy , 0.0 );
	TVector3 recoPtunit = recoPt; recoPtunit.SetMag(1.0);
	double sigmaPx2 = recoJet->getCovMatrix()[ 0 ];
	double sigmaPxPy = recoJet->getCovMatrix()[ 1 ];
	double sigmaPy2 = recoJet->getCovMatrix()[ 2 ];
	double sigmaPxPz = recoJet->getCovMatrix()[ 3 ];
	double sigmaPyPz = recoJet->getCovMatrix()[ 4 ];
	double sigmaPz2 = recoJet->getCovMatrix()[ 5 ];
	double sigmaPxE = recoJet->getCovMatrix()[ 6 ];
	double sigmaPyE = recoJet->getCovMatrix()[ 7 ];
	double sigmaPzE = recoJet->getCovMatrix()[ 8 ];
	double sigmaE2 = recoJet->getCovMatrix()[ 9 ];
	double dTheta_dPx = recoJetPx * recoJetPz / ( recoJetP2 * recoJetPt );
	double dTheta_dPy = recoJetPy * recoJetPz / ( recoJetP2 * recoJetPt );
	double dTheta_dPz = -recoJetPt / recoJetP2;
	double dPhi_dPx = -recoJetPy / recoJetPt2;
	double dPhi_dPy = -recoJetPx / recoJetPt2;
	double sigmaTheta = std::sqrt( std::fabs( sigmaPx2 * std::pow( dTheta_dPx , 2 ) + sigmaPy2 * std::pow( dTheta_dPy , 2 ) + sigmaPz2 * std::pow( dTheta_dPz , 2 ) + 2 * ( sigmaPxPy * dTheta_dPx * dTheta_dPy ) + 2 * ( sigmaPxPz * dTheta_dPx * dTheta_dPz ) + 2 * ( sigmaPyPz * dTheta_dPy * dTheta_dPz ) ) );
	double sigmaPhi = std::sqrt( std::fabs( sigmaPx2 * std::pow( dPhi_dPx , 2 ) + sigmaPy2 * std::pow( dPhi_dPy , 2 ) + 2 * ( sigmaPxPy * dPhi_dPx * dPhi_dPy ) ) );
	m_NormalizedResidualPx.push_back( ( recoJetPx - trueJetPx ) / std::sqrt( sigmaPx2 ) );
	h_NormalizedResidualPx->Fill( ( recoJetPx - trueJetPx ) / std::sqrt( sigmaPx2 ) ); ++n_NormalizedResidualPx;
	m_NormalizedResidualPy.push_back( ( recoJetPy - trueJetPy ) / std::sqrt( sigmaPy2 ) );
	h_NormalizedResidualPy->Fill( ( recoJetPy - trueJetPy ) / std::sqrt( sigmaPy2 ) ); ++n_NormalizedResidualPy;
	m_NormalizedResidualPz.push_back( ( recoJetPz - trueJetPz ) / std::sqrt( sigmaPz2 ) );
	h_NormalizedResidualPz->Fill( ( recoJetPz - trueJetPz ) / std::sqrt( sigmaPz2 ) ); ++n_NormalizedResidualPz;
	m_NormalizedResidualE.push_back( ( recoJetE - trueJetE ) / std::sqrt( sigmaE2 ) );
	h_NormalizedResidualE->Fill( ( recoJetE - trueJetE ) / std::sqrt( sigmaE2 ) ); ++n_NormalizedResidualE;
	double ThetaResidual = ( ( recoJetTheta - trueJetTheta ) > 0 ? acos( truePunit.Dot(recoProtated) ) : -1 * acos( truePunit.Dot(recoProtated) ) );
	m_NormalizedResidualTheta.push_back( ThetaResidual / sigmaTheta );
	h_NormalizedResidualTheta->Fill( ThetaResidual / sigmaTheta ); ++n_NormalizedResidualTheta;
	double PhiResidual = ( ( recoJetPhi - trueJetPhi ) > 0 ? acos( truePtunit.Dot(recoPtunit) ) : -1 * acos( truePtunit.Dot(recoPtunit) ) );
	m_NormalizedResidualPhi.push_back( PhiResidual / sigmaPhi );
	h_NormalizedResidualPhi->Fill( PhiResidual / sigmaPhi ); ++n_NormalizedResidualPhi;
}


void JetErrorAnalysis::getTrackInformation( LCEvent* pLCEvent , EVENT::ReconstructedParticle *testPFO , double &KaonTrackEnergyinJet , double &ProtonTrackEnergyinJet )
{
	LCCollection *MarlinTrkTracks{};
	LCCollection *MarlinTrkTracksKAON{};
	LCCollection *MarlinTrkTracksPROTON{};
	try
	{
		MarlinTrkTracks = pLCEvent->getCollection(m_MarlinTrkTracks);
		MarlinTrkTracksKAON = pLCEvent->getCollection(m_MarlinTrkTracksKAON);
		MarlinTrkTracksPROTON = pLCEvent->getCollection(m_MarlinTrkTracksPROTON);
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "	Could not find one of  " << MarlinTrkTracks << " / " << MarlinTrkTracksKAON << " / " << MarlinTrkTracksPROTON << " Collection" << std::endl;
	}
	const EVENT::TrackVec& inputPFOtrkvec = testPFO->getTracks();
	int nTRKsofPFO = inputPFOtrkvec.size();
	TLorentzVector pfoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	for ( int i_trk = 0 ; i_trk < nTRKsofPFO ; ++i_trk )
	{
		Track *pfoTrk = (Track*)inputPFOtrkvec.at( i_trk );
		float trackMass = 0.0;
		TLorentzVector trackFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
		if ( getTrackIndex( MarlinTrkTracksPROTON , pfoTrk ) != -1 )
		{
			trackMass = m_proton_mass;
		}
		else if ( getTrackIndex( MarlinTrkTracksKAON , pfoTrk ) != -1 )
		{
			trackMass = m_kaon_mass;
		}
		else
		{
			trackMass = m_pion_mass;
		}//( getTrackIndex( MarlinTrkTracks , pfoTrk ) != -1  )
		trackFourMomentum = getTrackFourMomentum( pfoTrk , trackMass );
		if ( trackMass == m_proton_mass )
		{
			m_protonTrackEnergy.push_back( trackFourMomentum.E() );
			m_protonTrackEnergyTotal += trackFourMomentum.E();
			ProtonTrackEnergyinJet += trackFourMomentum.E();
		}
		else if ( trackMass == m_kaon_mass )
		{
			m_kaonTrackEnergy.push_back( trackFourMomentum.E() );
			m_kaonTrackEnergyTotal += trackFourMomentum.E();
			KaonTrackEnergyinJet += trackFourMomentum.E();
		}
		else
		{
			m_pionTrackEnergy.push_back( trackFourMomentum.E() );
			m_pionTrackEnergyTotal += trackFourMomentum.E();
		}

	}
}

int JetErrorAnalysis::getTrackIndex( EVENT::LCCollection *TrackCollection , EVENT::Track* inputTrk )
{
	LCCollection *MarlinTrkTracks{};
	try
	{
		MarlinTrkTracks = TrackCollection;//pLCEvent->getCollection( TrackCollection );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "	Could not find the " << TrackCollection << " Collection" << std::endl;
		return -1;
	}
	unsigned int nTRKs = MarlinTrkTracks->getNumberOfElements();
	int track_index = -1;
	for (unsigned int i_trk = 0; i_trk < nTRKs;  ++i_trk )
	{
		Track* pionTrack = dynamic_cast<EVENT::Track*>( MarlinTrkTracks->getElementAt( i_trk ) );
		if ( pionTrack == inputTrk )
		{
			track_index = i_trk;
		}
	}
	if ( track_index == -1)
	{
		streamlog_out(DEBUG1) << "	Coudln't find track_index!!!  " << std::endl;
	}
	else
	{
		streamlog_out(DEBUG1) << "	Track index in " << m_MarlinTrkTracks << " collection is " << track_index << std::endl;
	}
	return track_index;
}


TLorentzVector JetErrorAnalysis::getTrackFourMomentum( EVENT::Track* inputTrk , double trackMass )
{
	streamlog_out(DEBUG1) << "	------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "	Calculating PFO 4-momentum from track parameters" << std::endl;
	streamlog_out(DEBUG1) << "	------------------------------------------------" << std::endl;
	double Phi = inputTrk->getPhi();
	double Omega = inputTrk->getOmega();
	double tanLambda = inputTrk->getTanLambda();
	streamlog_out(DEBUG0) << "	Track parameters obtained" << std::endl;
	double pT = eB / fabs( Omega );
	double px = pT * TMath::Cos( Phi );
	double py = pT * TMath::Sin( Phi );
	double pz = pT * tanLambda;
	double E = sqrt( pow( trackMass , 2 ) + px * px + py * py + pz * pz);
	streamlog_out(DEBUG0) << "	Track parameters is converted to (p,E)" << std::endl;
	TLorentzVector trackFourMomentum( px , py , pz , E );
	return trackFourMomentum;
}

void JetErrorAnalysis::check()
{

}

void JetErrorAnalysis::InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle )
{
	histogram->Scale( 1.0 / scale );
	histogram->SetLineColor( color );
	histogram->SetLineWidth( lineWidth );
	histogram->SetMarkerSize( markerSize );
	histogram->SetMarkerStyle( markerStyle );
	histogram->SetMarkerColor( color );
	float fit_range = 2.0;
	float fit_min = -2.0;
	float fit_max = 2.0;
	doProperGaussianFit( histogram , fit_min , fit_max , fit_range );
	histogram->GetFunction("gaus")->SetLineColor( color );
	float y_max = 1.2 * histogram->GetMaximum();
	histogram->GetYaxis()->SetRangeUser(0.0, y_max);
	histogram->GetXaxis()->SetTitleSize(0.06);
	histogram->GetXaxis()->SetTitleOffset(1.10);
	histogram->GetXaxis()->SetLabelSize(0.06);
	histogram->GetYaxis()->SetTitleSize(0.06);
	histogram->GetYaxis()->SetTitleOffset(1.30);
	histogram->GetYaxis()->SetLabelSize(0.06);
	histogram->Write();
	/*
	gPad->Update();
	TPaveStats *tps = (TPaveStats *)histogram->FindObject("stats");
	TPaveStats *tps = (TPaveStats *)histogram->FindObject("stats");
	tps->SetTextColor( color );
	tps->SetLineColor( color );
	tps->SetX1NDC( 0.64 );
	tps->SetX2NDC( 0.94 );
	tps->SetY1NDC( 0.5 );
	tps->SetY2NDC( 0.95 );
*/
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
	InitializeHistogram( h_NormalizedResidualPx , n_NormalizedResidualPx , m_histColour , 1 , 1.0 , 1 );
	InitializeHistogram( h_NormalizedResidualPy , n_NormalizedResidualPy , m_histColour , 1 , 1.0 , 1 );
	InitializeHistogram( h_NormalizedResidualPz , n_NormalizedResidualPz , m_histColour , 1 , 1.0 , 1 );
	InitializeHistogram( h_NormalizedResidualE , n_NormalizedResidualE , m_histColour , 1 , 1.0 , 1 );
	InitializeHistogram( h_NormalizedResidualTheta , n_NormalizedResidualTheta , m_histColour , 1 , 1.0 , 1 );
	InitializeHistogram( h_NormalizedResidualPhi , n_NormalizedResidualPhi , m_histColour , 1 , 1.0 , 1 );
	m_pTFile->Close();
	delete m_pTFile;


}
