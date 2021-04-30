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
	m_pTTree->Branch("protonTrackEnergyTotal",&m_protonTrackEnergyTotal,"protonTrackEnergyTotal/F") ;
	m_pTTree->Branch("kaonTrackEnergy",&m_kaonTrackEnergy) ;
	m_pTTree->Branch("kaonTrackEnergyTotal",&m_kaonTrackEnergyTotal,"kaonTrackEnergyTotal/F") ;
	m_pTTree->Branch("trueJetType",&m_trueJetType) ;
	m_pTTree->Branch("trueJetFlavour",&m_trueJetFlavour) ;

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
	m_pionTrackEnergy.clear();
	m_pionTrackEnergyTotal = 0.0;
	m_protonTrackEnergy.clear();
	m_protonTrackEnergyTotal = 0.0;
	m_kaonTrackEnergy.clear();
	m_kaonTrackEnergyTotal = 0.0;
	m_trueProtonEnergy.clear();
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
			streamlog_out(DEBUG0) << "	Type of True Jet[ " << i_jet << " ]: " << type_jet( i_jet ) << " ( " << trueJetType[ abs( type_jet( i_jet ) ) ] << " ) ; 	PDG of Initial Colour Neutral = " << pdg_icn_parent( initial_cn( i_jet ) ) << " ; 	Type of Final Colour Neutral : " << type_icn_parent( initial_cn( i_jet ) ) << "( " << icnType[ type_icn_parent( initial_cn( i_jet ) ) ] << " )" << std::endl;
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
				bool trueJetCriteria = false;
				bool recoJetCriteria = false;
				const EVENT::MCParticleVec& mcpVec =  true_partics( trueHadronicJetIndices[ i_jet ] );
				streamlog_out(DEBUG0) << "	Number of all MCParticles in trueJet [ " << trueHadronicJetIndices[ i_jet ] << " ] : " << mcpVec.size() << std::endl;
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
				ReconstructedParticle *recoJet = dynamic_cast<ReconstructedParticle*>( recoJetCol->getElementAt( recoJetIndices[ i_jet ] ) );
				ReconstructedParticleVec jetRecoPFOs  = recoJet->getParticles();
				streamlog_out(DEBUG0) << "	Number of all Reconstructed Particles in recoJet [ " << recoJetIndices[ i_jet ] << " ] : " << jetRecoPFOs.size() << std::endl;
				for ( unsigned int i_pfo = 0 ; i_pfo < jetRecoPFOs.size() ; ++i_pfo )
				{
					EVENT::ReconstructedParticle *testPFO = jetRecoPFOs.at( i_pfo );
					streamlog_out(DEBUG0) << "	PFO [ " << i_pfo << " ] : 	PFO Type = " << testPFO->getType() << std::endl;
					getTrackInformation( pLCEvent , testPFO );

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

void JetErrorAnalysis::getTrackInformation( LCEvent* pLCEvent , EVENT::ReconstructedParticle *testPFO )
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
		}
		else if ( trackMass == m_kaon_mass )
		{
			m_kaonTrackEnergy.push_back( trackFourMomentum.E() );
			m_kaonTrackEnergyTotal += trackFourMomentum.E();
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

void JetErrorAnalysis::end()
{
	m_pTFile->cd();
	m_pTTree->Write();
	m_pTFile->Close();
	delete m_pTFile;


}
