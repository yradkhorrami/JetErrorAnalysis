#ifndef JetErrorAnalysis_h
#define JetErrorAnalysis_h 1

//#include "TrueJet.h"
#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <iostream>
#include <string>
#include "TrueJet_Parser.h"
#include "TLorentzVector.h"
#include <TFile.h>
#include <TTree.h>
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;

using namespace lcio ;
using namespace marlin ;

class JetErrorAnalysis : public Processor , public TrueJet_Parser
{

	public:

		virtual Processor*  newProcessor()
		{
			return new JetErrorAnalysis;
		}

		JetErrorAnalysis();
		virtual ~JetErrorAnalysis() = default;
		JetErrorAnalysis(const JetErrorAnalysis&) = delete;
		JetErrorAnalysis& operator=(const JetErrorAnalysis&) = delete;

		/*
		* Called at the begin of the job before anything is read.
		* Use to initialize the processor, e.g. book histograms.
		*/
		virtual void init();
		virtual void Clear();

		/*
		* Called for every run.
		*/
		virtual void processRunHeader();

		/*
		* Called for every event - the working horse.
		*/
		virtual void processEvent( LCEvent * evt );

		/*
		* called for every pfo of reconstructed jet
		*/
		virtual void getTrackInformation( LCEvent* pLCEvent , EVENT::ReconstructedParticle *testPFO , double &KaonTrackEnergyinJet , double &ProtonTrackEnergyinJet );

		/*
		* called for every pair of true and reconstructed jets
		*/
		virtual void getJetResiduals( TLorentzVector trueJetFourMomentum , EVENT::ReconstructedParticle *recoJet );

		/*
		*
		*/
		int getTrackIndex( EVENT::LCCollection *TrackCollection , EVENT::Track* inputTrk );

		/*
		*
		*/
		TLorentzVector getTrackFourMomentum( EVENT::Track* inputTrk , double trackMass );


		virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
		virtual void doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange );


		virtual void check();

		/*
		* Called after data processing for clean up.
		*/
		virtual void end();

		std::string get_recoMCTruthLink()
		{
			return _recoMCTruthLink;
		};

	protected:

		/*
		* Input collection name.
		*/
		typedef	std::vector<int>		IntVector;
		typedef	std::vector<float>		floatVector;

		float					m_Bfield;
		double					c;
		double					mm2m;
		double					eV2GeV;
		double					eB;
		float					m_pion_mass;
		float					m_proton_mass;
		float					m_kaon_mass;
		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		int					m_nTrueJets;
		int					m_nTrueLeptons;
		int					m_nRecoJets;
		int					m_nRecoLeptons;
		int					m_HDecayMode;
		int					m_nSLDecayBHadron;
		int					m_nSLDecayCHadron;
		int					m_nSLDecayTotal;
		IntVector				m_trueJetType{};
		IntVector				m_trueJetFlavour{};
		floatVector				m_trueKaonEnergy{};
		float					m_trueKaonEnergyTotal;
		floatVector				m_trueProtonEnergy{};
		float					m_trueProtonEnergyTotal;
		floatVector				m_pionTrackEnergy{};
		float					m_pionTrackEnergyTotal;
		floatVector				m_protonTrackEnergy{};
		floatVector				m_ProtonTrackEnergyinJet{};
		float					m_protonTrackEnergyTotal;
		floatVector				m_kaonTrackEnergy{};
		floatVector				m_KaonTrackEnergyinJet{};
		float					m_kaonTrackEnergyTotal;
		floatVector				m_ResidualPx{};
		floatVector				m_ResidualPy{};
		floatVector				m_ResidualPz{};
		floatVector				m_ResidualE{};
		floatVector				m_ResidualTheta{};
		floatVector				m_ResidualPhi{};
		int					n_ResidualPx;
		int					n_ResidualPy;
		int					n_ResidualPz;
		int					n_ResidualE;
		int					n_ResidualTheta;
		int					n_ResidualPhi;
		floatVector				m_NormalizedResidualPx{};
		floatVector				m_NormalizedResidualPy{};
		floatVector				m_NormalizedResidualPz{};
		floatVector				m_NormalizedResidualE{};
		floatVector				m_NormalizedResidualTheta{};
		floatVector				m_NormalizedResidualPhi{};
		int					n_NormalizedResidualPx;
		int					n_NormalizedResidualPy;
		int					n_NormalizedResidualPz;
		int					n_NormalizedResidualE;
		int					n_NormalizedResidualTheta;
		int					n_NormalizedResidualPhi;
		TH1F					*h_ResidualPx{};
		TH1F					*h_ResidualPy{};
		TH1F					*h_ResidualPz{};
		TH1F					*h_ResidualE{};
		TH1F					*h_ResidualTheta{};
		TH1F					*h_ResidualPhi{};
		TH1F					*h_NormalizedResidualPx{};
		TH1F					*h_NormalizedResidualPy{};
		TH1F					*h_NormalizedResidualPz{};
		TH1F					*h_NormalizedResidualE{};
		TH1F					*h_NormalizedResidualTheta{};
		TH1F					*h_NormalizedResidualPhi{};

	private:

		std::string				m_referenceJetCollection{};
		std::string				m_recoJetCollectionName{};
		std::string				_MCParticleColllectionName{};
		std::string				m_MarlinTrkTracks{};
		std::string				m_MarlinTrkTracksKAON{};
		std::string				m_MarlinTrkTracksPROTON{};
		std::string				_recoParticleCollectionName{};
		std::string				_recoMCTruthLink{};
//		std::string				_trueJetCollectionName{};
		std::string				m_outputFile{};
		std::string				m_histName{};
		int					m_histColour{};
		float					m_minKaonTrackEnergy{};
		float					m_minProtonTrackEnergy{};
		TFile					*m_pTFile;
	        TTree					*m_pTTree;

};

#endif
