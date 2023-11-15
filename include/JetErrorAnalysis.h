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
#include "DDMarlinCED.h"
#include <iostream>
#include <string>
#include "TrueJet_Parser.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TStyle.h"
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

		typedef std::vector<EVENT::MCParticle*>			mcpVector;
		typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;

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
		* called for every pair of true and reconstructed jets
		*/
		virtual void	getResiduals( TVector3 jetTrueMomentum , double jetTrueEnergy , TVector3 jetRecoMomentum , double jetRecoEnergy , double &jetPxResidual , double &jetPyResidual , double &jetPzResidual , double &jetEnergyResidual , double &jetThetaResidual , double &jetPhiResidual );

		/*
		* called for every reconstructed jet
		*/
		virtual void	getResolutions( TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , double &sigmaE , double &sigmaTheta , double &sigmaPhi );
		EVENT::ReconstructedParticle* getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO );

		virtual void InitializeHistogram( TH1F *histogram , bool scale , int color , int lineWidth , int markerSize , int markerStyle , bool fitGaussian );
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

		int						m_nRun;
		int						m_nEvt;
		int						m_nRunSum;
		int						m_nEvtSum;
		int						m_nTrueJets;
		int						m_nTrueIsoLeps;
		int						m_nRecoJets;
		int						m_nIsoLeptons;
		IntVector				m_nSeenISRs{};
		int						m_nSLDecayBHadron{};
		int						m_nSLDecayCHadron{};
		int						m_nSLDecayTauLepton{};
		int						m_nSLDecayTotal{};
		int						m_nCorrectedSLD{};
		floatVector				m_quarkPx{};
		floatVector				m_quarkPy{};
		floatVector				m_quarkPz{};
		floatVector				m_quarkE{};
		floatVector				m_quarkTheta{};
		floatVector				m_quarkPhi{};
		floatVector				m_trueJetPx{};
		floatVector				m_trueJetPy{};
		floatVector				m_trueJetPz{};
		floatVector				m_trueJetE{};
		floatVector				m_trueJetTheta{};
		floatVector				m_trueJetPhi{};
		floatVector				m_trueSeenJetPx{};
		floatVector				m_trueSeenJetPy{};
		floatVector				m_trueSeenJetPz{};
		floatVector				m_trueSeenJetE{};
		floatVector				m_trueSeenJetTheta{};
		floatVector				m_trueSeenJetPhi{};
		floatVector				m_seenJetPx{};
		floatVector				m_seenJetPy{};
		floatVector				m_seenJetPz{};
		floatVector				m_seenJetE{};
		floatVector				m_seenJetTheta{};
		floatVector				m_seenJetPhi{};
		floatVector				m_recoJetPx{};
		floatVector				m_recoJetPy{};
		floatVector				m_recoJetPz{};
		floatVector				m_recoJetE{};
		floatVector				m_recoJetTheta{};
		floatVector				m_recoJetPhi{};
		floatVector				m_jetSigmaPx2{};
		floatVector				m_jetSigmaPxPy{};
		floatVector				m_jetSigmaPy2{};
		floatVector				m_jetSigmaPxPz{};
		floatVector				m_jetSigmaPyPz{};
		floatVector				m_jetSigmaPz2{};
		floatVector				m_jetSigmaPxE{};
		floatVector				m_jetSigmaPyE{};
		floatVector				m_jetSigmaPzE{};
		floatVector				m_jetSigmaE2{};
		floatVector				m_jetSigmaTheta2{};
		floatVector				m_jetSigmaPhi2{};
		floatVector				m_ResidualJetPx{};
		floatVector				m_ResidualJetPy{};
		floatVector				m_ResidualJetPz{};
		floatVector				m_ResidualJetE{};
		floatVector				m_ResidualJetTheta{};
		floatVector				m_ResidualJetPhi{};
		floatVector				m_NormalizedResidualJetPx{};
		floatVector				m_NormalizedResidualJetPy{};
		floatVector				m_NormalizedResidualJetPz{};
		floatVector				m_NormalizedResidualJetE{};
		floatVector				m_NormalizedResidualJetTheta{};
		floatVector				m_NormalizedResidualJetPhi{};
		floatVector				m_ResidualJetPxSeen{};
		floatVector				m_ResidualJetPySeen{};
		floatVector				m_ResidualJetPzSeen{};
		floatVector				m_ResidualJetESeen{};
		floatVector				m_ResidualJetThetaSeen{};
		floatVector				m_ResidualJetPhiSeen{};
		floatVector				m_NormalizedResidualJetPxSeen{};
		floatVector				m_NormalizedResidualJetPySeen{};
		floatVector				m_NormalizedResidualJetPzSeen{};
		floatVector				m_NormalizedResidualJetESeen{};
		floatVector				m_NormalizedResidualJetThetaSeen{};
		floatVector				m_NormalizedResidualJetPhiSeen{};
		IntVector				m_trueJetType{};
		IntVector				m_trueJetFlavour{};
		IntVector				m_recoJetFlavour{};
		TH1F					*h_ResidualJetPx{};
		TH1F					*h_ResidualJetPy{};
		TH1F					*h_ResidualJetPz{};
		TH1F					*h_ResidualJetE{};
		TH1F					*h_ResidualJetTheta{};
		TH1F					*h_ResidualJetPhi{};
		TH1F					*h_NormalizedResidualJetPx{};
		TH1F					*h_NormalizedResidualJetPy{};
		TH1F					*h_NormalizedResidualJetPz{};
		TH1F					*h_NormalizedResidualJetE{};
		TH1F					*h_NormalizedResidualJetTheta{};
		TH1F					*h_NormalizedResidualJetPhi{};
		TH1F					*h_ResidualJetPxSeen{};
		TH1F					*h_ResidualJetPySeen{};
		TH1F					*h_ResidualJetPzSeen{};
		TH1F					*h_ResidualJetESeen{};
		TH1F					*h_ResidualJetThetaSeen{};
		TH1F					*h_ResidualJetPhiSeen{};
		TH1F					*h_NormalizedResidualJetPxSeen{};
		TH1F					*h_NormalizedResidualJetPySeen{};
		TH1F					*h_NormalizedResidualJetPzSeen{};
		TH1F					*h_NormalizedResidualJetESeen{};
		TH1F					*h_NormalizedResidualJetThetaSeen{};
		TH1F					*h_NormalizedResidualJetPhiSeen{};

		floatVector				m_trueIsoLepPx{};
		floatVector				m_trueIsoLepPy{};
		floatVector				m_trueIsoLepPz{};
		floatVector				m_trueIsoLepE{};
		floatVector				m_trueIsoLepTheta{};
		floatVector				m_trueIsoLepPhi{};
		floatVector				m_seenIsoLepPx{};
		floatVector				m_seenIsoLepPy{};
		floatVector				m_seenIsoLepPz{};
		floatVector				m_seenIsoLepE{};
		floatVector				m_seenIsoLepTheta{};
		floatVector				m_seenIsoLepPhi{};
		floatVector				m_recoIsoLepPx{};
		floatVector				m_recoIsoLepPy{};
		floatVector				m_recoIsoLepPz{};
		floatVector				m_recoIsoLepE{};
		floatVector				m_recoIsoLepTheta{};
		floatVector				m_recoIsoLepPhi{};
		floatVector				m_isoLepSigmaPx2{};
		floatVector				m_isoLepSigmaPxPy{};
		floatVector				m_isoLepSigmaPy2{};
		floatVector				m_isoLepSigmaPxPz{};
		floatVector				m_isoLepSigmaPyPz{};
		floatVector				m_isoLepSigmaPz2{};
		floatVector				m_isoLepSigmaPxE{};
		floatVector				m_isoLepSigmaPyE{};
		floatVector				m_isoLepSigmaPzE{};
		floatVector				m_isoLepSigmaE2{};
		floatVector				m_isoLepSigmaTheta2{};
		floatVector				m_isoLepSigmaPhi2{};
		floatVector				m_ResidualIsoLepPx{};
		floatVector				m_ResidualIsoLepPy{};
		floatVector				m_ResidualIsoLepPz{};
		floatVector				m_ResidualIsoLepE{};
		floatVector				m_ResidualIsoLepTheta{};
		floatVector				m_ResidualIsoLepPhi{};
		floatVector				m_NormalizedResidualIsoLepPx{};
		floatVector				m_NormalizedResidualIsoLepPy{};
		floatVector				m_NormalizedResidualIsoLepPz{};
		floatVector				m_NormalizedResidualIsoLepE{};
		floatVector				m_NormalizedResidualIsoLepTheta{};
		floatVector				m_NormalizedResidualIsoLepPhi{};
		floatVector				m_ResidualIsoLepPxSeen{};
		floatVector				m_ResidualIsoLepPySeen{};
		floatVector				m_ResidualIsoLepPzSeen{};
		floatVector				m_ResidualIsoLepESeen{};
		floatVector				m_ResidualIsoLepThetaSeen{};
		floatVector				m_ResidualIsoLepPhiSeen{};
		floatVector				m_NormalizedResidualIsoLepPxSeen{};
		floatVector				m_NormalizedResidualIsoLepPySeen{};
		floatVector				m_NormalizedResidualIsoLepPzSeen{};
		floatVector				m_NormalizedResidualIsoLepESeen{};
		floatVector				m_NormalizedResidualIsoLepThetaSeen{};
		floatVector				m_NormalizedResidualIsoLepPhiSeen{};
		IntVector				m_trueIsoLepType{};
		TH1F					*h_ResidualIsoLepPx{};
		TH1F					*h_ResidualIsoLepPy{};
		TH1F					*h_ResidualIsoLepPz{};
		TH1F					*h_ResidualIsoLepE{};
		TH1F					*h_ResidualIsoLepTheta{};
		TH1F					*h_ResidualIsoLepPhi{};
		TH1F					*h_NormalizedResidualIsoLepPx{};
		TH1F					*h_NormalizedResidualIsoLepPy{};
		TH1F					*h_NormalizedResidualIsoLepPz{};
		TH1F					*h_NormalizedResidualIsoLepE{};
		TH1F					*h_NormalizedResidualIsoLepTheta{};
		TH1F					*h_NormalizedResidualIsoLepPhi{};
		TH1F					*h_ResidualIsoLepPxSeen{};
		TH1F					*h_ResidualIsoLepPySeen{};
		TH1F					*h_ResidualIsoLepPzSeen{};
		TH1F					*h_ResidualIsoLepESeen{};
		TH1F					*h_ResidualIsoLepThetaSeen{};
		TH1F					*h_ResidualIsoLepPhiSeen{};
		TH1F					*h_NormalizedResidualIsoLepPxSeen{};
		TH1F					*h_NormalizedResidualIsoLepPySeen{};
		TH1F					*h_NormalizedResidualIsoLepPzSeen{};
		TH1F					*h_NormalizedResidualIsoLepESeen{};
		TH1F					*h_NormalizedResidualIsoLepThetaSeen{};
		TH1F					*h_NormalizedResidualIsoLepPhiSeen{};

	private:

		std::string				m_inputIsolatedleptonCollection{};
		std::string				m_recoJetCollectionName{};
		std::string				_MCParticleColllectionName{};
		std::string				_recoParticleCollectionName{};
		std::string				_recoMCTruthLink{};
		std::string				_MCTruthRecoLink{};
		int						m_jetMatchingMethod{};
//		std::string				_trueJetCollectionName{};
		std::string				m_outputFile{};
		std::string				m_histName{};
		int						m_histColour{};
		bool					m_fillRootTree = true;
		bool					m_normalizeHistograms = true;
		float					m_minKaonTrackEnergy{};
		float					m_minProtonTrackEnergy{};
		TFile					*m_pTFile;
	    TTree					*m_pTTree;

};

#endif
