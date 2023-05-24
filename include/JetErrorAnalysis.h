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
		virtual void	getJetResiduals( TVector3 jetTrueMomentum , double jetTrueEnergy , TVector3 jetRecoMomentum , double jetRecoEnergy , double &jetPxResidual , double &jetPyResidual , double &jetPzResidual , double &jetEnergyResidual , double &jetThetaResidual , double &jetPhiResidual );

		/*
		* called for every reconstructed jet
		*/
		virtual void	getJetResolutions( TLorentzVector jetFourMomentum , std::vector<float> jetCovMat , double &sigmaE , double &sigmaTheta , double &sigmaPhi );


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

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		int					m_nTrueJets;
		int					m_nRecoJets;
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
		floatVector				m_ResidualPx{};
		floatVector				m_ResidualPy{};
		floatVector				m_ResidualPz{};
		floatVector				m_ResidualE{};
		floatVector				m_ResidualTheta{};
		floatVector				m_ResidualPhi{};
		floatVector				m_NormalizedResidualPx{};
		floatVector				m_NormalizedResidualPy{};
		floatVector				m_NormalizedResidualPz{};
		floatVector				m_NormalizedResidualE{};
		floatVector				m_NormalizedResidualTheta{};
		floatVector				m_NormalizedResidualPhi{};
		floatVector				m_ResidualPxSeen{};
		floatVector				m_ResidualPySeen{};
		floatVector				m_ResidualPzSeen{};
		floatVector				m_ResidualESeen{};
		floatVector				m_ResidualThetaSeen{};
		floatVector				m_ResidualPhiSeen{};
		floatVector				m_NormalizedResidualPxSeen{};
		floatVector				m_NormalizedResidualPySeen{};
		floatVector				m_NormalizedResidualPzSeen{};
		floatVector				m_NormalizedResidualESeen{};
		floatVector				m_NormalizedResidualThetaSeen{};
		floatVector				m_NormalizedResidualPhiSeen{};
		IntVector				m_trueJetType{};
		IntVector				m_trueJetFlavour{};
		IntVector				m_recoJetFlavour{};
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
		TH1F					*h_ResidualPxSeen{};
		TH1F					*h_ResidualPySeen{};
		TH1F					*h_ResidualPzSeen{};
		TH1F					*h_ResidualESeen{};
		TH1F					*h_ResidualThetaSeen{};
		TH1F					*h_ResidualPhiSeen{};
		TH1F					*h_NormalizedResidualPxSeen{};
		TH1F					*h_NormalizedResidualPySeen{};
		TH1F					*h_NormalizedResidualPzSeen{};
		TH1F					*h_NormalizedResidualESeen{};
		TH1F					*h_NormalizedResidualThetaSeen{};
		TH1F					*h_NormalizedResidualPhiSeen{};

	private:

		std::string				m_recoJetCollectionName{};
		std::string				_MCParticleColllectionName{};
		std::string				_recoParticleCollectionName{};
		std::string				_recoMCTruthLink{};
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
