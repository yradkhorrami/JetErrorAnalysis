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
#include <string>
#include "TrueJet_Parser.h"
#include <string>
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
		int					m_nTrueLeptons;
		int					m_nRecoJets;
		int					m_nRecoLeptons;
		int					m_HDecayMode;
		int					m_nSLDecayBHadron;
		int					m_nSLDecayCHadron;
		int					m_nSLDecayTotal;
		IntVector				m_trueJetType{};
		floatVector				m_trueKaonEnergy{};
		float					m_trueKaonEnergyTotal;
		floatVector				m_trueProtonEnergy{};
		float					m_trueProtonEnergyTotal;

	private:


		std::string				m_recoJetCollectionName{};
		std::string				_MCParticleColllectionName{};
		std::string				_recoParticleCollectionName{};
		std::string				_recoMCTruthLink{};
//		std::string				_trueJetCollectionName{};
		std::string				m_outputFile{};
		TFile					*m_pTFile;
	        TTree					*m_pTTree;

};

#endif
