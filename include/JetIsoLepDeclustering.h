#ifndef JetIsoLepDeclustering_h
#define JetIsoLepDeclustering_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "lcio.h"
#include "EVENT/LCStrVec.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include "TFile.h"
#include "TH1I.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include "TLorentzVector.h"
class TFile;
class TH1F;
class TH1I;
class TH2F;
class TTree;
using namespace lcio ;
using namespace marlin ;
class JetIsoLepDeclustering : public Processor
{
public:
	virtual Processor *newProcessor()
	{
		return new JetIsoLepDeclustering;
	}
	JetIsoLepDeclustering();
	virtual ~JetIsoLepDeclustering() = default;
	JetIsoLepDeclustering( const JetIsoLepDeclustering& ) = delete;
	JetIsoLepDeclustering &operator = ( const JetIsoLepDeclustering& ) = delete;
	virtual void init();
	virtual void Clear();
	virtual void processRunHeader();
	virtual void processEvent( EVENT::LCEvent *pLCEvent );
 	virtual void check( EVENT::LCEvent *pLCEvent );
	virtual void end();
private:

	typedef std::vector<int>		IntVector;
	typedef std::vector<double>		DoubleVector;
	typedef std::vector<float>		FloatVector;

	std::string				m_inputJetCollection{};
	std::string				m_inputIsoLepCollection{};
	std::string				m_outputPfoCollection{};
	std::string				m_outputIsolepCollection{};
	std::string				m_rootFile{};

	int					m_nRun;
	int					m_nEvt;
	int					m_nRunSum;
	int					m_nEvtSum;
	int					m_nInJets;
	int					m_nInJetPFOs;
	int					m_nOutJets;
	int					m_nOutJetPFOs;
	int					m_nInIsoLeps;
	int					m_nOutIsoLeps;
	float					m_IsoLepsInvMass;
	FloatVector				m_IsoLepPairsInvMass;
	int					n_nJetnIsolep = 0;
	bool					m_fillRootTree = true;

	int					m_nJets = 0;
	int					m_nIsoLeps = 0;
	double					m_diLepInvMass = 91.2;

	int					m_useEvent = 0;

	TFile					*m_pTFile{};
	TTree					*m_pTTree{};
	TH2F					*h_nJetnIsolep{};
};
#endif
