#include "JetIsoLepDeclustering.h"
#include <iostream>
#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TTree.h"
using namespace lcio ;
using namespace marlin ;

JetIsoLepDeclustering JetIsoLepDeclustering;

JetIsoLepDeclustering::JetIsoLepDeclustering():
	Processor("JetIsoLepDeclustering"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_nInJets(0),
	m_nInJetPFOs(0),
	m_nOutJets(0),
	m_nOutJetPFOs(0),
	m_nInIsoLeps(0),
	m_nOutIsoLeps(0),
	m_IsoLepsInvMass(0.f),
	m_IsoLepPairsInvMass{}
{
	_description = "JetIsoLepDeclustering checks the number of IsolatedLeptons and jets, removes unwanted leptons from IsolatedLeptons and adds the lepton to declustered jet particles. The declustered jet particles should be clustered again!";

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"JetCollection",
					"Name of input jet collection",
					m_inputJetCollection,
					std::string("Durham_nJets")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"isoLepCollection",
					"Name of input Isolated Lepton collection",
					m_inputIsoLepCollection,
					std::string("ISOLeptons")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"outputPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("declusteredJetIsoleps")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"outputIsolepCollection",
					"Name of output isolated lepton collection",
					m_outputIsolepCollection,
					std::string("IsolatedLeptons")
				);

	registerProcessorParameter(	"nJets",
					"Number of jet should be in the event",
					m_nJets,
					int(0)
				);

	registerProcessorParameter(	"nIsoLeps",
					"Number of Isolated Leptons should be in the event",
					m_nIsoLeps,
					int(0)
				);

	registerProcessorParameter(	"diLepInvMass",
					"Invariant mass of di-lepton system in Isolated Leptons [GeV]",
					m_diLepInvMass,
					double(91.2)
				);

	registerProcessorParameter(	"fillRootTree",
					"Fill root tree to check processor performance",
					m_fillRootTree,
					bool(true)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("Output.root")
				);

}

void JetIsoLepDeclustering::init()
{

	streamlog_out(DEBUG0) << "   init called  " << std::endl ;
	printParameters();
	if ( m_fillRootTree )
	{
		streamlog_out(DEBUG0) << "	Creating root file/tree/histograms" << std::endl ;
		m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
		m_pTTree = new TTree("SLDCorrection", "SLDCorrection");
		m_pTTree->SetDirectory(m_pTFile);
		m_pTTree->Branch("event", &m_nEvt, "event/I");
		m_pTTree->Branch("nInJets", &m_nInJets, "nInJets/I");
		m_pTTree->Branch("nInJetPFOs", &m_nInJetPFOs, "nInJetPFOs/I");
		m_pTTree->Branch("nOutJets", &m_nOutJets, "nOutJets/I");
		m_pTTree->Branch("nOutJetPFOs", &m_nOutJetPFOs, "nOutJetPFOs/I");
		m_pTTree->Branch("nInIsoLeps", &m_nInIsoLeps, "nInIsoLeps/I");
		m_pTTree->Branch("nOutIsoLeps", &m_nOutIsoLeps, "nOutIsoLeps/I");
		m_pTTree->Branch("IsoLepsInvMass", &m_IsoLepsInvMass, "IsoLepsInvMass/F");
		m_pTTree->Branch("IsoLepPairsInvMass", &m_IsoLepPairsInvMass);
		h_nJetnIsolep = new TH2F( "nJetnIsolep" , "; nJets ; nIsoLeps" , 5 , 0.0 , 5.0 , 5 , 0.0 , 5.0 );
		h_nJetnIsolep->GetXaxis()->SetBinLabel(1,"0");
		h_nJetnIsolep->GetXaxis()->SetBinLabel(2,"1");
		h_nJetnIsolep->GetXaxis()->SetBinLabel(3,"2");
		h_nJetnIsolep->GetXaxis()->SetBinLabel(4,"3");
		h_nJetnIsolep->GetXaxis()->SetBinLabel(5,"4");
		h_nJetnIsolep->GetYaxis()->SetBinLabel(1,"0");
		h_nJetnIsolep->GetYaxis()->SetBinLabel(2,"1");
		h_nJetnIsolep->GetYaxis()->SetBinLabel(3,"2");
		h_nJetnIsolep->GetYaxis()->SetBinLabel(4,"3");
		h_nJetnIsolep->GetYaxis()->SetBinLabel(5,"4");
		streamlog_out(DEBUG0) << "	Created root file/tree/histograms" << std::endl ;
	}
	this->Clear();
	streamlog_out(DEBUG0) << "   init finished successfully" << std::endl ;
}

void JetIsoLepDeclustering::Clear()
{
	streamlog_out(DEBUG0) << "   clear called" << std::endl ;
	m_nInJets = 0;
	m_nInJetPFOs = 0;
	m_nOutJets = 0;
	m_nOutJetPFOs = 0;
	m_nInIsoLeps = 0;
	m_nOutIsoLeps = 0;
	m_IsoLepsInvMass = 0;
	m_IsoLepPairsInvMass.clear();
	m_useEvent = 1;
	streamlog_out(DEBUG0) << "   clear finished successfully" << std::endl ;
}

void JetIsoLepDeclustering::processRunHeader()
{
	streamlog_out(DEBUG0) << "   processRunHeader called" << std::endl ;
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
	streamlog_out(DEBUG0) << "   processRunHeader finished successfully" << std::endl ;
}

void JetIsoLepDeclustering::processEvent( EVENT::LCEvent *pLCEvent )
{
	streamlog_out(DEBUG0) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	processEvent Called	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl ;
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	m_nEvtSum++;
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(DEBUG4) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(DEBUG4) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	this->Clear();

	bool trueNJets = false;
	bool trueNIsoLeps = false;

	const EVENT::LCCollection *JetCollection{};
	const EVENT::LCCollection *IsoleptonCollection{};
	IMPL::LCCollectionVec* m_IsolepCol(NULL);
	m_IsolepCol = new IMPL::LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	m_IsolepCol->setSubset( true );
	IMPL::LCCollectionVec* m_NewPFOsCol(NULL);
	m_NewPFOsCol = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	m_NewPFOsCol->setSubset( true );
	try
	{
		JetCollection = pLCEvent->getCollection( m_inputJetCollection );
		IsoleptonCollection = pLCEvent->getCollection( m_inputIsoLepCollection );

		m_nInJets = JetCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << m_nInJets << " jets in event, looking for " << m_nJets << " jets" << std::endl;
		if ( m_nInJets != m_nJets )
		{
			trueNJets = false;
			streamlog_out( DEBUG3 ) << "	Number of jets in the event mismatches the asked number of jets, --------EVENT REJECTED--------" << std::endl;
		}
		else
		{
			trueNJets = true;
			streamlog_out( DEBUG3 ) << "	Number of jets in the event matches the asked number of jets, --------EVENT ACCEPTED--------" << std::endl;
		}
		for ( int i_jet = 0 ; i_jet < m_nInJets ; ++i_jet )
		{
			ReconstructedParticle* jet = static_cast<ReconstructedParticle*>( JetCollection->getElementAt( i_jet ) );
			m_nInJetPFOs += jet->getParticles().size();
			streamlog_out( DEBUG3 ) << "	Jet " << i_jet << " has " << jet->getParticles().size() << " PFOs" << std::endl;
			for ( unsigned int i_pfo = 0 ; i_pfo < jet->getParticles().size() ; ++i_pfo )
			{
				m_NewPFOsCol->addElement( jet->getParticles()[ i_pfo ] );
				++m_nOutJetPFOs;
			}
		}

		m_nInIsoLeps = IsoleptonCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << m_nInIsoLeps << " isolated leptons in event, looking for " << m_nIsoLeps << " isolated leptons" << std::endl;
		if ( m_fillRootTree )
		{
			h_nJetnIsolep->Fill( m_nInJets + 0.5 , m_nInIsoLeps + 0.5 );
			++n_nJetnIsolep;
		}
		if ( m_nInIsoLeps < m_nIsoLeps )
		{
			trueNIsoLeps = false;
			streamlog_out( DEBUG3 ) << "	Number of isolated leptons in the event is less than the asked number of isolated leptons, --------EVENT REJECTED--------" << std::endl;
		}
		else if ( m_nInIsoLeps == m_nIsoLeps )
		{
			trueNIsoLeps = true;
			for ( int i_lep = 0 ; i_lep < m_nInIsoLeps ; ++i_lep )
			{
				m_IsolepCol->addElement( IsoleptonCollection->getElementAt( i_lep ) );
				streamlog_out( DEBUG3 ) << "	One Isolated lepton added to new isolated leptons" << std::endl;
				++m_nOutIsoLeps;
			}
			ReconstructedParticle* lepton1 = static_cast<ReconstructedParticle*>( IsoleptonCollection->getElementAt( 0 ) );
			ReconstructedParticle* lepton2 = static_cast<ReconstructedParticle*>( IsoleptonCollection->getElementAt( 1 ) );
			TLorentzVector lep1FourMomentum = TLorentzVector( lepton1->getMomentum()[ 0 ] , lepton1->getMomentum()[ 1 ] , lepton1->getMomentum()[ 2 ] , lepton1->getEnergy() );
			TLorentzVector lep2FourMomentum = TLorentzVector( lepton2->getMomentum()[ 0 ] , lepton2->getMomentum()[ 1 ] , lepton2->getMomentum()[ 2 ] , lepton2->getEnergy() );
			m_IsoLepsInvMass = ( lep1FourMomentum + lep2FourMomentum ).M();
			streamlog_out( DEBUG3 ) << "	Number of isolated leptons in the event matches the asked number of isolated leptons, --------EVENT ACCEPTED--------" << std::endl;
		}
		else
		{
			std::vector<double> massDiff{};
			std::vector<int> lep1Index{};
			std::vector<int> lep2Index{};
			for ( int i_lep1 = 0 ; i_lep1 < m_nInIsoLeps - 1 ; ++i_lep1 )
			{
				ReconstructedParticle* lepton1 = static_cast<ReconstructedParticle*>( IsoleptonCollection->getElementAt( i_lep1 ) );
				TLorentzVector lep1FourMomentum = TLorentzVector( lepton1->getMomentum()[ 0 ] , lepton1->getMomentum()[ 1 ] , lepton1->getMomentum()[ 2 ] , lepton1->getEnergy() );
				for ( int i_lep2 = i_lep1 + 1 ; i_lep2 < m_nInIsoLeps ; ++i_lep2 )
				{
					ReconstructedParticle* lepton2 = static_cast<ReconstructedParticle*>( IsoleptonCollection->getElementAt( i_lep2 ) );
					TLorentzVector lep2FourMomentum = TLorentzVector( lepton2->getMomentum()[ 0 ] , lepton2->getMomentum()[ 1 ] , lepton2->getMomentum()[ 2 ] , lepton2->getEnergy() );
					double diLepInvMass = ( lep1FourMomentum + lep2FourMomentum ).M();
					m_IsoLepPairsInvMass.push_back( diLepInvMass );
					massDiff.push_back( fabs( diLepInvMass - m_diLepInvMass ) );
					lep1Index.push_back( i_lep1 );
					lep2Index.push_back( i_lep2 );
				}
			}
			double smallestMassDiff = m_diLepInvMass;
			int iLepton1 = -1;
			int iLepton2 = -1;
			for ( unsigned int i_pair = 0 ; i_pair < massDiff.size() ; ++i_pair )
			{
				if ( massDiff[ i_pair ] < smallestMassDiff )
				{
					smallestMassDiff = massDiff[ i_pair ];
					iLepton1 = lep1Index[ i_pair ];
					iLepton2 = lep2Index[ i_pair ];
					m_IsoLepsInvMass = m_IsoLepPairsInvMass[ i_pair ];
				}
			}
			for ( int i_lep = 0 ; i_lep < m_nInIsoLeps ; ++i_lep )
			{
				if ( i_lep == iLepton1 || i_lep == iLepton2 )
				{
					m_IsolepCol->addElement( IsoleptonCollection->getElementAt( i_lep ) );
					streamlog_out( DEBUG3 ) << "	One Isolated lepton added to new isolated leptons" << std::endl;
					++m_nOutIsoLeps;
				}
				else
				{
					m_NewPFOsCol->addElement( IsoleptonCollection->getElementAt( i_lep ) );
					streamlog_out( DEBUG3 ) << "	One lepton removed from IsolatedLeptons and added to new PFO collection" << std::endl;
					++m_nOutJetPFOs;
				}
			}
			if ( iLepton1 != -1 && iLepton2 != -1 )
			{
				trueNIsoLeps = true;
				streamlog_out( DEBUG3 ) << "	Two leptons remained as isolated in the event, rest added to PFOsminusIsoLep to be clustered in jets, --------EVENT ACCEPTED--------" << std::endl;
			}
		}
		m_useEvent = ( trueNJets && trueNIsoLeps ? 1 : 0 );
		m_IsolepCol->parameters().setValue( "useEvent" , ( int )m_useEvent );

		streamlog_out( DEBUG3 ) << "	All New PFOs are " << m_NewPFOsCol->getNumberOfElements() << " PFOs in total" << std::endl;
		streamlog_out( DEBUG3 ) << "	All New Isolated Leptons are " << m_IsolepCol->getNumberOfElements() << " IsoLeps in total" << std::endl;

		pLCEvent->addCollection( m_IsolepCol , m_outputIsolepCollection.c_str() );
		pLCEvent->addCollection( m_NewPFOsCol , m_outputPfoCollection.c_str() );
		if ( m_fillRootTree )
		{
			m_pTTree->Fill();
		}
	}
	catch( DataNotAvailableException &e )
	{
		streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
	}
}

void JetIsoLepDeclustering::check( EVENT::LCEvent *pLCEvent )
{
	const EVENT::LCCollection *inJetCollection{};
	const EVENT::LCCollection *inIsoleptonCollection{};
	const EVENT::LCCollection *outPFOCollection{};
	const EVENT::LCCollection *outIsoleptonCollection{};
	try
	{
		inJetCollection = pLCEvent->getCollection( m_inputJetCollection );
		inIsoleptonCollection = pLCEvent->getCollection( m_inputIsoLepCollection );
		outPFOCollection = pLCEvent->getCollection( m_outputPfoCollection );
		outIsoleptonCollection = pLCEvent->getCollection( m_outputIsolepCollection );
		int nJets = inJetCollection->getNumberOfElements();
		int nInPFOs = 0;
		for ( int i_jet = 0 ; i_jet < nJets ; ++i_jet )
		{
			ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( inJetCollection->getElementAt( i_jet ) );
			nInPFOs += jet->getParticles().size();
		}
		int nInIsoLep = inIsoleptonCollection->getNumberOfElements();
		int nOutPFOs = outPFOCollection->getNumberOfElements();
		int nOutIsoLep = outIsoleptonCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << nJets << " jets with " << nInPFOs << " PFOs and " << nInIsoLep << " Isolated Leptons converted to " << nOutPFOs << " PFOs and " << nOutIsoLep << " Isolated Leptons" << std::endl;
	}
	catch( DataNotAvailableException &e )
        {
          streamlog_out( WARNING ) << "	Input/Output collections not found in event: " << m_nEvt << std::endl;
        }
}

void JetIsoLepDeclustering::end()
{
	if ( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree->Write();
		h_nJetnIsolep->Scale( 100.0 / n_nJetnIsolep );
		h_nJetnIsolep->Write();
		m_pTFile->Close();
		delete m_pTFile;
	}
}
