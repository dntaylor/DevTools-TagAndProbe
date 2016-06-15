//systemincludes
#include <memory>
#include <vector>
#include <iostream>

// CMS includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


class MiniAODElectronVIDEmbedder : public edm::EDProducer
{
public:
  explicit MiniAODElectronVIDEmbedder(const edm::ParameterSet&);
  ~MiniAODElectronVIDEmbedder() {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // Methods
  virtual void beginJob();
  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob();

  // Data
  edm::EDGetTokenT<edm::View<pat::Electron> > electronCollectionToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  std::auto_ptr<std::vector<pat::Electron> > out; // Collection we'll output at the end
};


// Constructors and destructors
// class member functions
MiniAODElectronVIDEmbedder::MiniAODElectronVIDEmbedder(const edm::ParameterSet& iConfig):
	electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.exists("src") ? 
			iConfig.getParameter<edm::InputTag>("src") :
			edm::InputTag("slimmedElectrons"))),
	vtxToken_(consumes<reco::VertexCollection>(iConfig.exists("vtxSrc") ? 
			iConfig.getParameter<edm::InputTag>("vtxSrc") :
			edm::InputTag("offlineSlimmedPrimaryVertices")))
{
	produces<std::vector<pat::Electron> >();
}


void MiniAODElectronVIDEmbedder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	out = std::auto_ptr<std::vector<pat::Electron> >(new std::vector<pat::Electron>);

	edm::Handle<edm::View<pat::Electron> > electronsIn;

	iEvent.getByToken(electronCollectionToken_, electronsIn);


	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	if (vertices->empty()) return;//skip event if no vertices found
        //std::cout<<"size of IDMap tokens: "<<idMapTokens_.size()<<std::endl;

	for(edm::View<pat::Electron>::const_iterator ei = electronsIn->begin();
			ei != electronsIn->end(); ei++) // loop over electrons
	{
		const edm::Ptr<pat::Electron> eptr(electronsIn, ei - electronsIn->begin());

		out->push_back(*ei); // copy electron to save correctly in event

		//Add electron isolation to the tree
		float eleIso03 = 999;
		eleIso03 = (eptr->pfIsolationVariables().sumChargedHadronPt + std::max(eptr->pfIsolationVariables().sumNeutralHadronEt +
					eptr->pfIsolationVariables().sumPhotonEt - 0.5 * eptr->pfIsolationVariables().sumPUPt, 0.0)) / ei->pt(); 

		int pc = 0;
		//electron conversion
		if ((ei->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS))<=1&&ei->passConversionVeto()){
			pc = 1;
		}

		//embed dxy and dz for electron

		float dZ = std::abs(ei->gsfTrack()->dz(vertices->at(0).position()));
		float dXY = std::abs(ei->gsfTrack()->dxy(vertices->at(0).position()));
 		
		int PASSID = 0;
		int PASSIDISO = 0;
		int PASSISO = 0;
                if ( pc==1 && dZ < 0.2 && dXY < 0.045) PASSID=1;
                if ( pc==1 && dZ < 0.2 && dXY < 0.045 && eleIso03 < 0.1) PASSIDISO=1;
                if ( eleIso03 < 0.1) PASSISO=1;

		out->back().addUserFloat("PassID", PASSID);
		out->back().addUserFloat("PassIDISO",PASSIDISO);
		out->back().addUserFloat("PassISO",PASSISO);

	}

	iEvent.put(out);
}


void MiniAODElectronVIDEmbedder::beginJob()
{}


void MiniAODElectronVIDEmbedder::endJob()
{}


void
MiniAODElectronVIDEmbedder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODElectronVIDEmbedder);
