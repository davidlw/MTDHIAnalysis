// -*- C++ -*-
//
// Package:    MTDHIAnalysis/GenPIDMatchingJets
// Class:      GenPIDMatchingJets
//
/**\class GenPIDMatchingJets GenPIDMatchingJets.cc MTDHIAnalysis/GenPIDMatchingJets/plugins/GenPIDMatchingJets.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wei Li
//         Created:  Wed, 03 Oct 2018 19:59:33 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>

#include <TH2.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

const double massPi = 0.13957018;
const double massK = 0.493677;
const double massP = 0.938272013;

const double clight = 0.3;
const double Bmag = 3.8;
const double R = 1.161;
const double Length = 3.0;

class GenPIDMatchingJets : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenPIDMatchingJets(const edm::ParameterSet&);
      ~GenPIDMatchingJets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      double GetFlightDistance(double eta, double pt, double q);
 
      // ----------member data ---------------------------
//      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    edm::EDGetTokenT<reco::GenJetCollection> tok_genJet_;

    edm::Service<TFileService> fs;

    TH2F* hBetavsP;
    TH2F* hBetaRatiovsP_pi;
    TH2F* hBetaRatiovsP_K;
    TH2F* hBetaRatiovsP_p;

    TH2F* hJetsEtaVsPt;

    TH2F* hJetDrVsPt;
    TH2F* hJetDrVsPt_pi;
    TH2F* hJetDrVsPt_k;
    TH2F* hJetDrVsPt_p;
    TH2F* hJetDrVsPt_notp;

    bool isETL;

    double sigmaT_;
    double nSigmaT_;
   
    double minEtaJet_;
    double maxEtaJet_;
    double minPtJet_;
    double maxPtJet_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenPIDMatchingJets::GenPIDMatchingJets(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));
  tok_genJet_ = consumes<reco::GenJetCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenJetCollection")));

  isETL = iConfig.getUntrackedParameter<bool>("isETL");

  sigmaT_ = iConfig.getUntrackedParameter<double>("sigmaT");
  nSigmaT_ = iConfig.getUntrackedParameter<double>("nSigmaT");

  minEtaJet_ = iConfig.getUntrackedParameter<double>("minEtaJet");
  maxEtaJet_ = iConfig.getUntrackedParameter<double>("maxEtaJet");
  minPtJet_ = iConfig.getUntrackedParameter<double>("minPtJet");
  maxPtJet_ = iConfig.getUntrackedParameter<double>("maxPtJet");
}


GenPIDMatchingJets::~GenPIDMatchingJets()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenPIDMatchingJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_,genpars);

    edm::Handle<reco::GenJetCollection> genjets;
    iEvent.getByToken(tok_genJet_,genjets);

    TF1* func_pi = new TF1("func_pi","sqrt([0]*[0]/x/x+1)",0,100); 
    func_pi->SetParameter(0,massPi);
    TF1* func_K = new TF1("func_K","sqrt([0]*[0]/x/x+1)",0,100);
    func_K->SetParameter(0,massK);
    TF1* func_p = new TF1("func_p","sqrt([0]*[0]/x/x+1)",0,100);
    func_p->SetParameter(0,massP);

    std::vector<TLorentzVector> pVect_pos;
    std::vector<TLorentzVector> pVect_neg;
    std::vector<bool> isPionVect_pos;
    std::vector<bool> isKaonVect_pos;
    std::vector<bool> isProtonVect_pos;
    std::vector<bool> isPionVect_neg;
    std::vector<bool> isKaonVect_neg;
    std::vector<bool> isProtonVect_neg;

    std::vector<TLorentzVector> pVect_jets;

    for(unsigned it=0; it<genjets->size(); ++it){

      const reco::GenJet & jet = (*genjets)[it];

      hJetsEtaVsPt->Fill(jet.eta(),jet.pt());

      if(fabs(jet.eta())<minEtaJet_ || fabs(jet.eta())>maxEtaJet_) continue;
      if(fabs(jet.pt())<minPtJet_ || fabs(jet.pt())>maxPtJet_) continue;

      TLorentzVector pjet(jet.px(),jet.py(),jet.pz(),jet.energy());
      pVect_jets.push_back(pjet);
    }

    for(unsigned it=0; it<genpars->size(); ++it){

      const reco::GenParticle & trk = (*genpars)[it];

      if(!trk.charge()) continue;
      if(trk.status()!=1) continue;
      if(trk.pt()<0.001) continue;
      if(fabs(trk.eta())>4) continue;

      if(trk.pt()<0.7) continue;

      int id = trk.pdgId();
      double mass = trk.mass();

      double eta = trk.eta();
      double pt = trk.pt();
      double p = trk.p();
      double q = trk.charge();
      double L = -9999.9;
      L = GetFlightDistance(eta,pt,q);
//      if(L<=0) continue;

      if(fabs(id)!=211 && fabs(id)!=321 && fabs(id)!=2212 && fabs(id)!=11 && fabs(id)!=13) continue;

      double beta_inv = sqrt(mass*mass/pt/pt/cosh(eta)/cosh(eta)+1);
      double tsmear = gRandom->Gaus(0,sigmaT_);
      double tof = beta_inv*L/clight-tsmear;

      double beta_measured = tof*clight/L;
      double mass2_measured = (beta_measured*beta_measured-1)*pt*pt*cosh(eta)*cosh(eta);
      double energy_measured = sqrt(mass2_measured+p*p);

      if(fabs(eta)<3)
      {
        hBetavsP->Fill(p,beta_measured);
        hBetaRatiovsP_pi->Fill(p,beta_measured/func_pi->Eval(p));
        hBetaRatiovsP_K->Fill(p,beta_measured/func_K->Eval(p));
        hBetaRatiovsP_p->Fill(p,beta_measured/func_p->Eval(p));
      }
      TLorentzVector particle(trk.px(),trk.py(),trk.pz(),energy_measured);

      bool isKaon = false;
      bool isPion = false;
      bool isProton = false;

      double sigma_beta = sigmaT_*clight/L;
/*
      if((fabs(eta)<3 && fabs(beta_measured-func_pi->Eval(p))/sigma_beta<nSigmaT_) || fabs(eta)>3) isPion = true;
      if((fabs(eta)<3 && fabs(beta_measured-func_K->Eval(p))/sigma_beta<nSigmaT_) || fabs(eta)>3) isKaon = true;
      if((fabs(eta)<3 && fabs(beta_measured-func_p->Eval(p))/sigma_beta<nSigmaT_) || fabs(eta)>3) isProton = true;
*/
/*
      if(fabs(eta)<3 && fabs(id)==211 && fabs(func_pi->Eval(p)-func_K->Eval(p))/sigma_beta>nSigmaT_ && fabs(func_pi->Eval(p)-func_p->Eval(p))/sigma_beta>nSigmaT_ ) isPion = true;
      if(fabs(eta)<3 && fabs(id)==321 && fabs(func_pi->Eval(p)-func_K->Eval(p))/sigma_beta>nSigmaT_ && fabs(func_K->Eval(p)-func_p->Eval(p))/sigma_beta>nSigmaT_ ) isKaon = true;
      if(fabs(eta)<3 && fabs(id)==2212 && fabs(func_p->Eval(p)-func_K->Eval(p))/sigma_beta>nSigmaT_ && fabs(func_pi->Eval(p)-func_p->Eval(p))/sigma_beta>nSigmaT_ ) isProton = true;
*/
      if(fabs(id)==211 ) isPion = true;
      if(fabs(id)==321 ) isKaon = true;
      if(fabs(id)==2212 ) isProton = true;

      if(q>0)
      {
        pVect_pos.push_back(particle);
        isPionVect_pos.push_back(isPion);
        isKaonVect_pos.push_back(isKaon);
        isProtonVect_pos.push_back(isProton);
      }

      if(q<0) 
      { 
        pVect_neg.push_back(particle);
        isPionVect_neg.push_back(isPion);
        isKaonVect_neg.push_back(isKaon);
        isProtonVect_neg.push_back(isProton);
      }
    }

    for(unsigned int i=0;i<pVect_jets.size();i++)
    {
      TLorentzVector pjet = pVect_jets[i];

      for(unsigned int i=0;i<pVect_pos.size();i++)
      {
        TLorentzVector p1 = pVect_pos[i];
        double deltaR = p1.DeltaR(pjet);
     
        hJetDrVsPt->Fill(deltaR,p1.Pt());

        bool isPion = isPionVect_pos[i];
        bool isKaon = isKaonVect_pos[i];
        bool isProton = isProtonVect_pos[i];

        if(isPion && !isKaon && !isProton) hJetDrVsPt_pi->Fill(deltaR,p1.Pt());
        if(!isPion && isKaon && !isProton) hJetDrVsPt_k->Fill(deltaR,p1.Pt());
        if(!isPion && !isKaon && isProton) hJetDrVsPt_p->Fill(deltaR,p1.Pt());
        if(!isProton) hJetDrVsPt_notp->Fill(deltaR,p1.Pt());
      }

      for(unsigned int i=0;i<pVect_neg.size();i++)
      {
        TLorentzVector p1 = pVect_neg[i];
        double deltaR = p1.DeltaR(pjet);

        hJetDrVsPt->Fill(deltaR,p1.Pt());

        bool isPion = isPionVect_neg[i];
        bool isKaon = isKaonVect_neg[i];
        bool isProton = isProtonVect_neg[i];

        if(isPion && !isKaon && !isProton) hJetDrVsPt_pi->Fill(deltaR,p1.Pt());
        if(!isPion && isKaon && !isProton) hJetDrVsPt_k->Fill(deltaR,p1.Pt());
        if(!isPion && !isKaon && isProton) hJetDrVsPt_p->Fill(deltaR,p1.Pt());
        if(!isProton) hJetDrVsPt_notp->Fill(deltaR,p1.Pt());

      }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
GenPIDMatchingJets::beginJob()
{
   hBetavsP = fs->make<TH2F>("BetavsP",";p (GeV/c);1/#beta;",1200,0,6,2500,0,5);
   hBetaRatiovsP_pi = fs->make<TH2F>("BetaRatiovsP_pi",";p (GeV/c);(1/#beta)_{ratio};",1200,0,6,2500,0,3);
   hBetaRatiovsP_K = fs->make<TH2F>("BetaRatiovsP_K",";p (GeV/c);(1/#beta)_{ratio};",1200,0,6,2500,0,3);
   hBetaRatiovsP_p = fs->make<TH2F>("BetaRatiovsP_p",";p (GeV/c);(1/#beta)_{ratio};",1200,0,6,2500,0,3);

   hJetsEtaVsPt = fs->make<TH2F>("JetsEtaVsPt",";#eta;p_{T} (GeV/c);",100,-5,5,200,0,200);

   hJetDrVsPt = fs->make<TH2F>("JetDrVsPt",";#Delta r;p_{T} (GeV/c);",100,0,2,100,0,20);
   hJetDrVsPt_pi = fs->make<TH2F>("JetDrVsPt_pi",";#Delta r;p_{T} (GeV/c);",100,0,2,100,0,20);
   hJetDrVsPt_k = fs->make<TH2F>("JetDrVsPt_k",";#Delta r;p_{T} (GeV/c);",100,0,2,100,0,20);
   hJetDrVsPt_p = fs->make<TH2F>("JetDrVsPt_p",";#Delta r;p_{T} (GeV/c);",100,0,2,100,0,20);
   hJetDrVsPt_notp = fs->make<TH2F>("JetDrVsPt_notp",";#Delta r;p_{T} (GeV/c);",100,0,2,100,0,20);
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenPIDMatchingJets::endJob()
{
}

double GenPIDMatchingJets::GetFlightDistance(double eta, double pt, double q)
{
      double L = -9999.0;
      if(fabs(eta)<1.6 && fabs(1-Bmag*Bmag*fabs(q)*fabs(q)*R*R/pt/pt*clight*clight/2.)<=1)
        L = pt/clight*cosh(eta)/Bmag/fabs(q)*acos(1-Bmag*Bmag*fabs(q)*fabs(q)*R*R/pt/pt*clight*clight/2.);
      else if(fabs(eta)>1.6 && fabs(eta)<3 && isETL)
        L = 3.04/tanh(fabs(eta));

      return L;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenPIDMatchingJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPIDMatchingJets);
