// -*- C++ -*-
//
// Package:    MTDHIAnalysis/GenPIDMatching
// Class:      GenPIDMatching
//
/**\class GenPIDMatching GenPIDMatching.cc MTDHIAnalysis/GenPIDMatching/plugins/GenPIDMatching.cc

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

class GenPIDMatching : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenPIDMatching(const edm::ParameterSet&);
      ~GenPIDMatching();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      double GetFlightDistance(double eta, double pt, double q);
 
      // ----------member data ---------------------------
//      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;

    edm::Service<TFileService> fs;

    TH2F* hBetavsP;
    TH2F* hBetaRatiovsP_pi;
    TH2F* hBetaRatiovsP_K;
    TH2F* hBetaRatiovsP_p;

    TH2F* hDMassVsPt;
    TH2F* hDMassVsY;
    TH2F* hDYVsPt;
    TH2F* hDMassVsPtPID;
    TH2F* hDMassVsYPID;
    TH2F* hDYVsPtPID;

    TH2F* hLamCMassVsPt;
    TH2F* hLamCMassVsY;
    TH2F* hLamCYVsPt;
    TH2F* hLamCMassVsPtPID;
    TH2F* hLamCMassVsYPID;
    TH2F* hLamCYVsPtPID;

    bool isETL;
    bool isLambdaC;

    double sigmaT_;
    double nSigmaT_;
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
GenPIDMatching::GenPIDMatching(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

  isETL = iConfig.getUntrackedParameter<bool>("isETL");
  isLambdaC = iConfig.getUntrackedParameter<bool>("isLambdaC");

  sigmaT_ = iConfig.getUntrackedParameter<double>("sigmaT");
  nSigmaT_ = iConfig.getUntrackedParameter<double>("nSigmaT");
}


GenPIDMatching::~GenPIDMatching()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenPIDMatching::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_,genpars);

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
      if((fabs(eta)<3 && fabs(beta_measured-func_pi->Eval(p))/sigma_beta<nSigmaT_) || fabs(eta)>3) isPion = true;
      if((fabs(eta)<3 && fabs(beta_measured-func_K->Eval(p))/sigma_beta<nSigmaT_) || fabs(eta)>3) isKaon = true;
      if((fabs(eta)<3 && fabs(beta_measured-func_p->Eval(p))/sigma_beta<nSigmaT_) || fabs(eta)>3) isProton = true;

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

    for(unsigned int i=0;i<pVect_pos.size();i++)
    {
      TLorentzVector p1 = pVect_pos[i];
      TLorentzVector p1_pi(p1.Px(),p1.Py(),p1.Pz(),sqrt(p1.P()*p1.P()+massPi*massPi));
      TLorentzVector p1_K(p1.Px(),p1.Py(),p1.Pz(),sqrt(p1.P()*p1.P()+massK*massK));
      TLorentzVector p1_p(p1.Px(),p1.Py(),p1.Pz(),sqrt(p1.P()*p1.P()+massP*massP));

      bool isPion1 = isPionVect_pos[i];
      bool isKaon1 = isKaonVect_pos[i];
      bool isProton1 = isProtonVect_pos[i];
      for(unsigned int j=0;j<pVect_neg.size();j++)
      {
        TLorentzVector p2 = pVect_neg[j];
        TLorentzVector p2_pi(p2.Px(),p2.Py(),p2.Pz(),sqrt(p2.P()*p2.P()+massPi*massPi));
        TLorentzVector p2_K(p2.Px(),p2.Py(),p2.Pz(),sqrt(p2.P()*p2.P()+massK*massK));
        TLorentzVector p2_p(p2.Px(),p2.Py(),p2.Pz(),sqrt(p2.P()*p2.P()+massP*massP));

        bool isPion2 = isPionVect_neg[j];
        bool isKaon2 = isKaonVect_neg[j];
        bool isProton2 = isProtonVect_neg[j];

        TLorentzVector dCand1 = p1_pi + p2_K;
        TLorentzVector dCand2 = p2_pi + p1_K;  

        double massdiff1 = fabs(dCand1.Mag()-1.86);
        double massdiff2 = fabs(dCand2.Mag()-1.86);
        double cand1y = dCand1.Rapidity();
        double cand2y = dCand2.Rapidity();
        double cand1pt = dCand1.Pt();
        double cand2pt = dCand2.Pt();

//        hDMassVsPt->Fill(dCand1.Mag(),dCand1.Pt());
//        hDMassVsPt->Fill(dCand2.Mag(),dCand2.Pt());
        if(massdiff1<0.2) hDYVsPt->Fill(cand1y,cand1pt);
        if(massdiff2<0.2) hDYVsPt->Fill(cand2y,cand2pt);
//        hDMassVsY->Fill(dCand1.Mag(),dCand1.Rapidity());
//        hDMassVsY->Fill(dCand2.Mag(),dCand2.Rapidity());
        
        if(isKaon1 && isPion2)
        { 
//          hDMassVsPtPID->Fill(dCand2.Mag(),dCand2.Pt());
//          hDMassVsYPID->Fill(dCand2.Mag(),dCand2.Rapidity());
          if(massdiff2<0.2) hDYVsPtPID->Fill(cand2y,cand2pt);
        }

        if(isKaon2 && isPion1)
        {
//          hDMassVsPtPID->Fill(dCand1.Mag(),dCand1.Pt());
//          hDMassVsYPID->Fill(dCand1.Mag(),dCand1.Rapidity());
          if(massdiff1<0.2) hDYVsPtPID->Fill(cand1y,cand1pt);
        } 

       if(!isLambdaC) continue;

        for(unsigned int k=j+1;k<pVect_neg.size();k++)
        {
          TLorentzVector p3 = pVect_neg[k];

//          TLorentzVector p3_pi(p3.Px(),p3.Py(),p3.Pz(),sqrt(p3.P()*p3.P()+massPi*massPi));
          TLorentzVector p3_K(p3.Px(),p3.Py(),p3.Pz(),sqrt(p3.P()*p3.P()+massK*massK));
          TLorentzVector p3_p(p3.Px(),p3.Py(),p3.Pz(),sqrt(p3.P()*p3.P()+massP*massP));

          TLorentzVector lamCCand1 = p1_pi + p2_K + p3_p;
          TLorentzVector lamCCand2 = p1_pi + p3_K + p2_p;
  
          double lmassdiff1 = fabs(lamCCand1.Mag()-2.286);
          double lmassdiff2 = fabs(lamCCand2.Mag()-2.286);

          if(lmassdiff1>0.3 && lmassdiff2>0.3) continue;

          bool isKaon3 = isKaonVect_neg[k];
          bool isProton3 = isProtonVect_neg[k];

          double lcand1y = lamCCand1.Rapidity();
          double lcand2y = lamCCand2.Rapidity();
          double lcand1pt = lamCCand1.Pt();
          double lcand2pt = lamCCand2.Pt();

//          hLamCMassVsPt->Fill(lamCCand1.Mag(),lamCCand1.Pt());
//          hLamCMassVsY->Fill(lamCCand1.Mag(),lamCCand1.Rapidity());
          if(lmassdiff1<0.3) hLamCYVsPt->Fill(lcand1y,lcand1pt);

//          hLamCMassVsPt->Fill(lamCCand2.Mag(),lamCCand2.Pt());
//          hLamCMassVsY->Fill(lamCCand2.Mag(),lamCCand2.Rapidity());
          if(lmassdiff2<0.3) hLamCYVsPt->Fill(lcand2y,lcand2pt);

          if(isKaon2 && isPion1 && isProton3)
          {
//            hLamCMassVsPtPID->Fill(lamCCand1.Mag(),lamCCand1.Pt());
//            hLamCMassVsYPID->Fill(lamCCand1.Mag(),lamCCand1.Rapidity());
            if(lmassdiff1<0.3) hLamCYVsPtPID->Fill(lcand1y,lcand1pt);
          }

          if(isKaon3 && isPion1 && isProton2)
          {
//            hLamCMassVsPtPID->Fill(lamCCand2.Mag(),lamCCand2.Pt());
//            hLamCMassVsYPID->Fill(lamCCand2.Mag(),lamCCand2.Rapidity());
            if(lmassdiff2<0.3) hLamCYVsPtPID->Fill(lcand2y,lcand2pt);
          }
        } 

        for(unsigned int k=i+1;k<pVect_pos.size();k++)
        {
          TLorentzVector p3 = pVect_pos[k];

//          TLorentzVector p3_pi(p3.Px(),p3.Py(),p3.Pz(),sqrt(p3.P()*p3.P()+massPi*massPi));
          TLorentzVector p3_K(p3.Px(),p3.Py(),p3.Pz(),sqrt(p3.P()*p3.P()+massK*massK));
          TLorentzVector p3_p(p3.Px(),p3.Py(),p3.Pz(),sqrt(p3.P()*p3.P()+massP*massP));

          TLorentzVector lamCCand1 = p2_pi + p1_K + p3_p;
          TLorentzVector lamCCand2 = p2_pi + p3_K + p1_p;

          double lmassdiff1 = fabs(lamCCand1.Mag()-2.286);
          double lmassdiff2 = fabs(lamCCand2.Mag()-2.286);

          if(lmassdiff1>0.3 && lmassdiff2>0.3) continue;

          bool isKaon3 = isKaonVect_pos[k];
          bool isProton3 = isProtonVect_pos[k];

          double lcand1y = lamCCand1.Rapidity();
          double lcand2y = lamCCand2.Rapidity();
          double lcand1pt = lamCCand1.Pt();
          double lcand2pt = lamCCand2.Pt();

//          hLamCMassVsPt->Fill(lamCCand1.Mag(),lamCCand1.Pt());
//          hLamCMassVsY->Fill(lamCCand1.Mag(),lamCCand1.Rapidity());
          if(lmassdiff1<0.3) hLamCYVsPt->Fill(lcand1y,lcand1pt);

//          hLamCMassVsPt->Fill(lamCCand2.Mag(),lamCCand2.Pt());
//          hLamCMassVsY->Fill(lamCCand2.Mag(),lamCCand2.Rapidity());
          if(lmassdiff2<0.3) hLamCYVsPt->Fill(lcand2y,lcand2pt);

          if(isKaon1 && isPion2 && isProton3)
          {
//            hLamCMassVsPtPID->Fill(lamCCand1.Mag(),lamCCand1.Pt());
//            hLamCMassVsYPID->Fill(lamCCand1.Mag(),lamCCand1.Rapidity());
            if(lmassdiff1<0.3) hLamCYVsPtPID->Fill(lcand1y,lcand1pt);
          }

          if(isKaon3 && isPion2 && isProton1)
          {
//            hLamCMassVsPtPID->Fill(lamCCand2.Mag(),lamCCand2.Pt());
//            hLamCMassVsYPID->Fill(lamCCand2.Mag(),lamCCand2.Rapidity());
            if(lmassdiff2<0.3) hLamCYVsPtPID->Fill(lcand2y,lcand2pt);
          }
        }

      }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
GenPIDMatching::beginJob()
{
   hBetavsP = fs->make<TH2F>("BetavsP",";p (GeV/c);1/#beta;",1200,0,6,2500,0,5);
   hBetaRatiovsP_pi = fs->make<TH2F>("BetaRatiovsP_pi",";p (GeV/c);(1/#beta)_{ratio};",1200,0,6,2500,0,3);
   hBetaRatiovsP_K = fs->make<TH2F>("BetaRatiovsP_K",";p (GeV/c);(1/#beta)_{ratio};",1200,0,6,2500,0,3);
   hBetaRatiovsP_p = fs->make<TH2F>("BetaRatiovsP_p",";p (GeV/c);(1/#beta)_{ratio};",1200,0,6,2500,0,3);

   hDMassVsPt = fs->make<TH2F>("DMassVsPt",";Mass (GeV);p_{T} (GeV/c);",100,1.66,2.06,200,0,20);
   hDMassVsY = fs->make<TH2F>("DMassVsY",";Mass (GeV);Y;",100,1.66,2.06,80,-4,4);
   hDYVsPt = fs->make<TH2F>("DEtaVsPt",";#eta;p_{T} (GeV/c);",80,-4,4,200,0,20);
   hDMassVsPtPID = fs->make<TH2F>("DMassVsPtPID",";Mass (GeV);p_{T} (GeV/c);",100,1.66,2.06,200,0,20);
   hDMassVsYPID = fs->make<TH2F>("DMassVsYPID",";Mass (GeV);Y;",100,1.66,2.06,80,-4,4);
   hDYVsPtPID = fs->make<TH2F>("DEtaVsPtPID",";#eta;p_{T} (GeV/c);",80,-4,4,200,0,20);

   hLamCMassVsPt = fs->make<TH2F>("LamCMassVsPt",";Mass (GeV);p_{T} (GeV/c);",100,2.286-0.3,2.286+0.3,200,0,20);
   hLamCMassVsY = fs->make<TH2F>("LamCMassVsY",";Mass (GeV);Y;",100,2.286-0.3,2.286+0.3,80,-4,4);
   hLamCYVsPt = fs->make<TH2F>("LamCEtaVsPt",";#eta;p_{T} (GeV/c);",80,-4,4,200,0,20);
   hLamCMassVsPtPID = fs->make<TH2F>("LamCMassVsPtPID",";Mass (GeV);p_{T} (GeV/c);",100,2.286-0.3,2.286+0.3,200,0,20);
   hLamCMassVsYPID = fs->make<TH2F>("LamCMassVsYPID",";Mass (GeV);Y;",100,2.286-0.3,2.286+0.3,80,-4,4);
   hLamCYVsPtPID = fs->make<TH2F>("LamCEtaVsPtPID",";#eta;p_{T} (GeV/c);",80,-4,4,200,0,20);
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenPIDMatching::endJob()
{
}

double GenPIDMatching::GetFlightDistance(double eta, double pt, double q)
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
GenPIDMatching::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(GenPIDMatching);
