//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file B4bSteppingAction.cc
/// \brief Implementation of the B4bSteppingAction class

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"

#include "B4bSteppingAction.hh"
#include "B4DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "CaloID.h"  // including CaloID, CaloHit, CaloTree
#include "CaloHit.h"  // including CaloID, CaloHit, CaloTree
#include "CaloTree.h"

#include "TH1D.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bSteppingAction::B4bSteppingAction(B4bEventAction* eventAction,CaloTree* histo)
  : G4UserSteppingAction(),
    fEventAction(eventAction),
    `(histo)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bSteppingAction::~B4bSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bSteppingAction::UserSteppingAction(const G4Step* step)
{

   G4Track* track = step ->GetTrack();

// Collect energy and track length step by step

   //   === begin of checking optical photon ===
   static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

    G4StepPoint* endPoint   = step->GetPostStepPoint();
    G4StepPoint* startPoint = step->GetPreStepPoint();

    const G4DynamicParticle* theParticle = track->GetDynamicParticle();
    const G4ParticleDefinition* particleDef =
    theParticle->GetParticleDefinition();

   if(particleDef == opticalphoton)
   {
      // cout<<"this is an optical photon.  Need to stop it"<<endl;
      track->SetTrackStatus(fStopAndKill); 
      return;
   }
   //   === end of checking optical photon ===


  // get volume of the current step
  // auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
 
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
      
  G4TrackStatus tkstatus=step->GetTrack()->GetTrackStatus();

  const G4DynamicParticle* dynamicParticle= track ->GetDynamicParticle();
  G4int pdgcode=dynamicParticle->GetPDGcode();
  // G4int absPdgCode=abs(pdgcode);
  G4ParticleDefinition* particle = dynamicParticle->GetDefinition();
  G4String particleName= particle ->GetParticleName();
  G4double kinEnergy=dynamicParticle->GetKineticEnergy();

  // Analysis of physics process...
  // if(abs(pdgcode)==2112) {
  //     track->SetTrackStatus(fStopAndKill);
  //     return;
  //  }




  // if(tkstatus==fStopAndKill && absPdgCode>100 && track->GetTrackID()==1) {
  // if(tkstatus==fStopAndKill && absPdgCode>100) {
  // save secondaries from hadronic interactions and photon/lepton DIS...
  if(tkstatus==fStopAndKill) {
     // fEventAction->FillSecondaries(step);
  }

  if(charge == 0.0) return;

  int skdebug=0;
  if(skdebug>0) std::cout<<"skdebug (step)  1..."<<std::endl;

  //  energy deposit in cell...
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto depth = touchable->GetHistory()->GetDepth();
  auto thisPhysical = touchable->GetVolume(); // mother
  auto thisCopyNo = thisPhysical->GetCopyNo();
  auto thisName = thisPhysical->GetName();
  G4ThreeVector posA = step->GetPreStepPoint()->GetPosition();

  double birks=1.0;
  vector<double> ncer;

  int caloType=0;
  int fiberNumber=0;
  int holeNumber=0;
  int rodNumber=0;
  int layerNumber=0;
  // int holeReplicaNumber=0;
  // int rodReplicaNumber=0;
  // int  layerReplicaNumber=0;
  if(thisName.compare(0,3,"Rod")==0) {
     caloType=1;
     fiberNumber=0;
     holeNumber=0;
     rodNumber=touchable->GetCopyNumber(0);
     layerNumber=touchable->GetCopyNumber(1);
     // holeReplicaNumber=touchable->GetReplicaNumber(2);
     // rodReplicaNumber=touchable->GetReplicaNumber(3);
     // layerReplicaNumber=touchable->GetReplicaNumber(4);
  }


  if(skdebug>0)  std::cout<<"skdebug (step)  2..."<<std::endl;

  if(thisName.compare(0,18,"fiberCoreScintPhys")==0)  {
     caloType=2;
     birks=getBirk(step);
  }
  if(thisName.compare(0,18,"fiberCoreCherePhys")==0)  {
     caloType=3;
     ncer=UserCerenkov(step);   // cerenkov photons;
  }
  if(skdebug>0)std::cout<<"skdebug (step)  3..."<<std::endl;

  if(caloType==2 || caloType==3) {  
     fiberNumber=touchable->GetCopyNumber(1);
     holeNumber=touchable->GetCopyNumber(2);
     rodNumber=touchable->GetCopyNumber(3);
     layerNumber=touchable->GetCopyNumber(4);
     // holeReplicaNumber=touchable->GetReplicaNumber(2);
     // rodReplicaNumber=touchable->GetReplicaNumber(3);
     // layerReplicaNumber=touchable->GetReplicaNumber(4);
  }
    if(skdebug>0) std::cout<<"skdebug (step)  4..."<<std::endl;
     // std::cout<<"Stepping Action: depth "<<depth<<"  volume "<<thisName<<"  copy no "<<thisCopyNo;
     // std::cout<<" calotype "<<caloType;
     // std::cout<<" f "<<fiberNumber;
     // std::cout<<" h "<<holeNumber;
     // std::cout<<" r "<<rodNumber<<" lyr "<<layerNumber;
     //  std::cout<<" "<<std::endl;
     // std::cout<<" replicas "<<holeReplicaNumber<<"  "<<rodReplicaNumber<<"  "<<layerReplicaNumber<<std::endl;

     CaloID caloid(caloType,fiberNumber,rodNumber,layerNumber,posA.z(),track->GetGlobalTime());

     CaloHit aHit;
     aHit.caloid=caloid;
     aHit.x=posA.x()/10.0;   // in cm
     aHit.y=posA.y()/10.0;
     aHit.z=posA.z()/10.0;
     aHit.pid=pdgcode;
     aHit.trackid=track->GetTrackID();
     aHit.globaltime=track->GetGlobalTime();
     aHit.steplength=track->GetTrackLength();
     aHit.edep=edep/1000.;   //  in GeV
     aHit.edepbirk=edep*birks/1000.;
     if(ncer.size()>0) {
        aHit.ncer=ncer[2];  // use captured cherenkov light, instead of all, ncer[0]
        aHit.ncercap=ncer[2];
     }
    if(skdebug>0) std::cout<<"skdebug (step)  5..."<<std::endl;

     // aHit.print();

     hh->accumulateHits(aHit);

     // hh->histo1D["edepX"]->Fill(aHit.x/10.0,edep);
     // hh->histo1D["edepY"]->Fill(aHit.y/10.0,edep);
     // hh->histo1D["cerX"]->Fill(aHit.x/10.0,ncer[0]);
     // hh->histo1D["cerY"]->Fill(aHit.y/10.0,ncer[0]);

     // fEventAction->AccumulateCaloHits(aHit);
     // fEventAction->StepAnalysisSensor(step,ncer);   // analysis for the sensor volume.

  // fEventAction->StepAnalysis(step,ncer[0],ncer[1]);   // in original sim. Moved to above.
    if(skdebug>0) std::cout<<"skdebug (step)  6..."<<std::endl;

}  // end of B4bSteppingAction::UserSteppingAction.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


vector<double> B4bSteppingAction::UserCerenkov(const G4Step* step)
{
    G4int n_scint = 0;
    G4int n_cer   = 0;
    G4int n_cerhad   = 0;
    int nCERlocal=0;
    int nCERlocalElec=0;
    int nCERlocalCap=0;
    int nCERlocalElecCap=0;

  // cout<<"B4bSteppingAction::UserCerenkov is called..."<<endl;

  //  Code from examples/extended/optical/OpNovice2/src/SteppingAction.cc
  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  G4Track* track          = step->GetTrack();
  G4StepPoint* endPoint   = step->GetPostStepPoint();
  G4StepPoint* startPoint = step->GetPreStepPoint();

  const G4DynamicParticle* theParticle = track->GetDynamicParticle();
  const G4ParticleDefinition* particleDef =
    theParticle->GetParticleDefinition();
  G4int pdgcode=abs(theParticle->GetPDGcode());

  if(particleDef == opticalphoton)
  {
      // do noting for now
  }
  else
  {  // particle != opticalphoton
    // print how many Cerenkov and scint photons produced this step
    // this demonstrates use of GetNumPhotons()
    auto proc_man =
      track->GetDynamicParticle()->GetParticleDefinition()->GetProcessManager();
    G4ProcessVector* proc_vec = proc_man->GetPostStepProcessVector(typeDoIt);
    G4int n_proc              = proc_vec->entries();

    // G4int n_scint = 0;
    // G4int n_cer   = 0;
    for(G4int i = 0; i < n_proc; ++i)
    {
      G4String proc_name = (*proc_vec)[i]->GetProcessName();
      if(proc_name.compare("Cerenkov") == 0)
      {
        auto cer = (G4Cerenkov*) (*proc_vec)[i];
        n_cer    = cer->GetNumPhotons();
        if(pdgcode==11){ n_cerhad=n_cer;}
      }
      else if(proc_name.compare("Scintillation") == 0)
      {
        auto scint = (G4Scintillation*) (*proc_vec)[i];
        n_scint    = scint->GetNumPhotons();
      }
    }
    int fVerbose=-1;    //  set a value here for now...
    if(fVerbose > 0)
    {
      if(n_cer > 0 || n_scint > 0)
      {
        G4cout << "In this step, " << n_cer << " Cerenkov and " << n_scint
               << " scintillation photons were produced." << G4endl;
      }
    }
  }

    // loop over secondaries, create statistics
    const std::vector<const G4Track*>* secondaries =
      step->GetSecondaryInCurrentStep();

    for(auto sec : *secondaries)
    {
      if(sec->GetDynamicParticle()->GetParticleDefinition() == opticalphoton)
      {
        G4String creator_process = sec->GetCreatorProcess()->GetProcessName();
        if(creator_process.compare("Cerenkov") == 0)
        {
          nCERlocal=nCERlocal+1;
          G4double en = sec->GetKineticEnergy();
          double wavelength=1239.8*eV/en;
          G4ThreeVector pvec=sec->GetMomentumDirection();

          int capture=0;
          if(abs(pvec.theta())<0.336) capture=1;   // NA=sin(theta)=0.33

          hh->histo1D["cerWL"]->Fill(wavelength);
          if(capture==1) { 
              nCERlocalCap=nCERlocalCap+1;
              hh->histo1D["cerWLcaptured"]->Fill(wavelength);
          }

          if(pdgcode==11) {
              nCERlocalElec=nCERlocalElec+1;
              if(capture==1) {
                  nCERlocalElecCap=nCERlocalElecCap+1;
                   hh->histo1D["cerWLcapturedELEC"]->Fill(wavelength);
              }
          }
          // cout<<"cerenkov phton  en="<<en<<endl;
          // run->AddCerenkovEnergy(en);
          // run->AddCerenkov();
          // analysisMan->FillH1(1, en / eV);
        }
        else if(creator_process.compare("Scintillation") == 0)
        {
          G4double en = sec->GetKineticEnergy();
          // run->AddScintillationEnergy(en);
          // run->AddScintillation();
          // analysisMan->FillH1(2, en / eV);

          // G4double time = sec->GetGlobalTime();
          // analysisMan->FillH1(3, time / ns);
        }
      }
    }  //  end of for(auto sec : *secondaries)

    // double NCER=double(n_cer)/10000.0;
    vector<double> NCER;
    NCER.push_back(double(nCERlocal));
    NCER.push_back(double(nCERlocalElec));
    NCER.push_back(double(nCERlocalCap));
    NCER.push_back(double(nCERlocalElecCap));
    return NCER;
}

// ========================================================================================    
double B4bSteppingAction::getBirk(const G4Step* step){
   double weight = 1.0;

     double edep = step->GetTotalEnergyDeposit();
     double steplength = step->GetStepLength() / CLHEP::cm ;   //  convert from mm to cm.
     double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
     double density = step->GetPreStepPoint()->GetMaterial()->GetDensity() / (CLHEP::g/CLHEP::cm3);
     string materialName=step->GetPreStepPoint()->GetMaterial()->GetName();
     // std::cout<<"B4bEventAction::getBirk:   name="<<materialName<<std::endl;

/*
     if(materialName.compare(0,8,"G4_PbWO4")==0) {
        weight=getBirkL3(edep,steplength,charge,density);
        // std::cout<<"B4bEventAction::getBirk:   name="<<materialName<<"   weight="<<weight<<std::endl;
     }

     if(materialName.compare(0,14,"H_Scintillator")==0 || materialName.compare(0,26,"G4_PLASTIC_SC_VINYLTOLUENE")==0) {
        weight=getBirkHC(edep,steplength,charge,density);
        // std::cout<<"B4bEventAction::getBirk:   name="<<materialName<<"   weight="<<weight<<std::endl;
     }
*/

     if(materialName.compare(0,11,"Polystyrene")==0) {
        weight=getBirkHC(edep,steplength,charge,density);
        // std::cout<<"B4bEventAction::getBirk:   name="<<materialName<<"   weight="<<weight<<std::endl;
     }

   return weight;;
}

double B4bSteppingAction::getBirkHC(double dEStep, double step, double charge, double density) {
  double weight = 1.;
  if (charge != 0. && step > 0.) {
      double birkC1HC_ = 0.0052;
      double birkC2HC_ = 0.142;
      double birkC3HC_ = 1.75;
    double dedx = dEStep / step;
    double rkb = birkC1HC_ / density;
    double c = birkC2HC_ * rkb * rkb;
    if (std::abs(charge) >= 2.)
      rkb /= birkC3HC_;
    weight = 1. / (1. + rkb * dedx + c * dedx * dedx);
  }
  return weight;
}

double B4bSteppingAction::getBirkL3(double dEStep, double step, double charge, double density) {
  double weight = 1.;
  if (charge != 0. && step > 0.) {
       double birkC1EC_ = 0.03333;
       double birkSlopeEC_ = 0.253694;
       double birkCutEC_ = 0.1;
    double dedx = dEStep / step;
    double rkb = birkC1EC_ / density;
    if (dedx > 0) {
      weight = 1. - birkSlopeEC_ * log(rkb * dedx);
      if (weight < birkCutEC_)
        weight = birkCutEC_;
      else if (weight > 1.)
        weight = 1.;
    }
  }
  return weight;
}

