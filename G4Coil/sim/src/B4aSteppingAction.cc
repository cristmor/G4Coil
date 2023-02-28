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
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4OpticalPhoton.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include <iostream>
#include <vector>

// #include "TH1D.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step ->GetTrack();

    //   === begin of checking optical photon ===
   static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

    const G4DynamicParticle* theParticle = track->GetDynamicParticle();
    const G4ParticleDefinition* particleDef =
    theParticle->GetParticleDefinition();

   if(particleDef == opticalphoton)
   {  
      // std::cout<<"this is an optical photon.  Need to stop it"<<std::endl;
      PrintStep(step);
      // track->SetTrackStatus(fStopAndKill); 
      return;
   }
   //   === end of checking optical photon ===

  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto thisPhysical = touchable->GetVolume(); // mother
  auto thisName = thisPhysical->GetName();

  if(thisName.compare("cylinderPV")==0){
    std::vector <double> cerenkov = UserCerenkov(step);
    // if(cerenkov[0]){
    //   std::cout<<"Cerenkov_: ";
    //   std::cout<<"nCERlocal: "<<cerenkov[0]
    //       <<", nCERlocalElec: "<<cerenkov[1]<<", nCERlocalCap: "<<cerenkov[2]
    //       <<",nCERlocalElecCap: "<<cerenkov[3]<<" ";

    //   std::cout<<"\n"<<std::endl;;
    // }
  }
}

// ============================================================================
void B4aSteppingAction::PrintStep(const G4Step* step) {
  //std::cout<<"Print step 1 ..."<<std::endl;
  if(step->GetTotalEnergyDeposit()<1.0E-10) return;
  // auto charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
  auto preStepPoint = step->GetPreStepPoint();
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto depth = touchable->GetHistory()->GetDepth();
  auto thisPhysical = touchable->GetVolume(); // mother
  auto thisCopyNo = thisPhysical->GetCopyNo();
  auto thisName = thisPhysical->GetName();
  auto worldPos = preStepPoint->GetPosition();
  auto localPos = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
  auto just_enterd=preStepPoint->GetStepStatus();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  G4StepPoint* point1 = step->GetPreStepPoint();
  G4ThreeVector pos1 = point1->GetPosition();
  double xx=pos1.x();
  double yy=pos1.y();
  double zz=pos1.z();
  G4Track* track = step ->GetTrack();
  G4int steps = track ->GetCurrentStepNumber();

  const G4DynamicParticle* dynamicParticle= track ->GetDynamicParticle();
  G4ParticleDefinition* particle = dynamicParticle->GetDefinition();
  G4String particleName= particle ->GetParticleName();
  G4double kinEnergy=dynamicParticle->GetKineticEnergy();

  std::cout<<" tk="<<track->GetTrackID();
  std::cout<<" step="<<steps;;
  std::cout<<" "<<particleName;
  std::cout<<" ke="<<kinEnergy;
  std::cout<<" edep="<<edep;
  std::cout<<" xyz=("<<xx<<","<<yy<<","<<zz<<")";
  std::cout<<" depth "<<depth;
  std::cout<<" "<<thisName;
  std::cout<<" "<<thisCopyNo;;
  //std::cout<<" ( "<<motherName;
  //std::cout<<"  "<<motherCopyNo<<" )";;
  std::cout<<" status="<<just_enterd;
  std::cout<<" fGeomBoundary="<<fGeomBoundary;
  std::cout<<" "<<std::endl;
  return ;
}

std::vector <double> B4aSteppingAction::UserCerenkov(const G4Step* step)
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
  static G4ParticleDefinition* opticalphoton = G4OpticalPhoton::OpticalPhotonDefinition();

  G4Track* track          = step->GetTrack();
  G4StepPoint* endPoint   = step->GetPostStepPoint();
  G4StepPoint* startPoint = step->GetPreStepPoint();

  const G4DynamicParticle* theParticle = track->GetDynamicParticle();
  const G4ParticleDefinition* particleDef =
    theParticle->GetParticleDefinition();
  G4int pdgcode=abs(theParticle->GetPDGcode());
  if(particleDef == opticalphoton)
  {
      std::cout<<particleDef<<std::endl;
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
        // std::cout<<proc_name<<std::endl;
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
      // std::cout<<sec<<std::endl;
      if(sec->GetDynamicParticle()->GetParticleDefinition() == opticalphoton)
      {
        // PrintStep(step);
        // PrintStepSecondary(sec);
        G4String creator_process = sec->GetCreatorProcess()->GetProcessName();
        // std::cout<<"creator_process:"<<creator_process<<std::endl;
        if(creator_process.compare("Cerenkov") == 0)
        {
          nCERlocal=nCERlocal+1;
          G4double en = sec->GetKineticEnergy();
          double wavelength=1239.8*eV/en;
          G4ThreeVector pvec=sec->GetMomentumDirection();

          int capture=0;
          if(abs(pvec.theta())<0.336) capture=1;   // NA=sin(theta)=0.33

          // hh->histo1D["cerWL"]->Fill(wavelength);
          if(capture==1) { 
              nCERlocalCap=nCERlocalCap+1;
              // std::cout<<"Capture ="<<wavelength<<std::endl;
              // hh->histo1D["cerWLcaptured"]->Fill(wavelength);
          }

          if(pdgcode==11) {
              nCERlocalElec=nCERlocalElec+1;
              if(capture==1) {
                  nCERlocalElecCap=nCERlocalElecCap+1;
                  // std::cout<<"Electron Capture ="<<wavelength<<std::endl;
                  //  hh->histo1D["cerWLcapturedELEC"]->Fill(wavelength);
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
    std::vector<double> NCER;
    NCER.push_back(double(nCERlocal));
    NCER.push_back(double(nCERlocalElec));
    NCER.push_back(double(nCERlocalCap));
    NCER.push_back(double(nCERlocalElecCap));
    return NCER;
}

void B4aSteppingAction::PrintStepSecondary(const G4Track* secondary){
  G4int TrackID = secondary->GetTrackID();
  G4int steps = secondary ->GetCurrentStepNumber();

  const G4DynamicParticle* dynamicParticle= secondary ->GetDynamicParticle();
  G4ParticleDefinition* particle = dynamicParticle->GetDefinition();
  G4String particleName= particle ->GetParticleName();

  G4String CreationProcess = secondary->GetCreatorProcess()->GetProcessName();
  G4double KineticEnergy = secondary->GetKineticEnergy();
  G4ThreeVector Positison = secondary->GetPosition();
  G4ThreeVector MomentumDirection = secondary->GetMomentumDirection();
  G4cout<<" tk=" << TrackID
        <<"  ,steps=" << steps
        <<"  ,particleName=" << particleName
        <<"  ,CreationProcess=" << CreationProcess
        <<"  ,KineticEnergy=" <<  KineticEnergy
        <<"  ,Positison=" << Positison
        <<"  ,MomentumDirection=" << MomentumDirection
        << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
