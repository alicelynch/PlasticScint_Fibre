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
// $Id: PSPrimaryGeneratorAction.cc,v 1.6 2006-06-29 17:54:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PSPrimaryGeneratorAction.hh"
#include "PSPrimaryGeneratorMessenger.hh"

#include "G4GeneralParticleSource.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PSPrimaryGeneratorAction::PSPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4GeneralParticleSource(); //G4ParticleGun(n_particle);

  //create a messenger for this class
  gunMessenger = new PSPrimaryGeneratorMessenger(this);

  //default kinematic
  //
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition*  particle = particleTable->FindParticle("e+");
   timeConstant = 0.;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PSPrimaryGeneratorAction::~PSPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   static G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton"); // Or Alpha or electron Or Ion
   static G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
   particleGun->SetParticleDefinition(particle);
   particleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.75*cm));

   // Give particle a random momentum direction
   G4double Phi = CLHEP::twopi * G4UniformRand();
   G4double Theta = 0.5 * CLHEP::pi * G4UniformRand();
   G4ThreeVector v(cos(Theta), sin(Theta)*cos(Phi), sin(Theta)*sin(Phi));
   particleGun->SetParticleMomentumDirection(v);

   // If the particle is an optical photon, set its polarisation.
   if(particleGun->GetParticleDefinition()->GetParticleName()=="opticalphoton"){
     G4double angle = G4UniformRand() * 360.0*deg;
     
     SetOptPhotonPolar(angle);


     G4double time = -std::log(G4UniformRand())*timeConstant;
     particleGun->SetParticleTime(time);
     }



  particleGun->SetParticleEnergy(1.0*MeV);//(2.88*eV);

  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSPrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (particleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }

 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton);
 G4double modul2       = product*product;

 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);

 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 particleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
