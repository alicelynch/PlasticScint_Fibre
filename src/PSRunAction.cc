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
// $Id: PSRunAction.cc,v 1.10 2006-06-29 17:54:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make this appear first!
#include "G4Timer.hh"

#include "PSRunAction.hh"

#include "G4Run.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PSRunAction::PSRunAction()
{
  timer = new G4Timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PSRunAction::~PSRunAction()
{
  delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSRunAction::BeginOfRunAction(const G4Run* aRun)
{
  tally = 0;
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  timer->Start();

  outFile = new std::ofstream("PlasticScint_Fibre.out");

	// determine time and date
	time_t ts;
	time(&ts);
	struct tm *date;
	date = gmtime(&ts);
	char date_out[25];
	strftime(date_out, 25, "%Y-%m-%dT%H:%M:%S", date);


	*outFile << "# ========================================================================" << std::endl;
	*outFile << "# Simulation started (GMT): " << date_out << std::endl;
	*outFile << "# ========================================================================" << std::endl;
	get detector parameters
	aRun->Print(outFile);
	*outFile << std::endl << std::endl << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PSRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();

  G4cout << "EndRunAction for " << aRun->GetNumberOfEvent()  << G4endl;
  G4cout << "Successful hist: " << tally << G4endl;
  G4cout<< "Time:" << *timer << G4endl;

  delete outFile;
  outFile = 0;
}

void PSRunAction::fillPerEvent(G4double energy)
{
  if (energy > 0) {
    ++tally;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
