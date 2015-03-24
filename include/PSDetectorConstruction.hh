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
// $Id: PSDetectorConstruction.hh,v 1.5 2006-06-29 17:53:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PSDetectorConstruction_h
#define PSDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
class G4Box;
class G4Tubs;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PSDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    PSDetectorConstruction();
   ~PSDetectorConstruction();
    
  public:
    G4VPhysicalVolume* Construct();

  private:
    void setDefaultDimensions();
    void defineMaterials();
    void defineSurfaces();
    G4VPhysicalVolume* ConstructDetector();
    
    G4double expHall_x;
    G4double expHall_y;
    G4double expHall_z;
    G4double fScintVol_x;
    G4double fScintVol_y;
    G4double fScintVol_z ;
    G4double fPMTDiameter;
    G4double fWindowThickness;
    G4double fEnvelopeThickness;
    G4double fGreaseThickness;
    G4double fAirGap ;
    G4double fPTFEThickness; 
    G4double fCathodeThickness;
   
    G4double           wlsfiberRX;
    G4double           wlsfiberRY;
    G4double           wlsfiberZ;

    G4double           cladRX;
    G4double           cladRY;
    G4double           cladZ;

    G4double           clrfiberHalfL;
    G4double           clrfiberZ;
    
    G4double           coupleRX;
    G4double           coupleRY;
    G4double           coupleZ;
 
    G4double           mirrorRmax;
    G4double           mirrorZ;
    G4bool             mirrorToggle;

    G4double wlsfiberOrigin;
    G4double coupleOrigin;

    G4double mirrorPolish;
    G4double mirrorReflectivity;
    G4double XYRatio;
    G4double holeRadius;
    G4double holeLength;
    G4double surfaceRoughness;
    G4double mirrorOrigin;

    G4Box*                   expHall_box;
    G4LogicalVolume*         expHall_log;
    G4VPhysicalVolume*       expHall_phys;
    G4Box*                   scint_vol;
    G4LogicalVolume*         scint_log;
    G4VPhysicalVolume*       scint_phys;
    G4Tubs*                  solidGrease;
    G4LogicalVolume*         logicGrease;
    G4VPhysicalVolume*       physiGrease;
    G4VSolid*                solidTeflon;
    G4LogicalVolume*         logicTeflon;
    G4VPhysicalVolume*       physiTeflon;
    G4Tubs*                  solidWindow;
    G4LogicalVolume*         logicWindow;
    G4VPhysicalVolume*       physiWindow;
    G4Tubs*                  solidCathode;
    G4LogicalVolume*         logicCathode;
    G4VPhysicalVolume*       physiCathode;
    G4Tubs*                  solidHole;
    G4LogicalVolume*         logicHole;
    G4VPhysicalVolume*       physiHole;
 
    void ConstructFiber();
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PSDetectorConstruction_h*/
