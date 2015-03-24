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
// $Id: PSDetectorConstruction.cc,v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PSDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "CathodeSD.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SDManager.hh"
#include "G4SurfaceProperty.hh"
#include "G4UImanager.hh"
#include "G4OpticalSurface.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4NistManager.hh"
#include "G4UImanager.hh"
#include "G4UserLimits.hh"
#include "G4String.hh"
#include "G4GeometryManager.hh"


PSDetectorConstruction::PSDetectorConstruction() : expHall_box(0),expHall_log(0),expHall_phys(0),  scint_vol(0),   scint_log(0),  scint_phys(0),  solidGrease(0),   logicGrease(0),  physiGrease(0),  solidTeflon(0),   logicTeflon(0),  physiTeflon(0),  solidWindow(0),  logicWindow(0),  physiWindow(0),  solidCathode(0),  logicCathode(0),  physiCathode(0)
{
  setDefaultDimensions();

  // set default lengths/sizes/ reflectivities
}


PSDetectorConstruction::~PSDetectorConstruction(){}

G4VPhysicalVolume* PSDetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(10.*m);
  setDefaultDimensions();
  defineMaterials();
  defineSurfaces();
  G4VPhysicalVolume* detector = ConstructDetector();
  return detector;
}


void PSDetectorConstruction::setDefaultDimensions()
{
  expHall_x = expHall_y = expHall_z = 0.10*m;
  fScintVol_x =  8.75*mm;
  fScintVol_y = 10.00*mm;
  fScintVol_z =  5.025*mm;  // HALF LENGTHS
  fPMTDiameter = 2.82*cm;
  fWindowThickness = 2.25*mm
  fGreaseThickness = 0.1*mm;
  fAirGap = 0.1*mm;
  fPTFEThickness = 0.2*mm;
  fCathodeThickness = 1.0*mm;


  //Fiber:

  surfaceRoughness = 1;

  mirrorToggle = true;
  mirrorPolish = 1.;
  mirrorReflectivity = 1.;

  XYRatio = 1.0;

  wlsfiberZ     = fScintVol_z + 1.0*cm;   // half length
  wlsfiberRY  = 0.5*mm;
  wlsfiberOrigin = fScintVol_z + 1.0*cm;

  //clrfiberZ  = mppcZ + 10.*nm ;  // clear fibre 10nm longer in z than the mppc. made of air.
  mirrorZ    = 0.1*mm;          // a mirror at end of fiber? or length 0.1mm

  holeRadius       = wlsfiberRY +0.1*mm;
  holeLength       = 2*fScintVol_z;       // full length

  wlsfiberRX  = XYRatio * wlsfiberRY;
  cladRX = wlsfiberRX + 0.03*wlsfiberRX;
  cladRY = wlsfiberRY + 0.03*wlsfiberRY;
  cladZ  = wlsfiberZ;

  mirrorOrigin = wlsfiberOrigin - wlsfiberZ - mirrorZ/2;
}


void PSDetectorConstruction::defineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();

  //------------------------------------
  //------------- Elements -------------
  //------------------------------------


  G4Element *N = man->FindOrBuildElement("N");
  G4Element *C = man->FindOrBuildElement("C");
  G4Element *H = man->FindOrBuildElement("H");
  G4Element *Si = man->FindOrBuildElement("Si");
  G4Element *O = man->FindOrBuildElement("O");
  G4Element *B = man->FindOrBuildElement("B");
  G4Element *Na = man->FindOrBuildElement("Na");
  G4Element *Al = man->FindOrBuildElement("Al");
  G4Element *K = man->FindOrBuildElement("K");
  G4Element* F = man->FindOrBuildElement("F");
  G4Element *Sb = man->FindOrBuildElement("Sb");
  G4Element *Rb = man->FindOrBuildElement("Rb");
  G4Element *Cs = man->FindOrBuildElement("Cs");

  G4double density;
  G4int nElem, nAtoms;

  const G4int nEntries = 50;

  // Photon Energies for Spectral Response of Materials
  G4double PhotonEnergy[nEntries] =
    {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
     2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
     2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
     2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
     2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
     2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
     2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
     3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
     3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
     3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

  // -----------------------------------------
  // --------------Aluminium------------------
  // -----------------------------------------

  density = 2.700*g/cm3;
  G4double z = 13.0;
  G4double a = 26.98*g/mole;
  G4Material* Alum = new G4Material( "Aluminium", z, a, density);

  const G4int num = 2;
  G4double Ephoton[num] = {1.65*eV, 12.4*eV};
  G4double  AlAbs[num]= {710*nm, 11.4*nm};
  G4double mirror_RIND[num] = {100,100};

  G4MaterialPropertiesTable *mirror_prop = new G4MaterialPropertiesTable();
  mirror_prop->AddProperty("RINDEX", Ephoton, mirror_RIND, num)->SetSpline(true);
  Alum->SetMaterialPropertiesTable(mirror_prop);


  //-------------------------------------------
  //---------------------Air-------------------
  //-------------------------------------------

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nElem=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  G4double RefractiveIndex[nEntries] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, nEntries);

  Air->SetMaterialPropertiesTable(myMPT2);


  //------------------------------------------------------
  //------------- Plastic Scintillator BC408--------------
  //------------------------------------------------------

  density =1.032*g/cm3;
  G4Material*  BC408 = new G4Material("BC408", density, nElem=2);
  BC408->AddElement(H, nAtoms=521);
  BC408->AddElement(C, nAtoms=474);

  const G4int nEntries2= 35;
  G4double PhotonEnergy2[nEntries2]= { 3.55, 3.47, 3.393, 3.315, 3.256, 3.195, 3.141, 3.103, 3.077, 3.056, 3.031, 3.011, 2.997, 2.983, 2.965, 2.956, 2.939, 2.921, 2.906, 2.890, 2.867, 2.852, 2.834, 2.83, 2.81, 2.796, 2.766, 2.72, 2.686, 2.653, 2.60, 2.55, 2.49, 2.38, 2.25  };
  G4double RINDEX_BC408[nEntries2] ={ 1.58, 1.58, 1.58, 1.58, 1.58,1.58, 1.58, 1.58, 1.58, 1.58,1.58, 1.58,1.58,1.58,1.58, 1.58, 1.58, 1.58, 1.58, 1.58,1.58, 1.58, 1.58, 1.58, 1.58,1.58, 1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58 };
  G4double ABSORPTION_BC408[nEntries2] ={ 210*cm, 210*cm, 210*cm, 210*cm, 210*cm,210*cm, 210*cm, 210*cm, 210*cm, 210*cm,210*cm, 210*cm , 210*cm, 210*cm , 210*cm ,210*cm, 210*cm, 210*cm, 210*cm, 210*cm,210*cm, 210*cm, 210*cm, 210*cm, 210*cm, 210*cm, 210*cm , 210*cm, 210*cm , 210*cm, 210*cm, 210*cm , 210*cm, 210*cm , 210*cm };
  G4double SCINTILLATION_BC408[nEntries2] = {0,0.017,0.03,0.049,0.071,0.101,0.152,0.220,0.288,0.374,0.478,0.585,0.668,0.735,0.804,0.846,0.897,0.947,0.984, 1.0, 0.981,0.945, 0.877,0.796,0.698,0.614,0.534,0.451,0.366,0.278,0.191, 0.133, 0.076,0.025,0};

  G4MaterialPropertiesTable *BC408_mt = new G4MaterialPropertiesTable();
  BC408_mt->AddProperty("RINDEX", PhotonEnergy2, RINDEX_BC408, nEntries2);
  BC408_mt->AddProperty("ABSLENGTH", PhotonEnergy2, ABSORPTION_BC408, nEntries2);
  BC408_mt->AddProperty("FASTCOMPONENT", PhotonEnergy2, SCINTILLATION_BC408,nEntries2);
  BC408_mt->AddConstProperty("SCINTILLATIONYIELD",500./MeV);
  BC408_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  BC408_mt->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  BC408_mt->AddConstProperty("YIELDRATIO",1.);
  // what about alpha quenching factor?
  BC408->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  BC408->SetMaterialPropertiesTable(BC408_mt);

  //--------------------------------------------------------------
  //------------------- Polydimethylsiloxane (Grease)-------------
  //--------------------------------------------------------------

  density = 0.97*g/cm3;
  G4Material* Polydimethylsiloxane = new G4Material("Polydimethylsiloxane", density, nElem=4, kStateLiquid);
  Polydimethylsiloxane->AddElement(Si,nAtoms=1);
  Polydimethylsiloxane->AddElement(O, nAtoms= 1);
  Polydimethylsiloxane->AddElement(C, nAtoms=2);
  Polydimethylsiloxane->AddElement(H, nAtoms=6);

  // ------------ Generate & Add Material Properties Table ------------
  G4MaterialPropertiesTable* polydimethylsiloxaneprop = new G4MaterialPropertiesTable();
  const G4int numentriespolydimethylsiloxane = 3;
  G4double polydimethylsiloxaneenergy[numentriespolydimethylsiloxane] = {1.2*eV, 3.1*eV, 6.5*eV};
  G4double polydimethylsiloxaneabsorp[numentriespolydimethylsiloxane] = {10.*cm, 10.*cm, 10.*cm};
  G4double polydimethylsiloxanerindex[numentriespolydimethylsiloxane] = {1.4, 1.4, 1.4};
  polydimethylsiloxaneprop->AddProperty("ABSLENGTH", polydimethylsiloxaneenergy, polydimethylsiloxaneabsorp, numentriespolydimethylsiloxane);
  polydimethylsiloxaneprop->AddProperty("RINDEX", polydimethylsiloxaneenergy, polydimethylsiloxanerindex, numentriespolydimethylsiloxane);
  Polydimethylsiloxane->SetMaterialPropertiesTable(polydimethylsiloxaneprop);

  //----------------------------------------------------------
  // ------------PMT Borosilicate glass-----------------------
  //----------------------------------------------------------

  density = 2.23*g/cm3;
  G4Material* BoroSiGlass = new G4Material("BoroSiGlass", density, nElem = 6, kStateSolid);
  BoroSiGlass->AddElement(B,0.040064);
  BoroSiGlass->AddElement(O,0.539562);
  BoroSiGlass->AddElement(Na,0.028191);
  BoroSiGlass->AddElement(Al,0.011644);
  BoroSiGlass->AddElement(Si,0.377220);
  BoroSiGlass->AddElement(K,0.003321);
  // by relative atomic abundance in material

  // ------------ Generate & Add Material Properties Table ------------
  G4MaterialPropertiesTable* BoroSiGlassProp = new G4MaterialPropertiesTable();
  const  G4int nGlass = 4;
  G4double BoroSiGlassEnergy[nGlass]= {3.1*eV, 3.87*eV, 4.13*eV, 4.96*eV} ;// {1.2*eV, 3.1*eV, 6.5*eV};
  G4double BoroSiGlassRindex[nGlass]= {1.49, 1.49, 1.49, 1.49};
  BoroSiGlassProp->AddProperty("RINDEX", BoroSiGlassEnergy, BoroSiGlassRindex, nGlass);
  BoroSiGlass->SetMaterialPropertiesTable(BoroSiGlassProp);

  //----------------------------------------------------------
  //------------- Bialkali Cathode --------------------
  //----------------------------------------------------------

  G4Material* BialkaliCathode = new G4Material("BialkaliCathode", 3*g/cm3, 3, kStateSolid);
  BialkaliCathode->AddElement(Sb, 1);
  BialkaliCathode->AddElement(Rb, 1);
  BialkaliCathode->AddElement(Cs, 1);
  G4MaterialPropertiesTable* BialkaliCathodeProp = new G4MaterialPropertiesTable();
  const G4int numentriesbialkalicath = 2;
  G4double bialkalicathodeenergy[numentriesbialkalicath] = {1.2*eV, 6.5*eV};
  G4double bialkalicathodeabsorp[numentriesbialkalicath] = {1.e-6*mm, 1.e-6*mm}; // absorb all
  BialkaliCathodeProp->AddProperty("ABSLENGTH", bialkalicathodeenergy, bialkalicathodeabsorp, numentriesbialkalicath);
  BialkaliCathodeProp->AddProperty("RINDEX", BoroSiGlassEnergy,BoroSiGlassRindex ,nGlass ); // use values from window to prevent refraction
  BialkaliCathode->SetMaterialPropertiesTable(BialkaliCathodeProp);


  //------------------------------------------------------
  //-------------PTFE Wrapping----------------------------
  //------------------------------------------------------

  G4Material *Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
  Teflon->AddElement(C, 0.240183);
  Teflon->AddElement(F, 0.759817);

  const G4int iNbEntries=2;
  G4double pdTeflonPhotonMomentum[iNbEntries]  = { 2.0*eV, 3.5*eV};
  G4double pdTeflonRefractiveIndex[iNbEntries] = {1.3, 1.3};
  G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

  // Optical Properties
  pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, iNbEntries);
  pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, iNbEntries);
  pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularLobe, iNbEntries);
  pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter, iNbEntries);
  pTeflonPropertiesTable->AddProperty("EFFICIENCY", pdTeflonPhotonMomentum, pdTeflonEfficiency, iNbEntries);
  Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

  //--------------------------------------------------------------------
  //-----Wavelength Shifting Optical Fiber Core-  WLSfiber PMMA
  //-------------------------------------------------------------------

  std::vector<G4int> natoms;
  std::vector<G4String> elements;

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  density = 1.190*g/cm3;

  G4Material*   PMMA = man->
          ConstructNewMaterial("PMMA", elements, natoms, density);

  elements.clear();
  natoms.clear();
 G4double RefractiveIndexWLSfiber[nEntries] =
  { 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
    1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60};

  G4double AbsWLSfiber[nEntries] =
  {5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
   5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
   5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
   1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
   1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};

  G4double EmissionFib[nEntries] =
  {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
   3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
   12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
   15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};               // WLS part of the material : emission

  // Add entries into properties table
  G4MaterialPropertiesTable* MPTWLSfiber = new G4MaterialPropertiesTable();
  MPTWLSfiber->
           AddProperty("RINDEX",PhotonEnergy,RefractiveIndexWLSfiber,nEntries);
  MPTWLSfiber->AddProperty("WLSABSLENGTH",PhotonEnergy,AbsWLSfiber,nEntries);
  MPTWLSfiber->AddProperty("WLSCOMPONENT",PhotonEnergy,EmissionFib,nEntries);
  MPTWLSfiber->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);

  PMMA->SetMaterialPropertiesTable(MPTWLSfiber); ////// core of fibre is WLS PMMA


  //----------------------------------------------------------------------
  //----- Wavelength Shifting Optical Fiber Cladding (polyethylene)
  //----------------------------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.200*g/cm3;

  G4Material* Pethylene = man->
          ConstructNewMaterial("Pethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

   G4double RefractiveIndexClad1[nEntries] =
  { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49};

  G4double AbsClad[nEntries] =
  {20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m};

  // Add entries into properties table
  G4MaterialPropertiesTable* MPTClad1 = new G4MaterialPropertiesTable();
  MPTClad1->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexClad1,nEntries);
  MPTClad1->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntries);

  Pethylene->SetMaterialPropertiesTable(MPTClad1);
  ////// cladding is Polyethylene: not a wls material


  //--------------------------------------------------------------------------------
  //---- Wavelength Shifting Optical Fiber Double Cladding (fluorinated polyethylene)
  //--------------------------------------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.400*g/cm3;

  G4Material* FPethylene = man->
          ConstructNewMaterial("FPethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  G4double RefractiveIndexClad2[nEntries] =
   { 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42};

  // Add entries into properties table
  G4MaterialPropertiesTable* MPTClad2 = new G4MaterialPropertiesTable();
  MPTClad2->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexClad2,nEntries);
  MPTClad2->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntries);

  FPethylene->SetMaterialPropertiesTable(MPTClad2);

  //--------------------------------------------------
  //-------------Polystyrene--------------------------
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(8);
  elements.push_back("H");     natoms.push_back(8);

  density = 1.050*g/cm3;

  G4Material*  Polystyrene = man->
          ConstructNewMaterial("Polystyrene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  G4double RefractiveIndexPS[nEntries] =
    { 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
      1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
      1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
      1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
      1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50};

  G4double AbsPS[nEntries] =
    {2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
     2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
     2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
     2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
     2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm};

  G4double ScintilFast[nEntries] =
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  // Add entries into properties table
  G4MaterialPropertiesTable* MPTPolystyrene = new G4MaterialPropertiesTable();
  MPTPolystyrene->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexPS,nEntries);
  MPTPolystyrene->AddProperty("ABSLENGTH",PhotonEnergy,AbsPS,nEntries);
  MPTPolystyrene->
               AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,nEntries);
  MPTPolystyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  MPTPolystyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPolystyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);

  Polystyrene->SetMaterialPropertiesTable(MPTPolystyrene);
  // Set the Birks Constant for the Polystyrene scintillator

  Polystyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);


}

void PSDetectorConstruction::defineSurfaces()
{
  //-----------------------------------------
  //------------- Surfaces ------------------
  //-----------------------------------------


  //--------------------------------------------
  //-----Plastic Scintillator & ExpHall---------
  //--------------------------------------------
  G4OpticalSurface* OpScintSurface = new G4OpticalSurface("ScintSurface");
  OpScintSurface->SetType(dielectric_dielectric);
  OpScintSurface->SetFinish(ground);
  OpScintSurface->SetModel(unified);

  G4LogicalBorderSurface * logicScintSurface = new G4LogicalBorderSurface("lScintSurface",scint_phys,expHall_phys,OpScintSurface);


  //---------------------------------------------------------
  //------------ Scint to Air-------------------------------
  //---------------------------------------------------------

  G4OpticalSurface* ScintAirSurface = new G4OpticalSurface("ScintAirSurface", unified,  ground,  dielectric_dielectric);
  G4double ScintRefl[2] = {1.0, 1.0};
  G4double energy[2] = {0.1*eV, 100.0*eV};
  G4MaterialPropertiesTable * ScintSurfaceProp = new G4MaterialPropertiesTable();
  ScintSurfaceProp->AddProperty("REFLECTIVITY",energy,ScintRefl,2);
  ScintAirSurface->SetMaterialPropertiesTable(ScintSurfaceProp);
  G4LogicalBorderSurface * logicScintAirSurface = new G4LogicalBorderSurface("lScintAirSurface",scint_phys, physiAir,ScintAirSurface);

  //---------------------------------------------------------
  //------------ Air to Teflon-------------------------------
  //---------------------------------------------------------

  G4OpticalSurface* AirTefSurface = new G4OpticalSurface("AirTefSurface" , unified,  polished,  dielectric_metal);
  G4double TefRefl[2] = {0.98, 0.98};
  G4double AirRINDEX[2] = {1.0,1.0};
  G4MaterialPropertiesTable* AirTefProp = new G4MaterialPropertiesTable();
  AirTefProp->AddProperty("REFLECTIVITY", energy, TefRefl, 2);
  AirTefSurface->SetMaterialPropertiesTable(AirTefProp);
  G4LogicalBorderSurface * logicAirTefSurface = new G4LogicalBorderSurface("lAirTefSurface",physiAir, physiTeflon, AirTefSurface);


  //-------------------------------------------
  //-----------Teflon & expHall Surface--------
  //-------------------------------------------
  G4OpticalSurface* OpTeflonSurface = new G4OpticalSurface("TeflonSurface");
  OpTeflonSurface->SetType(dielectric_LUT);
  OpTeflonSurface->SetModel(LUT);
  OpTeflonSurface->SetFinish(groundteflonair);

  //-------------------------------------------
  //-----------Scint & grease Surface--------
  //-------------------------------------------
  G4OpticalSurface* OpGreaseSurface = new G4OpticalSurface("GreaseSurface");
  OpGreaseSurface->SetModel(unified);
  OpGreaseSurface->SetType(dielectric_dielectric);
  OpGreaseSurface->SetFinish(ground);
  OpGreaseSurface->SetSigmaAlpha(0.2);
  G4double p_gr[2] = {2.00*eV, 3.47*eV};
  G4double refl_gr[2] = {1.0,1.0};
  G4double effi_gr[2] = {0, 0};
  G4MaterialPropertiesTable* OpGreaseSurfaceProp =  new G4MaterialPropertiesTable();
  OpGreaseSurfaceProp->AddProperty("REFLECTIVITY",p_gr,refl_gr,2);
  OpGreaseSurfaceProp->AddProperty("EFFICIENCY",p_gr,effi_gr,2);

  G4double grspecularlobe[2]={0,0}; //this gives the percentage of reflected photons that neighbor the pure specular reflection
  G4double grspecularspike[2]={0,0};//pure specular reflection percentage
  G4double grbackscatter[2]={0,0}; //backscatter percentage
  //pure Lambertian percentage is defined by previous 3 parameters as 1-(specularlobe+specularspike+backscatter)
  OpGreaseSurfaceProp->AddProperty("SPECULARLOBECONSTANT",p_gr ,grspecularlobe,2);
  OpGreaseSurfaceProp->AddProperty("SPECULARSPIKECONSTANT",p_gr ,grspecularspike,2);
  OpGreaseSurfaceProp->AddProperty("BACKSCATTERCONSTANT",p_gr ,grbackscatter,2);

  OpGreaseSurface -> SetMaterialPropertiesTable(OpGreaseSurfaceProp);
  G4LogicalBorderSurface * logicGreaseSurface = new G4LogicalBorderSurface("lGreaseSurface", scint_phys, physiGrease,  OpGreaseSurface);
  G4LogicalSurface * logicGreaseSurface2 = new G4LogicalBorderSurface("lGreaseSurface2",  physiGrease, scint_phys, OpGreaseSurface);
  G4LogicalSurface * logicGreaseSurface3 = new G4LogicalBorderSurface("lGreaseSurface3",  physiGrease, physiWindow, OpGreaseSurface);
  G4LogicalSurface * logicGreaseSurface4 = new G4LogicalBorderSurface("lGreaseSurface4",physiWindow,   physiGrease, OpGreaseSurface);

  //-------------------------------------------
  //-----------PMT Reflector Surface--------
  //-------------------------------------------
  G4OpticalSurface* PMTReflSurface = new G4OpticalSurface("PMTReflSurface",
							 glisur,
							 polished,
							  dielectric_metal, 1.0);

  G4MaterialPropertiesTable* PMTSurfaceProperty =  new G4MaterialPropertiesTable();
  G4double p_pmt[2] = {1.00*eV, 10*eV};
  G4double refl_pmt[2] = {1.0,1.0};
  G4double effi_pmt[2] = {0.0, 0.0};


  PMTSurfaceProperty->AddProperty("REFLECTIVITY",p_pmt,refl_pmt,2);//->SetSpline(true);
  PMTSurfaceProperty->AddProperty("EFFICIENCY",p_pmt,effi_pmt,2);//->SetSpline(true);

  PMTReflSurface->SetMaterialPropertiesTable(PMTSurfaceProperty);

  G4LogicalSkinSurface * SkinPMTReflSurface = new G4LogicalSkinSurface("SkinPMTReflSurface", logicRefl, PMTReflSurface);



}

G4VPhysicalVolume* PSDetectorConstruction::ConstructDetector()
{

  G4UImanager* UI = G4UImanager::GetUIpointer();


  //----------------------------------------------
  //-------------- Determine Positions -----------
  //----------------------------------------------

  G4ThreeVector positionScint = G4ThreeVector(0,0,fScintVol_z);
  G4ThreeVector positionGrease = G4ThreeVector(0,0,2*fScintVol_z +fGreaseThickness );//+fGreaseThickness);
  G4ThreeVector positionCathode = G4ThreeVector(0,0,2*fScintVol_z +2*fGreaseThickness +2*fWindowThickness +fCathodeThickness);
  G4ThreeVector positionWindow = G4ThreeVector(0,0,2*fScintVol_z +  2*fGreaseThickness + fWindowThickness  );
  G4ThreeVector positionTeflon = G4ThreeVector(0,0,fScintVol_z);//-fPTFEThickness/2 - fAirGap/2);
  G4ThreeVector positionRefl = G4ThreeVector(0,0,2*fScintVol_z +2*fGreaseThickness + fCathodeThickness + fWindowThickness);


  //----------------Construct Volumes ------------------

  //----------------------------------------------------
  //------------ The experimental Hall------------------
  //----------------------------------------------------

  expHall_box = new G4Box("sWorld",expHall_x,expHall_y,expHall_z);
  expHall_log = new G4LogicalVolume(expHall_box,G4Material::GetMaterial("Air"),"lWorld",0,0,0);
  expHall_phys= new G4PVPlacement(0,G4ThreeVector(),expHall_log,"pWorld",0,false,0);

  //-------------------------------------------------------
  //------------- The Plastic Scintillator-----------------
  //-------------------------------------------------------

  scint_vol = new G4Box("sScint",fScintVol_x,fScintVol_y,fScintVol_z);
  scint_log = new G4LogicalVolume(scint_vol,G4Material::GetMaterial("BC408"),"lScint",0,0,0);
  scint_phys= new G4PVPlacement(0,positionScint,scint_log,"pScint", expHall_log,false,0);


  //-------------------------------------------------------
  //-----------Hole in Plastic Scint for fiber-------------
  //-------------------------------------------------------

  solidHole = new G4Tubs("sHole",
			 0.0*cm,
			 holeRadius,
			 holeLength/2,
			 0.*deg,
			 360.*deg);

  logicHole = new G4LogicalVolume(solidHole,
				  G4Material::GetMaterial("Air"),
				  "lHole");

  physiHole = new G4PVPlacement(0,
				G4ThreeVector(0,0,0),
				logicHole,
				"pHole",
				scint_log,
				false,
				0);

  //------------------------------------------------------
  //-------------------Grease-----------------------------
  //------------------------------------------------------


  solidGrease = new G4Tubs("sGrease", 0.*cm, cladRX+ 0.5*mm, fGreaseThickness, 0.*deg, 360.*deg);
  logicGrease = new G4LogicalVolume(solidGrease, G4Material::GetMaterial("Polydimethylsiloxane"), "lGrease", 0,0,0);
  physiGrease = new G4PVPlacement(0, positionGrease, logicGrease, "pGrease", expHall_log, false, 0);

  //------------------------------------------------------
  //-------------------Teflon Wrapping -------------------
  //------------------------------------------------------

  // make a hollow box made of ptfe surrounding scintillator
  G4VSolid * t1= new G4Box("t1",fScintVol_x+fAirGap+fPTFEThickness, fScintVol_y+fAirGap+fPTFEThickness, fScintVol_z+fAirGap+fPTFEThickness);
  G4VSolid *t2 = new G4Box("t2",fScintVol_x+fAirGap, fScintVol_y+fAirGap, fScintVol_z+fAirGap);
  G4double z = fPTFEThickness/2+fAirGap/2;///2 +fAirGap/2;// +2*fAirGap  ;
  G4ThreeVector zTrans(0, 0,0 );

  G4SubtractionSolid* solidFullTeflon = new G4SubtractionSolid("hollowteflon", t1, t2, 0, zTrans);

  //Make hole for fiber exit: subtract solid tubs size of (0, hole, ptfethickness)
  G4VSolid *tHole = new G4Tubs("sTeflonHole",
			       0.0*cm,
			       holeRadius,
			       fPTFEThickness/2,
			       0.*deg,
			       360.*deg);
  G4ThreeVector zHole(0, 0, fScintVol_z +fAirGap + fPTFEThickness/2 );

  G4SubtractionSolid* solidHoleTeflon = new G4SubtractionSolid("hollowHoleteflon", solidFullTeflon, tHole, 0, zHole);
  logicTeflon = new G4LogicalVolume(solidHoleTeflon, G4Material::GetMaterial("Teflon"),"lTeflon",0,0,0);
  physiTeflon = new G4PVPlacement(0, positionTeflon, logicTeflon, "pTeflon", expHall_log, false, 0);

  //-----------------------------------------------------
  //----------------- PMT Window/Envelope----------------
  // ----------------------------------------------------
  // borosilicate glass n=1.49 spectral range 280-680nm
  // photocathode active size 25mm

  solidWindow = new G4Tubs("sWindow", 0.0*cm, fPMTDiameter/2, fWindowThickness, 0.0*deg, 360.0*deg);
  logicWindow = new G4LogicalVolume(solidWindow,G4Material::GetMaterial("BoroSiGlass"), "lWindow", 0,0,0);
  physiWindow = new G4PVPlacement(0, positionWindow, logicWindow, "pWindow", expHall_log, false, 0);


  // ----------------------------------------------------
  // -----------------PMT Cathode------------------------
  // ----------------------------------------------------
  // uncertain as to the position & size of cathode. Does it go behind the window?

  solidCathode = new G4Tubs("sCath", 0.0*cm, fPMTDiameter/2, fCathodeThickness, 0.0*deg, 360.0*deg);
  logicCathode = new G4LogicalVolume(solidCathode,G4Material::GetMaterial("BialkaliCathode") , "lCath", 0,0,0);
  physiCathode = new G4PVPlacement(0, positionCathode, logicCathode, "pCath", expHall_log, false, 0);


  //-----------------------------------------------------
  //---------------Define sensitive detector ------------
  // ----------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String sensitiveDetectorName = "/detector/sensitiveDetector";
  CathodeSD* theCathodeSD = new CathodeSD(sensitiveDetectorName, physiCathode); // include "CathodeSD.hh"

  CathodeSD::SaveData e =CathodeSD::Energy;
  theCathodeSD->setSaveOption( e , true);
  SDman->AddNewDetector(theCathodeSD);
  logicWindow->SetSensitiveDetector(theCathodeSD);


  //-----------------------------------------------------------
  //---------------Visualisation attributions------------------
  //-----------------------------------------------------------

  G4VisAttributes* ScintVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0, 0.6)); //magenta
  ScintVisAtt->SetForceSolid(true);
  G4VisAttributes* CathodeVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.5)); // cyan
  CathodeVisAtt->SetForceSolid(true);
  G4VisAttributes* WindowVisAtt = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.2)); // blue
  WindowVisAtt->SetForceSolid(true);
  G4VisAttributes* GreaseVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0, 0.5)); // yellow
  GreaseVisAtt->SetForceSolid(true);
  G4VisAttributes* TeflonVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5, 0.5)); //grey
  TeflonVisAtt->SetForceSolid(true);

  scint_log->SetVisAttributes(ScintVisAtt);
  logicCathode->SetVisAttributes(CathodeVisAtt);
  logicWindow->SetVisAttributes(WindowVisAtt);
  logicGrease->SetVisAttributes(GreaseVisAtt);
  logicTeflon->SetVisAttributes(TeflonVisAtt);

  //-----------------Trajectory Colours--------------------------------------
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID;
  G4TrajectoryDrawByParticleID* model2 = new G4TrajectoryDrawByParticleID("test");

  model->SetDefault("cyan");
  model->Set("opticalphoton", "cyan");
  model->Set("gamma", "green");
  model->Set("e+", "magenta");
  model->Set("e-", G4Colour(1.0, 0.5, 0.0));
  model->Set("alpha", "yellow");

  visManager->RegisterModel(model);

  visManager->SelectTrajectoryModel(model->Name());


  //--------------------------------------------------
  //---------------Fiber------------------------------
  //--------------------------------------------------

  ConstructFiber();

  //--------------------------------------------------
  //-------End of Construction------------------------
  //--------------------------------------------------

  return expHall_phys;
}


void PSDetectorConstruction::ConstructFiber()
{
  if (!(logicHole) || !(physiHole) ) {
     std::ostringstream o;
     o << "The Fiber Hole has not been constructed";
     G4Exception("WLSDetectorConstruction::ConstructFiber","",
                  FatalException,o.str().c_str());
  }

  // Pointers to the most recently constructed volume
  G4LogicalVolume* logicPlacement = logicHole;
  G4VPhysicalVolume* physiPlacement = physiHole;

  //--------------------------------------------------
  //----------Fiber Construction----------------------
  //--------------------------------------------------

  // Boundary Surface Properties
  G4OpticalSurface* OpSurface = NULL;

  if (surfaceRoughness < 1.)
     OpSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                      glisur,                  // SetModel
                                      ground,                  // SetFinish
                                      dielectric_dielectric,   // SetType
                                      surfaceRoughness);       // SetPolish

  G4LogicalVolume   *logicClad;
  G4VPhysicalVolume *physiClad;


  //--------------------------------------------------
  //----------------Cladding--------------------------
  //--------------------------------------------------

  G4VSolid* solidClad;

  if (XYRatio == 1.)
    solidClad = new G4Tubs("sClad",0.,cladRX,cladZ,0.0*rad,twopi*rad);
  else
    solidClad = new G4EllipticalTube("sClad",cladRX,cladRY,cladZ);

  logicClad = new G4LogicalVolume(solidClad,
				  G4Material::GetMaterial("Pethylene"),
				  "lClad");
  G4double  z = wlsfiberZ- holeLength/2 ;
  physiClad = new G4PVPlacement(0,
				G4ThreeVector(0.0,0.0,z),
				logicClad,
				"pClad",
				logicPlacement,
				false,
				0);



  // Place the rough surface only if needed
  if (OpSurface) {
    new G4LogicalBorderSurface("surfaceCladOut",
			       physiClad,
			       physiPlacement,
			       OpSurface);
    new G4LogicalBorderSurface("surfaceCladIn",
			       physiPlacement,
			       physiClad,
			       OpSurface);
  }

  logicPlacement = logicClad;
  physiPlacement = physiClad;


  //--------------------------------------------------
  //----------------WLS Fiber-------------------------
  //--------------------------------------------------

  G4VSolid* solidWLSfiber;

  if (XYRatio == 1.)
    solidWLSfiber =
      new G4Tubs("sWLSFiber",0.,wlsfiberRX,wlsfiberZ,0.0*rad,twopi*rad);
  else
    solidWLSfiber =
      new G4EllipticalTube("sWLSFiber",wlsfiberRX,wlsfiberRY,wlsfiberZ);

  G4LogicalVolume*   logicWLSfiber =
    new G4LogicalVolume(solidWLSfiber,
			G4Material::GetMaterial("PMMA"),
			"lWLSFiber");

  logicWLSfiber->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));           //   ?

  G4VPhysicalVolume* physiWLSfiber = new G4PVPlacement(0,
						       G4ThreeVector(0.0,0.0,0.0),
						       logicWLSfiber,
						       "pWLSFiber",
						       logicPlacement,
						       false,
						       0);


  G4VisAttributes* WLSVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.8)); // cyan
  WLSVisAtt->SetForceSolid(true);
  logicWLSfiber->SetVisAttributes(WLSVisAtt);

  // Place the rough surface only if needed
  if (OpSurface) {
    new G4LogicalBorderSurface("surfaceWLSOut",
			       physiWLSfiber,
			       physiPlacement,
			       OpSurface);
    new G4LogicalBorderSurface("surfaceWLSIn",
			       physiPlacement,
			       physiWLSfiber,
			       OpSurface);
  }


  //--------------------------------------------------
  // Mirror for reflection at one of the end
  //--------------------------------------------------

  // Place the mirror only if the user wants the mirror
  if (mirrorToggle) {

    G4VSolid* solidMirror = new G4Tubs("sMirror",0.,cladRX,mirrorZ/2,0.0*rad,twopi*rad);
    G4LogicalVolume* logicMirror = new G4LogicalVolume(solidMirror,
						       G4Material::GetMaterial("G4_Al"),
						       "lMirror");
    G4OpticalSurface* MirrorSurface = new G4OpticalSurface("MirrorSurface",
							   glisur,
							   ground,
							   dielectric_metal,
							   mirrorPolish);

    G4MaterialPropertiesTable* MirrorSurfaceProperty =
      new G4MaterialPropertiesTable();

    G4double p_mirror[2] = {2.00*eV, 3.47*eV};
    G4double refl_mirror[2] = {mirrorReflectivity,mirrorReflectivity};
    G4double effi_mirror[2] = {0, 0};

    MirrorSurfaceProperty->AddProperty("REFLECTIVITY",p_mirror,refl_mirror,2);
    MirrorSurfaceProperty->AddProperty("EFFICIENCY",p_mirror,effi_mirror,2);

    MirrorSurface -> SetMaterialPropertiesTable(MirrorSurfaceProperty);

    new G4PVPlacement(0,
		      G4ThreeVector(0.0,0.0,mirrorOrigin),
		      logicMirror,
		      "Mirror",
		      expHall_log,
		      false,
		      0);

    new G4LogicalSkinSurface("MirrorSurface",logicMirror,MirrorSurface);

    G4VisAttributes* MirrorVisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75, 0.5)); // silver
    MirrorVisAtt->SetForceSolid(true);
    logicMirror->SetVisAttributes(MirrorVisAtt);


  }

  G4VisAttributes* CladVisAtt = new G4VisAttributes(G4Colour(0.17, 0.74, 0.74, 0.5)); //darker cyan
  CladVisAtt->SetForceSolid(true);
  logicClad->SetVisAttributes(CladVisAtt);

}
