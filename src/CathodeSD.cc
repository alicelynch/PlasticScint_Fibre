#include "CathodeSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "PSStackingAction.hh"

//extern std::ofstream outFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CathodeSD::CathodeSD(G4String name, G4VPhysicalVolume *cathode)
  : G4VSensitiveDetector(name), savetime(true), saveposition(false), saveenergy(false), cathode(cathode)
{
  G4String HCname;
  collectionName.insert(HCname="pmtCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CathodeSD::~CathodeSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CathodeSD::Initialize(G4HCofThisEvent* HCE)
{
  photonCollection = new CathodeHitsCollection
                          (SensitiveDetectorName,collectionName[0]);
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, photonCollection );
  outFile.open("PlasticScint_Fibre.out", std::ofstream::app);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CathodeSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  if (aStep->GetTrack()->GetDefinition()->GetParticleType() != "opticalphoton")
    return false;

  if (aStep->GetPostStepPoint()->GetPhysicalVolume() != cathode)
    return false;

  G4double energy = aStep->GetPostStepPoint()->GetTotalEnergy();

  if(energy==0.)
    return false;

  CathodeHit* newHit = new CathodeHit;
  newHit->SetEnergy(energy);
  newHit->SetTime  (aStep->GetPostStepPoint()->GetGlobalTime());
  newHit->SetPos   (aStep->GetPostStepPoint()->GetPosition());
  photonCollection->insert( newHit );

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CathodeSD::EndOfEvent(G4HCofThisEvent*)
{
  // save results
  G4int NbHits = photonCollection->entries();

  if (outFile) {
    G4cout<< "Here1 "<<G4endl;
    outFile << "# " << NbHits << " Hits detected." << std::endl;
    outFile << "# ";
    if (savetime)
      outFile << "Time (ns)";
    if (saveposition) {
      if (savetime)
        outFile << "\t";
      outFile << "Position (x in mm)\tPosition (y in mm)";
    }
    if (saveenergy) {
      if (savetime || saveposition)
        outFile << "\t";
      outFile << "Energy (eV)";
    }

    outFile << std::endl;
    for (G4int i=0;i<NbHits;i++) (*photonCollection)[i]->Print(outFile, savetime, saveposition, saveenergy);
    outFile << std::endl << std::endl;
    outFile.close();
  }

  else {

    G4cout << "=============" << G4endl << "\nCathode: " << NbHits << " optical photons detected" << G4endl;
    
    outFile.open("PlasticScint_Fibre.out", std::ofstream::app);
    if (outFile.is_open()) G4cout<< " File Open" <<G4endl;
    else G4cout<< "File not open"<<G4endl;
    G4cout << "# " << NbHits << " Hits detected." <<G4endl;

    outFile << "# " << NbHits << " Hits detected." << std::endl;
    outFile << "# ";


    if (savetime)
      outFile << "Time (ns)";
    if (saveposition) {
      if (savetime)
        outFile << "\t";
      outFile << "Position (x in mm)\tPosition (y in mm)";
    }
    if (saveenergy) {
      if (savetime || saveposition)
        outFile << "\t";
      outFile << "Energy (eV)";
    }
    outFile << std::endl;
    for (G4int i=0;i<NbHits;i++) (*photonCollection)[i]->Print(outFile, savetime, saveposition, saveenergy);
    outFile << std::endl << std::endl;

    outFile.close();
  }

}

void CathodeSD::setSaveOption(CathodeSD::SaveData type, bool save)
{
  switch (type) {
    case Time:
      savetime = save;
      break;
    case Position:
      saveposition = save;
      break;
    case Energy:
      saveenergy = save;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
