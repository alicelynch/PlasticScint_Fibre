#ifndef CathodeSD_h
#define CathodeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "CathodeHit.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CathodeSD : public G4VSensitiveDetector
{
  public:
      CathodeSD(G4String, G4VPhysicalVolume *cathode);
     ~CathodeSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

      enum SaveData {
        Time,
        Position,
        Energy
      };

      void setSaveOption(SaveData type, bool save);

  private:
      CathodeHitsCollection* photonCollection;
      bool savetime;
      bool saveposition;
      bool saveenergy;
      const G4VPhysicalVolume *cathode;
      std::ofstream outFile;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

