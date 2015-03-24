#ifndef CathodeHit_h
#define CathodeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CathodeHit : public G4VHit
{
  public:

      CathodeHit();
     ~CathodeHit();
      CathodeHit(const CathodeHit&);
      const CathodeHit& operator=(const CathodeHit&);
      G4int operator==(const CathodeHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print(std::ostream &stream = G4cout, bool printtime = true, bool printposition = false, bool printenergy = false);

  public:

      void SetTime(G4double t)      { time = t; }
      void SetEnergy   (G4double en)      { energy = en; }
      void SetPos      (G4ThreeVector xyz){ pos = xyz; }

      G4double GetTime()  { return time; }
      G4double GetEnergy()    { return energy; }
      G4ThreeVector GetPos(){ return pos; }
      
  private:
  
      G4double      time;
      G4double      energy;
      G4ThreeVector pos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<CathodeHit> CathodeHitsCollection;

extern G4Allocator<CathodeHit> CathodeHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CathodeHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) CathodeHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void CathodeHit::operator delete(void *aHit)
{
  CathodeHitAllocator.FreeSingle((CathodeHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
