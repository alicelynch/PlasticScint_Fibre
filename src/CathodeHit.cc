#include "CathodeHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<CathodeHit> CathodeHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CathodeHit::CathodeHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CathodeHit::~CathodeHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CathodeHit::CathodeHit(const CathodeHit& right)
  : G4VHit()
{
  time = right.time;
  energy      = right.energy;
  pos       = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const CathodeHit& CathodeHit::operator=(const CathodeHit& right)
{
  time = right.time;
  energy      = right.energy;
  pos       = right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int CathodeHit::operator==(const CathodeHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CathodeHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    circle.SetScreenDiameter(4.0);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CathodeHit::Print(std::ostream &stream, bool printtime, bool printposition, bool printenergy)
{
  if (printtime)
    stream << time/ns;
  if (printposition) {
    if (printtime)
      stream << "\t";
    stream  << pos.x()/mm << "\t" << pos.y()/mm;
  }
  if (printenergy) {
    if (printtime || printposition)
      stream << "\t";
    stream << energy/eV;
  }
  stream << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

