
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
// $Id: B4DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class
#include <math.h> 
#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Sphere.hh" // included by rp for sphere
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4TessellatedSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//#define USE_CADMESH_TETGEN

//#define USE_CADMESH_TETGEN
#include "CADMesh.hh"
//#include "tetgen.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
	fSurfaceMPT = new G4MaterialPropertiesTable(); // Inits Object

	fSurface = new G4OpticalSurface("Surface"); // Inits Object
	fSurface->SetType(dielectric_dielectric); // Type of material on Each side of the surface
	fSurface->SetFinish(ground); // Finish of the Surface
	fSurface->SetModel(unified);
	fSurface->SetMaterialPropertiesTable(fSurfaceMPT); // Material Properties
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
	delete fSurface;
	delete fSurfaceMPT;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Some materials are defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density;
  G4int ncomponents, natoms; 
  G4Element* C  = nistManager->FindOrBuildElement(6);
  G4Element* H  = nistManager->FindOrBuildElement(1);
  G4Element* O  = nistManager->FindOrBuildElement(8);
  G4Element* F  = nistManager->FindOrBuildElement(9);

  G4Material* Scintillator = // Scintillator material
	new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
	Scintillator->AddElement(C, natoms=9);
	Scintillator->AddElement(H, natoms=10);
	Scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  nistManager->FindOrBuildMaterial("G4_AIR"); // Air material

  G4Material* polystyrene =
  new G4Material("Polystyrene",density= 1.05*g/cm3, ncomponents=2);
  polystyrene->AddElement(C, natoms=8);
  polystyrene->AddElement(H, natoms=8);

  G4Material* pmma =
  new G4Material("PMMA",density= 1.19*g/cm3, ncomponents=3);
  pmma->AddElement(C, natoms=5);
  pmma->AddElement(H, natoms=8);
  pmma->AddElement(O, natoms=2);

  G4Material* fluorinatedPolymer =
  new G4Material("Fluorinated_Polymer", density= 1.43*g/cm3, ncomponents=2);
  fluorinatedPolymer->AddElement(C,2);
  fluorinatedPolymer->AddElement(F,2);
  G4MaterialPropertiesTable* mpFS;

  //--- Generate and add material properties table ---
    G4double PhotonEnergy[] = {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
        2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
        2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
        2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
        2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
        2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
        2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
        3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
        3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
        3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

        const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

  //--- Fluorinated Polymer (FS) ---
    G4double RefractiveIndex_FluorinatedPolymer[nEntries] =
    {       
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42
    };
    mpFS = new G4MaterialPropertiesTable();
    mpFS->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_FluorinatedPolymer,nEntries);
    fluorinatedPolymer->SetMaterialPropertiesTable(mpFS);

  G4Material* core_S_Material = polystyrene;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  /*  G4 Gerometry Tree
     World
       - Coil
  */

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
  auto sCoilMaterial  = G4Material::GetMaterial("Fluorinated_Polymer");

  if ( ! defaultMaterial || ! sCoilMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  } 
   
  //     
  // World
  //

  worldSizeX = 25.0*cm ;  // half width
  worldSizeY = 25.0*cm ;  // half width
  worldSizeZ = 25.0*cm ;  // half width

  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,
                 fCheckOverlaps);		 // copy number

	//
	// G4Coil
	//
	int CreateCoil = 0;
	if(CreateCoil){
		auto G4Coil = CADMesh::TessellatedMesh::FromSTL("G4CoilSingleTurn.stl");
		G4Coil->SetScale(1.0);

		G4double xoffset=0.0;
		G4double yoffset=0.0;
		G4double zoffset=0.0;

		G4Coil->SetOffset(G4ThreeVector(xoffset, yoffset, zoffset));

		auto G4CoilS = G4Coil->GetSolid();

		auto G4CoilLV
				= new G4LogicalVolume(
						G4CoilS,
						sCoilMaterial,
						"G4CoilLV");

		for(int i=0;i<20;i++){
			new G4PVPlacement(
					0,
					G4ThreeVector(0*m,0*m,((-0.09)+i*(0.010))*m),
					G4CoilLV,
					"G4CoilPV",
					worldLV,
					false,
					0,
				fCheckOverlaps);
	}
	}

	//
	// Test Cylinder
	//

	G4double innerRadius = 0.0*cm;
	G4double outerRadius = 0.5*cm;
	G4double height = 20.0*cm;

	G4RotationMatrix* rotation = new G4RotationMatrix();
	rotation->rotateX(5.0*deg);

	// Create a G4Tubs object with the specified dimensions
	G4Tubs* cylinderS = new G4Tubs("cylinder", innerRadius, outerRadius, height/2, 0.0, 2.0*M_PI);

	auto cylinderLV
		= new G4LogicalVolume(
				cylinderS,
				sCoilMaterial,
				"cylinderLV");

	auto cylinderPV =
		new G4PVPlacement(
				rotation,
				G4ThreeVector(0.0, 0.0, 0.0),
				cylinderLV,
				"cylinderPV",
				worldLV,
				false,
				0,
				fCheckOverlaps);


	//
	// Test Code for Surface Component
	//

	/*Things to add:
	1. Add the extra properties that are not on the constructor of this file
	2.
	*/

	G4LogicalBorderSurface* surface =
	new G4LogicalBorderSurface("Surface", cylinderPV, worldPV, fSurface);

	G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
	surface->GetSurface(cylinderPV, worldPV)->GetSurfaceProperty());
	G4cout << "******  opticalSurface->DumpInfo:" << G4endl;
	if(opticalSurface)
	{
	opticalSurface->DumpInfo();
	}
	G4cout << "******  end of opticalSurface->DumpInfo" << G4endl;

	SetSurfaceSigmaAlpha(0.2);
	G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
	mpv->InsertValues(0.000002, 0.1);
	mpv->InsertValues(0.000008, 0.1);
	AddSurfaceMPV("SPECULARLOBECONSTANT",mpv);
	AddSurfaceMPV("SPECULARSPIKECONSTANT",mpv);
	AddSurfaceMPV("BACKSCATTERCONSTANT",mpv);

	AddSurfaceMPV("TRANSMITTANCE",mpv);
	AddSurfaceMPV("REFLECTIVITY",mpv);

	G4MaterialPropertyVector* mpv1 = new G4MaterialPropertyVector();
	mpv1->InsertValues(0.000002, 0.8);
	mpv1->InsertValues(0.000008, 0.8);

	AddSurfaceMPV("BACKSCATTERCONSTANT",mpv1);

	//                                       
	// Visualization attributes
	//

	worldLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0)));
	cylinderLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.627,0.125,0.941)));

	//
	// Always return the physical World
	//
	return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::AddSurfaceMPV(const G4String& prop,
                                         G4MaterialPropertyVector* mpv)
{
  fSurfaceMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::AddSurfaceMPC(const G4String& prop, G4double v)
{
  fSurfaceMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::SetSurfaceSigmaAlpha(G4double v)
{
  fSurface->SetSigmaAlpha(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface sigma alpha set to: " << fSurface->GetSigmaAlpha()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B4DetectorConstruction::SetSurfacePolish(G4double v)
{
  fSurface->SetPolish(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface polish set to: " << fSurface->GetPolish() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......