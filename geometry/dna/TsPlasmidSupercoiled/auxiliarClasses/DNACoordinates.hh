//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//
/*
 *  Developed by Nicholas Henthorn, PRECISE, University of Manchester
 *  nicholas.henthorn@manchester.ac.uk
 *  https://tinyurl.com/y7xuxw9d
 *  17/10/2018, LINK TO PUBLICATION
 */

#ifndef DNACoordinates_h
#define DNACoordinates_h
#include <vector>
#include <string>
#include <cmath>
#include "DNA.hh"
#include "TsVGeometryComponent.hh"
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

using namespace std;
using namespace CLHEP;

#include "G4SystemOfUnits.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class DNACoordinates{
public:
    DNACoordinates(){}
    ~DNACoordinates(){}
    void Generate(vector<G4ThreeVector> &path,
                  vector<DNA*> &DNAPts,
                  G4bool BuildHalfCyl,
                  G4bool BuildQuartCyl,
                  G4bool BuildSphere,
				  G4bool segment,
                  // Sphere-model placement. TsNucleus and this class historically used
                  // different hardcoded values for all of these, which made the two "Sphere"
                  // models different geometries under one name and stopped any scavenging
                  // probability calibrated in one from being valid in the other.
                  G4double helixRadius = 1.15*CLHEP::nanometer,
                  G4double backboneRadius = 0.29*CLHEP::nanometer,
                  G4double baseRadialSize = 0.3*CLHEP::nanometer,
                  G4bool nucleusConvention = false);
    

    
private:
    
    void Segment(vector<G4ThreeVector> &path, vector<G4ThreeVector> &newPath);
    void BuildDNA(vector<G4ThreeVector> &newPath, vector<DNA*> &DNAPts);
    void BuildSphereDNA(vector<G4ThreeVector> &newPath, vector<DNA*> &DNAPts,
                        G4double helixRadius, G4double backboneRadius,
                        G4double baseRadialSize, G4bool nucleusConvention);
    void ApplyRotation(G4ThreeVector &rotated, G4ThreeVector &vector, G4RotationMatrix*rot);
};

#endif /* DNACoordinates_h */
