#include <cmath>
#include <fstream>
#include <iostream>

#include "pmepot.hpp"

/** Test driver for Tim's Particle Mesh Ewald (PME) code.
 */
int main()
{
    // Conversion constants.
    const double amberChargeFactor = 18.2223;
    const double avogadrosNumber = 6.02214129E23;
    const double boltzmannConstant = 1.380658E-23;
    const double coulombConstant = 8.9875517873681764E9;
    const double coulombPerElementaryCharge = 1.602176565E-19;
    const double hydrogenCharge_cSPCE = 0.4238;
    const double temperature = 300; // in Kelvin
        
    // PME parameters for phosphonium (PH4+).
    // Input units are Angstroms and elementary charges.
    int numAtoms = 5;

    double positions[3*numAtoms]; 
    positions[0] = 15.000; positions[1] = 15.000; positions[2] = 15.000;
    positions[3] = 15.224; positions[4] = 14.910; positions[5] = 13.609;
    positions[6] = 15.563; positions[7] = 16.198; positions[8] = 15.491;
    positions[9] = 15.602; positions[10] = 13.898; positions[11] = 15.645;
    positions[12] = 13.611; positions[13] = 14.993; positions[14] = 15.254;

    double charges[numAtoms];
    charges[0] = 1.88345693E+01;
    charges[1] = -1.53067320E-01; charges[2] = -1.53067320E-01;
    charges[3] = -1.53067320E-01; charges[4] = -1.53067320E-01;
    for (int i = 0; i < numAtoms; i++)
    {
        // Convert Amber units to elementary charge.
        charges[i] /= amberChargeFactor;
    }
    
    double boxLenX = 30;
    double boxLenY = 30;
    double boxLenZ = 30;
    double gridSpacing = 0.5;
    int gridDimX = boxLenX / gridSpacing;
    int gridDimY = boxLenY / gridSpacing;
    int gridDimZ = boxLenZ / gridSpacing;
    // double smear = 0.11;
    double smear = 1;
    int splineOrder = 10;
    double cutoff = 100;

    const int numGridPoints = gridDimX * gridDimY * gridDimZ;
    
    double p[numGridPoints];
    
    // Run PME.
    pmepot_(numAtoms, positions, charges,
            gridDimX, gridDimY, gridDimZ, boxLenX, boxLenY, boxLenZ,
            smear, splineOrder, cutoff, p);

    // Write potential to OpenDX file for inspection.
    std::ofstream fout("pmepot.dx");

    if (!fout.is_open())
    {
        std::cout << "Failed to open file." << std::endl;
        return 1;
    }

    double sumPotential = 0;
    for (int i = 0; i < numGridPoints; i++)
    {
        // Convert potential (elementary charge per Angstrom) to Amber
        // comparable potential energy (kT).
        p[i] *= 1E10 * hydrogenCharge_cSPCE * coulombConstant
            * pow(coulombPerElementaryCharge, 2) / (boltzmannConstant * temperature);
        p[i] *= 1E-1;
        sumPotential += p[i];
    }

    fout << "# sumPotential: " << sumPotential << std::endl;

    fout << "object 1 class gridpositions counts "
         << gridDimX << " " << gridDimY << " " << gridDimZ << std::endl;
    fout << "origin 0 0 0" << std::endl;
    fout << "delta  " << gridSpacing << "  0.00000000  0.00000000" << std::endl;
    fout << "delta  0.0000000  " << gridSpacing << "  0.00000000" << std::endl;
    fout << "delta  0.0000000  0.00000000  " << gridSpacing << std::endl;
    fout << "object 2 class gridconnections counts  "
         << gridDimX << " " << gridDimY << " " << gridDimZ << std::endl;
    fout << "object 3 class array type double rank 0 items    "
         << numGridPoints << " data follows" << std::endl;
    
    for (int i = 0; i < numGridPoints; i += 3)
    {
        fout << p[i] << " ";
        fout << p[i + 1] << " ";
        fout << p[i + 2] << std::endl;
    }

    fout << "object \"Untitled\" class field" << std::endl;

    return 0;
}
