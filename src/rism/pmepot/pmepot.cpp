#include "pmepot.hpp"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "fftw3.h"

static const bool DEBUG = false;
std::ofstream fdebug, fdebug2, fdebug3, fdebug4, fdebug5, fdebug6, fdebug7, fdebug8, fdebug9;

double * const ccdl::cross_product(double const * const a, double const * const b)
{
    double * c = new double[3];
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}


double ccdl::dot_product(double const * const a, double const * const b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void ccdl::vector_product_test()
{
    double a[3] = {3, 4, 5};
    double b[3] = {4, 3, 5};
    double c[3] = {-5, -12, -13};
    if (ccdl::dot_product(a, b) != 49)
    {
        std::cout << "ERROR: dot_product(a, b) failed!" << std::endl;
    }
    double * const aXb = ccdl::cross_product(a, b);
    if (aXb[0] != 5 or aXb[1] != 5 or aXb[2] != -7)
    {
        std::cout << "ERROR: cross_product(a, b) failed!" << std::endl;
    }
    if (ccdl::dot_product(a, ccdl::cross_product(b, c)) != 6)
    {
        std::cout << "ERROR: dot_product(a, cross_product(b, c)) failed!" << std::endl;
    }
    double * const aXbXc = ccdl::cross_product(a, ccdl::cross_product(b, c));
    if (aXbXc[0] != -267 or aXbXc[1] != 204 or aXbXc[2] != -3)
    {
        std::cout << "ERROR: cross_product(a, cross_product(b, c)) failed!" << std::endl;
    }
}


double ccdl::anint(double const x)
{
    return (x > 0 ? std::floor(x + 0.5) : std::ceil(x - 0.5));
}


double ccdl::wrap(double const x, double const l)
{
    //TODO: Make this triclinic.
    return x - l * ccdl::anint(x / l);
}


void ccdl::bspline_recursion(double const offset, int const order, double * values)
{
    int orderM1 = order - 1;
    double const div = 1. / orderM1;

    values[orderM1] = div * offset * values[orderM1 - 1];

    for (int j = 1; j < orderM1; ++j)
    {
        values[orderM1 - j] = div * ((offset + j) * values[orderM1 - j - 1]
                                + (order - j - offset) * values[orderM1 - j]);
    }
    values[0] = div * (1 - offset) * values[0];
}


void ccdl::bspline(double const offset, int const order, double * values)
{
    values[1] = offset;
    values[0] = 1. - offset;
    if (order > 2)
    {
        // One pass to order 3:
        values[2] = 0.5 * offset * values[1];
        values[1] = 0.5 * ((offset + 1.) * values[0] + (2. - offset) * values[1]);
        values[0] = 0.5 * (1. - offset) * values[0];
        if (order > 3)
        {
            // One pass to order 4:
            double const div = 1. / 3.;
            values[3] = div * offset * values[2];
            values[2] = div * ((offset + 1.) * values[1] + (3. - offset) * values[2]);
            values[1] = div * ((offset + 2.) * values[0] + (2. - offset) * values[1]);
            values[0] = div * (1. - offset) * values[0];
            // And the rest.
            for (int k = 5; k <= order; ++k)
            {
                ccdl::bspline_recursion(offset, k, values);
            }
        }
    }
}


void ccdl::bspline_periodic_box
(double const x,
 double const boxLen,
 int const gridDim,
 int const order,
 double * values,
 int * nearestGridPoints)
{
    // Find the distance (in terms of fraction grid points) to the
    // closest grid point left of x.
    double xOnGrid = x / (boxLen / gridDim) - (order % 2) * 0.5;
    int nearbyGridPoint = std::floor(xOnGrid);
    double distance = xOnGrid - nearbyGridPoint;
    // Identify the grid points where the b-spline will be evaluated,
    // starting from the leftmost and considering periodicity such
    // that all grid points are kept within the reference box.
    nearbyGridPoint = ((nearbyGridPoint - (order / 2 - 1)) % gridDim + gridDim) % gridDim;
    for (int b = 0; b < order; ++b, ++nearbyGridPoint)
    {
        nearestGridPoints[b] = nearbyGridPoint % gridDim;
    }
    
    ccdl::bspline(distance, order, values);
}


void ccdl::nearby_grid_points_periodic
(double const x,
 double const boxLen,
 int const gridDim,
 int const order,
 int * nearestGridPoints)
{
    double xOnGrid = x / (boxLen / gridDim) - (order % 2) * 0.5;
    int nearbyGridPoint = std::floor(xOnGrid);
    nearbyGridPoint = ((nearbyGridPoint - (order / 2 - 1)) % gridDim + gridDim) % gridDim;
    for (int b = 0; b < order; ++b, ++nearbyGridPoint)
    {
        nearestGridPoints[b] = nearbyGridPoint % gridDim;
    }
}


void ccdl::bspline_Fourier_coefficients
(int const maxIndex,
 int const order,
 double const * bsplineValues,
 double * fourierCoefficients)
{
    for (int i = 0; i < maxIndex; ++i) fourierCoefficients[i] = 0.;

    // o = order / 2 - 1 for even order, order / 2 for odd order.
    int const o = (order - 1) / 2;
    int const Np = order / 2;
    int const Nh = maxIndex - o;
    
    double const tpion = ccdl::TWO_PI / maxIndex;
    for (int m = 0; m < maxIndex; ++m)
    {
        double const k = m * tpion;
        // Iterates from bsplineValues[(order - 1) / 2] to bsplineValues[order - 1].
        for (int n = 0; n <= Np; ++n)
        {
            fourierCoefficients[m] += bsplineValues[n + o] * std::cos(k * n);
        }
        // Iterates from bsplineValues[0] to bsplineValues[order / 2 - 2]
        // for even order, bsplineValues[order / 2 - 1] for odd order.
        for (int i = 0; i < o; ++i)
        {
            fourierCoefficients[m] += bsplineValues[i] * std::cos(k * (Nh + i));
        }
    }
}


void ccdl::bspline_half_Fourier_coefficients
(int const maxIndex,
 int const order,
 double const * bsplineValues,
 double * fourierCoefficients)
{
    for (int i = 0; i < maxIndex / 2 + 1; ++i) fourierCoefficients[i] = 0.;

    int const o = (order - 1) / 2;
    int const Np = order / 2;
    int const Nh = maxIndex - o;

    double const tpion = ccdl::TWO_PI / maxIndex;
    for (int m = 0; m <= maxIndex / 2; ++m)
    {
        double const k = m * tpion;
        for (int n = 0; n <= Np; ++n)
        {
            fourierCoefficients[m] += bsplineValues[n + o] * std::cos(k * n);
        }
        for (int i = 0; i < o; ++i)
        {
            fourierCoefficients[m] += bsplineValues[i] * std::cos(k * (Nh + i));
        }
    }
}


void pmepot_
(int const & numAtoms, double const * positions, double const * charges,
 int const * gridDim, double const * boxLen, double const & boxVolume,
 double const * unitCellUnitVectors, double const * unitCellVectors,
 double const * unitCellVectorsK,
 double const & smear, int const & splineOrder, double const & cutoff,
 double * potential)
{
    if (gridDim[0] % 2 > 0 or gridDim[1] % 2 > 0 or gridDim[2] % 2 > 0)
    {
        std::cerr << "pmepot: PME implementation only supports an even"
            " number of grid points on each axis." << std::endl;
        abort();
    }
    int const numGridPoints = gridDim[0] * gridDim[1] * gridDim[2];
    int const gridDimZ_k = (gridDim[2] / 2 + 1);
    int const numGridPointsK = gridDim[0] * gridDim[1] * gridDimZ_k;
    double const zeta = smear * smear;

    double gridSpacingX = boxLen[0] / gridDim[0];
    double gridSpacingY = boxLen[1] / gridDim[1];
    double gridSpacingZ = boxLen[2] / gridDim[2];

    // Angular wave numbers (assuming orthrorhombic unit cell).
    std::vector<double> kx(gridDim[0], 0.), ky(gridDim[1], 0.), kz(gridDimZ_k, 0.);
    for (int k = 0; k < gridDim[0]; ++k)
    {
        kx[k] = (ccdl::TWO_PI / boxLen[0])
            * (k < gridDim[0] / 2 + 1 ? k : -(gridDim[0] - k));
    }
    for (int k = 0; k < gridDim[1]; ++k)
    {
        ky[k] = (ccdl::TWO_PI / boxLen[1])
            * (k < gridDim[1] / 2 + 1 ? k : -(gridDim[1] - k));
    }
    for (int k = 0; k < gridDimZ_k; ++k)
    {
        kz[k] = (ccdl::TWO_PI / boxLen[2]) * k;
    }
    // double kxi[gridDim[0]][3], kyi[gridDim[1]][3], kzi[gridDimZ_k][3];
    // for (int k = 0; k < gridDim[0]; ++k)
    // {
    //     for (int index = 0; index < 3; ++index)
    //     {
    //         kxi[k][index] = ccdl::TWO_PI * unitCellVectorsK[0 + index * 3]
    //             * (k < gridDim[0] / 2 + 1 ? k : -(gridDim[0] - k));
    //     }    
    // }
    // for (int k = 0; k < gridDim[1]; ++k)
    // {
    //     for (int index = 0; index < 3; ++index)
    //     {
    //         kyi[k][index] = ccdl::TWO_PI * unitCellVectorsK[1 + index * 3]
    //             * (k < gridDim[1] / 2 + 1 ? k : -(gridDim[1] - k)); 
    //     }   
    // }
    // for (int k = 0; k < gridDimZ_k; ++k)
    // {
    //     for (int index = 0; index < 3; ++index)
    //     {
    //         kzi[k][index] = ccdl::TWO_PI * unitCellVectorsK[2 + index * 3] * k;
    //     }
    // }
    

    // Compute the discrete Fourier transform coefficients of the b-spline.
    std::vector<double> bsplineValuesX(splineOrder), bsplineValuesY(splineOrder), bsplineValuesZ(splineOrder);
    std::vector<double> bsplineFourierCoeffX(gridDim[0]), bsplineFourierCoeffY(gridDim[1]);
    std::vector<double> bsplineFourierCoeffZ(gridDimZ_k);

    ccdl::bspline((splineOrder % 2 == 0 ? 0.0 : 0.5), splineOrder, bsplineValuesX.data());
    ccdl::bspline_Fourier_coefficients(
        gridDim[0], splineOrder, bsplineValuesX.data(), bsplineFourierCoeffX.data());
    ccdl::bspline_Fourier_coefficients(
        gridDim[1], splineOrder, bsplineValuesX.data(), bsplineFourierCoeffY.data());
    //TODO: When switching to X_k_half, need to swap functions.
    ccdl::bspline_half_Fourier_coefficients(
        gridDim[2], splineOrder, bsplineValuesX.data(), bsplineFourierCoeffZ.data());


    if (DEBUG)
    {
        fdebug.precision(3);
        fdebug << std::fixed;
        fdebug.open("/home/jesse/rism/bspline_cppX", std::ios::out | std::ios::trunc);
        for (int i = 0; i < gridDim[0]; i++)
        {
            fdebug << bsplineFourierCoeffX[i] << "," ;
        }
        fdebug.close();
        fdebug.open("/home/jesse/rism/bspline_cppY", std::ios::out | std::ios::trunc);
        for (int i = 0; i < gridDim[1]; i++)
        {
            fdebug << bsplineFourierCoeffY[i] << "," ;
        }
        fdebug.close();
        fdebug.open("/home/jesse/rism/bspline_cppZ", std::ios::out | std::ios::trunc);
        for (int i = 0; i < gridDimZ_k; i++)
        {
            fdebug  << bsplineFourierCoeffZ[i] << ",";
        }
        fdebug.close();
    }    
    
    // Compute the discrete Fourier transform coefficients of a
    // Gaussian monopole.
    std::vector<double> gaussianFourierCoeffX(gridDim[0]), gaussianFourierCoeffY(gridDim[1]);
    std::vector<double> gaussianFourierCoeffZ(gridDimZ_k);
    double t = -0.25 / zeta;
    for (int i = 0; i < gridDim[0]; ++i)
    {
        gaussianFourierCoeffX[i] = std::exp(kx[i] * kx[i] * t);
    }
    for (int i = 0; i < gridDim[1]; ++i)
    {
        gaussianFourierCoeffY[i] = std::exp(ky[i] * ky[i] * t);
    }
    for (int i = 0; i < gridDimZ_k; ++i)
    {
        gaussianFourierCoeffZ[i] = std::exp(kz[i] * kz[i] * t);
    }

    if (DEBUG)
    {
        fdebug.open("/home/jesse/rism/gauss_cpp", std::ios::out | std::ios::trunc);
        for (int i = 0, iz = 0; iz < gridDimZ_k; ++iz)
        {
            for (int iy = 0; iy < gridDim[1]; ++iy)
            {
                for (int ix = 0; ix < gridDim[0]; ++ix)
                {
                    double k2 = kx[ix] * kx[ix] + ky[iy] * ky[iy] + kz[iz] * kz[iz];
                    double gxyz = gaussianFourierCoeffX[ix] * gaussianFourierCoeffY[iy] * gaussianFourierCoeffZ[iz];
                    fdebug  << gxyz << ",";
                }
            }
        }
        fdebug.close();
    }
    

    // Reciprocal space kernel with b-spline discrete Fourier
    // transform correction.
    std::vector<double>* kernel = new std::vector<double>(numGridPointsK);
    int ixyz = 0;
    if (DEBUG) fdebug.open("/home/jesse/rism/kernel_cpp", std::ios::out | std::ios::trunc);
    for (int ix = 0; ix < gridDim[0]; ++ix)
    {
        double kx2 = kx[ix] * kx[ix];
        for (int iy = 0; iy < gridDim[1]; ++iy)
        {
            double kx2ky2 = kx2 + ky[iy] * ky[iy];
            double gxy = gaussianFourierCoeffX[ix] * gaussianFourierCoeffY[iy];
            double bxy = bsplineFourierCoeffX[ix] * bsplineFourierCoeffY[iy];
            for (int iz = 0; iz < gridDimZ_k; ++iz, ++ixyz)
            {
                double kz2 = kz[iz] * kz[iz];
                double gxyz = gxy * gaussianFourierCoeffZ[iz];
                double bxyz = bxy * bsplineFourierCoeffZ[iz];
                double k2 = kx2ky2 + kz2;
                (*kernel)[ixyz] = (ccdl::FOUR_PI / k2) * gxyz / bxyz;
                // Note that bxyz is not squared here like it would normally be.

                if (DEBUG) fdebug  << (*kernel)[ixyz] << ",";
            }
        }
    }    
    // Remove the k = 0 term (tinfoil boundary conditions).
    (*kernel)[0] = 0.;
    if (DEBUG) fdebug.close();

    // Allocate the FFT.
    double * grid = fftw_alloc_real(numGridPoints);

    std::complex<double> * discreteFourierCoeff
        = reinterpret_cast< std::complex<double> * >
        (fftw_alloc_complex(numGridPointsK));

    fftw_plan * fplan = new fftw_plan
        (fftw_plan_dft_r2c_3d
          (gridDim[0], gridDim[1], gridDim[2],
            grid,
            reinterpret_cast< fftw_complex *>(discreteFourierCoeff),
            FFTW_ESTIMATE));

    fftw_plan * rplan = new fftw_plan
        (fftw_plan_dft_c2r_3d
          (gridDim[0], gridDim[1], gridDim[2],
            reinterpret_cast< fftw_complex *>(discreteFourierCoeff),
            grid,
            FFTW_ESTIMATE));

    
    // Initialize.
    for (int i = 0; i < numAtoms; ++i)
    {
        potential[i] = 0.;
    }
    
    for (int i = 0; i < numGridPoints; ++i)
    {
        grid[i] = 0.;
    }

    // std::cout << "unitCellVectorsK" << std::endl;
    // for (int i = 0; i < 3; i++)
    // {
    //     std::cout << i << ": "
    //               << unitCellVectorsK[i + 0 * 3]
    //               << " " << unitCellVectorsK[i + 1 * 3]
    //               << " " << unitCellVectorsK[i + 2 * 3]
    //               << std::endl;
    // }
    // std::cout << "positions" << std::endl;
    // for (int a = 0; a < numAtoms; a++)
    // {
    //     std::cout << a << ": "
    //               << positions[0 + a * 3]
    //               << " " << positions[1 + a * 3]
    //               << " " << positions[2 + a * 3]
    //               << std::endl;
    // }
        
    // Spread the charges to the grid.
    std::vector<int> gridPointsX(splineOrder), gridPointsY(splineOrder), gridPointsZ(splineOrder);
    // double atomCharge, totalCharge = 0;
    
    if (DEBUG)
    {
        fdebug.precision(6);
        fdebug2.precision(6);
        fdebug3.precision(6);
        fdebug4.precision(6);
        fdebug5.precision(6);
        fdebug6.precision(6);
        fdebug7.precision(6);
        fdebug8.precision(6);
        fdebug9.precision(6);
        fdebug.open("/home/jesse/rism/p_bsplineX_cpp", std::ios::out | std::ios::trunc);
        fdebug2.open("/home/jesse/rism/p_bsplineY_cpp", std::ios::out | std::ios::trunc);
        fdebug3.open("/home/jesse/rism/p_bsplineZ_cpp", std::ios::out | std::ios::trunc);
        fdebug4.open("/home/jesse/rism/pg_bsplineX_cpp", std::ios::out | std::ios::trunc);
        fdebug5.open("/home/jesse/rism/pg_bsplineY_cpp", std::ios::out | std::ios::trunc);
        fdebug6.open("/home/jesse/rism/pg_bsplineZ_cpp", std::ios::out | std::ios::trunc);
        fdebug7.open("/home/jesse/rism/p_bspline_cpp", std::ios::out | std::ios::trunc);
        fdebug8.open("/home/jesse/rism/pg_bspline_cpp", std::ios::out | std::ios::trunc);
        fdebug9.open("/home/jesse/rism/smear_cpp", std::ios::out | std::ios::trunc);
        fdebug << std::fixed;
    }

    for (int a = 0; a < numAtoms; ++a)
    {
#define pme_triclinic_long 1

#if pme_triclinic_long == 1
        // Project position to unit cell axes (which may be non-orthogonal).
        // std::cout << "unitCellProjections" << std::endl;
        // std::cout << "before: " << positions[0 + a * 3]
        //           << " " << positions[1 + a * 3]
        //           << " " << positions[2 + a * 3] << std::endl;
        // std::cout << "after:" << std::endl;
        
        double unitCellProjection[3];
        
        for (int index = 0; index < 3; ++index)
        {
            unitCellProjection[index] =
                unitCellVectorsK[index + 0 * 3] * positions[0 + a * 3]
                + unitCellVectorsK[index + 1 * 3] * positions[1 + a * 3]
                + unitCellVectorsK[index + 2 * 3] * positions[2 + a * 3];
            unitCellProjection[index] *= boxLen[index];
            // std::cout << unitCellProjection[index] << std::endl;
        }

        ccdl::bspline_periodic_box(unitCellProjection[0], boxLen[0], gridDim[0],
                                   splineOrder, bsplineValuesX.data(), gridPointsX.data());
        ccdl::bspline_periodic_box(unitCellProjection[1], boxLen[1], gridDim[1],
                                   splineOrder, bsplineValuesY.data(), gridPointsY.data());
        ccdl::bspline_periodic_box(unitCellProjection[2], boxLen[2], gridDim[2],
                                   splineOrder, bsplineValuesZ.data(), gridPointsZ.data());

        if (DEBUG)
        {
            for (int ib = 0; ib < splineOrder; ib++)
            {
                fdebug << bsplineValuesX[ib] << ",";
                fdebug2 << bsplineValuesY[ib] << ",";            
                fdebug3 << bsplineValuesZ[ib] << ",";
                fdebug4 << gridPointsX[ib] << ",";
                fdebug5 << gridPointsY[ib] << ",";
                fdebug6 << gridPointsZ[ib] << ",";
            }    
        }
#else        
        ccdl::bspline_periodic_box(positions[0 + a * 3], boxLen[0], gridDim[0],
                                   splineOrder, bsplineValuesX.data(), gridPointsX.data());
        ccdl::bspline_periodic_box(positions[1 + a * 3], boxLen[1], gridDim[1],
                                   splineOrder, bsplineValuesY.data(), gridPointsY.data());
        ccdl::bspline_periodic_box(positions[2 + a * 3], boxLen[2], gridDim[2],
                                   splineOrder, bsplineValuesZ.data(), gridPointsZ.data());
#endif
        for (int i = 0; i < splineOrder; ++i)
        {
            for (int j = 0; j < splineOrder; ++j)
            {
                int gridPointsXY = (gridPointsY[j] + gridPointsX[i] * gridDim[1]) * gridDim[2];
                double bxy = bsplineValuesX[i] * bsplineValuesY[j];
                for (int k = 0; k < splineOrder; ++k)
                {
                    grid[ gridPointsZ[k] + gridPointsXY ] += bxy * bsplineValuesZ[k] * charges[a];

                    if (DEBUG)
                    {
                        fdebug7 << bxy * bsplineValuesZ[k] << ",";
                        fdebug8 << gridPointsXY + gridPointsZ[k] << ",";    
                    }        
                }
            }
        }
    }
    if (DEBUG)
    {
        fdebug.close();
        fdebug2.close();
        fdebug3.close();
        fdebug4.close();
        fdebug5.close();
        fdebug6.close();
        fdebug7.close();
        fdebug8.close();
    }
    
    for (int i = 0; i < numGridPoints; ++i)
    {
        potential[i] = grid[i];
        if (DEBUG) fdebug9 << grid[i] << ",";
    }
    if (DEBUG) fdebug9.close();

    // // std::cout << "ROUND 1" << std::endl;
    // for (int iz = 0; iz < gridDim[2]; ++iz)
    // {
    //     for (int iy = 0; iy < gridDim[1]; ++iy)
    //     {
    //         for (int ix = 0; ix < gridDim[0]; ++ix)
    //         {        
    //             int ik_cpp = iz + (iy + ix * gridDim[1]) * gridDim[2];
    //             int ik_fort = ix + (iy + iz * gridDim[1]) * gridDim[0];
    //             // potential[ik_fort] = grid[ik_cpp];
    //             potential[ik_cpp] = grid[ik_cpp];
    //         }    
    //     }    
    // }
    
    // Compute the discrete Fourier transform coefficients of the
    // Gaussian charge density.
    fftw_execute(*fplan);

    // // std::cout << "ROUND 2" << std::endl;
    // for (int iz = 0; iz < gridDim[2]; ++iz)
    // {
    //     for (int iy = 0; iy < gridDim[1]; ++iy)
    //     {
    //         for (int ix = 0; ix < gridDim[0]; ++ix)
    //         {        
    //             int ik_cpp = iz + (iy + ix * gridDim[1]) * gridDim[2];
    //             int ik_fort = ix + (iy + iz * gridDim[1]) * gridDim[0];
    //             potential[ik_cpp] = discreteFourierCoeff[ik_cpp];
    //         }    
    //     }    
    // }

    // for (int k = 0; k < numGridPointsK; ++k)
    // {
    //     potential[2 * k] = discreteFourierCoeff[k].real();
    //     potential[2 * k + 1] = discreteFourierCoeff[k].imag();
    //     if (k < 10 or k > numGridPointsK - 10)
    //     {
    //         std::cout << potential[2 * k] << ", " << potential[2 * k + 1] << std::endl;
    //     }
    // }

    if (DEBUG) fdebug.open("/home/jesse/rism/r2c_cpp", std::ios::out | std::ios::trunc);
    for (int k = 0; k < numGridPointsK; ++k)
    {
        potential[2 * k] = discreteFourierCoeff[k].real();
        potential[2 * k + 1] = discreteFourierCoeff[k].imag();
        if (DEBUG)
        {
            fdebug << discreteFourierCoeff[k].real() <<  ", ";
            fdebug << discreteFourierCoeff[k].imag() <<  ", ";
        }
    }
    if (DEBUG) fdebug.close();

    // Convolute with the kernel.
    if (DEBUG) fdebug.open("/home/jesse/rism/convolution_cpp", std::ios::out | std::ios::trunc);
    for (int k = 0; k < numGridPointsK; ++k)
    {
        discreteFourierCoeff[k] *= (*kernel)[k];
        if (DEBUG)
        {
            fdebug << discreteFourierCoeff[k].real() <<  ", ";
            fdebug << discreteFourierCoeff[k].imag() <<  ", ";
        }
    }
    if (DEBUG) fdebug.close();
    
    // Evaluate the recip-space potential at the grid points.
    fftw_execute(*rplan);

    if (DEBUG)
    {
        fdebug.open("/home/jesse/rism/c2r_cpp", std::ios::out | std::ios::trunc);
        for (int i = 0; i < numGridPoints; ++i)
        {
            fdebug << grid[i] <<  ", ";
        }
        fdebug.close();
    }

    // Store in output array, scale by 1/boxVolume due to plane wave normalization
    // and apply the uniform background correction for charged systems.
    double chargeTotal = 0.;
    for (int a = 0; a < numAtoms; ++a)
    {
        chargeTotal += charges[a];
    }
    double chargeCorrection = - (ccdl::PI / zeta) * (chargeTotal / boxVolume);
    for (int i = 0; i < numGridPoints; ++i)
    {
        potential[i] = grid[i] / boxVolume + chargeCorrection;
    }

    // Deallocate the FFT.
    fftw_destroy_plan(*fplan);
    delete fplan;
    fplan = NULL;
    fftw_destroy_plan(*rplan);
    delete rplan;
    rplan = NULL;
    fftw_free(reinterpret_cast<fftw_complex * &> (discreteFourierCoeff));
    discreteFourierCoeff = NULL;
    fftw_free(grid);
    grid = NULL;

// #if 0
    
    /////////////////////////////////////////////////////////////////////////
    // For those grid points near the atoms, replace the Gaussian potential
    // with the point-charge potential.
    /////////////////////////////////////////////////////////////////////////
#define pme_short 1
#define pme_triclinic_short 0

#if pme_short == 1    
    int nlocx = std::min(2 * ((int)(cutoff / gridSpacingX) + 1), gridDim[0]);
    int nlocy = std::min(2 * ((int)(cutoff / gridSpacingY) + 1), gridDim[1]);
    int nlocz = std::min(2 * ((int)(cutoff / gridSpacingZ) + 1), gridDim[2]);
    gridPointsX.resize(nlocx);
    gridPointsY.resize(nlocy);
    gridPointsZ.resize(nlocz);
    std::vector<double> locx(nlocx);
    std::vector<double> locy(nlocy);
    std::vector<double> locz(nlocz);


    for (int a = 0; a < numAtoms; ++a)
    {
    
#if pme_triclinic_short == 1
        // Project positions to unit cell vectors.
        // int beginPosIndex = a * 3;
        // int endPosIndex = beginPosIndex + 3; // One past the z component.

        // double xProjection = std::inner_product(&positions[beginPosIndex],
        //                                         &positions[endPosIndex],
        //                                         &unitCellUnitVectors[0 * 3], 0.);
        // double yProjection = std::inner_product(&positions[beginPosIndex],
        //                                         &positions[endPosIndex],
        //                                         &unitCellUnitVectors[1 * 3], 0.);
        // double zProjection = std::inner_product(&positions[beginPosIndex],
        //                                         &positions[endPosIndex],
        //                                         &unitCellUnitVectors[2 * 3], 0.);

        double unitCellProjections[3];
        
        for (int index = 0; index < 3; ++index)
        {
            unitCellProjections[index] =
                unitCellUnitVectors[index + 0 * 3] * positions[0 + a * 3]
                + unitCellUnitVectors[index + 1 * 3] * positions[1 + a * 3]
                + unitCellUnitVectors[index + 2 * 3] * positions[2 + a * 3];
        }
        
        
        // Which grid points are nearby?
        // ccdl::nearby_grid_points_periodic(xProjection, boxLen[0], gridDim[0], nlocx, gridPointsX.data());
        // ccdl::nearby_grid_points_periodic(yProjection, boxLen[1], gridDim[1], nlocy, gridPointsY.data());
        // ccdl::nearby_grid_points_periodic(zProjection, boxLen[2], gridDim[2], nlocz, gridPointsZ.data());
        ccdl::nearby_grid_points_periodic(unitCellProjections[0], boxLen[0], gridDim[0], nlocx, gridPointsX.data());
        ccdl::nearby_grid_points_periodic(unitCellProjections[1], boxLen[1], gridDim[1], nlocy, gridPointsY.data());
        ccdl::nearby_grid_points_periodic(unitCellProjections[2], boxLen[2], gridDim[2], nlocz, gridPointsZ.data());
#else
        // Which grid points are nearby?
        ccdl::nearby_grid_points_periodic(positions[0 + a * 3], boxLen[0], gridDim[0], nlocx, gridPointsX.data());
        ccdl::nearby_grid_points_periodic(positions[1 + a * 3], boxLen[1], gridDim[1], nlocy, gridPointsY.data());
        ccdl::nearby_grid_points_periodic(positions[2 + a * 3], boxLen[2], gridDim[2], nlocz, gridPointsZ.data());
#endif
        
        // Where are those grid points located?
        // (There's a smarter way of doing this, but this is more clear)
#if 1
        // Natural grid.
        for (int i = 0; i < nlocx; ++i)
            locx[i] = gridSpacingX * gridPointsX[i];
        for (int i = 0; i < nlocy; ++i)
            locy[i] = gridSpacingY * gridPointsY[i];
        for (int i = 0; i < nlocz; ++i)
            locz[i] = gridSpacingZ * gridPointsZ[i];
#else
        // FFTW grid.
        for (int i = 0; i < nlocx; ++i)
            locx[i] = gridSpacingX * (gridPointsX[i] > nx / 2 ? gridPointsX[i] - nx : gridPointsX[i]);
        for (int i = 0; i < nlocy; ++i)
            locy[i] = gridSpacingY * (gridPointsY[i] > ny / 2 ? gridPointsY[i] - ny : gridPointsY[i]);
        for (int i = 0; i < nlocz; ++i)
            locz[i] = gridSpacingZ * (gridPointsZ[i] > nz / 2 ? gridPointsZ[i] - nz : gridPointsZ[i]);
#endif

        
        for (int i = 0; i < nlocx; ++i)
        {
#if pme_triclinic_short == 1            
            double x = ccdl::wrap(locx[i] - xProjection, boxLen[0]);
            double xVector[3] = {
                x * unitCellUnitVectors[0 + 0 * 3],
                x * unitCellUnitVectors[0 + 1 * 3],
                x * unitCellUnitVectors[0 + 2 * 3]
            }
            // std::valarray<double> xVector ({x * unitCellUnitVectors[0 + 0 * 3],
            //         x * unitCellUnitVectors[1 + 0 * 3],
            //         x * unitCellUnitVectors[2 + 0 * 3]
            //             }, 3);
#else
            double x = ccdl::wrap(locx[i] - positions[0 + a * 3], boxLen[0]);
            double x2 = x * x;
#endif
            for (int j = 0; j < nlocy; ++j)
            {
#if pme_triclinic_short == 1
                double y = ccdl::wrap(locy[j] - yProjection, boxLen[1]);
                double yVector[3] = {
                    y * unitCellUnitVectors[1 + 0 * 3],
                    y * unitCellUnitVectors[1 + 1 * 3],
                    y * unitCellUnitVectors[1 + 2 * 3]
                };
#else
                double y = ccdl::wrap(locy[j] - positions[1 + a * 3], boxLen[1]);
                double x2y2 = x2 + y * y;
#endif
                
                int gridPointsXY = (gridPointsY[j] + gridPointsX[i] * gridDim[1]) * gridDim[2];
                for (int k = 0; k < nlocz; ++k)
                {
#if pme_triclinic_short == 1
                    double z = ccdl::wrap(locz[k]  - zProjection, boxLen[2]);
                    double r2 = 0;
                    double zVector[3] = {
                        z * unitCellUnitVectors[2 + 0 * 3],
                        z * unitCellUnitVectors[2 + 1 * 3],
                        z * unitCellUnitVectors[2 + 2 * 3]
                    };
                    for (int axis = 0; axis < 3; axis++)
                    {
                        z = xVector[axis] + yVector[axis] + zVector[axis];
                        z *= z;
                        r2 += z;
                    }
#else
                    double z = ccdl::wrap(locz[k]  - positions[2 + a * 3], boxLen[2]);
                    double r2 = x2y2 + z * z;
#endif
                    
                    double e = 0.;

                    if (r2 > 1.e-8)
                    {
                        double r = std::sqrt(r2);
                        e = charges[a] * erfc(smear * r) / r;
                    }
                    else // Remove the Gaussian potential, only.
                    {
                        e = - charges[a] * std::sqrt(4 * zeta / ccdl::PI);
                    }

                    potential[ gridPointsZ[k] + gridPointsXY ] += e;

                }
            }
        }
    }
#endif // pme_short == 1
// #endif // debug
}

// void testcpp_
// (int const & numAtoms, double const * positions, double const * charges,
//  int const * gridDim, double const * boxLen, double const & boxVolume,
//  double const * unitCellUnitVectors, double const * unitCellVectors,
//  double const * unitCellVectorsK,
//  double const & smear, int const & splineOrder, double const & cutoff,
//  double * potential)
// {
//     std::vector<double> bsplineValuesX(splineOrder);
//     ccdl::bspline((splineOrder % 2 == 0 ? 0.0 : 0.5), splineOrder, bsplineValuesX.data());
//     for (auto value : bsplineValuesX)
//     {
//         std::cout << value << std::endl;
//     }
// }
