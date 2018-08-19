#ifndef _pmepot_hpp_
#define _pmepot_hpp_

#include <iostream>

// Depends on -lfftw3.

/** Evaluates the electrostatic potential of a bunch of points and
    returns the values at the grid points.

    pmepot_ is evaluating Eq. (44), (45), and line 1 of Eq. (46) 
    (at r=Rt) in the "pmepot.pdf" file.

    You probably don't want to use FFTW grid points, but rather natural ones.

    Caveats:
    (1) gridDimX,gridDimY,gridDimZ must all be EVEN numbers (only because I'm "lazy").
    (2) All input and output variables are assumed to be ATOMIC UNITS.
        (I have no idea what crazy units you use, so I won't even
         bother trying to use "correct" conversion constants.)
        The charge of a proton is +1.
    (3) This code only works for orthorhombic unit cells (again, "lazy").
    (4) Requires fftw3
    (5) This code was written in ~3 hours and has undergone 
        limited testing and optimization (it was written for clarity).
    (6) The location of the 'natural' grid points are:
           locx[i] = (boxLenX/gridDimX) * i
    (6*) The location of the FFTW grid points follows the convention used
        by FFTW (because FFTW outputs the potential from a DFT)
           locx[i] = (lx/nx) * (i > nx/2 ? i-nx : i)
        One would normally not need to worry about this as long as
        the grid was used consistently and the end result was
        interpolated back onto the atoms; however, you WILL need to
        worry about it because you are using the grid values directly.
    (7) Tinfoil boundary conditions are used.
    (8) A uniform background correction -(pi/beta^2)(Qtot/V)
        is applied to the potential to make it return the same answer
        you would if the grid points were extra atoms with zero charge.

    @param[in] numAtoms Number of atoms.
    @param[in] positions Atomic positions [3*numAtoms].
    @param[in] charges Atomic charges [numAtoms].
    @param[in] gridDimX, gridDimY, gridDimZ Number of grid points.
    @param[in] boxLenX, boxLenY, boxLenZ Unit cell side lengths.
    @param[in] boxAngleA, boxAngleB, boxAngleC Unit cell interior angles.
    @param[in] unitCellUnitVectors Unit axis vectors of the unit cell.
    @param[in] smear Gaussian beta (zeta=beta*beta).
    @param[in] splineOrder B-spline order.
    @param[in] cutoff Real-space cutoff.
    @param[out] potential Potential at each grid point [gridDimX*gridDimY*gridDimZ].
*/
extern "C"
{
    void pmepot_
    ( int const & numAtoms,
      double const * positions,
      double const * charges,
      int const * gridDim,
      double const * boxLen,
      double const & boxVolume,
      double const * unitCellUnitVectors,
      double const * unitCellVectors,
      double const * unitCellVectorsK,
      double const & smear,
      int const & splineOrder,
      double const & cutoff,
      double * potential );

    // void testcpp_
    // ( int const & numAtoms,
    //   double const * positions,
    //   double const * charges,
    //   int const * gridDim,
    //   double const * boxLen,
    //   double const & boxVolume,
    //   double const * unitCellUnitVectors,
    //   double const * unitCellVectors,
    //   double const * unitCellVectorsK,
    //   double const & smear,
    //   int const & splineOrder,
    //   double const & cutoff,
    //   double * potential );
}

// Utility functions for particle mesh Ewald.
namespace ccdl
{
    //TODO: Use PI from Amber.
    double const PI = 3.141592653589793238462643383279502884197;
    double const TWO_PI = 2. * PI;
    double const FOUR_PI = 4. * PI;

    /** Calculate the vector cross product of two 3-element arrays.
     */
    double * const cross_product(double const * const a, double const * const b);

    /// Calculate the vector dot product of two 3-element arrays.
    double dot_product(double const * const a, double const * const b);

    /// Unit test for the vector products.
    void vector_product_test();
    
    /// Round a real number towards nearest whole number.
    double anint(double const x);
    
    /// Wrap a 1D position at x based on a periodic line of length l.
    double wrap(double const x, double const l);

    /** Evaluate the recursive definition of the 1D normalized
        cardinal b-spline of given order at an offset within each knot
        interval. Values passed in are from evaluation of the b-spline
        at order - 1.

        @param[in] offset 
        Value between 0 and 1. See full description in ::bspline_eval.
        @param[in] order 
        The number of knot intervals in the evaluated b-spline.
        @param[in,out] values
        In: The values of the b-spline of order N - 1 at the N - 1
        evaluation points.
        Out: The values of the b-spline of order N at the N evaluation
        points. The first value in array is the leftmost point.
    */
    void bspline_recursion(double const offset, int const order, double * values);
    
    /** Evaluate a 1D normalized cardinal b-spline of given order at a
        fractional offset within each knot interval. Assumes order is
        2 or greater.

        A b-spline is evaluated on a grid, where each grid point is
        called a knot and the interval between neighboring knots is
        called a knot interval.

        A cardinal b-spline has a constant separation between
        neighboring knots (i.e., knot intervals have a fixed size). It
        is referred to as the the Irwin-Hall distribution in
        statistics.

        @param[in] offset
        Fractional offset of the b-spline evaluation point, valued
        between 0 and 1.

        If the order is even, then the offset is the percentage of the
        coordinate of the nearest knot to the right of the evaluation
        point, such that:
        
         offset * crd_nearest_right_knot + (1 - offset) * crd_nearest_left_knot
            = evaluation point position

        For example, for a pair of neighboring left and right knots
        with positions 1 and 2 respectively, an offset of 0.75 implies
        that the b-spline is being evaluated at 0.75 * 2 + 0.25 * 1 =
        1.75.

        If the order is odd, then the offset is a percentage of the
        coordinate of the knot interval midpoint to the right of the
        evaluation point. Thus for an offset of 0.5, the evaluation
        point is directly on a knot.

        In this case, for a set of three neighboring knots with
        positions 0, 1, and 2 (with implied knot interval midpoints
        0.5 and 1.5), an offset of 0.75 implies that the b-spline is
        being evaluated at 0.75 * 1.5 + 0.25 * 0.5 = 1.25.

        @param[in] order
        The number of knot intervals in the b-spline. Assumed to be 2
        or greater. For a cardinal b-spline, the order is one more
        than the degree and one less than the number of knots (when
        including the leftmost knot at 0).
        @param[out] values
        The values of the b-spline at the N evaluation points (where N =
        order). The first value in the array is the leftmost point.
    */
    void bspline(double const offset, int const order, double * values);
    
    /** Evaluates the 1D cardinal b-spline centered at x along a 1D
        periodic box with evenly spaced grid points. Evaluation points
        are at the nearest N grid points from x (where N = order).
        Also returns the grid indices of the evaluation points,
        considering periodicity such that all grid points are kept
        within the reference box.

        @param[in] x
        Cartesian center position of the b-spline.
        @param[in] boxLen
        Cartesian length of the 1D box.
        @param[in] gridDim
        Number of evenly spaced grid points along the 1D box.
        @param[in] order
        The number of knot intervals in the b-spline. Each interval
        contains one grid point.
        @param[out] values
        The values of the b-spline at the N evaluation points (where N =
        order). The first value in array is the leftmost point.
        @param[out] nearestGridPoints
        Grid indices of the evaluation points, considering periodicity.
    */
    void bspline_periodic_box(double const x, double const boxLen,
                              int const gridDim, int const order,
                              double * values, int * nearestGridPoints);
    
    /** Indexes the "order" nearest grid points to a particle located
        at x on a 1D periodic box.

        @param[in] x
        Cartesian particle position.
        @param[in] boxLen
        Length of the 1D periodic box.
        @param[in] gridDim
        Number of grid points in the 1D periodic box.
        @param[in] order
        Number of knot intervals in the b-spline.
        @param[out] nearestGridPoints 
        Array of nearst grid point indices, size equal to order.
    */
    void nearby_grid_points_periodic(double const x, double const boxLen, int const gridDim,
                                     int const order, int * nearestGridPoints);

    /** Computes the discrete Fourier transform coefficients of a 1D
        cardinal b-spline (all frequencies - the outer loops in FFTW).

        @param[in] maxIndex
        Compute all Fourier coefficients up to the max index.
        @param[in] order
        Number of knot intervals in the b-spline.
        @param[in] bsplineValues
        Real-space b-spline values, size equal to order.
        @param[out] fourierCoefficients
        Array of b-spline Fourier coefficients, size maxIndex.
    */
    void bspline_Fourier_coefficients(
        int const maxIndex, int const order,
        double const * bsplineValues, double * fourierCoefficients);
    
    /** Computes the discrete Fourier transform coefficients of a 1D
        cardinal b-spline (only half the frequencies - the inner loop
        in FFTW).

        For parameter descriptions, see ::bspline_Fourier_coefficients.
    */
    void bspline_half_Fourier_coefficients(
        int const maxIndex, int const order,
        double const * bsplineValues, double * fourierCoefficients);
}

#endif
