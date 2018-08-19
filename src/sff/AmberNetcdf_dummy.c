#include "AmberNetcdf.h"

// netcdfLoad()
/** Load the netcdf trajectory/restart file pointed to by filename and 
  * initialize the AmberNetcdf data structure. Allocate memory for coords 
  * (and velo if restart). Set currentFrame to 0.
  * \return 0 on success, 1 on error.
  */
int netcdfLoad(struct AmberNetcdf *A, char *filename) {
  return 0;
} 

// netcdfClose()
/** Close netcdf trajectory/restart file and free memory */
int netcdfClose(struct AmberNetcdf *A) {
  return 0;
}

// netcdfWriteRestart()
/** Create and write to an Amber netcdf restart file with name
  * specified by filename.
  * \param filename File name of netcdf restart
  * \param natom Number of atoms
  * \param X Coordinates
  * \param Velo If not NULL write velocity info
  * \param box If not NULL write box info
  * \param restart_time Restart time
  * \param restart_temp Restart temperature
  * \return 0 on successful write, 1 on error.
  */
int netcdfWriteRestart(char *filename, int natom, double *X, double *Velo, 
                       double *box, double restart_time, double restart_temp) 
{
  return 0;
}

// netcdfCreate()
/** Create an Amber netcdf trajectory for the specified #atoms with name 
  * filename. Also add box information if requested.
  */
int netcdfCreate(struct AmberNetcdf *A, char *filename, int natom, int isBox) {
  return 0;
}

// netcdfGetFrame()
/** Get the specified frame from amber netcdf trajectory/restart file.
  * Also get the box coords if present and box is not NULL.
  * Also get temperature T if TempVID defined and T is not NULL.
  * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
  * \return 0 on success, 1 on error.
  */
int netcdfGetFrame(struct AmberNetcdf *A, int set, double *X, double *box) {
  return 0;
}

// netcdfGetNextFrame()
/** Get the frame at currentFrame and increment. Return 0 if no more frames.
  */
int netcdfGetNextFrame(struct AmberNetcdf *A, double *X, double *box) {
  return 0;
}

// netcdfWriteFrame()
/** Write coords (and box if specified) to netcdf trajectory file. 
  * \return 0 on success, 1 on error. 
  */
int netcdfWriteFrame(struct AmberNetcdf *A, int set, double *X, double *box) {
  return 0;
}  

// netcdfWriteNextFrame()
/** Write coords to currentFrame. Increment currentFrame.
  * Return 1 on error.
  */
int netcdfWriteNextFrame(struct AmberNetcdf *A, double *X, double *box) {
  return 0;
}
  
// netcdfInfo()
/** Print general information about the netcdf file.
  */
int netcdfInfo(struct AmberNetcdf *A) {
  return 0;
}
