#AMBER BLAS CMake  buildfile
project(BLAS)

set(BLAS_SOURCES dasum.f daxpy.f dcabs1.f dcopy.f ddot.f 
         dgbmv.f dgemm.f dgemv.f dger.f dnrm2.f drot.f 
         drotg.f dsbmv.f dscal.f dspmv.f dspr.f dspr2.f 
         dswap.f dsymm.f dsymv.f dsyr.f dsyr2.f dsyr2k.f 
         dsyrk.f dtbmv.f dtbsv.f dtpmv.f dtpsv.f dtrmm.f 
         dtrmv.f dtrsm.f dtrsv.f dummy.f dzasum.f dznrm2.f 
         idamax.f isamax.f lsame.f xerbla.f 
         scopy.f sswap.f sscal.f sasum.f snrm2.f sdot.f ssyr2k.f
         saxpy.f sgemv.f ssyr2.f ssymv.f sgemm.f strmm.f strmv.f sger.f 
         zher2k.f zswap.f zgemv.f zhemv.f zscal.f zdotc.f zaxpy.f 
         zher2.f zdscal.f ztrmv.f zcopy.f ztrmm.f zgemm.f zgerc.f)

add_compile_options(${OPT_FFLAGS})

add_library(blas ${BLAS_SOURCES})
install_libraries(blas)

