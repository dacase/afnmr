c patterned after locmem in amber code, in order to integrate into that
c--------------------------------------------------------
      subroutine erfc_table_mem(startreal,endreal,
     $        mxerftab,lerf_arr,ltau,ew_coeff,cutoffnb,erftbdns)
c-----------------------------------------------------------------------------
      implicit none
      integer startreal,endreal
      integer mxerftab,lerf_arr,ltau
      double precision ew_coeff,cutoffnb,erftbdns

      integer mem_ptr

c erf table: assume all nonbond distances are less than 1.5 X cutoff
c between nonbond updates; i.e. no excess motion
      mxerftab = ew_coeff*erftbdns*cutoffnb*1.5
      write(6,*)'computed size of erftable = ',mxerftab

c do real array offsets
      mem_ptr = startreal

c permanent or heap REAL storage
      call adj_mem_ptr(mem_ptr,lerf_arr,4*mxerftab)

c stack REAL storage
      call adj_mem_ptr(mem_ptr,ltau,mxerftab)

      endreal =  mem_ptr

      return
      end
c---------------------------------------------------------------------------
