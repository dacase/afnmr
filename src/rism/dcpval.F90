      DOUBLE PRECISION             FUNCTION DCPVAL (P,NDEGP,X)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
! Copyright © 1996 California Institute of Technology, Pasadena,
! California. ALL RIGHTS RESERVED. Based on Government Sponsored Researc
! NAS7-03001. Redistribution and use in source and binary forms, with or
! without modification, are permitted provided that the following
! conditions are met: Redistributions of source code must retain this
! copyright notice, this list of conditions and the following
! disclaimer. Redistributions in binary form must reproduce the above
! copyright notice, this list of conditions and the following disclaimer
! in the documentation and/or other materials provided with the
! distribution. Neither the name of the California Institute of Technolo
! (Caltech) nor the names of its contributors may be used to endorse or
! promote products derived from this software without specific prior
! written permission. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILI
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIREC
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, B
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
! USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TOR
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  For
! those codes indicated with a Math a la Carte copyright, the same rules
! apply, except without the full force of the Caltech legal team.      
!>> 2018-06-13 DSI Luchko Change to double precision and converted to Fortran90
!>> 1994-10-20 SCPVAL Krogh  Changes to use M77CON
!>> 1994-04-20 SCPVAL CLL Edited to make DP & SP files similar.
!>> 1987-12-09 SCPVAL Lawson  Initial code.
!--S replaces "?": ?CPVAL
!
!     C.L.LAWSON,JPL, 1969 DEC 17   MODIFIED 1973 JULY 24
!
!     MODIFIED 1974 NOV 19
!
!     EVALUATE A POLYNOMIAL OF DEGREE NDEGP GIVEN TRANSFORMATION
!     PARAMETERS, P(1) AND P(2), AND COEFFICIENTS RELATIVE TO THE
!     CHEBYSHEV BASIS.
!
!     NDEGP               DEGREE OF POLYNOMIAL
!     (P(I),I=1,NDEGP+3)  PARAMETERS DEFINING THE POLYNOMIAL
!     X                  INPUT ARGUMENT
!             THE POLYNOMIAL'S VALUE AT X IS DEFINED AS FOLLOWS.
!
!                             S = ( X - P(1) ) / P(2)
!
!                SCPVAL=SUM OF P(I+3)*T(I,S) FOR I=0,1,...NDEGP
!
!                             WHERE T(I,S) DENOTES THE CHEBYSHEV
!                             POLYNOMIAL OF DEGREE I EVALUATED AT S .
!
      integer J, NDEGP
      DOUBLE PRECISION  P(*),W(3),S,S2,X
      W(1)=0.0D0
      W(2)=0.0D0
!                             TRANSFORM X TO S
      S=(X-P(1))/P(2)
      S2=S+S
!
!                             EVALUATE POLYNOMIAL USING RECURSION
!
      do 30 J = NDEGP+3, 4, -1
          W(3)=W(2)
          W(2)=W(1)
          W(1)=(S2*W(2)-W(3))+P(J)
   30 continue
      DCPVAL=(S*W(1)-W(2))+P(3)
      RETURN
      END
