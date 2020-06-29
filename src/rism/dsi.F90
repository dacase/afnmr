      DOUBLE PRECISION FUNCTION DSI (X)
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
!>> 1998-10-29 SSI Krogh  Moved external statement up for mangle.
!>> 1995-11-03 SSI Krogh  Removed blanks in numbers for C conversion.
!>> 1994-11-11 SSI Krogh  Declared all vars.
!>> 1994-10-20 SSI Krogh  Changes to use M77CON
!>> 1989-03-14 SSI Original W. V. Snyder at JPL
!
!     COMPUTE THE SINE INTEGRAL OF X =
!     INTEGRAL FROM 0 TO X OF (SIN(T)/T DT).
!
!     FOR ABS(X)<16, USE A CHEBYSHEV SERIES WITH ARGUMENT 2*Z*Z-1 WHERE
!     Z=X/16 TO EVALUATE SI(X)/Z, THEN MULTIPLY THE RESULT BY Z.  THIS
!     AVOIDS STORING ZERO COEFFICIENTS FOR EVEN ORDERS, AND PRESERVES
!     ACCURACY FOR SMALL Z.
!
!     FOR 16.LE.ABS(X).LE.100, USE CHEBYCHEV SERIES WITH ARGUMENT
!     2*Z*Z-1, WHERE Z=16/X ARE USED TO COMPUTE F(X)/X AND G(X)/(X*X).
!     THEN SI(X)=0.5*PI*SIGN(X)-F(X)/X*COS(X)-G(X)/(X*X)*SIN(X).
!
!     WHEN X.GT.100, USE ASYMPTOTIC APPROXIMATIONS FOR F(X)/X AND
!     G(X)/(X*X) AND COMPUTE SI(X) AS ABOVE.
!
!     THIS ALGORITHM YIELDS AT MOST 15 DIGITS OF PRECISION.
!
!--S replaces "?": ?SI, ?CPVAL
!
      INTEGER N
      DOUBLE PRECISION X
      EXTERNAL DCPVAL
      DOUBLE PRECISION PI2,Z,ZW,FZ,GZ,DCPVAL
      DOUBLE PRECISION FT,GT
      DOUBLE PRECISION S(23),F(13),G(13)
      DATA PI2/1.57079632679489662d0/
      DATA S/&
     & + 0.5D0,                 + 0.5D0,&
     & + 4.052926477680623D0, - 4.063980844911986D0,&
     & + 2.778756381742663D0, - 1.926565091150656D0,&
     & + 1.389308771171888D0, - 0.968322236987086D0,&
     & + 0.530148847916522D0, - 0.211263780976555D0,&
     & + 0.062033679432003D0, - 0.013867445589417D0,&
     & + 0.002436221404749D0, - 0.000345469155569D0,&
     & + 0.000040420271419D0, - 0.000003972908746D0,&
     & + 0.000000332988589D0, - 0.000000024100076D0,&
     & + 0.000000001522370D0, - 0.000000000084710D0,&
     & + 0.000000000004185D0, - 0.000000000000185D0,&
     & + 0.000000000000007D0/
      DATA F/&
     & + 0.5D0,                 + 0.5D0,&
     & + 0.062263729028927D0, - 0.000233756041393D0,&
     & + 0.000002453755677D0, - 0.000000058670317D0,&
     & + 0.000000002356196D0, - 0.000000000136096D0,&
     & + 0.000000000010308D0, - 0.000000000000964D0,&
     & + 0.000000000000107D0, - 0.000000000000014D0,&
     & + 0.000000000000002D0/
      DATA G/&
     & + 0.5D0,                 + 0.5D0,&
     & + 0.003862856096703D0, - 0.000042644182622D0,&
     & + 0.000000724995950D0, - 0.000000023468225D0,&
     & + 0.000000001169202D0, - 0.000000000079604D0,&
     & + 0.000000000006875D0, - 0.000000000000717D0,&
     & + 0.000000000000087D0, - 0.000000000000012D0,&
     & + 0.000000000000002D0/
!
      IF (ABS(X).LT.16.0) THEN
         Z = X/16.0d0
         ZW = Z*Z
         Z = Z*DCPVAL(S,20,ZW)
      ELSE
        IF (ABS(X).LE.100.0) THEN
!           16.LE.ABS(X).LE.100
            Z = 16.0d0/X
            ZW = Z*Z
            FZ = Z*DCPVAL(F,10,ZW)
            GZ = ZW*DCPVAL(G,10,ZW)
        ELSE
!           ABS(X).GT.100
            FZ = 1.0d0/X
            FT = FZ
            GZ = FZ/X
            GT = GZ
            Z = GZ
            DO 25 N = 2, 16, 2
               FT = -DBLE(N*(N-1))*Z*FT
               GT = -DBLE(N*(N+1))*Z*GT
               FZ = FZ + FT
               GZ = GZ + GT
25          CONTINUE
         END IF
         Z = SIGN(PI2,X) - FZ*COS(X) - GZ*SIN(X)
      END IF
!
      DSI = Z
      RETURN
      END
