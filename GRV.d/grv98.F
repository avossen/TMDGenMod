
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* Oct 2010, S. Gliske
*
* Originally downloaded form
* http://durpdg.dur.ac.uk/pdflib/grv/grv98/grv98.f
* Updates:
*    Package shortened to just parton routine
*    Changed calling routine, so pass init flag, and iok flag
*    Also only return requested flavor, in order
*    u = 0, d = 1, ub = 2, db = 3, s = sb = 4, gl = 5 
*    Changed name of routine to GRV98
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                   *
*     G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S      *
*                                                                   *
*                          1998 UPDATE                              *
*                                                                   *
*                  For a detailed explanation see                   *
*                   M. Glueck, E. Reya, A. Vogt :                   *
*        hep-ph/9806404  =  DO-TH 98/07  =  WUE-ITP-98-019          *
*                  (To appear in Eur. Phys. J. C)                   *
*                                                                   *
*   This package contains subroutines returning the light-parton    *
*   distributions in NLO (for the MSbar and DIS schemes) and LO;    * 
*   the respective light-parton, charm, and bottom contributions    *
*   to F2(electromagnetic); and the scale dependence of alpha_s.    *
*                                                                   *
*   The parton densities and F2 values are calculated from inter-   *
*   polation grids covering the regions                             *
*         Q^2/GeV^2  between   0.8   and  1.E6 ( 1.E4 for F2 )      *
*            x       between  1.E-9  and   1.                       *
*   Any call outside these regions stops the program execution.     *
*                                                                   *
*   At Q^2 = MZ^2, alpha_s reads  0.114 (0.125) in NLO (LO); the    *
*   heavy quark thresholds, QH^2 = mh^2, in the beta function are   *
*            mc = 1.4 GeV,  mb = 4.5 GeV,  mt = 175 GeV.            *
*   Note that the NLO alpha_s running is different from GRV(94).    * 
*                                                                   *
*    Questions, comments etc to:  avogt@physik.uni-wuerzburg.de     *
*                                                                   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*
*
      SUBROUTINE GRV98 (IDOINIT, ISET, X, Q2, IFLAV, VAL, IERR )
*********************************************************************
*                                                                   *
*   THE PARTON ROUTINE.                                             *
*                                     __                            *
*   INPUT:   IDOINIT = 0, don't initialize, just evaluate
*            IDOINIT != 0, only intialize, don't evaluate
*            ISET =  1 (LO),  2 (NLO, MS), or  3 (NLO, DIS)         *
*            X  =  Bjorken-x        (between  1.E-9 and 1.)         *
*            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)         *
*            IFLAV: u = 0, d = 1, ub = 2, db = 3,
*                   s = sb = 4, gl = 5 
*                                                                   *
*   OUTPUT:  VAL = value for given flavor
*            IERR = error code
*            Always x times the distribution is returned.           *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINIP / IINIP , and the     *
*            integer variable  IINIP  has always to be zero when    *
*            GRV98PA is called for the first time or when  ISET     *
*            has been changed.                                      *
*                                                                   *
*   GRIDS:   1. grv98lo.grid, 2. grv98nlm.grid, 3. grv98nld.grid,   *
*            (1+1809 lines with 6 columns, 4 significant figures)   *
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NPART=6, NX=68, NQ=27, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XDEF(NX,NQ), XUDF(NX,NQ),
     1     XSF(NX,NQ), XGF(NX,NQ), PARTON (NPART,NQ,NX-1), 
     2     QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      CHARACTER*80 LINE
      SAVE XUVF, XDVF, XDEF, XUDF, XSF, XGF, NA, ARRF
*     
*...  BJORKEN-X AND Q**2 VALUES OF THE GRID :
      DATA QS / 0.8E0, 
     1     1.0E0, 1.3E0, 1.8E0, 2.7E0, 4.0E0, 6.4E0,
     2     1.0E1, 1.6E1, 2.5E1, 4.0E1, 6.4E1,
     3     1.0E2, 1.8E2, 3.2E2, 5.7E2,
     4     1.0E3, 1.8E3, 3.2E3, 5.7E3,
     5     1.0E4, 2.2E4, 4.6E4,
     6     1.0E5, 2.2E5, 4.6E5, 
     7     1.E6 /
      DATA XB / 1.0E-9, 1.8E-9, 3.2E-9, 5.7E-9, 
     A     1.0E-8, 1.8E-8, 3.2E-8, 5.7E-8, 
     B     1.0E-7, 1.8E-7, 3.2E-7, 5.7E-7, 
     C     1.0E-6, 1.4E-6, 2.0E-6, 3.0E-6, 4.5E-6, 6.7E-6,
     1     1.0E-5, 1.4E-5, 2.0E-5, 3.0E-5, 4.5E-5, 6.7E-5,
     2     1.0E-4, 1.4E-4, 2.0E-4, 3.0E-4, 4.5E-4, 6.7E-4,
     3     1.0E-3, 1.4E-3, 2.0E-3, 3.0E-3, 4.5E-3, 6.7E-3,
     4     1.0E-2, 1.4E-2, 2.0E-2, 3.0E-2, 4.5E-2, 0.06, 0.08,
     5     0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     6     0.3, 0.325, 0.35, 0.375, 0.4,  0.45, 0.5, 0.55,
     7     0.6, 0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1. /
*     
*     check whether to inialize
*     
      IF (IDOINIT.EQ.0) GOTO 16
*     
*     assume error until proven otherwise
*     
      IERR = 1
*     
*     SELECTION AND READING OF THE GRID : 
*     (COMMENT: FIRST NUMBER IN THE FIRST LINE OF THE GRID)
      IF (ISET .EQ. 1) THEN
         OPEN (11,FILE='grv98lo.grid',STATUS='old',ERR=100) ! 7.332E-05
      ELSE IF (ISET .EQ. 2) THEN
         OPEN (11,FILE='grv98nlm.grid',STATUS='old',ERR=110) ! 1.015E-04
      ELSE IF (ISET .EQ. 3) THEN
         OPEN (11,FILE='grv98nld.grid',STATUS='old',ERR=120) ! 1.238E-04
      ELSE
         WRITE(6,93)
 93      FORMAT (2X,'NO OR INVALID PARTON SET CHOICE')
*     
*     goto return statement
*     
         GOTO 60
      END IF

      READ(11,89) LINE
 89   FORMAT(A80)
      DO 160 M = 1, NX-1 
         DO 15 N = 1, NQ
            READ(11,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1           PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M) 
 90         FORMAT (6(1PE10.3))
 15      CONTINUE
 160  CONTINUE
      CLOSE(11)
*     
*.... ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
         DO 20 IX = 1, NX-1
            XB0V = XB(IX)**0.5 
            XB0S = XB(IX)**(-0.2) 
            XB1 = 1.-XB(IX)
            XUVF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0V)
            XDVF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0V)
            XDEF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0V) 
            XUDF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0S)
            XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**7 * XB0S)
            XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0S)
 20      CONTINUE
         XUVF(NX,IQ) = 0.E0
         XDVF(NX,IQ) = 0.E0
         XDEF(NX,IQ) = 0.E0
         XUDF(NX,IQ) = 0.E0
         XSF(NX,IQ)  = 0.E0
         XGF(NX,IQ)  = 0.E0
 10   CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
         ARRF(IX) = DLOG(XB(IX))
 30   CONTINUE
      DO 40 IQ = 1, NQ
         ARRF(NX+IQ) = DLOG(QS(IQ))
 40   CONTINUE
*     
*     If this was an initialization run, then no need to proceed
*     
      IERR=0
      GOTO 60
*     
*...  if not initialization run, continue on from here
*
 16      CONTINUE
*
* assume all is OK
*
      IERR=0
*
*...CHECK OF X AND Q2 VALUES : 
*
      IF ( (X.LT.0.99D-9) .OR. (X.GT.1.D0) ) THEN
         WRITE(6,91) 
  91     FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
         IERR=3
         VAL=0
         GOTO 60
      ENDIF
      IF ( (Q2.LT.0.799) .OR. (Q2.GT.1.01E6) ) THEN
         WRITE(6,92) 
  92     FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
         IERR=4
         VAL=0
         GOTO 60
      ENDIF

*            UV = u - u(bar),  DV = d - d(bar),  US = u(bar),       *
*            DS = d(bar),  SS = s = s(bar),  GL = gluon.            *
*    Also only return requested flavor, in order
*    u = 0, d = 1, ub = 2, db = 3, s = sb = 4, gl = 5 

*
*...INTERPOLATION :
*
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)

      IF ( (IFLAV.EQ.0) .OR. (IFLAV.EQ.2) ) THEN
*
* u or ubar
*
         X1 = 1.- X
         XV = X**0.5
         UV = FINT(NARG,XT,NA,ARRF,XUVF) * X1**3 * XV
         DE = FINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
         UD = FINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
         US = 0.5 * (UD - DE)

         IERR = 0
         IF( IFLAV.EQ.0 ) THEN
            VAL=UV + US
         ELSE
            VAL=US
         END IF

      ELSE IF ( (IFLAV.EQ.1) .OR. (IFLAV.EQ.3) ) THEN
*
* d or dbar
*
         X1 = 1.- X
         XV = X**0.5
         DV = FINT(NARG,XT,NA,ARRF,XDVF) * X1**4 * XV
         DE = FINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
         UD = FINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
         DS = 0.5 * (UD + DE)

         IERR = 0
         IF( IFLAV.EQ.1 ) THEN
            VAL=DV + DS
         ELSE
            VAL=DS
         END IF

      ELSE IF (IFLAV.EQ.4) THEN
*
* s or sbar
*
         X1 = 1.- X
         XS = X**(-0.2)
         SS = FINT(NARG,XT,NA,ARRF,XSF)  * X1**7 * XS
         VAL = SS
      ELSE IF (IFLAV.EQ.5) THEN
*
* gluon
*
         XS = X**(-0.2)
         GL = FINT(NARG,XT,NA,ARRF,XGF)  * X1**5 * XS 
         VAL=GL
      ELSE
*
* all else
*
         VAL=0
         IERR=10
      ENDIF
*
* done with main body of code
*
 60   RETURN
*
*
* extra error messages
*
 100  WRITE(6,101)
 101  FORMAT (2X,'ERROR OPENING GRV FILE: grv98lo.grid' )
      GOTO 130

 110  WRITE(6,111)
 111  FORMAT (2X,'ERROR OPENING GRV FILE: grv98nlm.grid' )
      GOTO 130

 120  WRITE(6,121)
 121  FORMAT (2X,'ERROR OPENING GRV FILE: grv98nld.grid' )
      GOTO 130

 130  IERR=2
      RETURN

      END
