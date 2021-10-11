
        SUBROUTINE UMS2P1 (WHO, WHERE,
     $          N, NE, JOB, TRANS, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, LX, W, LW)
        INTEGER WHO, WHERE, N, NE, JOB, LVALUE, LINDEX, INDEX (LINDEX),
     $          KEEP (20), ICNTL (20), INFO (40), LX, LW
        REAL
     $          VALUE (LVALUE), CNTL (10), RINFO (20), B (LX),
     $          X (LX), W (LW)
        LOGICAL TRANS
        
C=== UMS2P1 ============================================================
C
C  Unsymmetric-pattern MultiFrontal Package (UMFPACK). Version 2.0s.
C  Copyright (C) 1995, Timothy A. Davis, University of Florida, USA.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  September 1995. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C***********************************************************************
C* NOTICE:  "The UMFPACK Package may be used SOLELY for educational,   *
C* research, and benchmarking purposes by non-profit organizations and *
C* the U.S. government.  Commericial and other organizations may make  *
C* use of UMFPACK SOLELY for benchmarking purposes only.  UMFPACK may  *
C* be modified by or on behalf of the User for such use but at no time *
C* shall UMFPACK or any such modified version of UMFPACK become the    *
C* property of the User.  UMFPACK is provided without warranty of any  *
C* kind, either expressed or implied.  Neither the Authors nor their   *
C* employers shall be liable for any direct or consequential loss or   *
C* damage whatsoever arising out of the use or misuse of UMFPACK by    *
C* the User.  UMFPACK must not be sold.  You may make copies of        *
C* UMFPACK, but this NOTICE and the Copyright notice must appear in    *
C* all copies.  Any other use of UMFPACK requires written permission.  *
C* Your use of UMFPACK is an implicit agreement to these conditions."  *
C*                                                                     *
C* The MA38 Package in the Harwell Subroutine Library (HSL) has        *
C* equivalent functionality (and identical calling interface) as       *
C* UMFPACK.  It is available for commercial use.   Technical reports,  *
C* information on HSL, and matrices are available via the World Wide   *
C* Web at http://www.cis.rl.ac.uk/struct/ARCD/NUM.html, or by          *
C* anonymous ftp at seamus.cc.rl.ac.uk/pub.  Also contact John         *
C* Harding, Harwell Subroutine Library, B 552, AEA Technology,         *
C* Harwell, Didcot, Oxon OX11 0RA, England.                            *
C* telephone (44) 1235 434573, fax (44) 1235 434340,                   *
C* email john.harding@aeat.co.uk, who will provide details of price    *
C* and conditions of use.                                              *
C***********************************************************************

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  print input/output arguments for UMS2FA, UMS2RF, and UMS2SO

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  This routine can be deleted on installation (replaced with a dummy
C  routine that just returns without printing) in order to completely
C  disable the printing of all input/output parameters.  To completely
C  disable all I/O, you can also replace the UMS2P2 routine with a
C  dummy subroutine.  If you make this modification, please do
C  not delete any original code - just comment it out instead.  Add a
C  comment and date to your modifications.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       who:            what routine called UMS2P1:
C                       1: UMS2FA, 2: UMS2RF, 3: UMS2SO
C       where:          called from where:
C                       1: entry of routine, else exit of routine
C       Icntl (3):      if < 3 then print nothing, if 3 then print
C                       terse output, if >= 4 print everything
C       Icntl (2):      I/O unit on which to print.  No printing
C                       occurs if < 0.
C
C       Parameters to print, see UMS2FA, UMS2RF, or UMS2SO for
C       descriptions:
C
C           n, ne, job, trans, lvalue, lindex, Value, Index, Keep,
C           Icntl, Info, Rinfo, B, X, lx, W, lw

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       on Icntl (2) I/O unit only

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  UMS2FA, UMS2RF, UMS2SO
C       functions called:       MIN
        INTRINSIC MIN

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        LOGICAL TRANSA, TRANSC, PRLU, BADLU, SGLTON, PRESRV, SYMBOL
        INTEGER IO, PRL, PRN, K, LUI1, LUI2, LUX1, LUX2, ROW, COL,
     $          FACNE, FACN, NZ, FACJOB, NBLKS, NZOFF, FACTRA, CPERMP,
     $          RPERMP, APP, AXP, AIP, OFFIP, OFFXP, LUBLPP, OFFPP,
     $          BLKPP, P1, P2, P, BLK, K1, K2, KN, LUIIP, LUXXP, NPIV,
     $          NLU, E, LUK, LUPP, LUIP, LUXP, LUDEGR, LUDEGC, LUNSON,
     $          LUSONP, LUCP, LURP, I, J, NZCOL, NZROW, UXP, SON,
     $          PRMAX, LUDIMR, LUDIMC, MAXDR, MAXDC, LUIR1, IP1, IP2,
     $          XP1
        PARAMETER (PRMAX = 10)

C  Printing control:
C  -----------------
C  io:      I/O unit for diagnostic messages
C  prl:     printing level
C  prn:     number of entries printed so far
C  prmax:   maximum number of entries to print if prl = 3
C  prlu:    true if printing LU factors
C
C  Location and status of LU factors:
C  ----------------------------------
C  transc:  TRANSC argument in UMS2SO
C  transa:  TRANSA argument in UMS2FA or UMS2RF when matrix factorized
C  badlu:   true if LU factors uncomputed or corrupted
C  presrv:  true if original matrix was preserved when factorized
C  symbol:  true if only symbolic part of LU factors needed on input
C  lui1:    integer part of LU factors start in Index (lui1...)
C  luir1:   Index (luir1 ... lui2) is needed for a call to UMS2RF
C  lui2:    integer part of LU factors end in Index (..lui2)
C  lux1:    real part of LU factors start in Value (lux1...)
C  lux2:    real part of LU factors end in Value (...lux1)
C  ip1:     pointer into leading part of LU factors in Index
C  ip2:     pointer into trailing part of LU factors in Index
C  xp1:     pointer into leading part of LU factors in Value
C
C  Arrays and scalars allocated in LU factors (in order):
C  ------------------------------------------------------
C  app:     Ap (1..n+1) array located in Index (app...app+n)
C  axp:     Ax (1..nz) array located in Value (axp...axp+nz-1)
C  aip:     Ai (1..nz) array located in Index (aip...aip+nz-1)
C  offip:   Offi (1..nzoff) array loc. in Index (offip...offip+nzoff-1)
C  offxp:   Offx (1..nzoff) array loc. in Value (offxp...offxp+nzoff-1)
C  ...      LU factors of each diagonal block located here
C  lublpp:  LUblkp (1..nblks) array in Index (lublpp..lublpp+nblks-1)
C  blkpp:   Blkp (1..nblks+1) array loc. in Index (blkpp...blkpp+nblks)
C  offpp:   Offp (1..n+1) array located in Index (offpp...offpp+n)
C  cpermp:  Cperm (1..n) array located in Index (cpermp...cpermp+n-1)
C  rpermp:  Rperm (1..n) array located in Index (rpermp...rpermp+n-1)
C  ...      seven scalars in Index (lui2-6...lui2):
C  factra:  0/1 if TRANSA argument was false/true in UMS2FA or UMS2RF
C  nzoff:   number of entries in off-diagonal part
C  nblks:   number of diagonal blocks
C  facjob:  JOB argument in UMS2FA or UMS2RF when matrix factorized 
C  nz:      entries in A
C  facn:    N argument in UMS2FA or UMS2RF when matrix factorized 
C  facne:   NE argument in UMS2FA or UMS2RF when matrix factorized 
C
C  A single diagonal block and its LU factors:
C  -------------------------------------------
C  blk:     current diagonal block
C  k1,k2:   current diagonal is A (k1..k2, k1..k2)
C  kn:      order of current diagonal block (= k2-k1+1)
C  sglton:  true if current diagonal block is 1-by-1 (a singleton)
C  luiip:   LU factors of a diagonal block start in Index (luiip...)
C  luxxp:   LU factors of a diagonal block start in Value (luxxp...)
C  npiv:    number of pivots in a diagonal block (0 <= npiv <= kn)
C  nlu:     number of elements in a diagonal block
C  lupp:    LUp (1..nlu) array located in Index (lupp...lupp+nlu-1)
C
C  An element in the LU factors of a single diagonal block:
C  --------------------------------------------------------
C  e:       element
C  luk:     number of pivots in element e
C  luip:    integer part of element is in Index (luip...)
C  luxp:    real part of element e is in Value (luxp...)
C  ludegr:  row degree (number of columns) of U2 block in element e
C  ludegc:  column degree (number of rows) of L2 block in element e
C  lunson:  number of sons of element e in the assembly DAG
C  lusonp:  list of sons of element e in Index(lusonp...lusonp+lunson-1)
C  lucp:    column pattern (row indices) of L2 block in Index (lucp..)
C  lurp:    row pattern (column indices) of U2 block in Index (lurp..)
C  nzcol:   entries in a column of L, including unit diagonal
C  nzrow:   entries in a row of U, including non-unit diagonal
C  uxp:     a row of the U2 block located in Value (uxp...)
C  son:     a son of the element e
C  ludimr:  row dimension (number of columns) in frontal matrix
C  ludimc:  column dimension (number of rows) in frontal matrix
C  maxdr:   largest ludimr for this block
C  maxdc:   largest ludimc for this block
C
C  Other:
C  ------
C  row:     row index
C  col:     column index
C  k:       kth pivot, and general loop index
C  i, j:    loop indices
C  p:       pointer
C  p1:      column of A starts Ai/Ax (p1...), or row Offi/x (p1...)
C  p2:      column of A ends in Ai/Ax (...p2), or row Offi/x (...p2)

C=======================================================================
C  EXECUTABLE STATEMENTS:
C       if (printing disabled on installation) return
C=======================================================================

C-----------------------------------------------------------------------
C  get printing control parameters
C-----------------------------------------------------------------------

        IO = ICNTL (2)
        PRL = ICNTL (3)
        IF (PRL .LT. 3 .OR. IO .LT. 0) THEN 
C          printing has not been requested
           RETURN
        ENDIF 

C-----------------------------------------------------------------------
C  who is this, and where.  Determine if LU factors are to be printed
C-----------------------------------------------------------------------

        IF (WHO .EQ. 1) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 200) 'UMS2FA input:       '
              PRLU = .FALSE.
           ELSE 
              WRITE (IO, 200) 'UMS2FA output:      '
              PRLU = .TRUE.
           ENDIF 
        ELSE IF (WHO .EQ. 2) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 200) 'UMS2RF input:       '
              PRLU = .TRUE.
           ELSE 
              WRITE (IO, 200) 'UMS2RF output:      '
              PRLU = .TRUE.
           ENDIF 
        ELSE IF (WHO .EQ. 3) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 200) 'UMS2SO input:       '
              PRLU = .TRUE.
           ELSE 
              WRITE (IO, 200) 'UMS2SO output:      '
              PRLU = .FALSE.
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print scalar input arguments: n, ne, job, trans, lvalue, lindex
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1) THEN 
           WRITE (IO, *)  '   Scalar arguments:'
           WRITE (IO, *)  '      N:      ', N, ' : order of matrix A'
           IF (WHO .EQ. 3) THEN 
C             UMS2SO:
C             was A or A^T factorized?
              LUI2 = KEEP (5)
              TRANSA = .FALSE.
              IF (LUI2-6 .GE. 1 .AND. LUI2-6 .LE. LINDEX) THEN 
                 TRANSA = INDEX (LUI2-6) .NE. 0
              ENDIF 
              TRANSC = TRANS
              IF (.NOT. TRANSC) THEN 
                 IF (JOB .EQ. 1) THEN 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve P''Lx=b'
                 ELSE IF (JOB .EQ. 2) THEN 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve UQ''x=b'
                 ELSE IF (.NOT. TRANSA) THEN 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve Ax=b',
     $                           ' (PAQ=LU was factorized)'
                 ELSE 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve A''x=b',
     $                           ' (PA''Q=LU was factorized)'
                 ENDIF 
              ELSE 
                 IF (JOB .EQ. 1) THEN 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve L''Px=b'
                 ELSE IF (JOB .EQ. 2) THEN 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve QU''x=b'
                 ELSE IF (.NOT. TRANSA) THEN 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve A''x=b',
     $                           ' (PAQ=LU was factorized)'
                 ELSE 
                    WRITE (IO, *)'      JOB:    ',JOB,' : solve Ax=b',
     $                           ' (PA''Q=LU was factorized)'
                 ENDIF 
              ENDIF 
              IF (TRANSC) THEN 
                 WRITE (IO, *)
     $                  '      TRANSC:   .true. : see JOB above '
              ELSE 
                 WRITE (IO, *)
     $                  '      TRANSC:   .false. : see JOB above '
              ENDIF 
           ELSE 
C             UMS2FA or UMS2RF:
              WRITE (IO, *)'      NE:     ', NE,' : entries in matrix A'
              IF (JOB .EQ. 1) THEN 
                 WRITE (IO, *)
     $                  '      JOB:    ',JOB,' : matrix A preserved'
              ELSE 
                 WRITE (IO, *)
     $                  '      JOB:    ',JOB,' : matrix A not preserved'
              ENDIF 
              TRANSA = TRANS
              IF (TRANSA) THEN 
                 WRITE (IO, *)
     $                  '      TRANSA:   .true. : factorize A transpose'
              ELSE 
                 WRITE (IO, *)
     $                  '      TRANSA:   .false. : factorize A'
              ENDIF 
           ENDIF 
           WRITE (IO, *)
     $  '      LVALUE: ',LVALUE,' : size of VALUE array'
           WRITE (IO, *)
     $  '      LINDEX: ',LINDEX,' : size of INDEX array'
        ENDIF 

C-----------------------------------------------------------------------
C  print control parameters:  Icntl, Cntl, and Keep (6..8)
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1) THEN 
           WRITE (IO, *)
     $  '   Control parameters, normally initialized by UMS2IN:'
           WRITE (IO, *)
     $  '      ICNTL (1...8): integer control parameters'
           WRITE (IO, *)
     $  '      ICNTL (1): ',ICNTL (1),' : I/O unit for error and',
     $                    ' warning messages'
           WRITE (IO, *)
     $  '      ICNTL (2): ',IO,' : I/O unit for diagnostics'
           WRITE (IO, *)
     $  '      ICNTL (3): ',PRL,' : printing control'
           IF (WHO .EQ. 1) THEN 
              IF (ICNTL (4) .EQ. 1) THEN 
                 WRITE (IO, *)
     $  '      ICNTL (4): ',ICNTL (4),' : use block triangular',
     $                    ' form (BTF)'
              ELSE 
                 WRITE (IO, *)
     $  '      ICNTL (4): ',ICNTL (4),' : do not permute to block',
     $                    ' triangular form (BTF)'
              ENDIF 
              WRITE (IO, *)
     $  '      ICNTL (5): ',ICNTL (5),' : columns examined during',
     $                    ' pivot search'
              IF (ICNTL (6) .NE. 0) THEN 
                 WRITE (IO, *)
     $  '      ICNTL (6): ',ICNTL (6),' : preserve symmetry'
              ELSE 
                 WRITE (IO, *)
     $  '      ICNTL (6): ',ICNTL (6),' : do not preserve symmetry'
              ENDIF 
           ENDIF 
           IF (WHO .NE. 3) THEN 
              WRITE (IO, *)
     $  '      ICNTL (7): ',ICNTL (7),' : block size for dense matrix',
     $                    ' multiply'
           ELSE 
              WRITE (IO, *)
     $  '      ICNTL (8): ',ICNTL (8),' : maximum number',
     $                    ' of iterative refinement steps'
           ENDIF 
           IF (WHO .EQ. 1) THEN 
              WRITE (IO, *)
     $  '      CNTL (1...3): real control parameters'
              WRITE (IO, *)
     $  '      CNTL (1):  ',CNTL (1),' : relative pivot tolerance'
              WRITE (IO, *)
     $  '      CNTL (2):  ',CNTL (2),' : frontal matrix',
     $                    ' growth factor'
              WRITE (IO, *)
     $  '      KEEP (6...8): integer control parameters',
     $                          ' not normally modified by user'
              WRITE (IO, *)
     $  '      KEEP (6):  ',KEEP(6),' : largest positive integer'
              WRITE (IO, *)
     $  '      KEEP (7):  ',KEEP(7),' : dense row/col control, d1'
              WRITE (IO, *)
     $  '      KEEP (8):  ',KEEP(8),' : dense row/col control, d2'
           ELSE IF (WHO .EQ. 3) THEN 
              WRITE (IO, *)
     $  '      CNTL (1...3): real control parameters'
              WRITE (IO, *)
     $  '      CNTL (3):  ',CNTL(3),' : machine epsilon'
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print the informational output
C-----------------------------------------------------------------------

        IF (WHERE .NE. 1) THEN 
           WRITE (IO, *)
     $  '   Output information:'
           WRITE (IO, *)
     $  '      INFO (1...24): integer output information'
           IF (INFO (1) .LT. 0) THEN 
              WRITE (IO, *)
     $  '      INFO (1):  ',INFO (1),' : error occurred!'
           ELSE IF (INFO (1) .GT. 0) THEN 
              WRITE (IO, *)
     $  '      INFO (1):  ',INFO (1),' : warning occurred'
           ELSE 
              WRITE (IO, *)
     $  '      INFO (1):  ',INFO (1),' : no error or warning',
     $                    ' occurred'
           ENDIF 
           IF (WHO .NE. 3) THEN 
              WRITE (IO, *)
     $  '      INFO (2):  ',INFO (2),' : duplicate entries in A'
              WRITE (IO, *)
     $  '      INFO (3):  ',INFO (3),' : invalid entries in A',
     $                          ' (indices not in 1..N)'
              WRITE (IO, *)
     $  '      INFO (4):  ',INFO (4),' : invalid entries in A',
     $                          ' (not in prior pattern)'
              WRITE (IO, *)
     $  '      INFO (5):  ',INFO (5),' : entries in A after adding'
              WRITE (IO, *)
     $  '                       duplicates and removing invalid entries'
              WRITE (IO, *)
     $  '      INFO (6):  ',INFO (6),' : entries in diagonal',
     $                    ' blocks of A'
              WRITE (IO, *)
     $  '      INFO (7):  ',INFO (7),' : entries in off-diagonal',
     $                    ' blocks of A'
              WRITE (IO, *)
     $  '      INFO (8):  ',INFO (8),' : 1-by-1 diagonal blocks',
     $                    ' in A'
              WRITE (IO, *)
     $  '      INFO (9):  ',INFO (9),' : diagonal blocks in A',
     $                    ' (>1 only if BTF used)'
              WRITE (IO, *)
     $  '      INFO (10): ',INFO (10),' : entries below diagonal in L'
              WRITE (IO, *)
     $  '      INFO (11): ',INFO (11),' : entries above diagonal in U'
              WRITE (IO, *)
     $  '      INFO (12): ',INFO (12),' : entries in L + U +',
     $                    ' offdiagonal blocks of A'
              WRITE (IO, *)
     $  '      INFO (13): ',INFO (13),' : frontal matrices'
              WRITE (IO, *)
     $  '      INFO (14): ',INFO (14),' : integer garbage',
     $                    ' collections'
              WRITE (IO, *)
     $  '      INFO (15): ',INFO (15),' : real garbage collections'
              WRITE (IO, *)
     $  '      INFO (16): ',INFO (16),' : diagonal pivots chosen'
              WRITE (IO, *)
     $  '      INFO (17): ',INFO (17),' : numerically valid pivots',
     $                    ' found in A'
              WRITE (IO, *)
     $  '      INFO (18): ',INFO (18),' : memory used in INDEX'
              WRITE (IO, *)
     $  '      INFO (19): ',INFO (19),' : minimum memory needed in',
     $                    ' INDEX'
              WRITE (IO, *)
     $  '      INFO (20): ',INFO (20),' : memory used in VALUE'
              WRITE (IO, *)
     $  '      INFO (21): ',INFO (21),' : minimum memory needed in',
     $                    ' VALUE'
              WRITE (IO, *)
     $  '      INFO (22): ',INFO (22),' : memory needed in',
     $                    ' INDEX for next call to UMS2RF'
              WRITE (IO, *)
     $  '      INFO (23): ',INFO (23),' : memory needed in',
     $                    ' VALUE for next call to UMS2RF'
           ELSE 
              WRITE (IO, *)
     $  '      INFO (24): ',INFO (24),' : steps of iterative',
     $                    ' refinement taken'
           ENDIF 
           IF (WHO .NE. 3) THEN 
              WRITE (IO, *)
     $  '      RINFO (1...8): real output information'
              WRITE (IO, *)
     $  '      RINFO (1): ',RINFO (1),' : total BLAS flop count'
              WRITE (IO, *)
     $  '      RINFO (2): ',RINFO (2),' : assembly flop count'
              WRITE (IO, *)
     $  '      RINFO (3): ',RINFO (3),' : pivot search flop count'
              WRITE (IO, *)
     $  '      RINFO (4): ',RINFO (4),' : Level-1 BLAS flop count'
              WRITE (IO, *)
     $  '      RINFO (5): ',RINFO (5),' : Level-2 BLAS flop count'
              WRITE (IO, *)
     $  '      RINFO (6): ',RINFO (6),' : Level-3 BLAS flop count'
           ELSE IF (LW .EQ. 4*N) THEN 
              WRITE (IO, *)
     $  '      RINFO (1...8): real output information'
              WRITE (IO, *)
     $  '      RINFO (7): ',RINFO (7),' : sparse error estimate',
     $                    ' omega1'
              WRITE (IO, *)
     $  '      RINFO (8): ',RINFO (8),' : sparse error estimate',
     $                    ' omega2'
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print input matrix A, in triplet form, for UMS2FA and UMS2RF
C-----------------------------------------------------------------------

        IF (WHERE .EQ. 1 .AND. WHO .NE. 3) THEN 

           IF (TRANSA) THEN 
              WRITE (IO, *) '   The input matrix A transpose:'
              WRITE (IO, *)
     $  '      VALUE (1 ... ',NE,' ): numerical values'
              WRITE (IO, *)
     $  '      INDEX (1 ... ',NE,' ): column indices'
              WRITE (IO, *)
     $  '      INDEX (',NE+1,' ... ',2*NE,' ): row indices'
              WRITE (IO, *) '   Entries in the matrix A transpose',
     $                      ' (entry number: row, column, value):'
           ELSE 
              WRITE (IO, *) '   The input matrix A:'
              WRITE (IO, *)
     $  '      VALUE (1 ... ',NE,' ): numerical values'
              WRITE (IO, *)
     $  '      INDEX (1 ... ',NE,' ): row indices'
              WRITE (IO, *)
     $  '      INDEX (',NE+1,' ... ',2*NE,' ): column indices'
              WRITE (IO, *) '   Entries in the matrix A',
     $                      ' (entry number: row, column, value):'
           ENDIF 

           PRN = MIN (PRMAX, NE)
           IF (PRL .GE. 4) THEN 
              PRN = NE
           ENDIF 
           DO 10 K = 1, PRN 
              IF (TRANSA) THEN 
                 ROW = INDEX (K+NE)
                 COL = INDEX (K)
              ELSE 
                 ROW = INDEX (K)
                 COL = INDEX (K+NE)
              ENDIF 
              WRITE (IO, *) '      ', K, ': ',ROW,' ',COL,' ', VALUE (K)
10         CONTINUE 
           IF (PRN .LT. NE) THEN 
              WRITE (IO, 220)
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  print the LU factors:  UMS2FA output, UMS2RF input/output,
C                         and UMS2SO input
C-----------------------------------------------------------------------

        IF (PRLU .AND. INFO (1) .LT. 0) THEN 
           WRITE (IO, *) '   LU factors not printed because of error',
     $                   ' flag, INFO (1) = ', INFO (1)
           PRLU = .FALSE.
        ENDIF 

        IF (PRLU) THEN 

C          -------------------------------------------------------------
C          description of what must be preserved between calls
C          -------------------------------------------------------------

           LUX1 = KEEP (1)
           LUX2 = KEEP (2)
           LUI1 = KEEP (3)
           LUIR1 = KEEP (4)
           LUI2 = KEEP (5)

           XP1 = LUX1
           IP1 = LUI1
           IP2 = LUI2

C          -------------------------------------------------------------
C          on input to UMS2RF, only the symbol information is used
C          -------------------------------------------------------------

           SYMBOL = WHO .EQ. 2 .AND. WHERE .EQ. 1

           IF (SYMBOL) THEN 
              WRITE (IO, *)
     $  '   KEEP (4...5) gives the location of LU factors'
              WRITE (IO, *)
     $  '      which must be preserved for calls to UMS2RF: '
           ELSE 
              WRITE (IO, *)
     $  '   KEEP (1...5) gives the location of LU factors'
              WRITE (IO, *)
     $  '      which must be preserved for calls to UMS2SO: '
              WRITE (IO, *)
     $  '         VALUE ( KEEP (1): ', LUX1,' ... KEEP (2): ', LUX2,' )'
              WRITE (IO, *)
     $  '         INDEX ( KEEP (3): ', LUI1,' ... KEEP (5): ', LUI2,' )'
              WRITE (IO, *)
     $  '      and for calls to UMS2RF: '
           ENDIF 
           WRITE (IO, *)
     $  '         INDEX ( KEEP (4): ',LUIR1,' ... KEEP (5): ', LUI2,' )'

           BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1 .OR.
     $        LUI2 .GT. LINDEX
           IF (.NOT. SYMBOL) THEN 
              BADLU = BADLU .OR. LUX1 .LE. 0 .OR.
     $        LUX1 .GT. LUX2 .OR. LUX2 .GT. LVALUE .OR. LUI1 .LE. 0 .OR.
     $        LUIR1 .LT. LUI1 .OR. LUIR1 .GT. LUI2
           ENDIF 

C          -------------------------------------------------------------
C          get the 7 scalars, and location of permutation vectors
C          -------------------------------------------------------------

           IF (BADLU) THEN 
C             pointers are bad, so these values cannot be obtained
              FACNE  = 0
              FACN   = 0
              NZ     = 0
              FACJOB = 0
              NBLKS  = 0
              NZOFF  = 0
              FACTRA = 0
           ELSE 
              FACNE  = INDEX (LUI2)
              FACN   = INDEX (LUI2-1)
              NZ     = INDEX (LUI2-2)
              FACJOB = INDEX (LUI2-3)
              NBLKS  = INDEX (LUI2-4)
              NZOFF  = INDEX (LUI2-5)
              FACTRA = INDEX (LUI2-6)
           ENDIF 

           PRESRV = FACJOB .NE. 0
           TRANSA = FACTRA .NE. 0
           RPERMP = (LUI2-6) - (FACN)
           CPERMP = RPERMP - (FACN)
           IP2 = CPERMP - 1

           IF (SYMBOL) THEN 
              WRITE (IO, *)'   Layout of LU factors in INDEX:'
           ELSE 
              WRITE (IO, *)'   Layout of LU factors in VALUE and INDEX:'
           ENDIF 

C          -------------------------------------------------------------
C          print location of preserved input matrix
C          -------------------------------------------------------------

           IF (PRESRV) THEN 
C             preserved column-form of original matrix
              APP = IP1
              AIP = APP + (FACN+1)
              IP1 = AIP + (NZ)
              AXP = XP1
              XP1 = XP1 + (NZ)
              IF (.NOT. SYMBOL) THEN 
                 WRITE (IO, *)'      preserved copy of original matrix:'
                 WRITE (IO, *)
     $  '         INDEX ( ',APP,' ... ', AIP-1,' ): column pointers'
                 WRITE (IO, *)
     $  '         INDEX ( ',AIP,' ... ',IP1-1,' ): row indices'
                 WRITE (IO, *)
     $  '         VALUE ( ',AXP,' ... ',XP1-1,' ): numerical values'
              ENDIF 
           ELSE 
              IF (.NOT. SYMBOL) THEN 
                 WRITE (IO, *) '      original matrix not preserved.'
              ENDIF 
           ENDIF 

           BADLU = BADLU .OR.
     $          N .NE. FACN .OR. NZ .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $          NBLKS .LE. 0 .OR. NBLKS .GT. N
           IF (.NOT. SYMBOL) THEN 
              BADLU = BADLU .OR. XP1 .GT. LUX2 .OR. NZOFF .LT. 0
           ENDIF 
           IF (BADLU) THEN 
              NBLKS = 0
           ENDIF 

           IF (NBLKS .LE. 1) THEN 

C             ----------------------------------------------------------
C             single block (or block triangular form not used),
C             or LU factors are corrupted
C             ----------------------------------------------------------

              WRITE (IO, *)
     $  '      collection of elements in LU factors:'
              WRITE (IO, *)
     $  '         (an "element" contains one or columns of L and'
              WRITE (IO, *)
     $  '         rows of U with similar nonzero pattern)'
              WRITE (IO, *)
     $  '         INDEX ( ',LUIR1,' ... ', IP2,
     $                                      ' ): integer data'
              IF (.NOT. SYMBOL) THEN 
                 WRITE (IO, *)
     $  '         VALUE ( ',XP1,' ... ', LUX2,' ): numerical values'
              ENDIF 

           ELSE 

C             ----------------------------------------------------------
C             block triangular form with more than one block
C             ----------------------------------------------------------

              OFFIP = IP1
              IP1 = IP1 + (NZOFF)
              OFFXP = XP1
              XP1 = XP1 + (NZOFF)
              OFFPP = CPERMP - (N+1)
              BLKPP = OFFPP - (NBLKS+1)
              LUBLPP = BLKPP - (NBLKS)
              IP2 = LUBLPP - 1
              BADLU = BADLU .OR. LUIR1 .GT. IP2
              IF (.NOT. SYMBOL) THEN 
                 BADLU = BADLU .OR. IP1 .GT. IP2 .OR.
     $           XP1 .GT. LUX2 .OR. LUIR1 .NE. IP1
              ENDIF 
              WRITE (IO, *)
     $  '      matrix permuted to upper block triangular form.'
              IF (NZOFF .NE. 0 .AND. .NOT. SYMBOL) THEN 
                 WRITE (IO, *) '      entries not in diagonal blocks:'
                 WRITE (IO, *)
     $  '         INDEX ( ',OFFIP,' ... ',LUIR1-1,' ): row indices'
                 WRITE (IO, *)
     $  '         VALUE ( ',OFFXP,' ... ',XP1-1,' ): numerical values'
              ENDIF 
              WRITE (IO, *)
     $  '      collection of elements in LU factors of diagonal blocks:'
              WRITE (IO, *)
     $  '         (an "element" contains one or columns of L and'
              WRITE (IO, *)
     $  '         rows of U with similar nonzero pattern)'
              IF (LUIR1 .LE. LUBLPP-1) THEN 
                 WRITE (IO, *)
     $  '         INDEX ( ',LUIR1,' ... ', IP2,
     $                                         ' ): integer data'
              ENDIF 
              IF (XP1 .LE. LUX2 .AND. .NOT. SYMBOL) THEN 
                 WRITE (IO, *)
     $  '         VALUE ( ',XP1,' ... ', LUX2,' ): numerical values'
              ENDIF 
              WRITE (IO, *) '      other block triangular data:'
              WRITE (IO, *)
     $  '         INDEX ( ',LUBLPP,' ... ',BLKPP-1,' ):',
     $                          ' pointers to block factors' 
              WRITE (IO, *)
     $  '         INDEX ( ', BLKPP,' ... ',OFFPP-1,' ):',
     $                          ' index range of blocks'
              IF (.NOT. SYMBOL) THEN 
                 WRITE (IO, *)
     $  '         INDEX ( ', OFFPP,' ... ',   LUI2-7,' ):',
     $          ' row pointers for off-diagonal part'
              ENDIF 
           ENDIF 

C          -------------------------------------------------------------
C          print location of permutation vectors and 7 scalars at tail
C          -------------------------------------------------------------

           WRITE (IO, *)
     $  '      permutation vectors (start at KEEP(4)-2*N-6):'
           WRITE (IO, *)
     $  '         INDEX ( ', CPERMP,' ... ',RPERMP-1,' ):',
     $                          ' column permutations'
           WRITE (IO, *)
     $  '         INDEX ( ', RPERMP,' ... ',LUI2-7,' ):',
     $                          ' row permutations'
           WRITE (IO, *) '      other data in INDEX: '
           WRITE (IO, *)
     $  '         INDEX ( ',LUI2-6,' ): ', FACTRA,' :',
     $                                  ' TRANSA UMS2FA/UMS2RF argument'
           WRITE (IO, *)
     $  '         INDEX ( ',LUI2-5,' ): ', NZOFF,' :',
     $                                  ' entries in off-diagonal part'
           WRITE (IO, *)
     $  '         INDEX ( ',LUI2-4,' ): ', NBLKS,' :',
     $                                  ' number of diagonal blocks'
           WRITE (IO, *)
     $  '         INDEX ( ',LUI2-3,' ): ', FACJOB,' :',
     $                                  ' JOB UMS2FA/UMS2RF argument'
           WRITE (IO, *)
     $  '         INDEX ( ',LUI2-2,' ): ', NZ,' :',
     $                                  ' entries in original matrix'
           WRITE (IO, *)
     $  '         INDEX ( ',LUI2-1,' ): ', FACN,' :',
     $                                  ' N UMS2FA/UMS2RF argument'
           WRITE (IO, *)
     $  '         INDEX ( ',LUI2  ,' ): ', FACNE,' :',
     $                                  ' NE UMS2FA/UMS2RF argument'

           IF (.NOT. SYMBOL) THEN 
              BADLU = BADLU .OR. IP1 .NE. LUIR1
           ENDIF 
           IP1 = LUIR1
           IF (BADLU) THEN 
              WRITE (IO, *) '   LU factors uncomputed or corrupted!'
              PRESRV = .FALSE.
              NBLKS = 0
           ENDIF 

C          -------------------------------------------------------------
C          copy of original matrix in column-oriented form
C          -------------------------------------------------------------

           IF (PRESRV .AND. .NOT. SYMBOL) THEN 
              WRITE (IO, 230)
              WRITE (IO, *) '   Preserved copy of original matrix:'
              DO 20 COL = 1, N 
                 P1 = INDEX (APP-1 + COL)
                 P2 = INDEX (APP-1 + COL+1) - 1
                 WRITE (IO, *) '      col: ', COL, ' nz: ', P2-P1+1
                 IF (PRL .EQ. 3) THEN 
                    P2 = MIN (PRMAX, P2)
                 ENDIF 
                 WRITE (IO, *) (INDEX (AIP-1 + P), P = P1, P2)
                 WRITE (IO, 210) (VALUE (AXP-1 + P), P = P1, P2)
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN 
C                   exit out of loop if done printing:
                    GO TO 30
                 ENDIF 
20            CONTINUE 
C             loop exit label:
30            CONTINUE
              IF (PRL .EQ. 3 .AND. NZ .GT. PRMAX) THEN 
                 WRITE (IO, 220)
              ENDIF 
           ENDIF 

C          -------------------------------------------------------------
C          entries in off-diagonal blocks, in row-oriented form
C          -------------------------------------------------------------

           IF (NBLKS .GT. 1 .AND. .NOT. SYMBOL) THEN 
              WRITE (IO, 230)
              WRITE (IO, *)'   Entries not in diagonal blocks:'
              IF (NZOFF .EQ. 0) THEN 
                 WRITE (IO, *)'      (none)'
              ENDIF 
              DO 40 ROW = 1, N 
                 P1 = INDEX (OFFPP-1 + ROW)
                 P2 = INDEX (OFFPP-1 + ROW+1) - 1
                 IF (P2 .GE. P1) THEN 
                    WRITE (IO, *) '      row: ', ROW, ' nz: ',P2-P1+1
                    IF (PRL .EQ. 3) THEN 
                       P2 = MIN (PRMAX, P2)
                    ENDIF 
                    WRITE (IO, *) (INDEX (OFFIP-1 + P), P = P1, P2)
                    WRITE (IO, 210) (VALUE (OFFXP-1 + P), P = P1, P2)
                 ENDIF 
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN 
C                   exit out of loop if done printing:
                    GO TO 50
                 ENDIF 
40            CONTINUE 
C             loop exit label:
50            CONTINUE
              IF (PRL .EQ. 3 .AND. NZ .GT. PRMAX) THEN 
                 WRITE (IO, 220)
              ENDIF 
           ENDIF 

C          -------------------------------------------------------------
C          LU factors of each diagonal block
C          -------------------------------------------------------------

           WRITE (IO, 230)
           IF (NBLKS .GT. 0) THEN 
              WRITE (IO, *) '   LU factors of each diagonal block:'
           ENDIF 
           PRN = 0
           DO 140 BLK = 1, NBLKS 

C             ----------------------------------------------------------
C             print the factors of a single diagonal block
C             ----------------------------------------------------------

              IF (NBLKS .GT. 1) THEN 
                 K1 = INDEX (BLKPP-1 + BLK)
                 K2 = INDEX (BLKPP-1 + BLK+1) - 1
                 KN = K2-K1+1
                 SGLTON = KN .EQ. 1
                 IF (SGLTON) THEN 
C                   this is a singleton
                    LUXXP = XP1-1 + INDEX (LUBLPP-1 + BLK)
                 ELSE 
                    LUIIP = IP1-1 + INDEX (LUBLPP-1 + BLK)
                 ENDIF 
              ELSE 
                 SGLTON = .FALSE.
                 K1 = 1
                 K2 = N
                 KN = N
                 LUIIP = IP1
              ENDIF 

              WRITE (IO, 240)
              IF (SGLTON) THEN 

C                -------------------------------------------------------
C                this is a singleton
C                -------------------------------------------------------

                 WRITE (IO, *)'   Singleton block: ', BLK,
     $                        ' at index : ', K1
                 IF (.NOT. SYMBOL) THEN 
                    WRITE (IO, *)
     $           '   located in VALUE ( ', LUXXP,' ): ', VALUE (LUXXP)
                 ENDIF 
                 IF (PRL .EQ. 3 .AND. PRN .GT. PRMAX) THEN 
C                   exit out of loop if done printing:
                    GO TO 150
                 ENDIF 
                 PRN = PRN + 1

              ELSE 

C                -------------------------------------------------------
C                this block is larger than 1-by-1
C                -------------------------------------------------------

                 LUXXP = XP1-1 + INDEX (LUIIP)
                 NLU = INDEX (LUIIP+1)
                 NPIV = INDEX (LUIIP+2)
                 MAXDC = INDEX (LUIIP+3)
                 MAXDR = INDEX (LUIIP+4)
                 LUPP = LUIIP+5
                 WRITE (IO, *) '   Block: ',BLK,' first index: ',K1,
     $                         ' last index: ',K2, '   order: ', KN
                 WRITE (IO, *) '   elements: ', NLU, '   pivots: ', NPIV
                 WRITE (IO, *) '   largest contribution block: ',
     $                         MAXDC, ' by ', MAXDR
                 WRITE (IO, *) '   located in INDEX ( ',LUIIP,' ... )'
                 IF (.NOT. SYMBOL) THEN 
                    WRITE (IO, *) '   and in VALUE ( ',LUXXP,' ... )'
                 ENDIF 
                 LUIIP = LUPP + NLU

C                Note: the indices of the LU factors of the block range
C                from 1 to kn, even though the kn-by-kn block resides in
C                A (k1 ... k2, k1 ... k2).
                 K = 0

                 DO 130 E = 1, NLU 

C                   ----------------------------------------------------
C                   print a single element
C                   ----------------------------------------------------

                    LUIP = LUIIP-1 + INDEX (LUPP-1 + E)
                    LUXP = LUXXP-1 + INDEX (LUIP)
                    LUK  = INDEX (LUIP+1)
                    LUDEGR = INDEX (LUIP+2)
                    LUDEGC = INDEX (LUIP+3)
                    LUNSON = INDEX (LUIP+4)
                    LUDIMR = INDEX (LUIP+5)
                    LUDIMC = INDEX (LUIP+6)
                    LUCP = LUIP + 7
                    LURP = LUCP + LUDEGC
                    LUSONP = LURP + LUDEGR
                    WRITE (IO, *) '      e: ',E, ' pivots: ', LUK,
     $                  ' children in dag: ', LUNSON,
     $                  ' frontal matrix: ', LUDIMR, ' by ', LUDIMC

C                   ----------------------------------------------------
C                   print the columns of L
C                   ----------------------------------------------------

                    P = LUXP
                    DO 80 J = 1, LUK 
                       COL = K+J
                       NZCOL = LUK-J+1+LUDEGC
                       WRITE (IO, *) '         col: ',COL,' nz: ',NZCOL
C                      L is unit diagonal
                       PRN = PRN + 1
                       ROW = COL
                       IF (SYMBOL) THEN 
                          WRITE (IO, *) '            ', ROW
                       ELSE 
                          WRITE (IO, *) '            ', ROW, '  1.0'
                       ENDIF 
                       P = P + 1
C                      pivot block
                       DO 60 I = J+1, LUK 
                          IF (PRL.EQ.3 .AND. PRN.GT.PRMAX) THEN 
C                            exit out of loop if done printing:
                             GO TO 150
                          ENDIF 
                          PRN = PRN + 1
                          ROW = K+I
                          IF (SYMBOL) THEN 
                             WRITE (IO, *)'            ', ROW
                          ELSE 
                             WRITE (IO, *)'            ', ROW, VALUE (P)
                          ENDIF 
                          P = P + 1
60                     CONTINUE 
C                      L block
                       DO 70 I = 1, LUDEGC 
                          IF (PRL.EQ.3 .AND. PRN.GT.PRMAX) THEN  
C                            exit out of loop if done printing:
                             GO TO 150
                          ENDIF 
                          PRN = PRN + 1
                          ROW = INDEX (LUCP-1+I)
                          IF (SYMBOL) THEN 
                             WRITE (IO, *)'            ', ROW
                          ELSE 
                             WRITE (IO, *)'            ', ROW, VALUE (P)
                          ENDIF 
                          P = P + 1
70                     CONTINUE 
                       P = P + J
80                  CONTINUE 

C                   ----------------------------------------------------
C                   print the rows of U
C                   ----------------------------------------------------

                    UXP = LUXP + LUK*(LUDEGC+LUK)
                    DO 110 I = 1, LUK 
                       ROW = K+I
                       NZROW = LUK-I+1+LUDEGR
                       WRITE (IO, *) '         row: ',ROW,' nz: ',NZROW
                       P = LUXP + (I-1) + (I-1) * (LUDEGC+LUK)
C                      pivot block
                       DO 90 J = I, LUK 
                          IF (PRL.EQ.3 .AND. PRN.GT.PRMAX) THEN 
C                            exit out of loop if done printing:
                             GO TO 150
                          ENDIF 
                          PRN = PRN + 1
                          COL = K+J
                          IF (SYMBOL) THEN 
                             WRITE (IO, *)'            ', COL
                          ELSE 
                             WRITE (IO, *)'            ', COL, VALUE (P)
                          ENDIF 
                          P = P + (LUDEGC+LUK)
90                     CONTINUE 
                       P = UXP
C                      U block
                       DO 100 J = 1, LUDEGR 
                          IF (PRL.EQ.3 .AND. PRN.GT.PRMAX) THEN 
C                            exit out of loop if done printing:
                             GO TO 150
                          ENDIF 
                          PRN = PRN + 1
                          COL = INDEX (LURP-1+J)
                          IF (SYMBOL) THEN 
                             WRITE (IO, *)'            ', COL
                          ELSE 
                             WRITE (IO, *)'            ', COL, VALUE (P)
                          ENDIF 
                          P = P + LUK
100                    CONTINUE 
                       UXP = UXP + 1
110                 CONTINUE 

C                   ----------------------------------------------------
C                   print the sons of the element in the assembly DAG
C                   ----------------------------------------------------

                    DO 120 I = 1, LUNSON 
                       PRN = PRN + 1
                       SON = INDEX (LUSONP-1+I)
                       IF (SON .LE. KN) THEN 
C                         an LUson
                          WRITE (IO, *) '         LUson: ', SON
                       ELSE IF (SON .LE. 2*KN) THEN 
C                         a Uson
                          WRITE (IO, *) '         Uson:  ', SON-KN
                       ELSE 
C                         an Lson
                          WRITE (IO, *) '         Lson:  ', SON-2*KN
                       ENDIF 
120                 CONTINUE 

C                   ----------------------------------------------------
C                   increment count of pivots within this block
C                   ----------------------------------------------------

                    K = K + LUK
130              CONTINUE 
              ENDIF 
140        CONTINUE 
C          if loop was not exited prematurely, do not print "..." :
           GO TO 160
C          loop exit label:
150        CONTINUE
           WRITE (IO, 220)
160        CONTINUE

C          -------------------------------------------------------------
C          row and column permutations
C          -------------------------------------------------------------

           IF (.NOT. BADLU) THEN 
              WRITE (IO, 230)
              WRITE (IO, *) '      Column permutations'
              IF (PRL .GE. 4 .OR. N .LE. PRMAX) THEN 
                 WRITE (IO, *) (INDEX (CPERMP+I-1), I = 1, N)
              ELSE 
                 WRITE (IO, *) (INDEX (CPERMP+I-1), I = 1, PRMAX)
                 WRITE (IO, 220)
              ENDIF 

              WRITE (IO, 230)
              WRITE (IO, *) '      Row permutations'
              IF (PRL .GE. 4 .OR. N .LE. PRMAX) THEN 
                 WRITE (IO, *) (INDEX (RPERMP+I-1), I = 1, N)
              ELSE 
                 WRITE (IO, *) (INDEX (RPERMP+I-1), I = 1, PRMAX)
                 WRITE (IO, 220)
              ENDIF 
           ENDIF 

        ENDIF 

C-----------------------------------------------------------------------
C  print B (on input) or W and X (on output) for UMS2SO
C-----------------------------------------------------------------------

        IF (WHO .EQ. 3) THEN 
           WRITE (IO, 230)
           PRN = MIN (PRMAX, N)
           IF (PRL .GE. 4) THEN 
C             print all of B, or W and X
              PRN = N
           ENDIF 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, *) '   W (1 ... ',LW,' ), work vector:',
     $                      ' not printed'
              WRITE (IO, *) '   B (1 ... ',N,' ), right-hand side: '
              DO 170 I = 1, PRN 
                 WRITE (IO, *) '      ', I, ': ', B (I)
170           CONTINUE 
              IF (PRN .LT. N) THEN 
                 WRITE (IO, 220)
              ENDIF 
           ELSE 
              IF (INFO (1) .LT. 0) THEN 
                 WRITE (IO, *) '   W (1 ... ',LW,' ), work vector, and'
                 WRITE (IO, *) '   X (1 ... ',N,' ), solution,'
                 WRITE (IO, *) '      not printed because of error',
     $                         ' flag, INFO (1) = ', INFO (1)
              ELSE 
                 IF (LW .EQ. 4*N) THEN 
C                   UMS2SO did iterative refinement
                    WRITE (IO, *) '   W (1 ... ',N,' ), residual: '
                    DO 180 I = 1, PRN 
                       WRITE (IO, *) '      ', I, ': ', W (I)
180                 CONTINUE 
                    IF (PRN .LT. N) THEN 
                       WRITE (IO, 220)
                    ENDIF 
                    WRITE (IO, *) '   W (',N+1,' ... ',LW,' )',
     $                            ', work vector: not printed'
                 ELSE 
C                   no iterative refinement
                    WRITE (IO, *) '   W (1 ... ',LW,' ),',
     $                            ' work vector: not printed'
                 ENDIF 
                 WRITE (IO, *) '   X (1 ... ',N,' ), solution: '
                 DO 190 I = 1, PRN 
                    WRITE (IO, *) '      ', I, ': ', X (I)
190              CONTINUE 
                 IF (PRN .LT. N) THEN 
                    WRITE (IO, 220)
                 ENDIF 
              ENDIF 
           ENDIF 
        ENDIF 

C-----------------------------------------------------------------------
C  who is this, and where:
C-----------------------------------------------------------------------

        IF (WHO .EQ. 1) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 200) 'end of UMS2FA input '
           ELSE 
              WRITE (IO, 200) 'end of UMS2FA output'
           ENDIF 
        ELSE IF (WHO .EQ. 2) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 200) 'end of UMS2RF input '
           ELSE 
              WRITE (IO, 200) 'end of UMS2RF output'
           ENDIF 
        ELSE IF (WHO .EQ. 3) THEN 
           IF (WHERE .EQ. 1) THEN 
              WRITE (IO, 200) 'end of UMS2SO input '
           ELSE 
              WRITE (IO, 200) 'end of UMS2SO output'
           ENDIF 
        ENDIF 

        RETURN

C-----------------------------------------------------------------------
C  format statments
C-----------------------------------------------------------------------

200     FORMAT (60('='), A20)
210     FORMAT (5E16.8)
220     FORMAT ('        ...')
230     FORMAT ('   ', 77 ('-'))
240     FORMAT ('   ', 77 ('.'))
        END 
