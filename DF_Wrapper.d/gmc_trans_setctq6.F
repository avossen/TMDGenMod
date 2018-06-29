C
C      Based on CTEQ6 routine SetCtq6
C      Modified for use in GMC_Trans
C
C      Changes include using an error flag rather than a stop statement
C      and creation of a common block to store the Isetold variable.
C      Also, sent error message to stderr
C
C      by Stephen Gliske
C      (sgliske@umich.edu)
C

      Subroutine GMC_Trans_SetCtq6 (Iset, Ierr)

      Implicit Double Precision (A-H,O-Z)
#include <pepadm.inc>
      Parameter (Isetmax0=3)
      Character Flnm(Isetmax0)*20, nn*3, Tablefile*40
      Data (Flnm(I), I=1,Isetmax0)
     > / './pdf/cteq6l.tbl', './pdf/cteq6d.tbl', './pdf/cteq6m.tbl' /

c  *  Allow common block to check for multiple initializations
      common /cteq_common_block/ Isetold

      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,101,140/
      save


      Ierr=0

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         Ierr = 1
         IU= NextUn6()
         If (Iset.ge.Isetmin0 .and. Iset.le.Isetmax0) Then
C            Tablefile=Flnm(Iset)//'.tbl'
            Tablefile=Flnm(Iset)
         Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
            write(nn,'(I3)') Iset
            Tablefile=Flnm(1)//nn//'.tbl'
         Else
            Print *, 'Invalid Iset number in SetCtq6 :', Iset
            Stop
         Endif
         Open(IU, File=Tablefile, Status='OLD', Err=100)
 21      Call ReadTbl6 (IU)
         Close (IU)
         Isetold=Iset
         Ierr=0
      Endif
      Return

C 100  Print 0, ' Data file ', Tablefile, ' cannot be opened '
C     >//'in SetCtq6!!'
C      RETURN

 100  WRITE(0,*) '\tERROR opening datafile ', Tablefile

C                             ********************
      End
