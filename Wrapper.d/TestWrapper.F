      program TestWrapper

      common/event/ x, y, z, P_hperp, M_h, cos_vartheta, phi_h,
     + phi_R, e2_rPhi, phi_S,
     + w_int, w_full, P_T

      character*80 filename
      integer IERR, P_T, narg;
      double precision x, y, z, P_hperp, M_h, cos_vartheta, phi_h,
     + phi_R, e2_rPhi, phi_S,
     + w_int, w_full

      save event_in
      save event_out

c     GET NAME OF FILE

      narg = IARGC ()

      if( narg .ne. 1 ) then
         print *, 'Error with usage'
         goto 100
      endif

      call getarg(1, filename)


c     INITIALIZE GMC_TRANS

      IERR = Init_TMDGen( filename )

      if( IERR .ne. 0 ) then
         print *, 'ERROR initializing TMDGen'
         goto 100
      end if


c     MAKE A FAKE EVENT

      x = 0.05
      y = 0.5
      z = 0.5
      P_hperp = 0.45
      M_h = 0.77
      cos_vartheta = 0.15
      phi_h = 2.
      phi_R = 2.
      e2_rPhi = 0.
      phi_S = 2.
      w_int = 123.
      w_full = 123.
      P_T = 1.


c     WEIGHT THE EVENT

      call Eval_TMDGen_XSec


c     REMEMBER TO AVERAGE w_int AND w_full BEFORE TAKING RATIO

      print *, ' F: Weights are ', w_int, ' and ', w_full


 100  continue

c     MUST CLEAR THE MEMORY

      call Clear_TMDGen



      end program
