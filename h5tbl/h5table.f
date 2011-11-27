!! Based on code by Timmes & Swesty.  
!! Adapted, by Neil Miller
!!  The purpose is to read in a tabular Free energy equation of state
!!  and provide the routine get_val 
!!     :  (Rho, T) -> P, E, S

!! h5table module
!!  provides: 
!!    type: h5tbl
!!    subroutine: setup_h5tbl(tbl, Xsize, Ysize, filename, name, ierr)
!!       - allocates memory and reads table from filename
!!    subroutine: free_h5tbl(tbl)
!!       - deallocate memory that was allocated in setup_h5tbl
!!    (Note that setup_h5tbl expects tbl to be allocated and free_h5tbl doesn't
!!     deallocate the tbl pointer itself.  This should be performed higher up)
!!    subroutine: get_val(tbl, Rho, T, Pvect, Evect, Svect, ierr)
!!       - get the value using the biquintic interpolation code

module h5table
  implicit none

  logical, parameter :: dbg = .false.
  integer, parameter :: ERR_FILE = 52

  !!! Indexes for Biquintic interpolation
  integer, parameter :: TBL_VAL = 1
  integer, parameter :: TBL_DX = 2
  integer, parameter :: TBL_DY = 3
  integer, parameter :: TBL_DX2 = 4
  integer, parameter :: TBL_DY2 = 5
  integer, parameter :: TBL_DXDY = 6
  integer, parameter :: TBL_DX2DY = 7
  integer, parameter :: TBL_DXDY2 = 8
  integer, parameter :: TBL_DX2DY2 = 9

  !!! What table are you interested in?  
  integer, parameter :: i_P = 1
  integer, parameter :: i_E = 2
  integer, parameter :: i_S = 3

  integer, parameter :: i_val = 1
  integer, parameter :: i_drho = 2
  integer, parameter :: i_dt = 3
  integer, parameter :: i_dRho2 = 4
  integer, parameter :: i_dt2 = 5
  integer, parameter :: i_dRhodt = 6
  
  double precision :: dblnan = -100000.
  
  type h5tbl
     character (len=256) :: name, datafilename
     
     ! sizes of the arrays
     integer :: imax
     integer :: jmax
     
     double precision :: minlnRho, maxlnRho, minlnT, maxlnT, dlnRho, dlnT
     
     !! let the dRho_v(i) = Rho_v(i+1) - Rho_v(i)
     !! Use the lnRho_v(i) as a hash
     !! Interpolation however occurs using the biquintic method which is linear
     double precision, dimension(:), pointer :: Rho_v, T_v, lnRho_v, lnT_v, &
          dRho_v, dT_v, dRho2_v, dT2_v, dRho_i, dRho2_i, dRho3_i, dT_i, dT2_i, dT3_i
     
     !!  (9, imax, jmax) 
     double precision, dimension(:,:,:), pointer :: Fij
     
     !! In some cases we may want these to already be allocated.  Here is where we will dump the data
     double precision, dimension(:), pointer :: Pvect, Evect, Svect
     
     double precision, dimension(:,:), pointer :: Zi_dv

     

     logical :: IsSetup = .false.               !! Has the table's members been allocated?
  end type h5tbl

  !! You may not need this concept of table handles.
  !! In some cases, you may want to be able to specify, which table only by an integer.
  !! This code allows you to do this by providing a map between int -> (tbl*)
  
  type h5_entry
     type (h5tbl), pointer :: tbl
     integer :: tbl_handle
  end type h5_entry

  integer, parameter :: REGSIZE = 20
  type (h5_entry) :: h5list(REGSIZE)
  
  contains
    
! This routine gets the table from a handle
    subroutine get_tbl(tbl_handle, tbl, ierr)
      integer, intent(in) :: tbl_handle
      type (h5tbl), pointer, intent(inout) :: tbl
      integer, intent(out) :: ierr
    
      ierr = -1
      if((tbl_handle .le. REGSIZE).and.(tbl_handle .gt.0)) then 
         if(h5list(tbl_handle)%tbl_handle .eq. tbl_handle) then
            ierr = 0
            tbl => h5list(tbl_handle)%tbl
         endif
      endif
    end subroutine get_tbl
  
  !! Associate an integer with the table pointer
    subroutine set_handle(tbl, tbl_handle, ierr)
      type (h5tbl), intent(in), pointer :: tbl
      integer, intent(out) :: tbl_handle, ierr
      integer :: index
    
      ierr = 0
      index = 1
      do while (h5list(index)%tbl_handle .eq. index)
         index = index + 1
      end do
    
      if((index .le. REGSIZE).and.(h5list(index)%tbl_handle .ne. index)) then
         tbl_handle = index
         h5list(index)%tbl_handle = index
         h5list(index)%tbl => tbl
      else
         ierr = -3         
      endif
    end subroutine set_handle
  
  !! The handle must not current be locked
    subroutine free_handle(tbl_handle, ierr)
      integer, intent(in) :: tbl_handle
      integer, intent(out) :: ierr
    
      ierr = 0
      if(h5list(tbl_handle)%tbl_handle.eq.tbl_handle) then
         h5list(tbl_handle)%tbl_handle = 0
         nullify(h5list(tbl_handle)%tbl)
      else 
         ierr = -6
      endif
    end subroutine free_handle
    
    !! setup_h5tbl performs two functions
    !!  (1) allocate members of tbl pointer
    !!  (2) load data from filename
    subroutine setup_h5tbl(tbl, Xsize, Ysize, &
         filename, name, ierr)
!      use alert_lib
      
      type (h5tbl), pointer, intent(in) :: tbl  !! This should already be allocated
      character (len=256), intent(in) :: filename, name
      integer, intent(in) :: Xsize, Ysize
      integer, intent(out) :: ierr ! 0 means OK

      integer :: ios, fid, i, j
      double precision :: Rhoval, Tval, lnRhoval, lnTval
      double precision, parameter :: epsilon = 1e-5
      double precision :: FvectIn(9)
      
      integer :: InterpSz !! Size of the interpolation dimension

      
      ierr = 0

      if(.not. associated(tbl)) then
         write(*,*) "You need to allocated the memory for tbl before passing to setup_h5tbl"
         stop
      endif

      tbl%imax = Xsize
      tbl%jmax = Ysize
      tbl%datafilename = filename
      tbl%name = name
      
      write(*,*) "Allocating table internal memory"
      call alloc_1d_array(tbl%Rho_v, Xsize)
      call alloc_1d_array(tbl%lnRho_v, Xsize)
      call alloc_1d_array(tbl%dRho_v, Xsize)
      call alloc_1d_array(tbl%dRho2_v, Xsize)
      call alloc_1d_array(tbl%dRho_i, Xsize)
      call alloc_1d_array(tbl%dRho2_i, Xsize)
      call alloc_1d_array(tbl%dRho3_i, Xsize)
      call alloc_1d_array(tbl%T_v, Ysize)      
      call alloc_1d_array(tbl%lnT_v, Ysize)
      call alloc_1d_array(tbl%dT_v, Ysize)
      call alloc_1d_array(tbl%dT2_v, Ysize)
      call alloc_1d_array(tbl%dT_i, Ysize)
      call alloc_1d_array(tbl%dT2_i, Ysize)
      call alloc_1d_array(tbl%dT3_i, Ysize)
      allocate(tbl%Fij(9, Xsize, Ysize), stat=ierr)
      allocate(tbl%Pvect(6), tbl%Evect(6), tbl%Svect(6), stat=ierr)
      allocate(tbl%Zi_dv(6,3))
      if(ierr /= 0) then
         write(*,*) 'failure to allocate memory for storing table'
      end if

      ios = 0
      ierr = 0
      fid = 20
      
      write(*,*) "Opening file - ", trim(filename)
      open(unit=fid, file=trim(filename), action='read', status ='old', iostat=ios)
      
      !! The format of the filename MUST be
      !!  Rho T Fij, dFdrho, dFdT, d2FdX2, d2FdY2, d2FdXdy, d3FdX2dY, d3FdXdY2, d4FdX2dY2.  
      !!  There should also be Xsize * Ysize lines to the file
      do i=1,Xsize
         do j=1,Ysize
            read(fid,*) lnRhoVal, lnTval, FvectIn
            tbl%Fij(1:9,i,j) = FvectIn

            RhoVal = exp(lnRhoVal)
            Tval = exp(lnTVal)
            if(isnan(RhoVal) .or. &
                 isnan(Tval) .or. &
                 isnan(tbl%Fij(TBL_VAL,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DX,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DY,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DX2,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DY2,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DXDY,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DX2DY,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DXDY2,i,j)) .or. &
                 isnan(tbl%Fij(TBL_DX2DY2,i,j))) then
               write(*,*) "NaN found - this is unnacceptable"  
               write(*,*) "Make sure that the table has NO NaN values in it before you pass it to the h5table"
               stop
            endif
            
            tbl%Rho_v(i) = RhoVal
            tbl%T_v(j) = Tval
            tbl%lnRho_v(i) = lnRhoVal
            tbl%lnT_v(j) = lnTval
         enddo
      enddo
      
      close(unit=fid)
      
      tbl%minlnRho = minval(tbl%lnRho_v)
      tbl%maxlnRho = maxval(tbl%lnRho_v)
      tbl%minlnT = minval(tbl%lnT_v)
      tbl%maxlnT = maxval(tbl%lnT_v)
      tbl%dlnRho = (tbl%maxlnRho - tbl%minlnRho) / (Xsize - 1)
      tbl%dlnT = (tbl%maxlnT - tbl%minlnT) / (Ysize - 1)

      !!! Set up the difference vectors and their inverses
      do i=1,Xsize-1
         tbl%dRho_v(i) = tbl%Rho_v(i+1) - tbl%Rho_v(i)
         tbl%dRho2_v(i) = tbl%dRho_v(i) * tbl%dRho_v(i)
         tbl%dRho_i(i) = 1d0 / tbl%dRho_v(i)
         tbl%dRho2_i(i) = tbl%dRho_i(i) * tbl%dRho_i(i)
         tbl%dRho3_i(i) = tbl%dRho_i(i) * tbl%dRho2_i(i)
      enddo      
      do j=1,Ysize-1
         tbl%dT_v(j) = tbl%T_v(j+1) - tbl%T_v(j)
         tbl%dT2_v(j) = tbl%dT_v(j) * tbl%dT_v(j)
         tbl%dT_i(j) = 1d0 / tbl%dT_v(j)
         tbl%dT2_i(j) = tbl%dT_i(j) * tbl%dT_i(j)
         tbl%dT3_i(j) = tbl%dT_i(j) * tbl%dT2_i(j)
      enddo

      
      !! Check spacing of Rho_v to make sure that it is logrithmically
      !! spaced with seperation dlnRho.  Note a common reason these might be mispaced
      !! would be if the user provides an inccorect Xsize or Ysize
      do i=2,Xsize  
         if(abs(log(tbl%Rho_v(i) / tbl%Rho_v(i-1)) - tbl%dlnRho) .gt. epsilon) then
            write(*,*) "There is something wrong with the grid or the dimensions of the table"
            write(*,*) "i : ", i
            write(*,*) "Rho_v(i) : ", tbl%Rho_v(i)
            write(*,*) "Rho_v(i-1) : ", tbl%Rho_v(i-1)
            write(*,*) "dlnRho : ", tbl%dlnRho
            stop
         endif
      enddo

      !! Check spacing of T_v that it is logrithmically spaced with seperation dlnT
      do j=2,Ysize
         if(abs(log(tbl%T_v(j) / tbl%T_v(j-1)) - tbl%dlnT) .gt. epsilon) then 
            write(*,*) " There is something wrong with the grid or the dimensions of the table"
            write(*,*) "j : ", j
            write(*,*) "T_v(j) : ", tbl%T_v(j)
            write(*,*) "T_v(j-1) : ", tbl%T_v(j-1)
            write(*,*) "dlnT : ", tbl%dlnT
            stop
         endif
      enddo

      !! IF we made it down to here, I guess we probably haven't had any problems
      tbl%IsSetup = .true.

    contains      
      subroutine alloc_1d_array(ptr,sz)
        double precision, dimension(:), pointer :: ptr
        integer, intent(in) :: sz        
        allocate(ptr(sz),stat=ierr)
        if (ierr /= 0) then 
           write(*,*) 'failure in attempt to allocate Aneos_Table storage'
        endif
      end subroutine alloc_1d_array
    end subroutine setup_h5tbl
    
    subroutine free_h5tbl(tbl)
      type (h5tbl), pointer, intent(in) :: tbl
      
      call free_1D(tbl%Rho_v)
      call free_1D(tbl%lnRho_v)
      call free_1D(tbl%dRho_v)
      call free_1D(tbl%dRho2_v)
      call free_1D(tbl%dRho_i)
      call free_1D(tbl%dRho2_i)
      call free_1D(tbl%dRho3_i)
      
      call free_1D(tbl%T_v)
      call free_1D(tbl%lnT_v)
      call free_1D(tbl%dT_v)
      call free_1D(tbl%dT2_v)
      call free_1D(tbl%dT_i)      
      call free_1D(tbl%dT2_i)
      call free_1D(tbl%dT3_i)
      if(associated(tbl%Fij)) deallocate(tbl%Fij)
      if(associated(tbl%Zi_dv)) deallocate(tbl%Zi_dv)
      
      call free_1D(tbl%Pvect)
      call free_1D(tbl%Evect)
      call free_1D(tbl%Svect)
      !! This doesn't deallocate tbl - that should be done higher up
      
    contains
      subroutine free_1D(vector)
        double precision, pointer, dimension(:) :: vector
        if(associated(vector)) deallocate(vector)
      end subroutine free_1D
    end subroutine free_h5tbl
    
    !! get_val : takes the (Rho, T) cordinates
    !!   (1) lookup the (Rho_index, T_index) using the logrithmic hash
    !!   (2) Use the Timmes biquintic iterpolation code - inline here

    !! TODO: Check that his code actually works
    subroutine h5_get_val(tbl, Rho, T, Pvect, Evect, Svect, ierr)
      type (h5tbl), pointer, intent(in) :: tbl
      double precision, intent(in) :: Rho, T
      double precision, intent(out), pointer, dimension(:) :: Pvect, Evect, Svect
      integer, intent(out) :: ierr
      
      integer :: i, dI, tI
      double precision :: lnRho, lnT, RhoSQ

      !! Define Variables for helm Holtz interpolation routine
      double precision :: fi(36) 
      double precision :: z
      double precision :: xt, xd, mxt, mxd!, &
!           w0t, w1t, w2t, w0mt, w1mt, w2mt, &
!           w0d, w1d, w2d, w0md, w1md, w2md
           
      double precision :: si0t, si1t, si2t, si0mt, si1mt, si2mt, &
           si0d, si1d, si2d, si0md, si1md, si2md, &
           dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, &
           dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md, &
           ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
           ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md, &
           dddsi0t, dddsi1t, dddsi2t, dddsi0mt, dddsi1mt, dddsi2mt, &
           dddsi0d, dddsi1d, dddsi2d, dddsi0md, dddsi1md, dddsi2md

      double precision :: free, df_d, df_t, df_dd, df_tt, df_dt, &
           df_ddd, df_ddt, df_dtt, df_ttt

      
      !! Define the types for the statement functions below
      double precision psi0, dpsi0, ddpsi0, dddpsi0,  &
           psi1, dpsi1, ddpsi1, dddpsi1, &
           psi2, dpsi2, ddpsi2, dddpsi2!, &
!           h5 

      !! The following statement functions are used in helm_interp.dek
      
      !..quintic hermite polynomial statement functions
      !..psi0 and its derivatives
      psi0(z)    = z**3 * ( z * (-6.0d0*z + 15.0d0) - 10.0d0) + 1.0d0
      dpsi0(z)   = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z)  = z* ( z*( -120.0d0*z + 180.0d0) - 60.0d0)
      dddpsi0(z) = z*( -360.0d0*z + 360.0d0) - 60.0d0


      !..psi1 and its derivatives
      psi1(z)    = z* (z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)   = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z)  = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
      dddpsi1(z) = z * (-180.0d0*z + 192.0d0) - 36.0d0


      !..psi2  and its derivatives
      psi2(z)    = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)   = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z)  = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) +2.0d0)
      dddpsi2(z) = 0.5d0*(z * (-60.0d0*z + 72.0d0) - 18.0d0)


      !..biquintic hermite polynomial statement function
!      h5(w0t, w1t, w2t, w0mt, w1mt, w2mt, w0d, w1d, w2d, w0md, w1md, w2md)= fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t  + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t  + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt + fi(13) *w1d*w0t   + fi(14) *w1md*w0t + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt + fi(17) *w2d*w0t   + fi(18) *w2md*w0t + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt + fi(21) *w1d*w1t   + fi(22) *w1md*w1t + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt + fi(25) *w2d*w1t   + fi(26) *w2md*w1t + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt + fi(29) *w1d*w2t   + fi(30) *w1md*w2t + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt + fi(33) *w2d*w2t   + fi(34) *w2md*w2t + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt 

      !! First try to make sure the function is being used appropriatedly
      if(.not.associated(tbl)) then
         write(*,*) "The table pointer is not associated!!!"
         stop
      endif

      if(.not.tbl%IsSetup) then
         write(*,*), "You need to call the procedure setup_h5tbl"
         write(*,*), " - this will allocate the member function memory"
         stop
      endif
      
      lnRho = log(Rho)
      lnT = log(T)
      RhoSQ = Rho*Rho

      !! Check that inputs are in the bounds of the table
      if((lnRho.gt.tbl%maxlnRho).or.(lnRho.lt.tbl%minlnRho).or.&
           (lnT.gt.tbl%maxlnT).or.(lnT.lt.tbl%minlnT)) then
         write(*,*) "Data table out of bounds ", tbl%name
         write(*,*) "lnRho range: ", tbl%minlnRho, tbl%maxlnRho
         write(*,*) "lnT range: ", tbl%minlnT, tbl%maxlnT
         write(*,*) "lnRho value: ", lnRho
         write(*,*) "lnT value: ", lnT
         stop
      endif

      dI = get_Rho_index(lnRho)
      tI = get_T_index(lnT)
      if(dbg) then
         write(*,*) "Rho range: ", tbl%Rho_v(dI), tbl%Rho_v(dI+1)
         write(*,*) "T range: ", tbl%T_v(tI), tbl%T_v(tI+1)
         write(*,*) "dRho_v: ", tbl%dRho_v(dI)
         write(*,*) "dRho2_v: ", tbl%dRho2_v(dI)
         write(*,*) "dRho_i: ", tbl%dRho_i(dI)
         write(*,*) "dRho2_i: ", tbl%dRho2_i(dI)
         write(*,*) "dRho3_i: ", tbl%dRho3_i(dI)
      endif
         
      include 'helm_interp.dek'
      !! Use Timmes Code HERE
    contains
      
      !! Hash functions.  The grid must actually be logrithmically spaced
      !!  we also check that this is in fact the case when the grid is
      !!  begin loaded.  
      integer function get_Rho_index(lnRho)
        double precision :: lnRho
        get_Rho_index = floor((lnRho - tbl%minlnRho) / tbl%dlnRho) + 1 !! Start at index 1
      end function get_Rho_index
      
      integer function get_T_index(lnT)
        double precision :: lnT      
        get_T_index = floor((lnT - tbl%minlnT) / tbl%dlnT) + 1 !! start at index 1
      end function get_T_index

      double precision function h5(w0t, w1t, w2t, w0mt, w1mt, w2mt, w0d, w1d, w2d, w0md, w1md, w2md)
        double precision, intent(in) :: w0t, w1t, w2t, w0mt, w1mt, w2mt, &
             w0d, w1d, w2d, w0md, w1md, w2md

        h5 = fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t  &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t  &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t  &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt 
      end function h5
    end subroutine h5_get_val
    
  end module h5table
