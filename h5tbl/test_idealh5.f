program test_idealh5
  
  !  PARTS:
  !  1. Compare the results of the h5table over a variety of points to the ideal equation of state for which it was generated
  !  2. Plot, a single segment along the grid in the temperature direction
  
  use ideal_water, only : ig_get_val
  use h5table
  
  implicit none  

  integer, parameter :: maxmat=21
  integer, parameter :: imax=128
  integer, parameter :: jmax=128
  integer :: imat,nntype(maxmat),klst,kinp,kpai
  double precision :: tti,ttti,rrhoi,ppi,uui,ssi,ccvi
  double precision :: ddpdti,ddpdri,ffkrosi,ccsi,ffve,ffva
  
  integer :: i, j, outfid, outfid2, ios
  double precision :: logtmin, logtmax, logrhomin, logrhomax, &
       logtstep, logrhostep
  double precision :: logrhovect(imax), rhovect(imax), lnrhovect(imax)
  double precision :: logtvect(jmax), tvect(jmax), lntvect(jmax)
  double precision :: logPmat(imax,jmax),logEmat(imax,jmax),logSmat(imax,jmax), &
       logPdiff(imax,jmax),logEdiff(imax,jmax),logSdiff(imax,jmax), &
       h5logP(imax,jmax), h5logE(imax,jmax), h5logS(imax,jmax), &
       free_energy(imax,jmax)
  double precision :: f_val, df_drho, df_dt, &
       d2f_drho2, d2f_dt2, d2f_drhodt, &
       d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, d4f_drho2dt2, rho, T, lnRho, lnT, RhoSQ
  
  
  
  type (h5tbl), pointer :: mytbl
  character (len=256) :: h5filename, h5name
  double precision, pointer, dimension(:) :: Pvect, Evect, Svect, IG_Pvect, IG_Evect, IG_Svect
  integer :: h5Xsz, h5Ysz
  integer :: ierr
  
  double precision :: free_cluster(3,3) !! take 9 points at every point to calculate numerical derivatives      
  
  
  double precision :: yD0, xD0, sigDX, sigDY, epsilon, mult, rad

  double precision :: delta = 1d-4
  double precision :: dX, dY, twodX, twodY, tmpvar
  double precision :: Xvect(3), Yvect(3)
  
  integer :: idx, idy
  double precision :: dfdx_num, dfdy_num, d2fdx2_num, d2fdy2_num, d2fdxdy_num, d3fdxdy2_num, d3fdx2dy_num, d4fdx2dy2_num
  
  !! Some variables for the 1D vector
  integer, parameter :: LEN = 1000
  
  character (len=256) :: filename
  integer, parameter :: VERBOSE = 0
  
  !! Initiate the h5table
  
  h5name = "water"
  ierr = 0
  h5filename = "idealwater.data"
  h5Xsz = 256
  h5Ysz = 256
  
  allocate(mytbl)
  allocate(Pvect(6), Evect(6), Svect(6), IG_Pvect(6), IG_Evect(6), IG_Svect(6))
  call setup_h5tbl(mytbl, h5Xsz, h5Ysz, h5filename, h5name, ierr)
  
  !! h5 Table initiated!!
  
  logtmin = 4.
  logtmax = 4.1
  logrhomin = -0.1
  logrhomax = 0.1
  
  logrhostep = (logrhomax - logrhomin) / (imax - 1)
  logtstep = (logtmax - logtmin) / (jmax - 1)
  
  do i=1,imax
     logrhovect(i) = logrhomin + logrhostep * (i-1)
     rhovect(i) = 10**logrhovect(i)
     lnrhovect(i) = log(rhovect(i))
  enddo
  
  do j=1,jmax
     logtvect(j) = logtmin + logtstep * (j-1)
     tvect(j) = 10**logtvect(j)
     lntvect(j) = log(tvect(j))
  enddo
  
  do i=1,imax
        
     Rho = rhovect(i)
     do j=1,jmax            
        T = tvect(j)
        
        !! Use the prewritten and packaged water ideal gas law equation of state
        lnRho = log(Rho)
        lnT = log(T)

        call ig_get_val(Rho, T, IG_Pvect, IG_Evect, IG_Svect, ierr)

!        write(*,*) IG_Pvect(i_val), IG_Svect(i_val), IG_Evect(i_val)
        logPMat(i,j) = log10(IG_Pvect(i_val))
        logSMat(i,j) = log10(IG_Svect(i_val))
        logEmat(i,j) = log10(IG_Evect(i_val))
 
        !! biquintic interpolation routine
        call h5_get_val(mytbl, Rho, T, Pvect, Evect, Svect, ierr)
        
        h5logP(i,j) = log10(Pvect(i_val))
        h5logE(i,j) = log10(Evect(i_val))
        h5logS(i,j) = log10(Svect(i_val))
        
        
        logPdiff(i,j) = logPmat(i,j) - h5logP(i,j)
        logEdiff(i,j) = logEmat(i,j) - h5logE(i,j)
        logSdiff(i,j) = logSmat(i,j) - h5logS(i,j)
        
           
     enddo
  enddo
     
     !! Write out the differences to a file
     
     
     
  filename = "h5tdata/h5t_logRhovect.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')
  write(19,'(e20.6)') logRhovect
  close(unit=19)
  
  filename = "h5tdata/h5t_logTvect.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')
  write(19,'(e20.6)') logtvect
  close(unit=19)
  
  filename = "h5tdata/h5t_logPdiff.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') logPdiff
  close(unit=19)
  
  filename = "h5tdata/h5t_logEdiff.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') logEdiff
  close(unit=19)
  
  filename = "h5tdata/h5t_logSdiff.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') logSdiff
  close(unit=19)
  
  filename = "h5tdata/an_logP.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') logPmat
  close(unit=19)
  
  filename = "h5tdata/an_logE.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') logEmat
  close(unit=19)
  
  filename = "h5tdata/an_logS.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') logSmat
  close(unit=19)
  
  filename = "h5tdata/h5t_logP.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') h5logP
  close(unit=19)
  
  filename = "h5tdata/h5t_logE.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') h5logE
  close(unit=19)
  
  filename = "h5tdata/h5t_logS.data"
  open(unit=19,file=trim(filename),status='replace', &
       iostat=ios,action='write',form='formatted')       
  write(19,'(e20.6)') h5logS
  close(unit=19)


  
  !! Clean up the h5 table
  
  call free_h5tbl(mytbl)
  deallocate(mytbl)
  deallocate(Pvect)
  deallocate(Evect)
  deallocate(Svect)
  deallocate(IG_Pvect)
  deallocate(IG_Evect)
  deallocate(IG_Svect)
  
contains
  subroutine donothing(temp)
    double precision, intent(inout) :: temp
  end subroutine donothing
end program test_idealh5

