program write_idealgas_tbl
  use ideal_water
  implicit none
  
  integer, parameter :: imax = 256
  integer, parameter :: jmax = 256
  double precision :: minlogRho, maxlogRho, minlogT, maxlogT, &
       dlogRho, dlogT, Rho, T, logRho, logT, lnRho, lnT, &
       f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
       d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, &
       d4f_drho2dt2, nden, N_A, mu, alpha, k, h, m
  integer :: i,j, fid
  character (len=256) :: filename
  
  minlogRho = -2.
  maxlogRho = 2.
  minlogT = 3.
  maxlogT = 6.

  dlogRho = (maxlogRho - minlogRho) / (imax - 1)
  dlogT = (maxlogT - minlogT) / (jmax - 1)
  
  fid = 20
  filename = "idealwater.data"
  open(fid, file=filename)
  
  do i=1,imax
     do j=1,jmax
        logRho = minlogRho + (i-1) * dlogRho
        logT = minlogT + (j-1) * dlogT
        
        Rho = 10.**logRho
        T = 10.**logT
        lnRho = log(Rho)
        lnT = log(T)
        
        !! Get the equation of state here

        call ig_get_f(lnRho, lnT, f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
             d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, d4f_drho2dt2)

        write(fid, *) lnRho, lnT, f_val, df_drho, df_dt, &
             d2f_drho2, d2f_dt2, d2f_drhodt, &
             d3f_drho2dt, d3f_drhodt2, d4f_drho2dt2
     enddo
  enddo

  close(fid)

end program write_idealgas_tbl
