module ideal_water
  implicit none

  integer, parameter :: i_P = 1
  integer, parameter :: i_E = 2
  integer, parameter :: i_S = 3

  integer, parameter :: i_val = 1
  integer, parameter :: i_drho = 2
  integer, parameter :: i_dt = 3
  integer, parameter :: i_dRho2 = 4
  integer, parameter :: i_dt2 = 5
  integer, parameter :: i_dRhodt = 6

contains
  
  !! This is the ideal gas law for water.  
  !! This does not include any effects such as dissociation or ionization.
  !! The purpose here though is to demonstrate how to use the h5table - not the physics
  
  subroutine ig_get_f(lnRho, lnT, f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
       d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, d4f_drho2dt2)
    double precision, intent(in) :: lnRho, lnT
    double precision, intent(out) :: f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
         d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, d4f_drho2dt2
    double precision :: N_A, mu, alpha, k, h, m, Rho, T
    
    Rho = exp(lnRho)
    T = exp(lnT)
    
    h = 6.6e-27 
    k = 1.38e-16 !! erg K^-1
    mu = 18. !! Water
    N_A = 6.e23 !! Number of particles per mole
    m = 18. / N_A !! g / particle
    alpha = (h / (2 * 3.14 * m * k)**0.5)**3.
    
    !! Get the equation of state here
    
    f_val = -(N_A / mu) * k * T * ( - log(rho) - log(N_A / mu) - log(alpha) + 1.5 * log(T) + 1.)
    df_drho =  (N_A / mu) * k * T / rho
    df_dt = -(N_A / mu) * k * ( - log(rho) - log(N_A / mu) - log(alpha) + 1.5 * log(T) + 2.5)
    d2f_drho2 = -(N_A / mu) * k * T / rho**2
    d2f_dt2 =  -(N_A / mu) * 1.5 * k  / T
    d2f_drhodt =  (N_A / mu) * k / rho
    d3f_drho3 = (N_A / mu) * 2. * k * T / rho**3
    d3f_drho2dt = -(N_A / mu) * k / rho**2
    d3f_drhodt2 = 0.
    d3f_dt3 = (N_A / mu) * 1.5 * k / T**2.
    d4f_drho2dt2 = 0.
    
  end subroutine ig_get_f
  
  subroutine convDF_thermo(rho, t, f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
       d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, &
       Pvect, Evect, Svect)
    
    double precision, intent(in) :: rho, t, f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
         d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3
    double precision, dimension(:), pointer :: Pvect, Evect, Svect  !! These should already be allocated 
    double precision :: RhoSQ
    
    RhoSQ = rho*rho
    Pvect(i_val) = RhoSQ * df_drho
    Pvect(i_drho) = RhoSQ * d2f_drho2 + 2 * rho * df_drho
    Pvect(i_dt) = RhoSQ * d2f_drhodt
    Pvect(i_drho2) = RhoSQ * d3f_drho3 + 4. * rho * d2f_drho2 + 2. * df_drho 
    Pvect(i_dt2) = RhoSQ * d3f_drhodt2
    Pvect(i_drhodt) = RhoSQ * d3f_drho2dt + 2. * rho * d2f_drhodt
    
    Svect(i_val) = -df_dt
    Svect(i_drho) = -d2f_dt2
    Svect(i_dt) = -d2f_drhodt
    Svect(i_drho2) = -d3f_drho2dt
    Svect(i_dt2) = - d3f_dt3
    Svect(i_drhodt) = - d3f_drhodt2
    
    Evect(i_val) = f_val - T * df_dt
    Evect(i_drho) = df_drho - T * d2f_drhodt
    Evect(i_dt) = -T * d2f_dt2
    Evect(i_drho2) = d2f_drho2 - T*d3f_drho2dt
    Evect(i_dt2) = -d2f_dt2 - T * d3f_dt3
    Evect(i_drhodt) = - T * d3f_drhodt2
  end subroutine convDF_thermo
  
  subroutine ig_get_val(Rho, T, Pvect, Evect, Svect, ierr)
    double precision, intent(in) :: Rho, T
    double precision, pointer, dimension(:) :: Pvect, Evect, Svect !! These should already be allocated
    integer, intent(out) :: ierr

    double precision :: lnRho, lnT, f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
         d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, d4f_drho2dt2
    
    lnRho = log(Rho)
    lnT   = log(T)
    
    call ig_get_f(lnRho, lnT, f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
         d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, d4f_drho2dt2)
    
    call convDF_thermo(Rho, T, f_val, df_drho, df_dt, d2f_drho2, d2f_dt2, d2f_drhodt, &
       d3f_drho3, d3f_drho2dt, d3f_drhodt2, d3f_dt3, &
       Pvect, Evect, Svect)    
  end subroutine ig_get_val
  
end module ideal_water
