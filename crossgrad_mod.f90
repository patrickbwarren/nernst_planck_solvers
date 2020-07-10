module crossgrad

  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: bc(8)
  integer :: nx, ny, ncmp, monitor = 0
  integer, allocatable :: z(:)
  real (kind=dp) :: dt, rmsdphi
  real (kind=dp) :: eps = 1.0e-14_dp
  real (kind=dp), allocatable :: d(:)
  real (kind=dp), allocatable :: rho(:, :, :), drho(:, :, :)
  real (kind=dp), allocatable :: jx(:, :, :), jy(:, :, :)
  real (kind=dp), allocatable :: sigma(:, :), g(:, :)
  real (kind=dp), allocatable :: phi(:, :)
  real (kind=dp), allocatable :: ix(:, :), iy(:, :)
  real (kind=dp), allocatable :: div(:, :), curl(:, :)
  character (len=100) :: run_name
  logical :: converged

contains

  subroutine read_rc(rc_file)
    implicit none
    character (len=*), intent(in) :: rc_file
    open (5, file=rc_file, action='read')
    read (5, '(A)') run_name
    read (5, *) ! skip ions
    read (5, *) nx, ny, ncmp
    read (5, *) bc ! specify boundary conditions
    call make_grid
    read (5, *) z
    read (5, *) d
    read (5, *) dt
    close (5)
    if (monitor.gt.1) then
       print *, 'read run controls from ', rc_file
       print *, 'run name ', run_name
       print *, 'nx, ny, ncmp = ', nx, ny, ncmp
       print *, 'z = ', z
       print '(A, 5F10.2)', 'd = ', d
       print '(A, F10.2)', 'dt = ', dt
    end if
  end subroutine read_rc

  subroutine make_grid
    implicit none
    if (monitor.gt.2) print *, 'crossgrad_mod: make_grid: entered'
    if (allocated(z)) deallocate(z)
    if (allocated(d)) deallocate(d)
    if (allocated(rho)) deallocate(rho)
    if (allocated(drho)) deallocate(drho)
    if (allocated(jx)) deallocate(jx)
    if (allocated(jy)) deallocate(jy)
    if (allocated(sigma)) deallocate(sigma)
    if (allocated(g)) deallocate(g)
    if (allocated(phi)) deallocate(phi)
    if (allocated(ix)) deallocate(ix)
    if (allocated(iy)) deallocate(iy)
    if (allocated(div)) deallocate(div)
    if (allocated(curl)) deallocate(curl)
    allocate(z(ncmp)) ; allocate(d(ncmp))
    allocate(rho(nx, ny, ncmp)) ; allocate(drho(nx, ny, ncmp))
    allocate(jx(nx, ny, ncmp)) ;  allocate(jy(nx, ny, ncmp))
    allocate(sigma(nx, ny)) ;  allocate(g(nx, ny)) ; allocate(phi(nx, ny))
    allocate(ix(nx, ny)) ;  allocate(iy(nx, ny))
    allocate(div(nx, ny)) ;  allocate(curl(nx, ny))
    phi = 0.0
    if (monitor.gt.1) print *, 'crossgrad_mod: make_grid: finished, nx, ny, ncmp', nx, ny, ncmp
  end subroutine make_grid

  subroutine set_sigma_g
    implicit none
    integer :: k
    if (monitor.gt.2) print *, 'crossgrad_mod: set_sigma_g: entered'
    g = 0.0
    sigma = 0.0
    do k = 1, ncmp
       g = g + z(k) * d(k) * rho(:, :, k)
       sigma = sigma + z(k)**2 * d(k) * rho(:, :, k)
    end do
    if (monitor.gt.1) print *, 'crossgrad_mod: set_sigma_g: finished'
  end subroutine set_sigma_g

  subroutine jacobi(niter)
    implicit none
    integer, intent(in) :: niter
    integer :: i, j, iter
    real (kind=dp) :: a, b, sigma_mid
    real (kind=dp) :: sumsq
    if (monitor.gt.2) print *, 'crossgrad_mod: jacobi: entered, niter =', niter
    sumsq = 0.0
    do iter = 1, niter
       do j = 1, ny
          do i = 1, nx
             a = 0.0
             b = 0.0
             if (i.gt.1) then
                sigma_mid = 0.5*(sigma(i, j) + sigma(i-1, j))
                a = a + (g(i-1, j) - g(i, j)) + sigma_mid * phi(i-1, j)
                b = b + sigma_mid
             end if
             if (i.lt.nx) then
                sigma_mid = 0.5*(sigma(i, j) + sigma(i+1, j))
                a = a + (g(i+1, j) - g(i, j)) + sigma_mid * phi(i+1, j)
                b = b + sigma_mid
             end if
             if (j.gt.1) then
                sigma_mid = 0.5*(sigma(i, j) + sigma(i, j-1))
                a = a + (g(i, j-1) - g(i, j)) + sigma_mid * phi(i, j-1)
                b = b + sigma_mid
             end if
             if (j.lt.ny) then
                sigma_mid = 0.5*(sigma(i, j) + sigma(i, j+1))
                a = a + (g(i, j+1) - g(i, j)) + sigma_mid * phi(i, j+1)
                b = b + sigma_mid
             end if
             if (iter.eq.niter) sumsq = sumsq + (phi(i, j) - a/b)**2
             phi(i, j) = a / b
          end do
       end do
    end do
    rmsdphi = sqrt(sumsq)
    if (monitor.gt.1) print *, 'crossgrad_mod: jacobi: finished, normdphi =', rmsdphi
  end subroutine jacobi

  subroutine converge(njacobi, niter, nmon)
    implicit none
    integer iter
    integer, intent(in) :: njacobi, niter, nmon
    real (kind=dp) :: normi
    call set_sigma_g
    do iter = 0, niter
       if (nmon.gt.0 .and. mod(iter, nmon).eq.0) then
          call testdiv(normi)
          print *, 'monitor: testdiv: ', iter*njacobi, '|ΣI|, |Δφ|/|φ| =', normi, rmsdphi
          call current_direct ; call divcurl
          print *, 'monitor: direct:  ', iter*njacobi, '|ΣI|, |Δφ|/|φ| =', sqrt(sum(div**2)), rmsdphi
          call fluxes ; call current ; call divcurl
          print *, 'monitor: fluxes:  ', iter*njacobi, '|ΣI|, |Δφ|/|φ| =', sqrt(sum(div**2)), rmsdphi
       end if
       if (iter.lt.niter) call jacobi(njacobi)
       if (rmsdphi .lt. eps) exit
    end do
    converged = rmsdphi .lt. eps
  end subroutine converge

  subroutine testdiv(normi)
    implicit none
    integer :: i, j
    real (kind=dp) :: local_div(nx, ny)
    real (kind=dp), intent(out) :: normi
    if (monitor.gt.2) print *, 'crossgrad_mod: testdiv: entered'
    local_div = 0.0
    do j = 1, ny
       do i = 1, nx
          if (i.gt.1) then
             local_div(i, j) = local_div(i, j) + (g(i-1, j) - g(i, j)) &
                  & + 0.5*(sigma(i, j) + sigma(i-1, j)) * (phi(i-1, j) - phi(i, j))
          end if
          if (i.lt.nx) then
             local_div(i, j) = local_div(i, j) + (g(i+1, j) - g(i, j)) &
                  & + 0.5*(sigma(i, j) + sigma(i+1, j)) * (phi(i+1, j) - phi(i, j))
          end if
          if (j.gt.1) then
             local_div(i, j) = local_div(i, j) + (g(i, j-1) - g(i, j)) &
                  & + 0.5*(sigma(i, j) + sigma(i, j-1)) * (phi(i, j-1) - phi(i, j))
          end if
          if (j.lt.ny) then
             local_div(i, j) = local_div(i, j) + (g(i, j+1) - g(i, j)) &
                  & + 0.5*(sigma(i, j) + sigma(i, j+1)) * (phi(i, j+1) - phi(i, j))
          end if
       end do
    end do
    normi = sqrt(sum(local_div**2))
    if (monitor.gt.1) print *, 'crossgrad_mod: testdiv: finished, |div I|    =', normi
  end subroutine testdiv

  subroutine testrhoz(sumrhoz)
    implicit none
    integer :: k
    real (kind=dp), intent(out) :: sumrhoz
    if (monitor.gt.2) print *, 'crossgrad_mod: testrhoz: entered'
    sumrhoz = 0.0
    do k = 1, ncmp
       sumrhoz = sumrhoz + z(k) * sum(rho(:, :, k))
    end do
    if (monitor.gt.1) print *, 'crossgrad_mod: testrhoz: finished, sumrhoz =', sumrhoz
  end subroutine testrhoz

  subroutine current_direct
    implicit none
    if (monitor.gt.2) print *, 'crossgrad_mod: current_direct: entered'
    ix(:, 1) = 0.0
    ix(:, 2:ny) = (g(:, 2:ny) - g(:, 1:ny-1)) &
         & + 0.5*(sigma(:, 1:ny-1) + sigma(:, 2:ny)) &
         &   * (phi(:, 2:ny) - phi(:, 1:ny-1))
    iy(1, :) = 0.0
    iy(2:nx, :) = (g(2:nx, :) - g(1:nx-1, :)) &
         & + 0.5*(sigma(1:nx-1, :) + sigma(2:nx, :)) &
         &    * (phi(2:nx, :) - phi(1:nx-1, :))
    if (monitor.gt.1) print *, 'crossgrad_mod: current_direct: finished'
  end subroutine current_direct

  subroutine fluxes
    implicit none
    integer :: k
    if (monitor.gt.2) print *, 'crossgrad_mod: fluxes: entered'
    do k = 1, ncmp
       jx(:, 1, k) = 0.0
       jx(:, 2:ny, k) = d(k) * (rho(:, 2:ny, k) - rho(:, 1:ny-1, k)) &
            & + 0.5 * z(k) * d(k) * (rho(:, 1:ny-1, k) + rho(:, 2:ny, k)) &
            &   * (phi(:, 2:ny) - phi(:, 1:ny-1))
       jy(1, :, k) = 0.0
       jy(2:nx, :, k) = d(k) * (rho(2:nx, :, k) - rho(1:nx-1, :, k)) &
            & + 0.5 * z(k) * d(k) * (rho(1:nx-1, :, k) + rho(2:nx, :, k)) &
            &   * (phi(2:nx, :) - phi(1:nx-1, :))
    end do
    if (monitor.gt.1) print *, 'crossgrad_mod: fluxes: finished'
  end subroutine fluxes

  subroutine current
    implicit none
    integer :: k
    if (monitor.gt.2) print *, 'crossgrad_mod: current_from_fluxes: entered'
    ix = 0.0
    iy = 0.0
    do k = 1, ncmp
       ix = ix + z(k) * jx(:, :, k)
       iy = iy + z(k) * jy(:, :, k)
    end do
    if (monitor.gt.1) print *, 'crossgrad_mod: current_from_fluxes: finished'
  end subroutine current
  
  subroutine divjdt_calc
    implicit none
    integer :: i, j, k
    if (monitor.gt.2) print *, 'crossgrad_mod: divjdt_calc: entered'
    drho = 0.0
    do k = 1, ncmp
       do j = 1, ny
          do i = 1, nx
             if (j.gt.1)  drho(i, j, k) = drho(i, j, k) - jx(i, j, k)
             if (j.lt.ny) drho(i, j, k) = drho(i, j, k) + jx(i, j+1, k)
             if (i.gt.1)  drho(i, j, k) = drho(i, j, k) - jy(i, j, k)
             if (i.lt.nx) drho(i, j, k) = drho(i, j, k) + jy(i+1, j, k)
          end do
       end do
    end do
    drho = drho * dt
    if (monitor.gt.1) print *, 'crossgrad_mod: divjdt_calc: finished'
  end subroutine divjdt_calc

  subroutine correct_for_bc
    implicit none
    integer :: i, j, k1, k2
    real (kind=dp) :: rho1, rho2, del1, del2, b, c, j_ex
    if (monitor.gt.2) print *, 'crossgrad_mod: correct_for_bc: entered'
    if (bc(1) .gt. 0) then ! bc(1:2) sets rho(1, :) 
       k1 = bc(1)
       k2 = bc(2)
       do j = 1, ny
          rho1 = rho(1, j, k1)
          rho2 = rho(1, j, k2)
          ! print *, 1, j, k1, k2, ':: rho1*rho2 = ', rho1, ' * ', rho2, ' = ', rho1*rho2
          del1 = drho(1, j, k1)
          del2 = drho(1, j, k2)
          b = rho1 + rho2 + del1 + del2
          c = rho1*del2 + rho2*del1 + del1*del2
          j_ex = 0.5*(sqrt(b**2 - 4*c) - b)
          drho(1, j, k1) = drho(1, j, k1) + j_ex
          drho(1, j, k2) = drho(1, j, k2) + j_ex
       end do
    end if
    if (bc(3) .gt. 0) then ! bc(3:4) sets rho(nx, :) 
       k1 = bc(3)
       k2 = bc(4)
       do j = 1, ny
          rho1 = rho(nx, j, k1)
          rho2 = rho(nx, j, k2)
          ! print *, nx, j, k1, k2, ':: rho1*rho2 = ', rho1, ' * ', rho2, ' = ', rho1*rho2
          del1 = drho(nx, j, k1)
          del2 = drho(nx, j, k2)
          b = rho1 + rho2 + del1 + del2
          c = rho1*del2 + rho2*del1 + del1*del2
          j_ex = 0.5*(sqrt(b**2 - 4*c) - b)
          drho(nx, j, k1) = drho(nx, j, k1) + j_ex
          drho(nx, j, k2) = drho(nx, j, k2) + j_ex
       end do
    end if
    if (bc(5) .gt. 0) then ! bc(5:6) sets rho(:, 1) 
       k1 = bc(5)
       k2 = bc(6)
       do i = 1, nx
          rho1 = rho(i, 1, k1)
          rho2 = rho(i, 1, k2)
          ! print *, i, 1, k1, k2, ':: rho1*rho2 = ', rho1, ' * ', rho2, ' = ', rho1*rho2
          del1 = drho(i, 1, k1)
          del2 = drho(i, 1, k2)
          b = rho1 + rho2 + del1 + del2
          c = rho1*del2 + rho2*del1 + del1*del2
          j_ex = 0.5*(sqrt(b**2 - 4*c) - b)
          drho(i, 1, k1) = drho(i, 1, k1) + j_ex
          drho(i, 1, k2) = drho(i, 1, k2) + j_ex
       end do
    end if
    if (bc(7) .gt. 0) then ! bc(7:8) sets rho(:, ny) 
       k1 = bc(7)
       k2 = bc(8)
       do i = 1, nx
          rho1 = rho(i, ny, k1)
          rho2 = rho(i, ny, k2)
          ! print *, i, ny, k1, k2, ':: rho1*rho2 = ', rho1, ' * ', rho2, ' = ', rho1*rho2
          del1 = drho(i, ny, k1)
          del2 = drho(i, ny, k2)
          b = rho1 + rho2 + del1 + del2
          c = rho1*del2 + rho2*del1 + del1*del2
          j_ex = 0.5*(sqrt(b**2 - 4*c) - b)
          drho(i, ny, k1) = drho(i, ny, k1) + j_ex
          drho(i, ny, k2) = drho(i, ny, k2) + j_ex
       end do
    end if
    if (monitor.gt.1) print *, 'crossgrad_mod: correct_for_bc: finished'
  end subroutine correct_for_bc
  
  subroutine neutralize(kk)
    implicit none
    integer, intent(in) :: kk
    integer :: k
    if (monitor.gt.2) print *, 'crossgrad_mod: neutralize: entered, kk = ', kk
    drho(:, :, kk) = 0.0
    do k = 1, ncmp
       if (k.ne.kk) drho(:, :, kk) = drho(:, :, kk) + z(k) * rho(:, :, k)
    end do
    rho(:, :, kk) = - z(kk) * drho(:, :, kk)
    if (monitor.gt.1) print *, 'crossgrad_mod: neutralize: finished'
  end subroutine neutralize
  
  subroutine divcurl
    implicit none
    integer :: i, j
    if (monitor.gt.2) print *, 'crossgrad_mod: divcurl: entered'
    do j = 1, ny
       do i = 1, nx
          div(i, j) = - ix(i, j) - iy(i, j)
          if (j.lt.ny) div(i, j) = div(i, j) + ix(i, j+1)
          if (i.lt.nx) div(i, j) = div(i, j) + iy(i+1, j)
       end do
    end do
    curl = 0.0
    do j = 2, ny
       do i = 2, nx
          curl(i, j) = iy(i, j) - ix(i, j) - iy(i, j-1) + ix(i-1, j)
       end do
    end do
    curl(1, :) = curl(2, :)
    curl(:, 1) = curl(:, 2)
    curl(1, 1) = 0.5*(curl(1, 2) + curl(2, 1))
    if (monitor.gt.1) print *, 'crossgrad_mod: divcurl: finished'
  end subroutine divcurl

  subroutine diagnostics
    implicit none
    real (kind=dp) :: sumcurl, lineint, sumrhoz
    call current_direct ; call divcurl
    print *, 'diagnostic: (direct) |I|, |∇ ⨯ I|, |∇·I| =', &
         & sqrt(sum(ix**2 + iy**2)), sqrt(sum(curl**2)), sqrt(sum(div**2))
    call fluxes ; call current ; call divcurl
    print *, 'diagnostic: (fluxes) |I|, |∇ ⨯ I|, |∇·I| =', &
         & sqrt(sum(ix**2 + iy**2)), sqrt(sum(curl**2)), sqrt(sum(div**2))
    sumcurl = sum(curl)
    lineint = sum(iy(:, ny)) - sum(ix(nx, :)) - sum(iy(:, 1)) + sum(ix(1, :))
    print *, 'diagnostic: Σ ∇ ⨯ I - Σ I·dl = ',sumcurl, -lineint, &
         & ' = ', sumcurl-lineint
    call testrhoz(sumrhoz)
    print *, 'diagnostic: Σi zi ρi = ', sumrhoz
  end subroutine diagnostics

  subroutine write_array(array_file)
    implicit none
    character (len=*), intent(in) :: array_file
    open (10, file=array_file, form='unformatted', action='write')
    write (10) rho, phi, ix, iy, div, curl
    close (10)
    print *, 'written array data to ', array_file
  end subroutine write_array

  subroutine read_array(array_file)
    implicit none
    character (len=*), intent(in) :: array_file
    open (11, file=array_file, form='unformatted', action='read')
    read (11) rho, phi
    close (11)
    print *, 'read array data from ', array_file
  end subroutine read_array

  subroutine write_rho(filename)
    implicit none
    integer :: i, j, k
    character (len=*), intent(in) :: filename
    open (10, file=filename, action='write')
    do k = 1, ncmp
       do j = 1, ny
          do i = 1, nx
             write (10, '(3I3, " : ", E20.5)') i, j, k, rho(i, j, k)
          end do
       end do
    end do
    close (10)
    print *, 'written density data to ', filename
  end subroutine write_rho
  
  subroutine write_fingerprint(finger_file, i0, j0, j1)
    implicit none
    integer, intent(in) :: i0, j0, j1
    character (len=*), intent(in) :: finger_file
    open (10, file=finger_file, action='write')
    write (10, '(3I3)') i0-1, j0-1, j1-1
    write (10, *) 'sigma = ', sigma(i0, j0:j1)
    write (10, *) '    g = ', g(i0, j0:j1)
    write (10, *) '  phi = ', phi(i0, j0:j1)
    close (10)
    print *, 'written fingerprint to ', finger_file
  end subroutine write_fingerprint

  subroutine print_status(step)
    implicit none
    integer, intent(in) :: step
    real (kind=dp) :: sumrhoz, normi, normdiv, normcurl
    call testrhoz(sumrhoz)
    normi = sqrt(sum(ix**2 + iy**2))
    normdiv = sqrt(sum(div**2))
    normcurl = sqrt(sum(curl**2))
    print '(I5, F10.3, 2F15.6, 2E10.2)', step, step*dt, normi, normcurl, normdiv, sumrhoz
  end subroutine print_status
    
end module crossgrad


  
