program crossgrad_solver

  use crossgrad
  
  implicit none
  integer :: skip, step
  integer :: njacobi, niter, nmon, nstep, ndiag, nsave
  character (len=32) :: rc_file
  character (len=48) :: array_file
  
  call get_command_argument(1, rc_file)

  if (len_trim(rc_file) == 0) then
     print *, 'missing run control file'
     call exit(1)
  end if

  call read_rc(rc_file)

  ! re-read the run control to extract the array data file etc

  open (5, file=rc_file, action='read')
  do skip = 1, 7 ! skip to past dt setting
     read (5, *)
  end do
  read (5, *) njacobi, niter, nmon
  read (5, *) nstep, ndiag, nsave
  read (5, '(A)') array_file
  close (5)

  print *, 'njacobi, niter, nmon = ', njacobi, niter, nmon
  print *, 'nstep, ndiag, nsave = ', nstep, ndiag, nsave
  print *, 'array file = ', trim(array_file)
  
  if (all(bc.eq.0)) print *, 'all boundary conditions zero flux'
  if (any(bc.gt.0)) print *, 'some chemical potential boundary conditions'

  call read_array(trim(array_file))

!  call write_rho('test_rho.dat')
  
  do step = 0, nstep
     
     call converge(njacobi, niter, nmon)
     if (.not. converged) print *, 'WARNING, not converged'

     if (mod(step, ndiag).eq.0) call diagnostics

     if (mod(step, nsave).eq.0) then
        write (array_file, '(A, I5.5, A)') trim(run_name)//'_f', step/nsave, '.dat'
        call write_array(array_file)
     end if

     call fluxes ; call current ; call divcurl ;
     call print_status(step)
     
     if (step.lt.nstep) then
        call divjdt_calc
        if (any(bc.gt.0)) call correct_for_bc
        rho = rho + drho
     end if

  end do
  
end program crossgrad_solver
