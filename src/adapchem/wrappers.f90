! wrappers.f90: A set of wrappers around Fortran intrinsics and
! Chemkin/AdapChem subroutines to make them easier to call from C/C++.

subroutine f90_open_unformatted(unit, name)
    implicit none
    character*(*) name
    integer unit

    open(unit, status='old', form='unformatted', file=name)

end subroutine f90_open_unformatted

subroutine f90_open_formatted(unit, name)
    implicit none
    character*(*) name
    integer unit

    open(unit, status='old', form='formatted', file=name)

end subroutine f90_open_formatted

subroutine f90_close(unit)
    implicit none
    integer unit

    close(unit)

end subroutine f90_close


subroutine f90_ckinit(in_unit, out_unit, len_iwork, iwork, len_rwork, rwork, len_cwork, cwork)
    implicit none

    integer in_unit, out_unit
    integer len_iwork, len_rwork, len_cwork

    real*8:: rwork(len_rwork)
    integer iwork(len_iwork)
    character cwork(len_cwork)*16

    call ckinit(len_iwork, len_rwork, len_cwork, in_unit, out_unit, iwork, rwork, cwork)

end subroutine f90_ckinit

! Proof-of-concept for accessing data in common blocks
subroutine f90_setrtol(rtol)
    implicit none
    double precision rtol, resid_tol

    common /adpchm1/ resid_tol
    resid_tol = rtol

end subroutine f90_setrtol
