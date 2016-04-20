! Module   : gyre_model_factory
! Purpose  : factory procedures for model_t
!
! Copyright 2016 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

$include 'core.inc'

module gyre_model_factory

  ! Uses

  use core_kinds

  use gyre_evol_model
  use gyre_hom_model
  use gyre_model
  use gyre_model_par
  use gyre_poly_model

  use gyre_amdl_file
  use gyre_famdl_file
  use gyre_fgong_file
  use gyre_losc_file
  use gyre_mesa_file
  use gyre_osc_file
  $if ($HDF5)
  use gyre_b3_file
  use gyre_gsm_file
  use gyre_poly_file
  $endif

  use ISO_FORTRAN_ENV
  
  ! No implicit typing

  implicit none

  ! Interfaces

  interface model_t
     module procedure model_t_
  end interface model_t

  ! Access specifiers

  private

  public :: model_t

  ! Procedures

contains

  function model_t_ (ml_p) result (ml)


    type(model_par_t), intent(in) :: ml_p
    class(model_t), pointer       :: ml

    ! Construct the model_t

    select case (ml_p%model_type)
    case ('EVOL')

       select case (ml_p%file_format)
       case ('AMDL')
          call read_amdl_model(ml_p, ml)
       case ('B3')
          $if($HDF5) 
          call read_b3_model(ml_p, ml)
          $else
          $ABORT(No HDF5 support, therefore cannot read B3-format files)
          $endif
       case ('FAMDL')
          call read_famdl_model(ml_p, ml)
       case ('FGONG')
          call read_fgong_model(ml_p, ml)
       case ('GSM')
          $if($HDF5)
          call read_gsm_model(ml_p, ml)
          $else
          $ABORT(No HDF5 support, therefore cannot read GSM-format files)
          $endif
       case ('LOSC')
          call read_losc_model(ml_p, ml)
       case ('MESA')
          call read_mesa_model(ml_p, ml)
       case ('OSC')
          call read_osc_model(ml_p, ml)
       case default
          $ABORT(Invalid file_format)
       end select

    case ('POLY')

       $if($HDF5)
       call read_poly_model(ml_p, ml)
       $else
       $ABORT(No HDF5 support, therefore cannot read POLY files)
       $endif

    case ('HOM')

       allocate(ml, SOURCE=hom_model_t(ml_p))

    case default

       $ABORT(Invalid model_type)

    end select

    ! Finish

    return

  end function model_t_

  !****

  

end module gyre_model_factory
