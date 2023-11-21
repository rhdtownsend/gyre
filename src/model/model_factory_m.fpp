! Module  : model_factory_m
! Purpose : factory procedures for model_t
!
! Copyright 2016-2023 Rich Townsend & The GYRE Team
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

module model_factory_m

  ! Uses

  use kinds_m

  use hom_model_m
  use model_m
  use model_par_m

  use amdl_file_m
  use famdl_file_m
  use fgong_file_m
  use losc_file_m
  use mesa_file_m
  use parfait_file_m
  use parfaitd_file_m
  use osc_file_m
  use wdec_file_m
  $if ($HDF5)
  use b3_file_m
  use gsm_file_m
  use poly_file_m
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
       case ('WDEC')
          call read_wdec_model(ml_p, ml)
       case default
          $ABORT(Invalid file_format)
       end select

    case ('POLY')

       $if($HDF5)
       call read_poly_model(ml_p, ml)
       $else
       $ABORT(No HDF5 support, therefore cannot read POLY files)
       $endif

    case ('PARFAIT')

       $if($HDF5)
       call read_parfait_model(ml_p, ml)
       $else
       $ABORT(No HDF5 support, therefore cannot read PARFAIT files)
       $endif

    case ('PARFAITD')

       $if($HDF5)
       call read_parfaitd_model(ml_p, ml)
       $else
       $ABORT(No HDF5 support, therefore cannot read PARFAITD files)
       $endif

    case ('HOM')

       allocate(ml, SOURCE=hom_model_t(ml_p))

    case default

       $ABORT(Invalid model_type)

    end select

    ! Finish

    return

  end function model_t_

end module model_factory_m
