module MMeasurementInterface
   use MDiagramMeasurements
   use MBathMeasurements
   use MLocalMeasurements
   use MWormMeas
   use MSamplerD
   use MSamplerZ
   use csamplers
   implicit none

   type, extends(TMeasurement) :: TAllTargetsCollection
      type(TSignM) :: mSign
      type(TGreenM) :: mGreen
      type(TG2iw) :: mG2iw
      type(TG4iw_full) :: mG4iw_full
      type(TDensityMatrixM) :: mDM
      type(TExpansionOrderM) :: mExpOrd
      !!!TODOcompl: reactivate n(tau)n(0)
      !type(TNtauN0M) :: mNtauN0

      ! worm measurements
      type(TWormP2phM) :: mWormP2ph
      type(TWormP2ppM) :: mWormP2pp
      type(TWormQQM) :: mWormQQ
      type(TWormRamanM) :: mWormRaman
      class(TMeasurement), allocatable :: mWormCustom
   contains
      procedure, pass :: measure => targetcoll_measure_all
      procedure, pass :: decouple => targetcoll_decouple_all
   end type TAllTargetsCollection

   ! FIXME: as soon as possible, change from module variable to
   ! function parameter or sth like that
   type(TAllTargetsCollection), target :: meastargets
contains

   ! FIXME: unusable as only some actually need to be measured in a
   ! given calculation which is currently checked by ctqmc_measure
   subroutine targetcoll_measure_all(meas, diag)
      class(TAllTargetsCollection), intent(inout) :: meas
      type(TDiagram), intent(in)                  :: diag

      call meas%mSign%measure(diag)
      call meas%mGreen%measure(diag)
      call meas%mG2iw%measure(diag)
      call meas%mG4iw_full%measure(diag)
      call meas%mDM%measure(diag)
      call meas%mExpOrd%measure(diag)
      !call meas%mNtauN0%measure(diag)
   end subroutine targetcoll_measure_all

   subroutine targetcoll_decouple_all(meas)
      class(TAllTargetsCollection), intent(inout) :: meas

      call meas%mSign%decouple()
      call meas%mGreen%decouple()
      call meas%mG2iw%decouple()
      call meas%mG4iw_full%decouple()
      call meas%mDM%decouple()
      call meas%mExpOrd%decouple()
      !call meas%mNtauN0%decouple()

      call meas%mWormP2ph%decouple()
      call meas%mWormP2pp%decouple()
      call meas%mWormQQ%decouple()
      call meas%mWormRaman%decouple()
      if (allocated(meas%mWormCustom)) then
         call meas%mWormCustom%decouple()
         deallocate(meas%mWormCustom)
      end if
   end subroutine targetcoll_decouple_all

   ! Sets the accumulator referred to by cacc as target for sign
   ! measurement of the target collection referred to by ctargets.
   subroutine set_accumulator_sign(ctargets, cacc) bind(C, name='set_accumulator_sign')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mSign%c_set_accumulator(cacc)
   end subroutine set_accumulator_sign

   ! Add a sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for Green's function
   ! measurement of the target collection referred to by ctargets.
   subroutine add_sampler_G(ctargets, csamp, cacc) bind(C, name='add_sampler_G')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mGreen%c_add_sampler(csamp, cacc)
   end subroutine add_sampler_G

   ! Add a sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for the 2-frequency Green's function
   ! measurement of the target collection referred to by ctargets.
   subroutine add_sampler_G2iw(ctargets, csamp, cacc) bind(C, name='add_sampler_G2iw')
      type(c_ptr), value                    :: ctargets      
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mG2iw%c_add_sampler(csamp, cacc)
   end subroutine add_sampler_G2iw

   ! Sets the accumulator referred to by cacc as target for the four-fermionic-frequency
   ! 2P - Green's function measurement of the target collection referred to by
   ! ctargets.
   subroutine set_accumulator_G4iw_full(ctargets, cacc, cngrid, niw, nbands, beta) bind(C, name='set_accumulator_G4iw_full')
      type(c_ptr), value                    :: ctargets      
      type(c_ptr), value                    :: cacc
      type(TAllTargetsCollection), pointer  :: targets
      integer(c_int64_t), value     , intent(in)             :: niw, nbands
      real(c_double), value, intent(in) :: beta
      integer(c_int64_t), intent(inout)           :: cngrid(niw)
      class(ZSampler2D), allocatable              :: g2samp
      integer(c_int64_t), pointer           :: ngrid(:)
      real(c_double) ::  coeff(2,2)
      integer :: niws(2), statsign(2), n, ncomp, nmax

      call c_f_pointer(ctargets, targets)

      if (.not. allocated(targets%mG4iw_full%sparse_freq)) then
         allocate(targets%mG4iw_full%sparse_freq(niw))
      end if
      targets%mG4iw_full%sparse_freq = cngrid(1:niw)

      nmax = maxval(targets%mG4iw_full%sparse_freq)
      n = (nmax - 1) / 2 + 1

      ! init sampler for G2(v,v')
      coeff(:,:) = 0
      coeff(1,1) = 1
      coeff(2,2) = 1
      statsign(1) = -1   
      statsign(2) = -1   
      niws(1) = n
      niws(2) = n
      ncomp = 4*nbands*nbands

      call init_fast_matsubara_sampler(targets%mG4iw_full%g2_sampler, ncomp, beta, statsign, niws, coeff, .true.)
      
      call targets%mG4iw_full%c_set_accumulator(cacc)
   end subroutine set_accumulator_G4iw_full

   ! Sets the accumulator referred to by cacc as target for density
   ! matrix measurement of the target collection referred to by
   ! ctargets.
   subroutine set_accumulator_dm(ctargets, cacc) bind(C, name='set_accumulator_dm')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mDM%c_set_accumulator(cacc)
   end subroutine set_accumulator_dm

   ! Sets the accumulator referred to by cacc as target for expansion
   ! order measurement of the target collection referred to by
   ! ctargets.
   subroutine set_accumulator_expord(ctargets, cacc) bind(C, name='set_accumulator_expord')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mExpOrd%c_set_accumulator(cacc)
   end subroutine set_accumulator_expord

   ! Sets the accumulator referred to by cacc as target for
   ! susceptibility measurement of the target collection referred to
   ! by ctargets.
   subroutine set_accumulator_ntaun0(ctargets, cacc) bind(C, name='set_accumulator_ntaun0')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      !call targets%mNtauN0%c_set_accumulator(cacc)
   end subroutine set_accumulator_ntaun0

   ! Add a sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for susceptibility
   ! measurement of the target collection referred to by ctargets.
   subroutine add_sampler_wormP2ph(ctargets, csamp, cacc) bind(C, name='add_sampler_wormP2ph')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mWormP2ph%c_add_sampler(csamp, cacc)
   end subroutine add_sampler_wormP2ph

   ! Add a sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for susceptibility
   ! measurement of the target collection referred to by ctargets.
   subroutine add_sampler_wormP2pp(ctargets, csamp, cacc) bind(C, name='add_sampler_wormP2pp')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mWormP2pp%c_add_sampler(csamp, cacc)
   end subroutine add_sampler_wormP2pp

   ! Add a sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for symm. impr. IE QQ
   ! measurement of the target collection referred to by ctargets.
   subroutine add_sampler_wormQQ(ctargets, csamp, cacc) bind(C, name='add_sampler_wormQQ')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mWormQQ%c_add_sampler(csamp, cacc)
   end subroutine add_sampler_wormQQ

   ! Add a sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for Raman diagram
   ! measurement of the target collection referred to by ctargets.
   subroutine add_sampler2d_wormRaman(ctargets, csamp, cacc) bind(C, name='add_sampler2d_wormRaman')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      call targets%mWormRaman%c_add_sampler(csamp, cacc)
   end subroutine add_sampler2d_wormRaman

   ! FIXME: check whether sampler dimensions are compatible with the
   ! amount of custom worm taus that have been set

   ! Add a 1d sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for a registered
   ! custom worm estimator of the target collection referred to by
   ! ctargets.
   subroutine add_sampler1d_wormCustom(ctargets, csamp, cacc) &
        bind(C, name='add_sampler1d_wormCustom')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      if (.not. allocated(targets%mWormCustom)) then
         allocate(TWormCustom1DM :: targets%mWormCustom)
      end if
      select type (s1d => targets%mWormCustom)
      type is (TWormCustom1DM)
         call s1d%c_add_sampler(csamp, cacc)
      class default
         error stop "add_sampler1d_wormCustom: wrong sampler dimensions"
      end select
   end subroutine add_sampler1d_wormCustom

   ! Add a 2d sampler referred to by csamp and its target
   ! accumulator referred to by cacc to the sampler register for a
   ! registered custom worm estimator of the target collection
   ! referred to by ctargets.
   subroutine add_sampler2d_wormCustom(ctargets, csamp, cacc) &
        bind(C, name='add_sampler2d_wormCustom')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      if (.not. allocated(targets%mWormCustom)) then
         allocate(TWormCustom2DM :: targets%mWormCustom)
      end if
      select type (s2d => targets%mWormCustom)
      type is (TWormCustom2DM)
         call s2d%c_add_sampler(csamp, cacc)
      class default
         error stop "add_sampler2d_wormCustom: wrong sampler dimensions"
      end select
   end subroutine add_sampler2d_wormCustom

   ! Add a 3d sampler referred to by csamp and its target accumulator
   ! referred to by cacc to the sampler register for a registered
   ! custom worm estimator of the target collection referred to by
   ! ctargets.
   subroutine add_sampler3d_wormCustom(ctargets, csamp, cacc) &
        bind(C, name='add_sampler3d_wormCustom')
      type(c_ptr), value                    :: ctargets
      type(c_ptr), value                    :: csamp, cacc
      type(TAllTargetsCollection), pointer  :: targets

      call c_f_pointer(ctargets, targets)
      if (.not. allocated(targets%mWormCustom)) then
         allocate(TWormCustom3DM :: targets%mWormCustom)
      end if
      select type (s3d => targets%mWormCustom)
      type is (TWormCustom3DM)
         call s3d%c_add_sampler(csamp, cacc)
      class default
         error stop "add_sampler3d_wormCustom: wrong sampler dimensions"
      end select
   end subroutine add_sampler3d_wormCustom

   ! Get a C pointer to the module variable containing all measurement
   ! targets used by CTQMC
   function get_global_meastargets_FIXME() result(modvar) bind(C, name='get_global_meastargets_FIXME')
      type(c_ptr) :: modvar

      modvar = c_loc(meastargets)
   end function get_global_meastargets_FIXME
end module MMeasurementInterface
