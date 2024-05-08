! radiation_config.F90 - Derived type to configure the radiation scheme
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2017-07-22  R. Hogan  Added Yi et al. ice optics model
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-03-15  R. Hogan  Added logicals controlling surface spectral treatment
!   2018-08-29  R. Hogan  Added monochromatic single-scattering albedo / asymmetry factor
!   2018-09-03  R. Hogan  Added min_cloud_effective_size
!   2018-09-04  R. Hogan  Added encroachment_scaling
!   2018-09-13  R. Hogan  Added IEncroachmentFractal
!   2019-01-02  R. Hogan  Added Cloudless solvers
!   2019-01-14  R. Hogan  Added out_of_bounds_[1,2,3]d for checker routines
!   2019-01-18  R. Hogan  Added albedo weighting
!   2019-02-03  R. Hogan  Added ability to fix out-of-physical-bounds inputs
!   2019-02-10  R. Hogan  Renamed "encroachment" to "entrapment"
!   2020-05-18  R. Hogan  Moved out_of_bounds_* to radiation_check.F90
!   2021-07-04  R. Hogan  Numerous changes for ecCKD and general cloud/aerosol optics
!   2022-01-18  P. Ukkonen Added support for RRTMGP gas optics, with optional neural networks
!
! Note: The aim is for ecRad in the IFS to be as similar as possible
! to the offline version, so if you make any changes to this or any
! files in this directory, please inform Robin Hogan.
!

#include "ecrad_config.h"

module radiation_config

  use parkind1

  use radiation_cloud_optics_data
  use radiation_general_cloud_optics_data
  use radiation_aerosol_optics_data
  use radiation_pdf_sampler
  use radiation_cloud_cover
  use radiation_ecckd
  use mo_gas_optics_rrtmgp
  use mod_network_rrtmgp

  implicit none
  

  ! Configuration codes: use C-style enumerators to avoid having to
  ! remember the numbers

  ! Solvers: can be specified for longwave and shortwave
  ! independently, except for "Homogeneous", which must be the same
  ! for both
  enum, bind(c)
     enumerator ISolverCloudless, ISolverHomogeneous, ISolverMcICA, &
          &     ISolverSpartacus, ISolverTripleclouds, ISolverMcICAACC
  end enum
  character(len=*), parameter :: SolverName(0:5) = (/ 'Cloudless   ', &
       &                                              'Homogeneous ', &
       &                                              'McICA       ', &
       &                                              'SPARTACUS   ', &
       &                                              'Tripleclouds', &
       &                                              'McICA ACC   ' /)

  ! SPARTACUS shortwave solver can treat the reflection of radiation
  ! back up into different regions in various ways
  enum, bind(c)
     enumerator &
       & IEntrapmentZero, &     ! No entrapment, as Tripleclouds
       & IEntrapmentEdgeOnly, & ! Only radiation passed through cloud edge is horizontally homogenized
       & IEntrapmentExplicit, & ! Estimate horiz migration dist, account for fractal clouds
       & IEntrapmentExplicitNonFractal, & ! As above but ignore fractal nature of clouds
       & IEntrapmentMaximum ! Complete horizontal homogenization within regions (old SPARTACUS assumption)
  end enum

  ! Names available in the radiation namelist for variable
  ! sw_entrapment_name
  character(len=*), parameter :: EntrapmentName(0:4)   = [ 'Zero       ', &
       &                                                   'Edge-only  ', &
       &                                                   'Explicit   ', &
       &                                                   'Non-fractal', &
       &                                                   'Maximum    ' ]
  ! For backwards compatibility, the radiation namelist also supports
  ! the equivalent variable sw_encroachment_name with the following
  ! names
  character(len=*), parameter :: EncroachmentName(0:4) = [ 'Zero    ', &
       &                                                   'Minimum ', &
       &                                                   'Fractal ', &
       &                                                   'Computed', &
       &                                                   'Maximum ' ]

  ! Two-stream models
  ! This is not configurable at run-time

  ! Gas models
  enum, bind(c)
     enumerator IGasModelMonochromatic, IGasModelIFSRRTMG, IGasModelECCKD, IGasModelRRTMGP, IGasModelRRTMGP_NN
  end enum
  character(len=*), parameter :: GasModelName(0:4) = (/ 'Monochromatic', &
       &                                                'RRTMG-IFS    ', &
       &                                                'ECCKD        ', &
       &                                                'RRTMGP       ', & !/)
       &                                                'RRTMGP-NN    ' /)

  ! Liquid cloud optics models for use with RRTMG gas optics
  enum, bind(c)
     enumerator ILiquidModelMonochromatic, &
          &     ILiquidModelSOCRATES, ILiquidModelSlingo
  end enum
  character(len=*), parameter :: LiquidModelName(0:2) = (/ 'Monochromatic', &
       &                                                   'SOCRATES     ', &
       &                                                   'Slingo       ' /)

  ! Ice optics models for use with RRTMG gas optics. Note that of the
  ! "Baran" parameterizations, only Baran2016 is published (Baran,
  ! J. Climate, 2016) - the others are experimental and not
  ! recommended.
  enum, bind(c)
     enumerator IIceModelMonochromatic, IIceModelFu, &
          &  IIceModelBaran, IIceModelBaran2016, IIceModelBaran2017,   &
          &  IIceModelYi
  end enum
  character(len=*), parameter :: IceModelName(0:5) = (/ 'Monochromatic         ', &
       &                                                'Fu-IFS                ', &
       &                                                'Baran-EXPERIMENTAL    ', &
       &                                                'Baran2016             ', &
       &                                                'Baran2017-EXPERIMENTAL', &
       &                                                'Yi                    ' /)

  ! Cloud PDF distribution shapes
  enum, bind(c)
    enumerator IPdfShapeLognormal, IPdfShapeGamma
  end enum
  character(len=*), parameter :: PdfShapeName(0:1) = (/ 'Lognormal', &
       &                                                'Gamma    ' /)

  ! Maximum number of different aerosol types that can be provided
  integer, parameter :: NMaxAerosolTypes = 256

  ! Maximum number of different cloud types that can be provided
  integer, parameter :: NMaxCloudTypes = 12

  ! Maximum number of shortwave albedo and longwave emissivity
  ! intervals
  integer, parameter :: NMaxAlbedoIntervals = 256

  ! Length of string buffer for printing config information
  integer, parameter :: NPrintStringLen = 60

  !---------------------------------------------------------------------
  ! Derived type containing all the configuration information needed
  ! to run the radiation scheme.  The intention is that this is fixed
  ! for a given model run.  The parameters are to list first those
  ! quantities that can be set directly by the user, for example using a
  ! namelist, and second those quantities that are computed afterwards
  ! from the user-supplied numbers, especially the details of the gas
  ! optics model.
  type config_type
    ! USER-CONFIGURABLE PARAMETERS

    ! Scale the solar spectrum per band (or g-point if
    ! do_cloud_aerosol_per_sw_g_point=true) via vector
    ! single_level%spectral_solar_scaling
    logical :: use_spectral_solar_scaling = .false.

    ! Modify the solar spectrum per g-point to account for the current
    ! phase of the solar cycle, via scalar
    ! single_level%spectral_solar_cycle_multiplier
    logical :: use_spectral_solar_cycle = .false.

    ! Directory in which gas, cloud and aerosol data files are to be
    ! found
    character(len=511) :: directory_name = '.'

    ! If this is true then support arbitrary hydrometeor types (not
    ! just ice and liquid) and arbitrary spectral discretization (not
    ! just RRTMG). It is required that this is true if the ecCKD gas
    ! optics model is selected. General cloud optics has only been
    ! available from ecRad version 1.5.
    logical :: use_general_cloud_optics = .true.

    ! If this is true then support aerosol properties at an arbitrary
    ! spectral discretization (not just RRTMG). It is required that
    ! this is true if the ecCKD gas optics model is selected.
    logical :: use_general_aerosol_optics = .true.

    ! Cloud is deemed to be present in a layer if cloud fraction
    ! exceeds this value
    real(jprb) :: cloud_fraction_threshold = 1.0e-6_jprb
    ! ...and total cloud water mixing ratio exceeds this value
    real(jprb) :: cloud_mixing_ratio_threshold = 1.0e-9_jprb

    ! Overlap scheme
    integer :: i_overlap_scheme = IOverlapExponentialRandom

    ! Use the Shonk et al. (2010) "beta" overlap parameter, rather
    ! than the "alpha" overlap parameter of Hogan and Illingworth
    ! (2000)?
    logical :: use_beta_overlap = .false.

    ! Use a more vectorizable McICA cloud generator, at the expense of
    ! more random numbers being generated?  This is the default on NEC
    ! SX.
    logical :: use_vectorizable_generator = .false.

    ! Shape of sub-grid cloud water PDF
    integer :: i_cloud_pdf_shape = IPdfShapeGamma

    ! The ratio of the overlap decorrelation length for cloud
    ! inhomogeneities to the overlap decorrelation length for cloud
    ! boundaries.  Observations suggest this has a value of 0.5
    ! (e.g. from the decorrelation lengths of Hogan and Illingworth
    ! 2003 and Hogan and Illingworth 2000).
    real(jprb) :: cloud_inhom_decorr_scaling = 0.5_jprb

    ! Factor controlling how much of the cloud edge length interfaces
    ! directly between the clear-sky region (region a) and the
    ! optically thick cloudy region (region c).  If Lxy is the length
    ! of the interfaces between regions x and y, and Lab and Lbc have
    ! been computed already, then
    !   Lac=clear_to_thick_fraction*min(Lab,Lbc).
    real(jprb) :: clear_to_thick_fraction = 0.0_jprb

    ! Factor allowing lateral transport when the sun is close to
    ! overhead; consider atand(overhead_sun_factor) to be the number
    ! of degrees that the sun angle is perturbed from zenith for the
    ! purposes of computing lateral transport.  A value of up to 0.1
    ! seems to be necessary to account for the fact that some forward
    ! scattered radiation is treated as unscattered by delta-Eddington
    ! scaling; therefore it ought to have the chance to escape.
    real(jprb) :: overhead_sun_factor = 0.0_jprb

    ! Minimum gas optical depth in a single layer at any wavelength,
    ! for stability
    real(jprb) :: min_gas_od_lw = 1.0e-15_jprb
    real(jprb) :: min_gas_od_sw = 0.0_jprb

    ! Maximum gas optical depth in a layer before that g-point will
    ! not be considered for 3D treatment: a limit is required to avoid
    ! expensive computation of matrix exponentials on matrices with
    ! large elements
    real(jprb) :: max_gas_od_3d = 8.0_jprb

    ! Maximum total optical depth of a cloudy region for stability:
    ! optical depth will be capped at this value in the SPARTACUS
    ! solvers
    real(jprb) :: max_cloud_od_lw = merge (10.3_jprb, 16.0_jprb, jprb /= jprd)
    real(jprb) :: max_cloud_od_sw = merge (11.5_jprb, 16.0_jprb, jprb /= jprd)
    ! How much longwave scattering is included?
    logical :: do_lw_cloud_scattering = .true.
    logical :: do_lw_aerosol_scattering = .true.

    ! Number of regions used to describe clouds and clear skies. A
    ! value of 2 means one clear and one cloudy region, so clouds are
    ! horizontally homogeneous, while a value of 3 means two cloudy
    ! regions with different optical depth, thereby representing
    ! inhomogeneity via the Shonk & Hogan (2008) "Tripleclouds"
    ! method.
    integer :: nregions = 3

    ! Code specifying the solver to be used: use the enumerations
    ! defined above
    integer :: i_solver_sw = ISolverMcICA
    integer :: i_solver_lw = ISolverMcICA

    ! Do shortwave delta-Eddington scaling on the cloud-aerosol-gas
    ! mixture (as in the original IFS scheme), rather than the more
    ! correct approach of separately scaling the cloud and aerosol
    ! scattering properties before merging with gases.  Note that
    ! .true. is not compatible with the SPARTACUS solver.
    logical :: do_sw_delta_scaling_with_gases = .false.

    ! Codes describing the gas model
    integer :: i_gas_model_sw = IGasModelIFSRRTMG
    integer :: i_gas_model_lw = IGasModelIFSRRTMG

    ! Optics if i_gas_model==IGasModelMonochromatic.
    ! The wavelength to use for the Planck function in metres. If this
    ! is positive then the output longwave fluxes will be in units of
    ! W m-2 um-1.  If this is zero or negative (the default) then
    ! sigma*T^4 will be used and the output longwave fluxes will be in
    ! W m-2.
    real(jprb) :: mono_lw_wavelength = -1.0_jprb
    ! Total zenith optical depth of the atmosphere in the longwave and
    ! shortwave, distributed vertically according to the pressure.
    ! Default is zero.
    real(jprb) :: mono_lw_total_od = 0.0_jprb
    real(jprb) :: mono_sw_total_od = 0.0_jprb
    ! Single-scattering albedo and asymmetry factor: values typical
    ! for liquid clouds with effective radius of 10 microns, at (SW)
    ! 0.55 micron wavelength and (LW) 10.7 microns wavelength
    real(jprb) :: mono_sw_single_scattering_albedo = 0.999999_jprb
    real(jprb) :: mono_sw_asymmetry_factor = 0.86_jprb
    real(jprb) :: mono_lw_single_scattering_albedo = 0.538_jprb
    real(jprb) :: mono_lw_asymmetry_factor = 0.925_jprb

    ! Codes describing particle scattering models
    integer :: i_liq_model = ILiquidModelSOCRATES
    integer :: i_ice_model = IIceModelBaran

    ! The mapping from albedo/emissivity intervals to SW/LW bands can
    ! either be done by finding the interval containing the central
    ! wavenumber of the band (nearest neighbour), or by a weighting
    ! according to the spectral overlap of each interval with each
    ! band
    logical :: do_nearest_spectral_sw_albedo = .false.
    logical :: do_nearest_spectral_lw_emiss  = .false.

    ! User-defined monotonically increasing wavelength bounds (m)
    ! between input surface albedo/emissivity intervals. Implicitly
    ! the first interval starts at zero and the last ends at
    ! infinity. These must be set with define_sw_albedo_intervals and
    ! define_lw_emiss_intervals.
    real(jprb) :: sw_albedo_wavelength_bound(NMaxAlbedoIntervals-1) = -1.0_jprb
    real(jprb) :: lw_emiss_wavelength_bound( NMaxAlbedoIntervals-1) = -1.0_jprb

    ! The index to the surface albedo/emissivity intervals for each of
    ! the wavelength bounds specified in sw_albedo_wavelength_bound
    ! and lw_emiss_wavelength_bound
    integer :: i_sw_albedo_index(NMaxAlbedoIntervals) = 0
    integer :: i_lw_emiss_index (NMaxAlbedoIntervals) = 0

    ! Do we compute longwave and/or shortwave radiation?
    logical :: do_lw = .true.
    logical :: do_sw = .true.

    ! Do we compute clear-sky fluxes and/or solar direct fluxes?
    logical :: do_clear = .true.
    logical :: do_sw_direct = .true.

    ! Do we include 3D effects?
    logical :: do_3d_effects = .true.

    character(len=511) :: cloud_type_name(NMaxCloudTypes) = ["","","","","","","","","","","",""]
! &
!         &   = ["mie_droplet                   ", &
!         &      "baum-general-habit-mixture_ice"]

    ! Spectral averaging method to use with generalized cloud optics;
    ! see Edwards & Slingo (1996) for definition.  Experimentation
    ! with ecRad suggests that "thick" averaging is more accurate for
    ! both liquid and ice clouds.
    logical :: use_thick_cloud_spectral_averaging(NMaxCloudTypes) &
         &  = [.true.,.true.,.true.,.true.,.true.,.true., &
         &     .true.,.true.,.true.,.true.,.true.,.true.]

    ! To what extent do we include "entrapment" effects in the
    ! SPARTACUS solver? This essentially means that in a situation
    ! like this
    !
    ! 000111
    ! 222222
    !
    ! Radiation downwelling from region 1 may be reflected back into
    ! region 0 due to some degree of homogenization of the radiation
    ! in region 2.  Hogan and Shonk (2013) referred to this as
    ! "anomalous horizontal transport" for a 1D model, although for 3D
    ! calculations it is desirable to include at least some of it. The
    ! options are described by the IEntrapment* parameters above.
    integer :: i_3d_sw_entrapment = IEntrapmentExplicit

    ! In the longwave, the equivalent process it either "on" (like
    ! maximum entrapment) or "off" (like zero entrapment):
    logical :: do_3d_lw_multilayer_effects = .false.

    ! Do we account for the effective emissivity of the side of
    ! clouds?
    logical :: do_lw_side_emissivity = .true.

    ! The 3D transfer rate "X" is such that if transport out of a
    ! region was the only process occurring then by the base of a
    ! layer only exp(-X) of the original flux would remain in that
    ! region. The transfer rate computed geometrically can be very
    ! high for the clear-sky regions in layers with high cloud
    ! fraction.  For stability reasons it is necessary to provide a
    ! maximum possible 3D transfer rate.
    real(jprb) :: max_3d_transfer_rate = merge (2.0_jprb, 10.0_jprb, jprd /= jprb)
    ! It has also sometimes been found necessary to set a minimum
    ! cloud effective size for stability (metres)
    real(jprb) :: min_cloud_effective_size = 500.0_jprb
! from IFS, previously 100. 500 seems necessary for stability

    ! Given a horizontal migration distance, there is still
    ! uncertainty about how much entrapment occurs associated with how
    ! one assumes cloud boundaries line up in adjacent layers. This
    ! factor can be varied between 0.0 (the boundaries line up to the
    ! greatest extent possible given the overlap parameter) and 1.0
    ! (the boundaries line up to the minimum extent possible). In the
    ! Hogan et al. entrapment paper it is referred to as the overhang
    ! factor zeta, and a value of 0 matches the Monte Carlo
    ! calculations best.
    real(jprb) :: overhang_factor = 0.0_jprb

    ! By default, the Meador & Weaver (1980) expressions are used
    ! instead of the matrix exponential whenever 3D effects can be
    ! neglected (e.g. cloud-free layers or clouds with infinitely
    ! large effective cloud size), but setting the following to true
    ! uses the matrix exponential everywhere, enabling the two
    ! methods to be compared. Note that Meador & Weaver will still be
    ! used for very optically thick g points where the matrix
    ! exponential can produce incorrect results.
    logical :: use_expm_everywhere = .false.

    ! Aerosol descriptors: aerosol_type_mapping must be of length
    ! n_aerosol_types, and contains 0 if that type is to be ignored,
    ! positive numbers to map on to the indices of hydrophobic
    ! aerosols in the aerosol optics configuration file, and negative
    ! numbers to map on to (the negative of) the indices of
    ! hydrophilic aerosols in the configuration file.
    logical :: use_aerosols = .false.
    integer :: n_aerosol_types = 0
    integer :: i_aerosol_type_map(NMaxAerosolTypes)

    ! Save the gas and cloud optical properties for each g point in
    ! "radiative_properties.nc"?
    logical :: do_save_radiative_properties = .false.

    ! Save the flux profiles in each band?
    logical :: do_save_spectral_flux = .false.

    ! Save the surface downwelling shortwave fluxes in each band?
    logical :: do_surface_sw_spectral_flux = .true.

    ! Save the TOA fluxes in each band?
    logical :: do_toa_spectral_flux = .false.

    ! Compute the longwave derivatives needed to apply the approximate
    ! radiation updates of Hogan and Bozzo (2015)
    logical :: do_lw_derivatives = .false.

    ! Save the flux profiles in each g-point (overrides
    ! do_save_spectral_flux if TRUE)?
    logical :: do_save_gpoint_flux = .false.

    ! In the IFS environment, setting up RRTM has already been done
    ! so not needed to do it again
    logical :: do_setup_ifsrrtm = .true.

    ! In the IFS environment the old scheme has a bug in the Fu
    ! longwave ice optics whereby the single scattering albedo is one
    ! minus what it should be.  Unfortunately fixing it makes
    ! forecasts worse. Setting the following to true reproduces the
    ! bug.
    logical :: do_fu_lw_ice_optics_bug = .false.

    ! Control verbosity: 0=none (no output to standard output; write
    ! to standard error only if an error occurs), 1=warning, 2=info,
    ! 3=progress, 4=detailed, 5=debug.  Separate settings for the
    ! setup of the scheme and the execution of it.
    integer :: iverbosesetup = 2
    integer :: iverbose = 1

    ! Are we doing radiative transfer in complex surface canopies
    ! (streets/vegetation), in which case tailored downward fluxes are
    ! needed at the top of the canopy?
    logical :: do_canopy_fluxes_sw = .false.
    logical :: do_canopy_fluxes_lw = .false.
    ! If so, do we use the full spectrum as in the atmosphere, or just
    ! the reduced spectrum in which the shortwave albedo and longwave
    ! emissivity are provided?
    logical :: use_canopy_full_spectrum_sw = .false.
    logical :: use_canopy_full_spectrum_lw = .false.
    ! Do we treat gas radiative transfer in streets/vegetation?
    logical :: do_canopy_gases_sw = .false.
    logical :: do_canopy_gases_lw = .false.

    ! Optics file names for overriding the ones generated from the
    ! other options. If these remain empty then the generated names
    ! will be used (see the "consolidate_config" routine below). If
    ! the user assigns one of these and it starts with a '/' character
    ! then that will be used instead. If the user assigns one and it
    ! doesn't start with a '/' character then it will be prepended by
    ! the contents of directory_name.
    character(len=511) :: ice_optics_override_file_name     = ''
    character(len=511) :: liq_optics_override_file_name     = ''
    character(len=511) :: aerosol_optics_override_file_name = ''
    character(len=511) :: gas_optics_sw_override_file_name  = ''
    character(len=511) :: gas_optics_lw_override_file_name  = ''

    ! Optionally override the default file describing variations in
    ! the spectral solar irradiance associated with the solar cycle
    character(len=511) :: ssi_override_file_name = ''

    ! Do we use the solar spectral irradiance file to update the solar
    ! irradiance in each g point? Only possible if
    ! use_spectral_solar_cycle==true.
    logical :: use_updated_solar_spectrum = .false.

    ! Optionally override the look-up table file for the cloud-water
    ! PDF used by the McICA solver
    character(len=511) :: cloud_pdf_override_file_name = ''

    ! Do we compute cloud, aerosol and surface optical properties per
    ! g point?  Not available with RRTMG gas optics model.
    logical :: do_cloud_aerosol_per_sw_g_point = .true.
    logical :: do_cloud_aerosol_per_lw_g_point = .true.

    ! Do we weight the mapping from surface emissivity/albedo to
    ! g-point/band weighting by a reference Planck function (more
    ! accurate) or weighting each wavenumber equally (less accurate
    ! but consistent with IFS Cycle 48r1 and earlier)?
    logical :: do_weighted_surface_mapping = .true.

    ! COMPUTED PARAMETERS

    ! Users of this library should not edit these parameters directly;
    ! they are set by the "consolidate" routine

    ! Has "consolidate" been called?
    logical :: is_consolidated = .false.

    ! Fraction of each g point in each wavenumber interval,
    ! dimensioned (n_wav_frac_[l|s]w, n_g_[l|s]w)
    real(jprb), allocatable, dimension(:,:) :: g_frac_sw, g_frac_lw

    ! If the nearest surface albedo/emissivity interval is to be used
    ! for each SW/LW band then the following arrays will be allocated
    ! to the length of the number of bands and contain the index to
    ! the relevant interval
    integer, allocatable, dimension(:) :: i_albedo_from_band_sw
    integer, allocatable, dimension(:) :: i_emiss_from_band_lw

    ! ...alternatively, this matrix dimensioned
    ! (n_albedo_intervals,n_bands_sw) providing the weights needed for
    ! computing the albedo in each ecRad band from the albedo in each
    ! native albedo band - see radiation_single_level.F90
    real(jprb), allocatable, dimension(:,:) :: sw_albedo_weights
    ! ...and similarly in the longwave, dimensioned
    ! (n_emiss_intervals,n_bands_lw)
    real(jprb), allocatable, dimension(:,:) :: lw_emiss_weights

    ! Arrays of length the number of g-points that convert from
    ! g-point to the band index
    integer, allocatable, dimension(:) :: i_band_from_g_lw
    integer, allocatable, dimension(:) :: i_band_from_g_sw

    ! We allow for the possibility for g-points to be ordered in terms
    ! of likely absorption (weakest to strongest) across the shortwave
    ! or longwave spectrum, in order that in SPARTACUS we select only
    ! the first n g-points that will not have too large an absorption,
    ! and therefore matrix exponentials that are both finite and not
    ! too expensive to compute.  The following two arrays map the
    ! reordered g-points to the original ones.
    integer, allocatable, dimension(:) :: i_g_from_reordered_g_lw
    integer, allocatable, dimension(:) :: i_g_from_reordered_g_sw

    ! The following map the reordered g-points to the bands
    integer, allocatable, dimension(:) :: i_band_from_reordered_g_lw
    integer, allocatable, dimension(:) :: i_band_from_reordered_g_sw

    ! The following map the reordered g-points to the spectral
    ! information being saved: if do_save_gpoint_flux==TRUE then this
    ! will map on to the original g points, but if only
    ! do_save_spectral_flux==TRUE then this will map on to the bands
    integer, pointer, dimension(:) :: i_spec_from_reordered_g_lw
    integer, pointer, dimension(:) :: i_spec_from_reordered_g_sw

    ! Number of spectral intervals used for the canopy radiative
    ! transfer calculation; they are either equal to
    ! n_albedo_intervals/n_emiss_intervals or n_g_sw/n_g_lw
    integer :: n_canopy_bands_sw = 1
    integer :: n_canopy_bands_lw = 1

    ! Data structures containing gas optics description in the case of
    ! ecCKD
    type(ckd_model_type)         :: gas_optics_sw, gas_optics_lw

    ! ----- RRTMGP settings ------
    ! Data structures containing gas optics description in the case of
    ! RRTMGP
    type(ty_gas_optics_rrtmgp)   :: k_dist_sw, k_dist_lw

    ! Data structures containing gas optics neural network models in the case of
    ! RRTMGP-NN
    type(rrtmgp_network_type), dimension(3) :: rrtmgp_neural_nets
    !! 1:SW-tau, 2:SW-ray, 3:LW-tau, 4:LW-planck   OR
    ! 1:SW-tau, 2:SW-ray, 3:LW-both

    ! Array of strings describing which RRTMGP gases are available in the host model
    ! in the case of RRTMGP / RRTMGP-NN
    character*32, dimension(:), allocatable :: rrtmgp_gas_names
    ! ----- RRTMGP settings ------

    ! Data structure containing cloud scattering data
    type(cloud_optics_type)      :: cloud_optics

    ! Number of general cloud types, default liquid and ice
    integer :: n_cloud_types = 2

    ! List of data structures (one per cloud type) containing cloud
    ! scattering data
    type(general_cloud_optics_type), allocatable :: cloud_optics_sw(:)
    type(general_cloud_optics_type), allocatable :: cloud_optics_lw(:)

    ! Data structure containing aerosol scattering data
    type(aerosol_optics_type)    :: aerosol_optics

    ! Object for sampling from a gamma or lognormal distribution
    type(pdf_sampler_type)       :: pdf_sampler

    ! Optics file names
    character(len=511) :: ice_optics_file_name, &
         &                liq_optics_file_name, &
         &                aerosol_optics_file_name, &
         &                gas_optics_sw_file_name, &
         &                gas_optics_lw_file_name, &
         &                rrtmgp_gas_optics_file_name_lw, &
         &                rrtmgp_gas_optics_file_name_sw, &
         &                rrtmgp_neural_net_sw_tau, &
         &                rrtmgp_neural_net_sw_ray, &
         &                rrtmgp_neural_net_lw


    ! Solar spectral irradiance file name
    character(len=511) :: ssi_file_name

    ! McICA PDF look-up table file name
    character(len=511) :: cloud_pdf_file_name

    ! Number of gpoints and bands in the shortwave and longwave - set
    ! to zero as will be set properly later
    integer :: n_g_sw = 0, n_g_lw = 0
    integer :: n_bands_sw = 0, n_bands_lw = 0

    ! Number of spectral points to save (equal either to the number of
    ! g points or the number of bands
    integer :: n_spec_sw = 0, n_spec_lw = 0

    ! Number of wavenumber intervals used to describe the mapping from
    ! g-points to wavenumber space
    integer :: n_wav_frac_sw = 0, n_wav_frac_lw = 0

    ! Dimensions to store variables that are only needed if longwave
    ! scattering is included. "n_g_lw_if_scattering" is equal to
    ! "n_g_lw" if aerosols are allowed to scatter in the longwave,
    ! and zero otherwise. "n_bands_lw_if_scattering" is equal to
    ! "n_bands_lw" if clouds are allowed to scatter in the longwave,
    ! and zero otherwise.
    integer :: n_g_lw_if_scattering = 0, n_bands_lw_if_scattering = 0

    ! Treat clouds as horizontally homogeneous within the gribox
    logical :: is_homogeneous = .false.

    ! If the solvers are both "Cloudless" then we don't need to do any
    ! cloud processing
    logical :: do_clouds = .true.

   contains
     
     
     
     
     
     
     
     
     
     
     

    
    
    
    

  procedure :: read_GPU => read_config_from_namelist_GPU

  procedure :: consolidate_GPU => consolidate_config_GPU

  procedure :: set_GPU  => set_config_GPU

  procedure :: print_GPU => print_config_GPU

  procedure :: get_sw_weights_GPU

  procedure :: get_sw_mapping_GPU

  procedure :: define_sw_albedo_intervals_GPU

  procedure :: define_lw_emiss_intervals_GPU

  procedure :: set_aerosol_wavelength_mono_GPU

  procedure :: consolidate_sw_albedo_intervals_GPU

  procedure :: consolidate_lw_emiss_intervals_GPU

  procedure :: create_device_GPU => radiation_config_create_device_GPU

  procedure :: update_host_GPU => radiation_config_update_host_GPU

  procedure :: update_device_GPU => radiation_config_update_device_GPU

  procedure :: delete_device_GPU => radiation_config_delete_device_GPU

  procedure :: read_CPU => read_config_from_namelist_CPU

  procedure :: consolidate_CPU => consolidate_config_CPU

  procedure :: set_CPU  => set_config_CPU

  procedure :: print_CPU => print_config_CPU

  procedure :: get_sw_weights_CPU

  procedure :: get_sw_mapping_CPU

  procedure :: define_sw_albedo_intervals_CPU

  procedure :: define_lw_emiss_intervals_CPU

  procedure :: set_aerosol_wavelength_mono_CPU

  procedure :: consolidate_sw_albedo_intervals_CPU

  procedure :: consolidate_lw_emiss_intervals_CPU

  end type config_type

!  procedure, private :: print_logical, print_real, print_int

contains


  !---------------------------------------------------------------------
  ! This subroutine reads configuration data from a namelist file, and
  ! anything that is not in the namelists will be set to default
  ! values. If optional output argument "is_success" is present, then
  ! on error (e.g. missing file) it will be set to .false.; if this
  ! argument is missing then on error the program will be aborted. You
  ! may either specify the file_name or the unit of an open file to
  ! read, but not both.
  


  !---------------------------------------------------------------------
  ! This routine is called by radiation_interface:setup_radiation and
  ! it converts the user specified options into some more specific
  ! data such as data file names
  


  !---------------------------------------------------------------------
  ! This subroutine sets members of the configuration object via
  ! optional arguments, and any member not specified is left
  ! untouched. Therefore, this should be called after taking data from
  ! the namelist.
  


  !---------------------------------------------------------------------
  ! Print configuration information to standard output
  


  !---------------------------------------------------------------------
  ! In order to estimate UV and photosynthetically active radiation,
  ! we need weighted sum of fluxes considering wavelength range
  ! required.  This routine returns information for how to correctly
  ! weight output spectral fluxes for a range of input wavelengths.
  ! Note that this is approximate; internally it may be assumed that
  ! the energy is uniformly distributed in wavenumber space, for
  ! example.  If the character string "weighting_name" is present, and
  ! iverbose>=2, then information on the weighting will be provided on
  ! nulout.
  


  !---------------------------------------------------------------------
  ! As get_sw_weights but suitable for a larger number of spectral
  ! diagnostics at once: a set of monotonically increasing wavelength
  ! bounds are provided (m), and a mapping matrix is allocated and
  ! returned such that y=matmul(mapping,x), where x is a set of
  ! band-wise fluxes after calling ecRad, e.g. flux%sw_dn_surf_band,
  ! and y is the resulting fluxes in each of the wavenumber
  ! intervals. If the character string "weighting_name" is present,
  ! and iverbose>=2, then information on the weighting will be
  ! provided on nulout.
  


  !---------------------------------------------------------------------
  ! The input shortwave surface albedo coming in is likely to be in
  ! different spectral intervals to the gas model in the radiation
  ! scheme. We assume that the input albedo is defined within
  ! "ninterval" spectral intervals covering the wavelength range 0 to
  ! infinity, but allow for the possibility that two intervals may be
  ! indexed back to the same albedo band.
  


  !---------------------------------------------------------------------
  ! As define_sw_albedo_intervals but for longwave emissivity
  


  !---------------------------------------------------------------------
  ! Set the wavelengths (m) at which monochromatic aerosol properties
  ! are required. This routine must be called before consolidation of
  ! settings.
  


  !---------------------------------------------------------------------
  ! Consolidate the surface shortwave albedo intervals with the
  ! band/g-point intervals
  


  !---------------------------------------------------------------------
  ! Consolidate the surface longwave emissivity intervals with the
  ! band/g-point intervals
  


  !---------------------------------------------------------------------
  ! Return the 0-based index for str in enum_str, or abort if it is
  ! not found
  


  !---------------------------------------------------------------------
  ! Print one line of information: logical
  


  !---------------------------------------------------------------------
  ! Print one line of information: integer
  


  !---------------------------------------------------------------------
  ! Print one line of information: real
  


  !---------------------------------------------------------------------
  ! Print one line of information: enum
  


  

  

  

  

  subroutine read_config_from_namelist_GPU(this, file_name, unit, is_success, lacc)
use yomhook
use radiation_io
class(config_type), intent(inout)         :: this
character(*),       intent(in),  optional :: file_name
integer,            intent(in),  optional :: unit
logical,            intent(out), optional :: is_success
 







































 







 


logical, intent (in) :: lacc

























































 
 
 
 
 
 
 





































































































































end subroutine read_config_from_namelist_GPU

  subroutine consolidate_config_GPU(this, lacc)
use parkind1
use yomhook
use radiation_io
class(config_type), intent(inout)         :: this


logical, intent (in) :: lacc












































end subroutine consolidate_config_GPU

  subroutine set_config_GPU(config, directory_name, &
&  do_lw, do_sw, &
&  do_lw_aerosol_scattering, do_lw_cloud_scattering, &
&  do_sw_direct, lacc)
class(config_type), intent(inout):: config
character(len=*), intent(in), optional  :: directory_name
logical, intent(in), optional           :: do_lw, do_sw
logical, intent(in), optional           :: do_lw_aerosol_scattering
logical, intent(in), optional           :: do_lw_cloud_scattering
logical, intent(in), optional           :: do_sw_direct
logical, intent (in) :: lacc






end subroutine set_config_GPU

  subroutine print_config_GPU(this, iverbose, lacc)
use radiation_io
class(config_type), intent(in) :: this
integer, optional,  intent(in) :: iverbose

logical, intent (in) :: lacc


end subroutine print_config_GPU

  subroutine get_sw_weights_GPU(this, wavelength1, wavelength2, &
&                    nweights, iband, weight, weighting_name, lacc)
use parkind1
use radiation_io
use radiation_spectral_definition
class(config_type), intent(in) :: this

real(jprb), intent(in) :: wavelength1, wavelength2

integer,    intent(out)   :: nweights



integer,    intent(out)   :: iband(:)
real(jprb), intent(out)   :: weight(:)
character(len=*), optional, intent(in) :: weighting_name


 
 

logical, intent (in) :: lacc 










end subroutine get_sw_weights_GPU

  subroutine get_sw_mapping_GPU(this, wavelength_bound, mapping, weighting_name, lacc)
use parkind1
use radiation_io
use radiation_spectral_definition
class(config_type), intent(in) :: this

real(jprb), intent(in)  :: wavelength_bound(:)
real(jprb), intent(out), allocatable :: mapping(:,:)
character(len=*), optional, intent(in) :: weighting_name




logical, intent (in) :: lacc  











end subroutine get_sw_mapping_GPU

  subroutine define_sw_albedo_intervals_GPU(this, ninterval, wavelength_bound, &
&                                i_intervals, do_nearest, lacc)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this

integer,              intent(in)    :: ninterval



real(jprb),           intent(in)    :: wavelength_bound(ninterval-1)

integer,              intent(in)    :: i_intervals(ninterval)
logical,    optional, intent(in)    :: do_nearest
logical, intent (in) :: lacc











end subroutine define_sw_albedo_intervals_GPU

  subroutine define_lw_emiss_intervals_GPU(this, ninterval, wavelength_bound, &
&                                i_intervals, do_nearest, lacc)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this

integer,              intent(in)    :: ninterval



real(jprb),           intent(in)    :: wavelength_bound(ninterval-1)

integer,              intent(in)    :: i_intervals(ninterval)
logical,    optional, intent(in)    :: do_nearest
logical, intent (in) :: lacc







end subroutine define_lw_emiss_intervals_GPU

  subroutine set_aerosol_wavelength_mono_GPU(this, wavelength_mono, lacc)
use radiation_io
class(config_type), intent(inout) :: this
real(jprb),         intent(in)    :: wavelength_mono(:)
logical, intent (in) :: lacc




end subroutine set_aerosol_wavelength_mono_GPU

  subroutine consolidate_sw_albedo_intervals_GPU(this, lacc)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this

logical, intent (in) :: lacc








end subroutine consolidate_sw_albedo_intervals_GPU

  subroutine consolidate_lw_emiss_intervals_GPU(this, lacc)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this

logical, intent (in) :: lacc








end subroutine consolidate_lw_emiss_intervals_GPU

  subroutine get_enum_code_GPU(str, enum_str, var_name, icode, lacc)
use radiation_io
character(len=*), intent(in)  :: str
character(len=*), intent(in)  :: enum_str(0:)
character(len=*), intent(in)  :: var_name
integer,          intent(out) :: icode


logical, intent (in) :: lacc



end subroutine get_enum_code_GPU

  subroutine print_logical_GPU(message_str, name, val, lacc)
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: name
logical,            intent(in) :: val


logical, intent (in) :: lacc



end subroutine print_logical_GPU

  subroutine print_integer_GPU(message_str, name, val, lacc)
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: name
integer,            intent(in) :: val

logical, intent (in) :: lacc


end subroutine print_integer_GPU

  subroutine print_real_GPU(message_str, name, val, lacc)
use parkind1
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: name
real(jprb),         intent(in) :: val

logical, intent (in) :: lacc


end subroutine print_real_GPU

  subroutine print_enum_GPU(message_str, enum_str, name, val, lacc)
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: enum_str(0:)
character(len=*),   intent(in) :: name
integer,            intent(in) :: val

logical, intent (in) :: lacc


end subroutine print_enum_GPU

  subroutine radiation_config_create_device_GPU(this, lacc)
class(config_type), intent(inout) :: this
logical, intent (in) :: lacc




















end subroutine radiation_config_create_device_GPU

  subroutine radiation_config_update_host_GPU(this, lacc)
class(config_type), intent(inout) :: this
logical, intent (in) :: lacc




















end subroutine radiation_config_update_host_GPU

  subroutine radiation_config_update_device_GPU(this, lacc)
class(config_type), intent(inout) :: this
logical, intent (in) :: lacc




















end subroutine radiation_config_update_device_GPU

  subroutine radiation_config_delete_device_GPU(this, lacc)
class(config_type), intent(inout) :: this
logical, intent (in) :: lacc




















end subroutine radiation_config_delete_device_GPU

  subroutine read_config_from_namelist_CPU(this, file_name, unit, is_success)
use yomhook
use radiation_io
class(config_type), intent(inout)         :: this
character(*),       intent(in),  optional :: file_name
integer,            intent(in),  optional :: unit
logical,            intent(out), optional :: is_success
 







































 







 



























































 
 
 
 
 
 
 





































































































































end subroutine read_config_from_namelist_CPU

  subroutine consolidate_config_CPU(this)
use parkind1
use yomhook
use radiation_io
class(config_type), intent(inout)         :: this














































end subroutine consolidate_config_CPU

  subroutine set_config_CPU(config, directory_name, &
&  do_lw, do_sw, &
&  do_lw_aerosol_scattering, do_lw_cloud_scattering, &
&  do_sw_direct)
class(config_type), intent(inout):: config
character(len=*), intent(in), optional  :: directory_name
logical, intent(in), optional           :: do_lw, do_sw
logical, intent(in), optional           :: do_lw_aerosol_scattering
logical, intent(in), optional           :: do_lw_cloud_scattering
logical, intent(in), optional           :: do_sw_direct






end subroutine set_config_CPU

  subroutine print_config_CPU(this, iverbose)
use radiation_io
class(config_type), intent(in) :: this
integer, optional,  intent(in) :: iverbose



end subroutine print_config_CPU

  subroutine get_sw_weights_CPU(this, wavelength1, wavelength2, &
&                    nweights, iband, weight, weighting_name)
use parkind1
use radiation_io
use radiation_spectral_definition
class(config_type), intent(in) :: this

real(jprb), intent(in) :: wavelength1, wavelength2

integer,    intent(out)   :: nweights



integer,    intent(out)   :: iband(:)
real(jprb), intent(out)   :: weight(:)
character(len=*), optional, intent(in) :: weighting_name


 
 
 










end subroutine get_sw_weights_CPU

  subroutine get_sw_mapping_CPU(this, wavelength_bound, mapping, weighting_name)
use parkind1
use radiation_io
use radiation_spectral_definition
class(config_type), intent(in) :: this

real(jprb), intent(in)  :: wavelength_bound(:)
real(jprb), intent(out), allocatable :: mapping(:,:)
character(len=*), optional, intent(in) :: weighting_name



  











end subroutine get_sw_mapping_CPU

  subroutine define_sw_albedo_intervals_CPU(this, ninterval, wavelength_bound, &
&                                i_intervals, do_nearest)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this

integer,              intent(in)    :: ninterval



real(jprb),           intent(in)    :: wavelength_bound(ninterval-1)

integer,              intent(in)    :: i_intervals(ninterval)
logical,    optional, intent(in)    :: do_nearest











end subroutine define_sw_albedo_intervals_CPU

  subroutine define_lw_emiss_intervals_CPU(this, ninterval, wavelength_bound, &
&                                i_intervals, do_nearest)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this

integer,              intent(in)    :: ninterval



real(jprb),           intent(in)    :: wavelength_bound(ninterval-1)

integer,              intent(in)    :: i_intervals(ninterval)
logical,    optional, intent(in)    :: do_nearest







end subroutine define_lw_emiss_intervals_CPU

  subroutine set_aerosol_wavelength_mono_CPU(this, wavelength_mono)
use radiation_io
class(config_type), intent(inout) :: this
real(jprb),         intent(in)    :: wavelength_mono(:)




end subroutine set_aerosol_wavelength_mono_CPU

  subroutine consolidate_sw_albedo_intervals_CPU(this)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this









end subroutine consolidate_sw_albedo_intervals_CPU

  subroutine consolidate_lw_emiss_intervals_CPU(this)
use radiation_io
use radiation_spectral_definition
class(config_type),   intent(inout) :: this









end subroutine consolidate_lw_emiss_intervals_CPU

  subroutine get_enum_code_CPU(str, enum_str, var_name, icode)
use radiation_io
character(len=*), intent(in)  :: str
character(len=*), intent(in)  :: enum_str(0:)
character(len=*), intent(in)  :: var_name
integer,          intent(out) :: icode





end subroutine get_enum_code_CPU

  subroutine print_logical_CPU(message_str, name, val)
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: name
logical,            intent(in) :: val





end subroutine print_logical_CPU

  subroutine print_integer_CPU(message_str, name, val)
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: name
integer,            intent(in) :: val



end subroutine print_integer_CPU

  subroutine print_real_CPU(message_str, name, val)
use parkind1
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: name
real(jprb),         intent(in) :: val



end subroutine print_real_CPU

  subroutine print_enum_CPU(message_str, enum_str, name, val)
use radiation_io
character(len=*),   intent(in) :: message_str
character(len=*),   intent(in) :: enum_str(0:)
character(len=*),   intent(in) :: name
integer,            intent(in) :: val



end subroutine print_enum_CPU

end module radiation_config

