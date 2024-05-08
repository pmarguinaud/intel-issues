! -- Field definitions - prototype AJG 24/10/2012 --
!
! This module contains all the field-specific information, e.g.
! names, dimensions, ID number, attributes, field classes. 
! There are a number of tedious but undemanding steps 
! required to add a new field:
!
!   (F1) Create a pointer in the field_access type
!   (F2) Add a field ID name...
!   (F3) ... and number (in sequence)
!   (F4) Put the field in any appropriate "classes" (field_class.F90)
!   (F5) Fill in other details (e.g. name as a string, dimensions, GRIB code?)
!   (F6) Set up the mapping from a field ID to the named pointer in field_access

! If a new model field is required, this should be the only place 
! that needs to be changed. Hopefully the tedious stuff can be eliminated
! by auto-generating this file in the future (e.g. from a central database of 
! model fields in xml)

! Problem 1: what about multi-dimension GOM fields, i.e. GHG, GRG, AERO, PHYS?
! Ideally these should be replaced by individual fields, which can be grouped in a "class"

! Problem 2: what about GOM-only "fake" fields like ul, vl, tl, ql? Ideally, GOM should do
! vertical interpolation or sub-sampling to a smaller number of model levels, thus removing the
! need for ul, vl, etc. (but it should be possible for GOM to create and use fields for itself
! if necessary).

! Problem 3: 2D and 3D canari model error fields used by GOM

module field_definitions

use field_definitions_base, only: set_fvar, type_fvar, field_access_base, field_metadata_base, &
  & jp_name_max_len, jp_necv_2d_max, jp_necv_3d_max, jp_name_max_len, jp_comments_max_len

use parkind1, only: jpim, jprb
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RV, RCPV, RCW, RCS
USE YOMLUN   , ONLY : NULOUT
USE YOM_GRIB_CODES

implicit none
public

! This type is to allow "named" and "shaped" user access to fields. (F1)
type, extends(field_access_base) :: field_access

  ! GMV
  real(kind=jprb), pointer :: u(:,:) => null()
  real(kind=jprb), pointer :: v(:,:) => null()
  real(kind=jprb), pointer :: div(:,:) => null()
  real(kind=jprb), pointer :: vor(:,:) => null()
  real(kind=jprb), pointer :: t(:,:) => null()
  real(kind=jprb), pointer :: pd(:,:) => null()
  real(kind=jprb), pointer :: vd(:,:) => null()   
  real(kind=jprb), pointer :: nhx(:,:) => null()
  real(kind=jprb), pointer :: edot(:,:) => null()
  real(kind=jprb), pointer :: vcasrs(:,:) => null() ! ky: obsolete field (deep-layer models pruned)
  real(kind=jprb), pointer :: rdlasr(:,:) => null() ! ky: obsolete field (deep-layer models pruned)
  real(kind=jprb), pointer :: sgrtl(:,:) => null()
  real(kind=jprb), pointer :: sgrtm(:,:) => null()
  real(kind=jprb), pointer :: unl(:,:) => null()
  real(kind=jprb), pointer :: unl_si(:,:) => null()
  real(kind=jprb), pointer :: vnl(:,:) => null()
  real(kind=jprb), pointer :: vnl_si(:,:) => null()
  real(kind=jprb), pointer :: tnl(:,:) => null()
  real(kind=jprb), pointer :: tnl_si(:,:) => null()
  real(kind=jprb), pointer :: spnl(:,:)  => null()
  real(kind=jprb), pointer :: spnl_si(:,:)  => null()
  real(kind=jprb), pointer :: vwvnl(:,:) => null()
  real(kind=jprb), pointer :: vdnl_si(:,:) => null()
  real(kind=jprb), pointer :: pdnl(:,:) => null()
  real(kind=jprb), pointer :: pdnl_si(:,:) => null()
  real(kind=jprb), pointer :: gw(:,:) => null()
  real(kind=jprb), pointer :: nhy(:,:) => null()
  real(kind=jprb), pointer :: curhs(:,:) => null()
  real(kind=jprb), pointer :: cvrhs(:,:) => null()
  real(kind=jprb), pointer :: ctrhs(:,:) => null()
  real(kind=jprb), pointer :: csprhs(:) => null()
  real(kind=jprb), pointer :: cspdrhs(:,:) => null()
  real(kind=jprb), pointer :: csvdrhs(:,:) => null()
  real(kind=jprb), pointer :: cunl(:,:) => null()
  real(kind=jprb), pointer :: cvnl(:,:) => null()
  real(kind=jprb), pointer :: ctnl(:,:) => null()
  real(kind=jprb), pointer :: cspnl(:,:) => null()
  real(kind=jprb), pointer :: cvwvnl(:,:) => null()
  real(kind=jprb), pointer :: cspdnl(:,:) => null()
  real(kind=jprb), pointer :: cpdnl(:,:) => null()
  real(kind=jprb), pointer :: cupt(:,:) => null()
  real(kind=jprb), pointer :: cvpt(:,:) => null()
  real(kind=jprb), pointer :: ctpt(:,:) => null()
  real(kind=jprb), pointer :: cpdpt(:,:) => null()
  real(kind=jprb), pointer :: cvdpt(:,:) => null()
  real(kind=jprb), pointer :: dphi(:,:) => null()
  real(kind=jprb), pointer :: nhxnl(:,:) => null()
  real(kind=jprb), pointer :: cnhxnl(:,:) => null()
  ! GFL
  real(kind=jprb), pointer :: q(:,:) => null()
  real(kind=jprb), pointer :: o3(:,:) => null()
  real(kind=jprb), pointer :: l(:,:) => null()
  real(kind=jprb), pointer :: i(:,:) => null()
  real(kind=jprb), pointer :: a(:,:) => null()
  real(kind=jprb), pointer :: s(:,:) => null()
  real(kind=jprb), pointer :: r(:,:) => null()
  real(kind=jprb), pointer :: g(:,:) => null()
  real(kind=jprb), pointer :: h(:,:) => null()
  real(kind=jprb), pointer :: lconv(:,:) => null()
  real(kind=jprb), pointer :: iconv(:,:) => null()
  real(kind=jprb), pointer :: rconv(:,:) => null()
  real(kind=jprb), pointer :: sconv(:,:) => null()
  real(kind=jprb), pointer :: tke(:,:) => null()
  real(kind=jprb), pointer :: tte(:,:) => null()
  real(kind=jprb), pointer :: efb1(:,:) => null()
  real(kind=jprb), pointer :: efb2(:,:) => null()
  real(kind=jprb), pointer :: efb3(:,:) => null()
  real(kind=jprb), pointer :: mxl(:,:) => null()
  real(kind=jprb), pointer :: src(:,:) => null()
  real(kind=jprb), pointer :: lmf(:,:) => null()
  real(kind=jprb), pointer :: imf(:,:) => null()
  real(kind=jprb), pointer :: amf(:,:) => null()
  real(kind=jprb), pointer :: shtur(:,:) => null()
  real(kind=jprb), pointer :: fqtur(:,:) => null()
  real(kind=jprb), pointer :: fstur(:,:) => null()
  real(kind=jprb), pointer :: cvv(:,:) => null()
  real(kind=jprb), pointer :: rkth(:,:) => null()
  real(kind=jprb), pointer :: rktqv(:,:) => null()
  real(kind=jprb), pointer :: rktqc(:,:) => null()
  real(kind=jprb), pointer :: cpf(:,:) => null()
  real(kind=jprb), pointer :: uom(:,:) => null()
  real(kind=jprb), pointer :: ual(:,:) => null()
  real(kind=jprb), pointer :: dom(:,:) => null()
  real(kind=jprb), pointer :: dal(:,:) => null()
  real(kind=jprb), pointer :: uen(:,:) => null()
  real(kind=jprb), pointer :: unebh(:,:) => null()
  real(kind=jprb), pointer :: spf(:,:) => null()
  real(kind=jprb), pointer :: cvgq(:,:) => null()
  real(kind=jprb), pointer :: lrad(:,:) => null()
  real(kind=jprb), pointer :: irad(:,:) => null()
  real(kind=jprb), pointer :: rspec(:,:) => null()
  real(kind=jprb), pointer :: lrch4(:,:) => null()
  real(kind=jprb), pointer :: extra3d(:,:,:) => null()
  real(kind=jprb), pointer :: ezdiag(:,:,:) => null()
  
  ! pressure
  real(kind=jprb), pointer :: pre(:,:) => null()
  real(kind=jprb), pointer :: pref(:,:) => null()
  real(kind=jprb), pointer :: delp(:,:) => null()
  ! GMVS
  real(kind=jprb), pointer :: sp(:)  => null()
  real(kind=jprb), pointer :: csppt(:)  => null()
  real(kind=jprb), pointer :: dbbc(:)  => null()
  real(kind=jprb), pointer :: gws(:)  => null()
  real(kind=jprb), pointer :: spnl2(:)  => null()
  real(kind=jprb), pointer :: cspnl2(:)  => null()
  real(kind=jprb), pointer :: prehyds(:)  => null()
  real(kind=jprb), pointer :: extcv2d(:,:)  => null()
  real(kind=jprb), pointer :: extcv3d(:,:)  => null()

  ! For LELAM
  real(kind=jprb), pointer :: u_mean(:) => null()
  real(kind=jprb), pointer :: v_mean(:) => null()
  
  ! SLB2
  real(kind=jprb), pointer :: vvel(:,:) => null()
  ! Constants
  real(kind=jprb), pointer :: orog(:) => null()
  real(kind=jprb), pointer :: cori(:) => null()
  real(kind=jprb), pointer :: gnordl(:) => null()
  real(kind=jprb), pointer :: gnordm(:) => null()
  ! Clouds
  real(kind=jprb), pointer :: ccc(:) => null()
  real(kind=jprb), pointer :: lcc(:) => null()
  real(kind=jprb), pointer :: mcc(:) => null()
  real(kind=jprb), pointer :: hcc(:) => null()
  real(kind=jprb), pointer :: tcc(:) => null()
  ! Surface
  ! Group SB
  real(kind=jprb), pointer :: sb_t(:,:) => null()
  real(kind=jprb), pointer :: sb_q(:,:) => null()
  real(kind=jprb), pointer :: sb_tl(:,:) => null()
  ! Group SG
  real(kind=jprb), pointer :: sg_f(:,:) => null()
  real(kind=jprb), pointer :: sg_a(:,:) => null()
  real(kind=jprb), pointer :: sg_r(:,:) => null()
  real(kind=jprb), pointer :: sg_t(:,:) => null()
  real(kind=jprb), pointer :: sg_w(:,:) => null()
  ! Group SL
  real(kind=jprb), pointer :: sl_lict(:) => null()
  real(kind=jprb), pointer :: sl_lmlt(:) => null()
  real(kind=jprb), pointer :: sl_ltlt(:) => null()
  real(kind=jprb), pointer :: sl_lblt(:) => null()
  real(kind=jprb), pointer :: sl_lshf(:) => null()
  real(kind=jprb), pointer :: sl_licd(:) => null()
  real(kind=jprb), pointer :: sl_lmld(:) => null()
  ! Group RR
  real(kind=jprb), pointer :: rr_t(:) => null()
  real(kind=jprb), pointer :: rr_tir(:) => null()
  real(kind=jprb), pointer :: rr_tmw(:) => null()
  real(kind=jprb), pointer :: rr_w(:) => null()
  real(kind=jprb), pointer :: rr_fc(:) => null()
  real(kind=jprb), pointer :: rr_ic(:) => null()
  real(kind=jprb), pointer :: rr_fp1(:) => null()
  ! Group CL
  real(kind=jprb), pointer :: cl_tcls(:) => null() 
  real(kind=jprb), pointer :: cl_hucls(:) => null() 
  ! Group OM=OML: prognostic quantities for ocean mixed layer model (KPP/TKE)
  real(kind=jprb), pointer :: om_to(:,:) => null() 
  real(kind=jprb), pointer :: om_so(:,:) => null() 
  real(kind=jprb), pointer :: om_uo(:,:) => null() 
  real(kind=jprb), pointer :: om_vo(:,:) => null() 
  ! Group EP=EXTRP: extra 3-d prognostic fields
  real(kind=jprb), pointer :: ep_ep(:,:) => null()
  ! Group X2=XTRP2: extra 2-d prognostic fields
  real(kind=jprb), pointer :: x2_x2(:,:) => null() 
  ! Group CI=CANRI: 2-d prognostic fields for CANARI
  real(kind=jprb), pointer :: ci_ci(:,:) => null() 
  ! Group VF
  real(kind=jprb), pointer :: vf_z0f(:) => null() 
  real(kind=jprb), pointer :: vf_albf(:) => null() 
  real(kind=jprb), pointer :: vf_emisf(:) => null() 
  real(kind=jprb), pointer :: vf_getrl(:) => null() 
  real(kind=jprb), pointer :: vf_lsm(:) => null() 
  real(kind=jprb), pointer :: vf_veg(:) => null() 
  real(kind=jprb), pointer :: vf_vrlan(:) => null() 
  real(kind=jprb), pointer :: vf_vrldi(:) => null() 
  real(kind=jprb), pointer :: vf_sig(:) => null() 
  real(kind=jprb), pointer :: vf_albsf(:) => null() 
  real(kind=jprb), pointer :: vf_lan(:) => null() 
  real(kind=jprb), pointer :: vf_sst(:) => null() 
  real(kind=jprb), pointer :: vf_sss(:) => null() 
  real(kind=jprb), pointer :: vf_lz0h(:) => null() 
  real(kind=jprb), pointer :: vf_cvl(:) => null() 
  real(kind=jprb), pointer :: vf_co2typ(:) => null() 
  real(kind=jprb), pointer :: vf_cvh(:) => null()
  real(kind=jprb), pointer :: vf_fwet(:) => null()
  real(kind=jprb), pointer :: vf_cur(:) => null()
  real(kind=jprb), pointer :: vf_tvl(:) => null() 
  real(kind=jprb), pointer :: vf_tvh(:) => null() 
  real(kind=jprb), pointer :: vf_lail(:) => null() 
  real(kind=jprb), pointer :: vf_laih(:) => null() 
  real(kind=jprb), pointer :: vf_soty(:) => null() 
  real(kind=jprb), pointer :: vf_clk(:) => null() 
  real(kind=jprb), pointer :: vf_dl(:) => null() 
  real(kind=jprb), pointer :: vf_ci(:) => null() 
  real(kind=jprb), pointer :: vf_ucur(:) => null() 
  real(kind=jprb), pointer :: vf_vcur(:) => null()
  real(kind=jprb), pointer :: vf_z0rlf(:) => null()
  real(kind=jprb), pointer :: vf_cgpp(:) => null()
  real(kind=jprb), pointer :: vf_crec(:) => null()
  real(kind=jprb), pointer :: vf_sdfor(:) => null()
  real(kind=jprb), pointer :: vf_aluvp(:) => null() 
  real(kind=jprb), pointer :: vf_aluvd(:) => null() 
  real(kind=jprb), pointer :: vf_alnip(:) => null() 
  real(kind=jprb), pointer :: vf_alnid(:) => null() 
  real(kind=jprb), pointer :: vf_aluvi(:) => null() 
  real(kind=jprb), pointer :: vf_aluvv(:) => null() 
  real(kind=jprb), pointer :: vf_aluvg(:) => null() 
  real(kind=jprb), pointer :: vf_alnii(:) => null() 
  real(kind=jprb), pointer :: vf_alniv(:) => null() 
  real(kind=jprb), pointer :: vf_alnig(:) => null() 
  real(kind=jprb), pointer :: vf_fp1(:) => null()     ! surface orography in the 2nd part of FULLPOS-927
  real(kind=jprb), pointer :: vf_so2dd(:) => null()   ! sulphate dry dep velocity
  real(kind=jprb), pointer :: vf_dmso(:) => null()    ! oceanic DMS
  real(kind=jprb), pointer :: vf_urbf(:) => null()   ! SOA from CO
  real(kind=jprb), pointer :: vf_fca1(:) => null()    ! fraction of calcite over dust 1st bin
  real(kind=jprb), pointer :: vf_fca2(:) => null()    ! fraction of calcite over dust 2st bin
  real(kind=jprb), pointer :: vf_aerdep(:) => null()  ! dust emission potential 
  real(kind=jprb), pointer :: vf_aerlts(:) => null()  ! dust lifting threshold speed 
  real(kind=jprb), pointer :: vf_aerscc(:) => null()  ! dust soil clay content
  real(kind=jprb), pointer :: vf_dsf(:) => null()     ! dust soil clay content
  real(kind=jprb), pointer :: vf_dsz(:) => null()     ! Dust size dist modulate
  real(kind=jprb), pointer :: vf_chemflxo(:,:) => null() !  total chemistry flux (emissions + deposition)
  real(kind=jprb), pointer :: vf_chemwdflx(:,:) => null() !  wet deposition chemistry flux 
  real(kind=jprb), pointer :: vf_chemddflx(:,:) => null() !  dry deposition chemistry flux 
  real(kind=jprb), pointer :: vf_chemdv(:,:) => null()  ! chemistry deposition velocity
  real(kind=jprb), pointer :: vf_nudm(:) => null()    ! nudging mask
  real(kind=jprb), pointer :: vf_emis2d(:,:) => null()  ! 2D emission fields for composition
  real(kind=jprb), pointer :: vf_emis2daux(:,:) => null()  ! 2D emission auxiliary fields for composition
  ! Group VP=VCLIP: deep soil diagnostic fields
  real(kind=jprb), pointer :: vp_tpc(:) => null()
  real(kind=jprb), pointer :: vp_wpc(:) => null()
  ! Group VV=VCLIV: vegetation diagnostic fields
  real(kind=jprb), pointer :: vv_arg(:) => null() 
  real(kind=jprb), pointer :: vv_sab(:) => null() 
  real(kind=jprb), pointer :: vv_d2(:) => null() 
  real(kind=jprb), pointer :: vv_iveg(:) => null() 
  real(kind=jprb), pointer :: vv_rsmin(:) => null() 
  real(kind=jprb), pointer :: vv_lai(:) => null() 
  real(kind=jprb), pointer :: vv_hv(:) => null() 
  real(kind=jprb), pointer :: vv_z0h(:) => null() 
  real(kind=jprb), pointer :: vv_als(:) => null() 
  real(kind=jprb), pointer :: vv_alv(:) => null() 
  ! Group VN=VCLIN: cloudiness diagnostic predictors:
  real(kind=jprb), pointer :: vn_top(:) => null() 
  real(kind=jprb), pointer :: vn_bas(:) => null() 
  real(kind=jprb), pointer :: vn_acpr(:) => null() 
  real(kind=jprb), pointer :: vn_accpr(:) => null() 
  real(kind=jprb), pointer :: vn_accpr5(:) => null() 
  ! Group VH=VCLIH: convective cloud diagnostic fields
   real(kind=jprb), pointer :: vh_tcch(:) => null()
   real(kind=jprb), pointer :: vh_scch(:) => null()
   real(kind=jprb), pointer :: vh_bcch(:) => null()
   real(kind=jprb), pointer :: vh_pblh(:) => null()
   real(kind=jprb), pointer :: vh_spsh(:) => null()
   real(kind=jprb), pointer :: vh_qsh (:) => null()
  ! Group  VA=VCLIA: aerosol diagnostic fields
   real(kind=jprb), pointer :: va_sea(:) => null()
   real(kind=jprb), pointer :: va_lan(:) => null()
   real(kind=jprb), pointer :: va_soo(:) => null()
   real(kind=jprb), pointer :: va_des(:) => null()
   real(kind=jprb), pointer :: va_sul(:) => null()
   real(kind=jprb), pointer :: va_vol(:) => null()
  ! Group  VG=VCLIG: ice-coupler diagnostic fields
   real(kind=jprb), pointer :: vg_icfr(:) => null()
   real(kind=jprb), pointer :: vg_soup(:) => null()
   real(kind=jprb), pointer :: vg_irup(:) => null()
   real(kind=jprb), pointer :: vg_chss(:) => null()
   real(kind=jprb), pointer :: vg_evap(:) => null()
   real(kind=jprb), pointer :: vg_taux(:) => null()
   real(kind=jprb), pointer :: vg_tauy(:) => null()
  ! Group VC=VO3ABC: A,B and C (Climatological ozone profiles) diagnostic fields
   real(kind=jprb), pointer :: vc_a(:) => null()
   real(kind=jprb), pointer :: vc_b(:) => null()
   real(kind=jprb), pointer :: vc_c(:) => null()
  ! Group V2=VDIAGO2: 2-D climatological/diagnostic fields for an ocean mixed layer model (KPP)
   real(kind=jprb), pointer :: v2_ocdep(:) => null()
   real(kind=jprb), pointer :: v2_ustrc(:) => null()
   real(kind=jprb), pointer :: v2_vstrc(:) => null()
  ! * Group V3=VDIAGO3: 3-D climatological/diagnostic fields for an ocean mixed layer model (KPP):
   real(kind=jprb), pointer :: v3_difm(:,:)  => null()
   real(kind=jprb), pointer :: v3_dift(:,:)  => null()
   real(kind=jprb), pointer :: v3_difs(:,:)  => null()
   real(kind=jprb), pointer :: v3_advt(:,:)  => null()
   real(kind=jprb), pointer :: v3_advs(:,:)  => null()
   real(kind=jprb), pointer :: v3_tri0(:,:)  => null()
   real(kind=jprb), pointer :: v3_tri1(:,:)  => null()
   real(kind=jprb), pointer :: v3_swdk(:,:)  => null()
   real(kind=jprb), pointer :: v3_zo(:,:)  => null()
   real(kind=jprb), pointer :: v3_ho(:,:)  => null()
   real(kind=jprb), pointer :: v3_do(:,:)  => null()
   real(kind=jprb), pointer :: v3_ho_inv(:,:)  => null()
   real(kind=jprb), pointer :: v3_uoc(:,:)  => null()
   real(kind=jprb), pointer :: v3_voc(:,:)  => null()
   real(kind=jprb), pointer :: v3_otke(:,:)  => null()
   ! Group VD
   real(kind=jprb), pointer :: vd_lsp(:) => null() 
   real(kind=jprb), pointer :: vd_cp(:) => null() 
   real(kind=jprb), pointer :: vd_sf(:) => null() 
   real(kind=jprb), pointer :: vd_fzra(:) => null() 
   real(kind=jprb), pointer :: vd_bld(:) => null() 
   real(kind=jprb), pointer :: vd_sshf(:) => null() 
   real(kind=jprb), pointer :: vd_slhf(:) => null() 
   real(kind=jprb), pointer :: vd_nee(:) => null() 
   real(kind=jprb), pointer :: vd_gpp(:) => null() 
   real(kind=jprb), pointer :: vd_rec(:) => null() 
   real(kind=jprb), pointer :: vd_msl(:) => null() 
   real(kind=jprb), pointer :: vd_sp(:) => null() 
   real(kind=jprb), pointer :: vd_tcc(:) => null() 
   real(kind=jprb), pointer :: vd_10u(:) => null() 
   real(kind=jprb), pointer :: vd_10v(:) => null() 
   real(kind=jprb), pointer :: vd_2t(:) => null() 
   real(kind=jprb), pointer :: vd_2d(:) => null() 
   real(kind=jprb), pointer :: vd_2sh(:) => null() 
   real(kind=jprb), pointer :: vd_ssr(:) => null() 
   real(kind=jprb), pointer :: vd_str(:) => null() 
   real(kind=jprb), pointer :: vd_tsr(:) => null() 
   real(kind=jprb), pointer :: vd_ttr(:) => null() 
   real(kind=jprb), pointer :: vd_ewss(:) => null() 
   real(kind=jprb), pointer :: vd_nsss(:) => null() 
   real(kind=jprb), pointer :: vd_e(:) => null() 
   real(kind=jprb), pointer :: vd_pev(:) => null() 
   real(kind=jprb), pointer :: vd_ccc(:) => null() 
   real(kind=jprb), pointer :: vd_lcc(:) => null() 
   real(kind=jprb), pointer :: vd_mcc(:) => null() 
   real(kind=jprb), pointer :: vd_hcc(:) => null() 
   real(kind=jprb), pointer :: vd_lgws(:) => null() 
   real(kind=jprb), pointer :: vd_mgws(:) => null() 
   real(kind=jprb), pointer :: vd_gwd(:) => null() 
   real(kind=jprb), pointer :: vd_mx2t(:) => null() 
   real(kind=jprb), pointer :: vd_mn2t(:) => null() 
   real(kind=jprb), pointer :: vd_mx2t6(:) => null() 
   real(kind=jprb), pointer :: vd_mn2t6(:) => null() 
   real(kind=jprb), pointer :: vd_ro(:) => null() 
   real(kind=jprb), pointer :: vd_sro(:) => null() 
   real(kind=jprb), pointer :: vd_ssro(:) => null() 
   real(kind=jprb), pointer :: vd_alb(:) => null() 
   real(kind=jprb), pointer :: vd_iewss(:) => null() 
   real(kind=jprb), pointer :: vd_insss(:) => null() 
   real(kind=jprb), pointer :: vd_isshf(:) => null() 
   real(kind=jprb), pointer :: vd_ie(:) => null() 
   real(kind=jprb), pointer :: vd_inee(:) => null() 
   real(kind=jprb), pointer :: vd_igpp(:) => null() 
   real(kind=jprb), pointer :: vd_irec(:) => null()
   real(kind=jprb), pointer :: vd_ich4(:) => null()
   real(kind=jprb), pointer :: vd_csf(:) => null() 
   real(kind=jprb), pointer :: vd_lssf(:) => null() 
   real(kind=jprb), pointer :: vd_mxtpr(:) => null() 
   real(kind=jprb), pointer :: vd_mntpr(:) => null() 
   real(kind=jprb), pointer :: vd_mxtpr6(:) => null() 
   real(kind=jprb), pointer :: vd_mntpr6(:) => null() 
   real(kind=jprb), pointer :: vd_tpr(:) => null() 
   real(kind=jprb), pointer :: vd_lsrr(:) => null() 
   real(kind=jprb), pointer :: vd_crr(:) => null() 
   real(kind=jprb), pointer :: vd_lssfr(:) => null() 
   real(kind=jprb), pointer :: vd_csfr(:) => null() 
   real(kind=jprb), pointer :: vd_ptype(:) => null() 
   real(kind=jprb), pointer :: vd_ilspf(:) => null() 
   real(kind=jprb), pointer :: vd_z0f(:) => null() 
   real(kind=jprb), pointer :: vd_lz0h(:) => null() 
   real(kind=jprb), pointer :: vd_tcw(:) => null() 
   real(kind=jprb), pointer :: vd_tcwv(:) => null() 
   real(kind=jprb), pointer :: vd_tclw(:) => null() 
   real(kind=jprb), pointer :: vd_tciw(:) => null() 
   real(kind=jprb), pointer :: vd_tcrw(:) => null() 
   real(kind=jprb), pointer :: vd_tcsw(:) => null() 
   real(kind=jprb), pointer :: vd_tcslw(:) => null() 
   real(kind=jprb), pointer :: vd_ssrd(:) => null() 
   real(kind=jprb), pointer :: vd_strd(:) => null() 
   real(kind=jprb), pointer :: vd_ssrdc(:) => null() 
   real(kind=jprb), pointer :: vd_strdc(:) => null() 
   real(kind=jprb), pointer :: vd_blh(:) => null() 
   real(kind=jprb), pointer :: vd_sund(:) => null() 
   real(kind=jprb), pointer :: vd_spar(:) => null() 
   real(kind=jprb), pointer :: vd_suvb(:) => null() 
   real(kind=jprb), pointer :: vd_sfdir(:) => null() 
   real(kind=jprb), pointer :: vd_scdir(:) => null() 
   real(kind=jprb), pointer :: vd_sdsrp(:) => null() 
   real(kind=jprb), pointer :: vd_cape(:) => null() 
   real(kind=jprb), pointer :: vd_capes(:) => null() 
   real(kind=jprb), pointer :: vd_mucape(:) => null() 
   real(kind=jprb), pointer :: vd_pdepl(:) => null() 
   real(kind=jprb), pointer :: vd_mlcape50(:) => null() 
   real(kind=jprb), pointer :: vd_mlcape100(:) => null() 
   real(kind=jprb), pointer :: vd_mlcin50(:) => null() 
   real(kind=jprb), pointer :: vd_mlcin100(:) => null() 
   real(kind=jprb), pointer :: vd_tropotp(:) => null() 
   real(kind=jprb), pointer :: vd_tsrc(:) => null() 
   real(kind=jprb), pointer :: vd_ttrc(:) => null() 
   real(kind=jprb), pointer :: vd_ssrc(:) => null() 
   real(kind=jprb), pointer :: vd_strc(:) => null() 
   real(kind=jprb), pointer :: vd_es(:) => null() 
   real(kind=jprb), pointer :: vd_smlt(:) => null() 
   real(kind=jprb), pointer :: vd_10fg(:) => null() 
   real(kind=jprb), pointer :: vd_10fg6(:) => null() 
   real(kind=jprb), pointer :: vd_10fgcv(:) => null() 
   real(kind=jprb), pointer :: vd_i10fg(:) => null() 
   real(kind=jprb), pointer :: vd_lspf(:) => null() 
   real(kind=jprb), pointer :: vd_tco3(:) => null() 
   real(kind=jprb), pointer :: vd_vimd(:) => null() 
   real(kind=jprb), pointer :: vd_sparc(:) => null() 
   real(kind=jprb), pointer :: vd_stinc(:) => null() 
   real(kind=jprb), pointer :: vd_cbase(:) => null() 
   real(kind=jprb), pointer :: vd_0degl(:) => null() 
   real(kind=jprb), pointer :: vd_m10degl(:) => null() 
   real(kind=jprb), pointer :: vd_visih(:) => null() 
   real(kind=jprb), pointer :: vd_cin(:) => null() 
   real(kind=jprb), pointer :: vd_kindex(:) => null() 
   real(kind=jprb), pointer :: vd_ttindex(:) => null() 
   real(kind=jprb), pointer :: vd_cbasea(:) => null() 
   real(kind=jprb), pointer :: vd_ctopc(:) => null() 
   real(kind=jprb), pointer :: vd_ztwetb0(:) => null() 
   real(kind=jprb), pointer :: vd_ztwetb1(:) => null() 
   real(kind=jprb), pointer :: vd_tcghg(:) => null() 
   real(kind=jprb), pointer :: vd_tcchem(:) => null() 
   real(kind=jprb), pointer :: vd_aerodiag(:,:) => null() 
   real(kind=jprb), pointer :: vd_aero_wvl_diag(:,:) => null() 
   real(kind=jprb), pointer :: vd_100u(:) => null() 
   real(kind=jprb), pointer :: vd_100v(:) => null() 
   real(kind=jprb), pointer :: vd_zust(:) => null() 
   real(kind=jprb), pointer :: vd_10nu(:) => null() 
   real(kind=jprb), pointer :: vd_10nv(:) => null() 
   real(kind=jprb), pointer :: vd_dndzn(:) => null() 
   real(kind=jprb), pointer :: vd_dndza(:) => null() 
   real(kind=jprb), pointer :: vd_dctb(:) => null() 
   real(kind=jprb), pointer :: vd_tplb(:) => null() 
   real(kind=jprb), pointer :: vd_tplt(:) => null() 
   real(kind=jprb), pointer :: vd_odss(:) => null() 
   real(kind=jprb), pointer :: vd_oddu(:) => null() 
   real(kind=jprb), pointer :: vd_odom(:) => null() 
   real(kind=jprb), pointer :: vd_odbc(:) => null() 
   real(kind=jprb), pointer :: vd_odsu(:) => null() 
   real(kind=jprb), pointer :: vd_odni(:) => null() 
   real(kind=jprb), pointer :: vd_odam(:) => null() 
   real(kind=jprb), pointer :: vd_odsoa(:) => null() 
   real(kind=jprb), pointer :: vd_odvfa(:) => null() 
   real(kind=jprb), pointer :: vd_odvsu(:) => null() 
   real(kind=jprb), pointer :: vd_odtoacc(:) => null() 
   real(kind=jprb), pointer :: vd_aepm1(:) => null() 
   real(kind=jprb), pointer :: vd_aepm25(:) => null() 
   real(kind=jprb), pointer :: vd_aepm10(:) => null() 
   real(kind=jprb), pointer :: vd_uvbed(:) => null() 
   real(kind=jprb), pointer :: vd_uvbedcs(:) => null() 
   real(kind=jprb), pointer :: vd_litoti(:) => null() 
   real(kind=jprb), pointer :: vd_licgi(:) => null() 
   real(kind=jprb), pointer :: vd_litota6(:) => null() 
   real(kind=jprb), pointer :: vd_licga6(:) => null()
   real(kind=jprb), pointer :: vd_ptypeocc6(:,:) => null()
   real(kind=jprb), pointer :: vd_200u(:) => null() 
   real(kind=jprb), pointer :: vd_200v(:) => null() 

   real(kind=jprb), pointer :: vd_sdsl(:) => null() 

! * Group SM=SATSIM: (ECMWF) simulated satellite images:
   real(kind=jprb), pointer :: sm_clbt(:,:) => null()
   real(kind=jprb), pointer :: sm_csbt(:,:) => null()
   
  ! Group WS
  real(kind=jprb), pointer :: ws_char(:) => null() 
  real(kind=jprb), pointer :: ws_charhq(:) => null() 
  real(kind=jprb), pointer :: ws_ustokes(:) => null() 
  real(kind=jprb), pointer :: ws_vstokes(:) => null() 
  real(kind=jprb), pointer :: ws_tauocx(:) => null() 
  real(kind=jprb), pointer :: ws_tauocy(:) => null() 
  real(kind=jprb), pointer :: ws_phioc(:) => null() 
  real(kind=jprb), pointer :: ws_wsemean(:) => null() 
  real(kind=jprb), pointer :: ws_wsfmean(:) => null()
  !  Group WW 
  real(kind=jprb), pointer :: ww_u10n(:) => null() 
  real(kind=jprb), pointer :: ww_v10n(:) => null() 
  real(kind=jprb), pointer :: ww_rho(:) => null() 
  real(kind=jprb), pointer :: ww_zil(:) => null() 
  real(kind=jprb), pointer :: ww_cif(:) => null() 
  real(kind=jprb), pointer :: ww_clk(:) => null() 
  real(kind=jprb), pointer :: ww_ucurw(:) => null() 
  real(kind=jprb), pointer :: ww_vcurw(:) => null() 

  ! Group VX
  real(kind=jprb), pointer :: vx_oro(:) => null() 
  real(kind=jprb), pointer :: vx_tsc(:) => null() 
  real(kind=jprb), pointer :: vx_pws(:) => null() 
  real(kind=jprb), pointer :: vx_pwp(:) => null() 
  real(kind=jprb), pointer :: vx_sno(:) => null() 
  real(kind=jprb), pointer :: vx_tpc(:) => null() 
  real(kind=jprb), pointer :: vx_sab(:) => null() 
  real(kind=jprb), pointer :: vx_xd2(:) => null() 
  real(kind=jprb), pointer :: vx_lsm(:) => null() 
  real(kind=jprb), pointer :: vx_iveg(:) => null() 
  real(kind=jprb), pointer :: vx_arg(:) => null() 
  real(kind=jprb), pointer :: vx_rsmin(:) => null() 
  real(kind=jprb), pointer :: vx_lai(:) => null() 
  real(kind=jprb), pointer :: vx_veg(:) => null() 

  ! * Group VK=VCLIK: Convective cloud pseudo-historic fields:
   real(kind=jprb), pointer :: vk_udgro(:) => null()

  class(field_access),pointer :: dm => null()
  class(field_access),pointer :: dl => null()

  contains

  procedure :: field_map_storage => main_field_map_storage

end type

type, extends(field_metadata_base) :: field_metadata

  contains

  procedure :: field_set_metadata => main_field_set_metadata
  procedure :: field_get_clevtype => main_field_get_clevtype

end type

! We will always work internally with lists of generic fields, identified 
! by their field ID. 
integer(kind=jpim),parameter :: jpnumfids=4200
! (F2)
type type_field_id
!!$  integer(kind=jpim) :: u, v, div,vor,t, pd, vd, nhx, edot,vcasrs,rdlasr,&
!!$   & sgrtl,sgrtm,unl,unl_si,vnl,vnl_si,tnl,tnl_si,spnl,spnl_si,vwvnl,&
!!$   & vdnl_si,pdnl,pdnl_si,gw,nhy,curhs,cvrhs,ctrhs,csprhs,cspdrhs,csvdrhs,&
!!$   & cunl,cvnl,ctnl,cspnl,cvwvnl,cspdnl,cupt,cvpt,ctpt,cpdpt,cvdpt,dphi,nhxnl,cnhxnl,&
!!$    & q, o3, l, i, a, ghg, grg, aero, phys, s, r, g, h, lconv,iconv,rconv,sconv,&
!!$    & tke,tte,efb1,efb2,efb3,mxl,src,shtur,fqtur,fstur,cvv,cpf,spf,cvgq,&
!!$    & uom,ual,dom,dal,uen,unebh,rkth,rktqv,rktqc,&
!!$    & lrad,irad,rspec,pre, pref, delp, extra3d,ezdiag,&
!!$    & sp, csppt,dbbc,gws,spnl2,cspnl2,prehyds,u_mean,v_mean,orog, cori, gnordl, gnordm, ccc, lcc, mcc, hcc, tcc, &
!!$    & sb_q, sg_f, sg_a, sg_r, rr_t, rr_w, rr_ic, cl_tcls, &
!!$    & cl_hucls, x2_prwa, x2_prsn, vf_z0f, vf_albf, vf_emisf, vf_getrl, vf_lsm, vf_veg, vf_cvl, vf_co2typ, &
!!$    & vf_cvh, vf_tvl, vf_tvh, vf_lail, vf_laih, vf_soty, vf_ci, vf_ucur, vf_vcur, &
!!$    & vf_aluvp, vf_aluvd, vf_alnip, vf_alnid, &
!!$    & vf_aluvi, vf_aluvv, vf_aluvg, vf_alnii, vf_alniv, vf_alnig, &
!!$    & vv_arg, vv_sab, vv_hv, vv_z0h, vn_top, vn_bas, vn_acpr, vn_accpr, &
!!$    & ws_char, vd_10nu, vd_10nv, vd_upd, vd_z0f, vx_oro, vx_tsc, vx_sno, est, esn, et2, eh2, &
!!$    & ev1, ez, eh, vvel
! From GMV
integer(kind=jpim) :: u=1
integer(kind=jpim) :: v=2
integer(kind=jpim) :: div=3
integer(kind=jpim) :: vor=4
integer(kind=jpim) :: t=5
integer(kind=jpim) :: pd=6
integer(kind=jpim) :: vd=7
integer(kind=jpim) :: nhx=8
integer(kind=jpim) :: edot=9
! ky: GMV fields VCASRS and RDLASR (containing a/r) have been pruned, but how to renumber fields?
integer(kind=jpim) :: vcasrs=10 ! ky: obsolete field (deep-layer models pruned)
integer(kind=jpim) :: rdlasr=11 ! ky: obsolete field (deep-layer models pruned)
! -----------------------------------------------------------------------------------------------
integer(kind=jpim) :: sgrtl=12
integer(kind=jpim) :: sgrtm=13
integer(kind=jpim) :: unl=14
integer(kind=jpim) :: unl_si=15
integer(kind=jpim) :: vnl=16
integer(kind=jpim) :: vnl_si=17
integer(kind=jpim) :: tnl=18
integer(kind=jpim) :: tnl_si=19
integer(kind=jpim) :: spnl=20
integer(kind=jpim) :: spnl_si=21
integer(kind=jpim) :: vwvnl=22
integer(kind=jpim) :: vdnl_si=23
integer(kind=jpim) :: pdnl=24
integer(kind=jpim) :: pdnl_si=25
integer(kind=jpim) :: gw=26
integer(kind=jpim) :: curhs=27
integer(kind=jpim) :: cvrhs=28
integer(kind=jpim) :: ctrhs=29
integer(kind=jpim) :: csprhs=30
integer(kind=jpim) :: cspdrhs=31
integer(kind=jpim) :: csvdrhs=32
integer(kind=jpim) :: cunl=33
integer(kind=jpim) :: cvnl=34
integer(kind=jpim) :: ctnl=35
integer(kind=jpim) :: cspnl=36
integer(kind=jpim) :: cvwvnl=37
integer(kind=jpim) :: cspdnl=38
integer(kind=jpim) :: cupt=39
integer(kind=jpim) :: cvpt=40
integer(kind=jpim) :: ctpt=41
integer(kind=jpim) :: cpdpt=42
integer(kind=jpim) :: cvdpt=43
integer(kind=jpim) :: dphi=44
integer(kind=jpim) :: nhxnl=45
integer(kind=jpim) :: cnhxnl=46
! From gfl
integer(kind=jpim) :: q=47
integer(kind=jpim) :: o3=48 
integer(kind=jpim) :: l=49 
integer(kind=jpim) :: i=50 
integer(kind=jpim) :: a=51 
!integer(kind=jpim) :: ghg=52   ! }
!integer(kind=jpim) :: chem=53  ! } Split into 4xxx
!integer(kind=jpim) :: aero=54  ! }
integer(kind=jpim) :: lrch4=52
integer(kind=jpim) :: phys=55 
integer(kind=jpim) :: s=56
integer(kind=jpim) :: r=57
integer(kind=jpim) :: g=58 
integer(kind=jpim) :: h=59 
integer(kind=jpim) :: lconv=60
integer(kind=jpim) :: iconv=61
integer(kind=jpim) :: rconv=62
integer(kind=jpim) :: sconv=63
integer(kind=jpim) :: tke=64
integer(kind=jpim) :: tte=65
integer(kind=jpim) :: efb1=66
integer(kind=jpim) :: efb2=67
integer(kind=jpim) :: efb3=68
integer(kind=jpim) :: mxl=69
integer(kind=jpim) :: src=70
integer(kind=jpim) :: shtur=71
integer(kind=jpim) :: fqtur=72
integer(kind=jpim) :: fstur=73
integer(kind=jpim) :: cvv=74
integer(kind=jpim) :: cpf=75
integer(kind=jpim) :: spf=76
integer(kind=jpim) :: cvgq=77
integer(kind=jpim) :: uom=78
integer(kind=jpim) :: ual=79
integer(kind=jpim) :: dom=80
integer(kind=jpim) :: dal=81
integer(kind=jpim) :: uen=82
integer(kind=jpim) :: unebh=83
integer(kind=jpim) :: rkth=84
integer(kind=jpim) :: rktqv=85
integer(kind=jpim) :: rktqc=86
integer(kind=jpim) :: lrad=87
integer(kind=jpim) :: irad=88
integer(kind=jpim) :: rspec=89
integer(kind=jpim) :: pre=90
integer(kind=jpim) :: pref=91 
integer(kind=jpim) :: delp=92 
integer(kind=jpim) :: extra3d=93
integer(kind=jpim) :: ezdiag=94
integer(kind=jpim) :: sp=95 
integer(kind=jpim) :: csppt=96
integer(kind=jpim) :: dbbc=97
integer(kind=jpim) :: gws=98
integer(kind=jpim) :: spnl2=99
integer(kind=jpim) :: cspnl2=100
integer(kind=jpim) :: u_mean=101
integer(kind=jpim) :: v_mean=102
integer(kind=jpim) :: orog=103
integer(kind=jpim) :: cori=104
integer(kind=jpim) :: gnordl=105 
integer(kind=jpim) :: gnordm=106
integer(kind=jpim) :: ccc=107 
integer(kind=jpim) :: lcc=108 
integer(kind=jpim) :: mcc=109 
integer(kind=jpim) :: hcc=110 
integer(kind=jpim) :: tcc=111
integer(kind=jpim) :: est=159
integer(kind=jpim) :: esn=160 
integer(kind=jpim) :: et2=161 
integer(kind=jpim) :: eh2=162 
integer(kind=jpim) :: ev1=163 
integer(kind=jpim) :: ez=164
integer(kind=jpim) :: eh=165 
integer(kind=jpim) :: vvel=166
integer(kind=jpim) :: nhy=167
integer(kind=jpim) :: lmf=168
integer(kind=jpim) :: imf=169
integer(kind=jpim) :: amf=170
integer(kind=jpim) :: extcv2d=171
integer(kind=jpim) :: extcv3d=171+jp_necv_2d_max

! ky: sorry, I don't know the rule of numbering this
integer(kind=jpim) :: prehyds=999

  ! Surface (see surface_fields_mix.F90)
  ! Group SB
integer(kind=jpim) :: sb_q=1101
integer(kind=jpim) :: sb_t=1102
integer(kind=jpim) :: sb_tl=1103
  ! Group SG
integer(kind=jpim) :: sg_f=1201 
integer(kind=jpim) :: sg_a=1202
integer(kind=jpim) :: sg_r=1203
integer(kind=jpim) :: sg_t=1204
integer(kind=jpim) :: sg_w=1205
! Group SL
integer(kind=jpim) :: sl_lict=1301
integer(kind=jpim) :: sl_lmlt=1302
integer(kind=jpim) :: sl_ltlt=1303
integer(kind=jpim) :: sl_lblt=1304
integer(kind=jpim) :: sl_lshf=1305
integer(kind=jpim) :: sl_licd=1306
integer(kind=jpim) :: sl_lmld=1307
! Group RR
integer(kind=jpim) :: rr_t=1401
integer(kind=jpim) :: rr_w=1402
integer(kind=jpim) :: rr_fc=1403 
integer(kind=jpim) :: rr_ic=1404
integer(kind=jpim) :: rr_fp1=1405
integer(kind=jpim) :: rr_tir=1406
integer(kind=jpim) :: rr_tmw=1407
!Group CL
integer(kind=jpim) :: cl_tcls=1501
integer(kind=jpim) :: cl_hucls=1502 
!Group OM
integer(kind=jpim) :: om_to=1601
integer(kind=jpim) :: om_so=1602
integer(kind=jpim) :: om_uo=1603
integer(kind=jpim) :: om_vo=1604
! Group EP
integer(kind=jpim) :: ep_ep=1701
! Group X2
integer(kind=jpim) :: x2_x2=1801
! Group CI
integer(kind=jpim) :: ci_ci=1901
! Group VF
integer(kind=jpim) :: vf_z0f=2001 
integer(kind=jpim) :: vf_albf=2002 
integer(kind=jpim) :: vf_emisf=2003
integer(kind=jpim) :: vf_getrl=2004 
integer(kind=jpim) :: vf_lsm=2005
integer(kind=jpim) :: vf_veg=2006
integer(kind=jpim) :: vf_vrlan=2007 
integer(kind=jpim) :: vf_vrldi=2008 
integer(kind=jpim) :: vf_sig=2009 
integer(kind=jpim) :: vf_albsf=2010
integer(kind=jpim) :: vf_lan=2011
integer(kind=jpim) :: vf_sst=2012 
integer(kind=jpim) :: vf_sss=2013 
integer(kind=jpim) :: vf_lz0h=2014 
integer(kind=jpim) :: vf_cvl=2015
integer(kind=jpim) :: vf_cvh=2016
integer(kind=jpim) :: vf_tvl=2017 
integer(kind=jpim) :: vf_tvh=2018 
integer(kind=jpim) :: vf_lail=2019
integer(kind=jpim) :: vf_laih=2020 
integer(kind=jpim) :: vf_soty=2021
integer(kind=jpim) :: vf_clk=2022
integer(kind=jpim) :: vf_dl=2023
integer(kind=jpim) :: vf_ci=2024 
integer(kind=jpim) :: vf_ucur=2025 
integer(kind=jpim) :: vf_vcur=2026
integer(kind=jpim) :: vf_z0rlf=2027
integer(kind=jpim) :: vf_cgpp=2032
integer(kind=jpim) :: vf_crec=2033
integer(kind=jpim) :: vf_sdfor=2036
integer(kind=jpim) :: vf_aluvp=2037 
integer(kind=jpim) :: vf_aluvd=2038 
integer(kind=jpim) :: vf_alnip=2039
integer(kind=jpim) :: vf_alnid=2040
integer(kind=jpim) :: vf_fp1=2041     ! surface orography in the 2nd part of FULLPOS-927
integer(kind=jpim) :: vf_so2dd=2051   ! sulphate dry dep velocity
integer(kind=jpim) :: vf_dmso=2056    ! oceanic dms
integer(kind=jpim) :: vf_urbf=2058   ! soa from co
integer(kind=jpim) :: vf_fca1=2060    ! fraction of calcite over dust 1st bin
integer(kind=jpim) :: vf_fca2=2061    ! fraction of calcite over dust 2st bin
integer(kind=jpim) :: vf_aerdep=2062  ! dust emission potential 
integer(kind=jpim) :: vf_aerlts=2063  ! dust lifting threshold speed 
integer(kind=jpim) :: vf_aerscc=2064  ! dust soil clay content
integer(kind=jpim) :: vf_dsf=2065     ! dust source function
integer(kind=jpim) :: vf_chemflxo=2067 ! total chemistry flux (emissions + deposition) 
integer(kind=jpim) :: vf_chemdv=2068  ! chemistry deposition velocity
integer(kind=jpim) :: vf_nudm=2069    ! nudging mask

! MODIS 6-component albedo coefficients
integer(kind=jpim) :: vf_aluvi=2070
integer(kind=jpim) :: vf_aluvv=2071
integer(kind=jpim) :: vf_aluvg=2072
integer(kind=jpim) :: vf_alnii=2073
integer(kind=jpim) :: vf_alniv=2074
integer(kind=jpim) :: vf_alnig=2075

integer(kind=jpim) :: vf_cur=2076 ! Urban

integer(kind=jpim) :: vf_emis2d=2077 ! 2D emission fields for composition
integer(kind=jpim) :: vf_emis2daux=2078 ! 2D emission auxiliary fields for composition

integer(kind=jpim) :: vf_dsz=2079     ! dust size variation

! Wet and dry deposition fluxes  
integer(kind=jpim) :: vf_chemwdflx=2080 ! wet deposition chemistry flux  
integer(kind=jpim) :: vf_chemddflx=2081 ! dry deposition chemistry flux  

integer(kind=jpim) :: vf_co2typ=2082 ! C3/C4 CO2 photosynthesis type for low vegetation

integer(kind=jpim) :: vf_fwet=2083 ! Wetland

! * Group VP=VCLIP: deep soil diagnostic fields
integer(kind=jpim) :: vp_tpc=2101 ! climatological deep layer temperature
integer(kind=jpim) :: vp_wpc=2102 ! climatological deep layer moisture
! * Group VV=VCLIV: vegetation diagnostic fields:
integer(kind=jpim) :: vv_arg=2201
integer(kind=jpim) :: vv_sab=2202
integer(kind=jpim) :: vv_d2=2203
integer(kind=jpim) :: vv_iveg=2204
integer(kind=jpim) :: vv_rsmin=2205
integer(kind=jpim) :: vv_lai=2206
integer(kind=jpim) :: vv_hv=2207 
integer(kind=jpim) :: vv_z0h=2208 
integer(kind=jpim) :: vv_als=2209
integer(kind=jpim) :: vv_alv=2210
! * Group VN=VCLIN: cloudiness diagnostic predictors:
integer(kind=jpim) :: vn_top=2301
integer(kind=jpim) :: vn_bas=2302
integer(kind=jpim) :: vn_acpr=2303
integer(kind=jpim) :: vn_accpr=2304
integer(kind=jpim) :: vn_accpr5=2305
! * Group VH=VCLIH: convective cloud diagnostic fields:
integer(kind=jpim) :: vh_tcch=2401
integer(kind=jpim) :: vh_scch=2402
integer(kind=jpim) :: vh_bcch=2403
integer(kind=jpim) :: vh_pblh=2404
integer(kind=jpim) :: vh_spsh=2405
integer(kind=jpim) :: vh_qsh =2406
! Group VA=VCLIA: aerosol diagnostic fields:
integer(kind=jpim) :: va_sea=2501
integer(kind=jpim) :: va_lan=2502
integer(kind=jpim) :: va_soo=2503
integer(kind=jpim) :: va_des=2504
integer(kind=jpim) :: va_sul=2505
integer(kind=jpim) :: va_vol=2506
! Group  VG=VCLIG: ice-coupler diagnostic fields
integer(kind=jpim) :: vg_icfr=2601
integer(kind=jpim) :: vg_soup=2602
integer(kind=jpim) :: vg_irup=2603
integer(kind=jpim) :: vg_chss=2604
integer(kind=jpim) :: vg_evap=2605
integer(kind=jpim) :: vg_taux=2606
integer(kind=jpim) :: vg_tauy=2607
!  Group VC=VO3ABC: A,B and C (Climatological ozone profiles) diagnostic fields
integer(kind=jpim) :: vc_a=2701
integer(kind=jpim) :: vc_b=2702
integer(kind=jpim) :: vc_c=2703
! Group V2=VDIAGO2: 2-D climatological/diagnostic fields for an ocean mixed layer model (KPP)
integer(kind=jpim) :: v2_ocdep=2801
integer(kind=jpim) :: v2_ustrc=2802
integer(kind=jpim) :: v2_vstrc=2803
! Group V3=VDIAGO3: 3-D climatological/diagnostic fields for an ocean mixed layer model (KPP):
integer(kind=jpim) :: v3_difm=2901
integer(kind=jpim) :: v3_dift=2902
integer(kind=jpim) :: v3_difs=2903
integer(kind=jpim) :: v3_advt=2904
integer(kind=jpim) :: v3_advs=2905
integer(kind=jpim) :: v3_tri0=2906
integer(kind=jpim) :: v3_tri1=2907
integer(kind=jpim) :: v3_swdk=2908
integer(kind=jpim) :: v3_zo=2909
integer(kind=jpim) :: v3_ho=2910
integer(kind=jpim) :: v3_do=2911
integer(kind=jpim) :: v3_ho_inv=2912
integer(kind=jpim) :: v3_uoc=2913
integer(kind=jpim) :: v3_voc=2914
integer(kind=jpim) :: v3_otke=2915

! * Group VD=VDIAG: (ECMWF) diagnostic fields:
integer(kind=jpim) :: vd_lsp=3001
integer(kind=jpim) :: vd_cp=3002 
integer(kind=jpim) :: vd_sf=3003 
integer(kind=jpim) :: vd_fzra=3004 
integer(kind=jpim) :: vd_bld=3005 
integer(kind=jpim) :: vd_sshf=3006 
integer(kind=jpim) :: vd_slhf=3007 
integer(kind=jpim) :: vd_nee=3008 
integer(kind=jpim) :: vd_gpp=3009 
integer(kind=jpim) :: vd_rec=3010 
integer(kind=jpim) :: vd_msl=3011 
integer(kind=jpim) :: vd_sp=3012 
integer(kind=jpim) :: vd_tcc=3013 
integer(kind=jpim) :: vd_10u=3014 
integer(kind=jpim) :: vd_10v=3015 
integer(kind=jpim) :: vd_2t=3016 
integer(kind=jpim) :: vd_2d=3017 
integer(kind=jpim) :: vd_ssr=3018 
integer(kind=jpim) :: vd_str=3019 
integer(kind=jpim) :: vd_tsr=3020 
integer(kind=jpim) :: vd_ttr=3021
integer(kind=jpim) :: vd_ewss=3022 
integer(kind=jpim) :: vd_nsss=3023 
integer(kind=jpim) :: vd_e=3024 
integer(kind=jpim) :: vd_pev=3025 
integer(kind=jpim) :: vd_ccc=3026 
integer(kind=jpim) :: vd_lcc=3027 
integer(kind=jpim) :: vd_mcc=3028 
integer(kind=jpim) :: vd_hcc=3029 
integer(kind=jpim) :: vd_lgws=3030 
integer(kind=jpim) :: vd_mgws=3031 
integer(kind=jpim) :: vd_gwd=3032 
integer(kind=jpim) :: vd_mx2t=3033 
integer(kind=jpim) :: vd_mn2t=3034 
integer(kind=jpim) :: vd_mx2t6=3035 
integer(kind=jpim) :: vd_mn2t6=3036 
integer(kind=jpim) :: vd_ro=3037 
integer(kind=jpim) :: vd_sro=3038 
integer(kind=jpim) :: vd_ssro=3039 
integer(kind=jpim) :: vd_alb=3040 
integer(kind=jpim) :: vd_iewss=3041 
integer(kind=jpim) :: vd_insss=3042 
integer(kind=jpim) :: vd_isshf=3043 
integer(kind=jpim) :: vd_ie=3044 
integer(kind=jpim) :: vd_inee=3045 
integer(kind=jpim) :: vd_igpp=3046 
integer(kind=jpim) :: vd_irec=3047 
integer(kind=jpim) :: vd_csf=3048 
integer(kind=jpim) :: vd_lssf=3049 
integer(kind=jpim) :: vd_mxtpr=3050 
integer(kind=jpim) :: vd_mntpr=3051 
integer(kind=jpim) :: vd_mxtpr6=3052 
integer(kind=jpim) :: vd_mntpr6=3053 
integer(kind=jpim) :: vd_tpr=3054 
integer(kind=jpim) :: vd_lsrr=3055 
integer(kind=jpim) :: vd_crr=3056
integer(kind=jpim) :: vd_lssfr=3057 
integer(kind=jpim) :: vd_csfr=3058
integer(kind=jpim) :: vd_ptype=3059 
integer(kind=jpim) :: vd_ilspf=3060 
integer(kind=jpim) :: vd_z0f=3061 
integer(kind=jpim) :: vd_lz0h=3062 
integer(kind=jpim) :: vd_tcw=3063 
integer(kind=jpim) :: vd_tcwv=3064 
integer(kind=jpim) :: vd_tclw=3065 
integer(kind=jpim) :: vd_tciw=3066 
integer(kind=jpim) :: vd_tcrw=3067 
integer(kind=jpim) :: vd_tcsw=3068 
integer(kind=jpim) :: vd_tcslw=3069 
integer(kind=jpim) :: vd_ssrd=3070
integer(kind=jpim) :: vd_strd=3071 
integer(kind=jpim) :: vd_ssrdc=3072 
integer(kind=jpim) :: vd_strdc=3073 
integer(kind=jpim) :: vd_blh=3074 
integer(kind=jpim) :: vd_sund=3075 
integer(kind=jpim) :: vd_spar=3076 
integer(kind=jpim) :: vd_suvb=3077 
integer(kind=jpim) :: vd_sfdir=3078 
integer(kind=jpim) :: vd_scdir=3079 
integer(kind=jpim) :: vd_sdsrp=3080 
integer(kind=jpim) :: vd_cape=3081
integer(kind=jpim) :: vd_capes=3082
integer(kind=jpim) :: vd_tsrc=3083
integer(kind=jpim) :: vd_ttrc=3084 
integer(kind=jpim) :: vd_ssrc=3085 
integer(kind=jpim) :: vd_strc=3086 
integer(kind=jpim) :: vd_es=3087 
integer(kind=jpim) :: vd_smlt=3088 
integer(kind=jpim) :: vd_10fg=3089 
integer(kind=jpim) :: vd_10fg6=3090 
integer(kind=jpim) :: vd_10fgcv=3091 
integer(kind=jpim) :: vd_i10fg=3092
integer(kind=jpim) :: vd_lspf=3093
integer(kind=jpim) :: vd_tco3=3094 
integer(kind=jpim) :: vd_vimd=3095 
integer(kind=jpim) :: vd_sparc=3096 
integer(kind=jpim) :: vd_stinc=3097 
integer(kind=jpim) :: vd_cbase=3098 
integer(kind=jpim) :: vd_0degl=3099 
integer(kind=jpim) :: vd_visih=3100 
integer(kind=jpim) :: vd_cin=3101 
integer(kind=jpim) :: vd_kindex=3102 
integer(kind=jpim) :: vd_ttindex=3103 
integer(kind=jpim) :: vd_cbasea=3104 
integer(kind=jpim) :: vd_ctopc=3105 
integer(kind=jpim) :: vd_ztwetb0=3106 
integer(kind=jpim) :: vd_ztwetb1=3107 
integer(kind=jpim) :: vd_tcghg=3108
integer(kind=jpim) :: vd_tcchem=3109
integer(kind=jpim) :: vd_aerodiag=3110
integer(kind=jpim) :: vd_aero_wvl_diag=3111
integer(kind=jpim) :: vd_100u=3112 
integer(kind=jpim) :: vd_100v=3113 
integer(kind=jpim) :: vd_zust=3114 
integer(kind=jpim) :: vd_10nu=3115 
integer(kind=jpim) :: vd_10nv=3116 
integer(kind=jpim) :: vd_dndzn=3117 
integer(kind=jpim) :: vd_dndza=3118 
integer(kind=jpim) :: vd_dctb=3119
integer(kind=jpim) :: vd_tplb=3120 
integer(kind=jpim) :: vd_tplt=3121 
integer(kind=jpim) :: vd_odss=3122 
integer(kind=jpim) :: vd_oddu=3123 
integer(kind=jpim) :: vd_odom=3124 
integer(kind=jpim) :: vd_odbc=3125 
integer(kind=jpim) :: vd_odsu=3126
integer(kind=jpim) :: vd_odni=3127
integer(kind=jpim) :: vd_odam=3128 
integer(kind=jpim) :: vd_odsoa=3129
integer(kind=jpim) :: vd_odvfa=3130 
integer(kind=jpim) :: vd_odvsu=3131 
integer(kind=jpim) :: vd_aepm1=3132 
integer(kind=jpim) :: vd_aepm25=3133 
integer(kind=jpim) :: vd_aepm10=3134 
integer(kind=jpim) :: vd_uvbed=3135
integer(kind=jpim) :: vd_uvbedcs=3136 
integer(kind=jpim) :: vd_litoti=3137
integer(kind=jpim) :: vd_licgi=3138
integer(kind=jpim) :: vd_litota6=3139
integer(kind=jpim) :: vd_licga6=3140
integer(kind=jpim) :: vd_200u=3141 
integer(kind=jpim) :: vd_200v=3142 
integer(kind=jpim) :: vd_2sh=3143
integer(kind=jpim) :: vd_odtoacc=3144
integer(kind=jpim) :: vd_m10degl=3145
integer(kind=jpim) :: vd_mucape=3146
integer(kind=jpim) :: vd_pdepl=3147
integer(kind=jpim) :: vd_mlcape50=3148
integer(kind=jpim) :: vd_mlcape100=3149
integer(kind=jpim) :: vd_mlcin50=3150
integer(kind=jpim) :: vd_mlcin100=3151
integer(kind=jpim) :: vd_tropotp=3152
integer(kind=jpim) :: vd_sdsl=3153
integer(kind=jpim) :: vd_ich4=3154
integer(kind=jpim) :: vd_ptypeocc6=3155

! Group SM=SATSIM: (ECMWF) simulated satellite images
integer(kind=jpim) :: sm_clbt=3201
integer(kind=jpim) :: sm_csbt=3202

! Group WS=WAVES: surface prognostic quantities over sea (used by IFS)
integer(kind=jpim) :: ws_char=3301
integer(kind=jpim) :: ws_charhq=3302
integer(kind=jpim) :: ws_ustokes=3303
integer(kind=jpim) :: ws_vstokes=3304
integer(kind=jpim) :: ws_tauocx=3305
integer(kind=jpim) :: ws_tauocy=3306
integer(kind=jpim) :: ws_phioc=3307
integer(kind=jpim) :: ws_wsemean=3308
integer(kind=jpim) :: ws_wsfmean=3309

! * Group VX=VCLIX: auxilary climatological diagnostic fields:
integer(kind=jpim) :: vx_oro=3401
integer(kind=jpim) :: vx_tsc=3402
integer(kind=jpim) :: vx_pws=3403
integer(kind=jpim) :: vx_pwp=3404
integer(kind=jpim) :: vx_sno=3405
integer(kind=jpim) :: vx_tpc=3406
integer(kind=jpim) :: vx_sab=3407
integer(kind=jpim) :: vx_xd2=3408
integer(kind=jpim) :: vx_lsm=3409
integer(kind=jpim) :: vx_iveg=3410
integer(kind=jpim) :: vx_arg=3411
integer(kind=jpim) :: vx_rsmin=3412
integer(kind=jpim) :: vx_lai=3413
integer(kind=jpim) :: vx_veg=3414
! * Group VK=VCLIK: Convective cloud pseudo-historic fields:
integer(kind=jpim) :: vk_udgro=3501

! * Dynamically-assigned COMPO 3D fields
!integer(kind=jpim) :: compo_3d_first = 3601
!integer(kind=jpim) :: compo_3d_last = 3800
  !  Group WW 

integer(kind=jpim) :: ww_u10n=3801 
integer(kind=jpim) :: ww_v10n=3802 
integer(kind=jpim) :: ww_rho= 3803
integer(kind=jpim) :: ww_zil= 3804
integer(kind=jpim) :: ww_cif= 3805
integer(kind=jpim) :: ww_clk= 3806
integer(kind=jpim) :: ww_ucurw=3807 
integer(kind=jpim) :: ww_vcurw= 3808

end type type_field_id

! (F3)
type type_dynamic_names
  character(len=JP_NAME_MAX_LEN) :: cname=''
  character(len=JP_COMMENTS_MAX_LEN) :: clongname=''
  integer(kind=jpim) :: igribcode
  integer(kind=jpim) :: ifid
end type type_dynamic_names

! Dynamic namespace (currently for COMPO)
integer(kind=jpim),parameter :: compo_3d_first = 3601
integer(kind=jpim),parameter :: compo_3d_last = 3800
integer(kind=jpim)           :: ndyn_names=0
type(type_dynamic_names) ::  dyn_nam(compo_3d_last-compo_3d_first+1)

type(type_field_id), parameter :: fid=type_field_id()
type(field_metadata) :: main_field_metadata

! 2nd dimension types
integer(kind=jpim), parameter :: d2none      = 0 ! no second dimension
integer(kind=jpim), parameter :: d2full      = 1 ! full model levels
integer(kind=jpim), parameter :: d2half      = 2 ! half model levels
integer(kind=jpim), parameter :: d2soil      = 3 ! soil model levels
integer(kind=jpim), parameter :: d2snow      = 4 ! snow model levels
integer(kind=jpim), parameter :: d2om        = 5 ! Ocean model "levels"
integer(kind=jpim), parameter :: d2chemflxo  = 6
integer(kind=jpim), parameter :: d2chemwdflx = 7
integer(kind=jpim), parameter :: d2chemddflx = 8
integer(kind=jpim), parameter :: d2chemdv    = 9
integer(kind=jpim), parameter :: d2emis2d    = 10
integer(kind=jpim), parameter :: d2emis2daux = 11

integer(kind=jpim), parameter :: nleveltypes_main = 11

! 3rd dimension info - AJGDB where does this come from when it's used for real?

integer(kind=jpim),parameter :: ndim3types_main=2
integer(kind=jpim), parameter :: d3extra3d = 1
integer(kind=jpim), parameter :: d3ezdiag  = 2

!integer(kind=jpim), parameter :: nlev_3d = 2

#ifndef FIELD_MOD_TEST
#include "abor1.intfb.h"
#endif

contains

! ---------------------------------------------------------
! (Essentially a constructor) Initialises names, dimensions 
! and other attributes
! ---------------------------------------------------------
subroutine main_field_set_metadata(this,metadata,kleveltypes,kdim3types)
class(field_metadata),    intent(in)    :: this
type(type_fvar), ALLOCATABLE, intent(inout) :: metadata(:)
integer(kind=jpim),       intent(  out) :: kleveltypes
integer(kind=jpim),       intent(  out) :: kdim3types



















 
 













































 





























































































































































































































































































































































































































































end subroutine main_field_set_metadata

function main_field_get_clevtype(this, kleveltype) result(clevtype)
class(field_metadata), intent(in) :: this
integer(kind=jpim), intent(in) :: kleveltype
character(len=3) :: clevtype

end function main_field_get_clevtype

! -------------------------------------------------
!
! Perhaps clunky. This kind of technique was not
! used for the GOMs, because it could be quite slow 
! in that context. It should be more appropriate 
! here when the number of calls should not be too 
! large.
! -------------------------------------------------
subroutine main_field_map_storage(this, kid, storage_1d, storage_2d, storage_3d, ld_nullify)
class(field_access),               intent(inout) :: this
integer(kind=jpim),                intent(in)    :: kid
real(kind=jprb), optional, target, intent(in)    :: storage_1d(:), storage_2d(:,:),storage_3d(:,:,:)
logical, optional,                 intent(in)    :: ld_nullify










end subroutine main_field_map_storage

subroutine read_dynamic_namespace












end subroutine read_dynamic_namespace
end module field_definitions

