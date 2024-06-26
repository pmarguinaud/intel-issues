MODULE YOMVAR

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Switches for variational assimilation
! LTEST    : .T. = TEST OF THE GRADIENT BEFORE THE MINIMIZATION
! LREPRO4DVAR: .T. = Runs 4D-Var in bit reproducible mode (slower)
! LSLADREP   : CONTROL BIT REPRODUCIBILITY
!            : T - SLAD IS BIT REPRODUCIBLE WHEN NUMBER OF PROCESSORS
!            :     OR PARTITITIONING IS CHANGED
!            : F - SLAD IS NOT BIT REPRODUCIBLE WHEN NUMBER OF PROCESSORS
!            :     OR PARTITITIONING IS CHANGED
! L_CHECK_CONVERGENCE : .T. = ABORT IF THE MINIMIZATION FAILS
! R_NORM_REDUCTION_ABORT_LEVEL: = ABORT IF THE RATIO OF FINAL/INITIAL
!                            GRADIENT NORM SQUARED EXCEEDS THIS FACTOR
! L_CHECK_GRADIENT: .T. => CONGRAD checks its gradient each iteration.
!                          (This performs an additional call to SIM4D each iteration.
!                           If a problem is found, the gradient is written out and
!                           the minimization aborts.)
! RTOL_CHECK_GRADIENT: Gradient error norm (relative to gradient norm) required
!                      to trigger an abort if L_CHECK_GRADIENT==.T..
! LFCOBS  : .F. = Switch to compute the forecast sensitivity with respect to
! LFCOBSTEST  : .F. = Switch to compute the spectral and grid point norms during the temporal loop in cnt4tl/ad
! LCLDSINK : .T. = SWITCH FOR CLOUD SINK VARIABLE
! CTOPBGE : 5.= BG ERROR   CLOUD SINK VARIABLE
! CAMTBGE : 0.01 = BG ERROR   CLOUD SINK VARIABLE
! LTOVSCV : .T. = SWITCH FOR TOVS CONTROL VARIABLE
! LTOVSREP: .T. = Reproducible TOVS scalar product
! LJC     : .T. = ADD A JC TERM IN THE COST FUNCTION
! LJCDFI  : .T. = compute a digital formulation of a Jc term
! LUSEJCDFI: .T. = a Jc-DFI term is added to the cost function when
! LINITCV : .T. = include initial condition in 4D-Var control variable.
! LMODERR : .T. = include model error term in 4D-Var.
! LVARBC  : .T. = include bias parameters in control variable.
! LPERTMRPAS : .T. = blacklist passive observations in perturbed EDA members
! LTRREF  : Controls what is read/write in array SPA5 or SPA7
!          .T.: reference trajectory  ---> SPA7
!          .F.: current trajectory    ---> SPA5
! LREFINC : Informs on the nature of the reference SPA7 content (with LFCOBS).
!          .T. SPA7 contains an increment
!              (differences between 2 historic fieldsets)
!          .F. SPA7 contains an historic fieldset.
! LZOWA   : .T. = DIAGNOSTIC OF COST FUNCTION BY ZONAL WAVE NUMBER
!           .F. = DIAGNOSTIC OF COST FUNCTION BY TOTAL WAVE NUMBER
! LTWANA  : .T. = WRITE ON FA FILES THE CURRENT ANALYSIS
! LTWGRA  : .T. = WRITE ON FA FILES THE CURRENT GRADIENT
! LTWLCZ  : .T. = WRITE ON FA FILES THE CURRENT SV
! LTWINC  : .T. = WRITE ON FA FILES THE CURRENT INCREMENT
! LTWBGV  : .T. = WRITE ON FA FILES A RANDOM BACKGROUND VECTOR
! LTWCGL  : .T. = WRITE ON FA FILES AN EIGENVECTOR OF THE HESSIAN
! LGRASCAL: .T. = WRITE GRADIENT WITH RESPECT TO SCALP INNER PRODUCT.
! LSKIPMIN: .T. = FORCE SKIPPING OF MINIMIZATION (For diagnostics)
! L_ABS_CONVERGENCE: .T. => convergence criterion is the absolute
!                           value of the final gradient.
! L_INFO_CONVERGENCE: .T. => convergence criterion is relative increase
!                            in information (DFI) per iteration.
! NITER   : MAX NUMBER OF ITERATIONS
! NSIMU   : MAX NUMBER OF SIMULATIONS
! NITER_MIN : MINIMUM NUMBER OF ITERATIONS (NB: CONGRAD ONLY)
! NINFRA  : EXPECTED FRACTION OF THE DIMINUTION OF THE COST FUNCTION
! RCVGE   : Convergence of the minimization is judged to have been achieved
!           when:
!             IF (L_ABS_CONVERGENCE) THEN
!               The norm of the gradient is less than RCVGE.
!              (NB: The value specified by RCVGE is not the true norm of the
!               proper cost function. It is the square of twice the norm.
!               I.e. it is the value reported as FINAL GRADIENT.)
!             ELSEIF (L_INFO_CONVERGENCE) THEN
!               The increase in Jb compared with the preceding iteration
!               is less than RCVGE. (Note that Jb represents "degrees
!               of freedom for signal", DFI, so the change in Jb per
!               iteration is a measure of the amount of information added
!               per iteration.
!             ELSE
!               The norm of the gradient is reduced by a factor of RCVGE.
!             ENDIF
! NMIMP   : CONTROLS THE IMPRESSIONS OF THE OPTIMIZER
! NSIM4D  : COUNTER (number of calls to SIM4D)
! NITER4D : COUNTER (number of iterations of the minimisation)
! NSIM4DL : VALUE OF NSIM4D AT LAST SIMULATION
! NDIAG   : Diagnostic management (switched on only for M1GCX)
! RDX     : starting point for GRTEST (generally 1.E-15)
! RXMIN   : minimum distance in the sup-norm distinguishable by the optimizer
! ALPHAG  : weighting coefficient for the weak constraint term Jc.
! ALPHAV  : weighting coefficient to scale the Jr cost function & gradient
! NOANEF  : Number of analysis error fields
! NBGVECS : Number of random background vectors to be generated.
! LBGTRUNC: .T. ===> filter the randomization estimate of sigma_b.
! NBGTRUNC: wavenumber above which sigma_b coefficients are zero.
! LBGOBS  : .T. ===> Compute Bg errors for observed quantities.
! LBGM    : .T. ===> Propagate bg-errors using TL model.
! LANOBS  : .T. ===> Compute and Store in ODB the analysis sensitivity to observations
! LWREINI : .T. ===> write perturbed initial file for LELAM and not LTEST
! LWRIBVEC: .T. ===> write random background vectors to files
! LWRIBVEC_FULL: .T. ===> write random background vectors to gribfull files
! LREABVEC: .T. ===> read  random background vectors from files
! LWRIEVEC: .T. ===> write eigenvectors of the Hessian to files
! NWRIEVEC: Maximum number of Hessian eigenvectors to write
! LEVECCNTL : .T. ===> write eigenvectors in control space (before CHAVARIN)
! LEVECGRIB : .T. ===> write eigenvectors in GRIB format
! N_DIAGS_CONVERGENCE: Convergence diagnostics level
!                      0 - No diagnostics
!                      1 - Print diagnostics if bad convergence
!                      2 - Print diagnostics
! N_DIAGS_EIGENVECS  : Save eigenvectors in GRIB (could replace LEVECGRIB)
!                      0 - Don't save eigenvectors
!                      1 - Save eigenvectors if bad convergence
!                      2 - Save eigenvectors
! LWRISIGB: .T. ===> write standard deviation of background error
! LWRISIGA: .T. ===> write standard deviation of analysis   error
! LWRISIGF: .T. ===> write standard deviation of forecast   error
! CFNSIGB : file to which standard deviation of background error is written
! CFNSIGA : file to which standard deviation of analysis   error is written
! CFNSIGF : file to which standard deviation of forecast   error is written
! MBGVEC  : COUNTER (number of random vector being constructed).
! NFGFCLEN: Length of first guess forecast (hours)
! LN1CG1  : .T. ===> Minimize using N1CG1
! LCONGRAD: .T. ===> Minimize using CONGRAD (cost function must be quadratic)
! L3DFGAT : .T. ===> Use 3d FGAT (4dVar with TL model=identity in minimization)
! NPRECO  : 0 ===> N1CG1 : unpreconditioned CG
!         : 2 ===> N1CG1 : L_BFGS preconditioning
! NBFGSB  : 0 ===> N1CG1 : don't build any preconditioner
!         : 2 ===> N1CG1 : build a L_BFGS preconditioner
! N1IMP   : Printing level
! NSELECT : Selection of the pairs to build the L_BFGS preconditioner
!         : 0 ===> FIFO strategy (last pairs saved)
!         : 1 ===> Uniform selection (pairs are distributed uniformly through CG run)
!         : 2 ===> Selection by the Rayleigh quotient
! ZEPSNEG : N1CG1 : control the positivity of the Hessian during the minimization
! LAMV_REASSIGN_PASSIVE : Passive calculation of predictors for AMV bias correction
! LAMV_HEIGHT_ADJUST: Reassign heights of AMV diagnosed above the model cloud
! LREO3_BCOR: Switch for ozone bias correction
! LCH4_BCOR: Switch for methane bias correction

! L_GUESS_RUNTIME     : .T. => estimate when the minimization will end.
! NITER_GUESS_RUNTIME : Estimate when the job will end at iteration NITER_GUESS_RUNTIME...
! NFREQ_GUESS_RUNTIME : ...and every NFREQ_GUESS_RUNTIME iterations thereafter.
! NMEM_GUESS_RUNTIME  : Base the estimate on the NMEM_GUESS_RUNTIME most recent iterations.

! FILTERFACTOR : Filtering factor applied to trajectory fields
! FILTEREXPO: Filtering exponent applied to trajectory fields
! FILTERRESOL: Filtering of trajectory fields only applied for this resolution and above

! LJBIMPACT: Obs impact on Jb+Jq
! LJCIMPACT: Obs impact on Jc
! LMONITOR_FCDEPAR: configuration for forecast departure monitoring
! ======== FORECAST DEPARTURE MONITORING============================
! NUPTRA_RANGE : index of forecat timerange
!========= SIMULATED OBSERVATION EVENTS MANAGEMENT==================

! NFRREF  : frequency of reference observation events
! NREFTS  : array containing observation events steps
!     EXPLANATION :
!     1) IF NREFTS(0)=0 ACTION IF MOD(JSTEP,NFRREF)=0
!     2) IF NREFTS(0)>0 NREFTS(0) SIGNIFICANT NUMBERS IN
! NREFTS ARE THEN CONSIDERED AND :
!       ACTION FOR JSTEP=NREFTS(.)*NFRREF

!========= WRITES CURRENT ANALYSIS ON FILES =================

! NFRANA  : frequency of writes, relative to the number of simulations
! NANATS  : array containing write events steps
!     EXPLANATION :
!     1) IF NANATS(0)=0 ACTION IF MOD(JSTEP,NFRANA)=0
!     2) IF NANATS(0)>0 NANATS(0) SIGNIFICANT NUMBERS IN
! NANATS ARE THEN CONSIDERED AND :
!       ACTION FOR NSIM4D=NANATS(.)*NFRANA
! MSIME   : number of the current emsemble member

!========= WRITES CURRENT GRADIENT ON FILES =================

! NFRGRA  : frequency of writes, relative to the number of simulations
! NGRATS  : array containing write events steps
!     EXPLANATION : As for NFRANA and NANATS above

!========= ASSIMILATION WITH THE T.L. MODEL (nconf 131) ==========
! NUPTRA  : number of updates of the trajectory during the minimisation.
! MUPTRA  : Maximum number of updates of the trajectory during the minimisation

!========= Combined conjugate-gradient and Lanczos algorithm =====
! LAVCGL  : .T. ====> use combined conjugate-gradient / Lanczos scheme
! LMPCGL  : .T. ====> precondition conjugate-gradient minimization
! R_MAX_CNUM_PC : Maximum allowed condition number for the preconditioner
! CFNPCV  : filename prefix for files containing preconditioner vectors
! NPCVECS : number of vectors which make up the preconditioner
! EVBCGL  : max relative error in eigenvalues of written-out vectors
! MCGLVEC : COUNTER identifies the vector being written
! LGCV    : .T. => Calculate the Generalized Cross Validation function
! NITERGCV: Number of iterations for trace calculation
! LGCVJO  : (Internal) calculate Jo gradient for random departues.
! GCVJO   : (Internal) Jo at analysis point
! LFDBERR :  .T. => write bg and an error estimates to FDB

! L_INFO_CONTENT      : .T. => Calculate information content (degrees of freedom for signal, etc.)
! N_INFO_CONTENT_METHOD : 1 => Use Bai et al's algorthm
!                         2 => Use analysis-difference method
! N_INFO_CONTENT_SEED : Random number seed for randomized trace estimate.
! LPROPTL: .T. => 3DFGAT through 4DVAR looping with backward TL propagation
! LAEOLUSAMD  : .F. = include Aeolus production of auxiliary met data
! LBACKGE    : .T. = compute background error variance from an ensemble for full variables
! LBACKGECV  : .T. = compute background error variance from an ensemble for control variables
! LBACKGERENORM : .F. = compute renormalisation coefficients induced by C matrix in wavelet space
! LCONSTANTZFCE : .F. = specify that the variance is constant (and equal to 1) when doing a randomisation of B matrix
! LUSEWAVRENORM : .F. = use renormalisation coefficients to renormalise the correlation matrix
! LWRISB_VPROF:.T. = write vertical profile of sigmab from writesd/fltbgerr (in bgvecs and bgevecs)
! LFAMEMBERS : .T. = read ensemble members in FA format (no need for femars/readvec anymore)
! LTOY42     : Small resolution 4DVAR
! LSUSPQLIM  : Disable SUSPQLIM calculation
! LSPINT     : .T. = the low resolution fg and analyses will be directly read into
!                    the high resolution run and padded with zeros (ECMWF only)
! NDATE_TIME_WINDOW_END  integer date and time of 4dvar window end
! LBGPERT  : .T. = RANDBG perturbations
! DELTA : scaling factor of the size of random perturbations
! LCHRESINCR : .F. : Read control vector, change resolution in control vector
!                  space, transform to analysis incremnt, add to background
!                  and write out (HIRLAM)
! LINC_TOVSCV : update skin temperature estimate starting from estimate in previous outer loop
! LUSE_EDA_SKT: use skin temperature background error calculated from the EDA
! LECV : extended control variable
! L_FGNOTBG : .T. = First minimization starts with first guess not equal to background
! L_BGREC : .T. = Re-centre background (in EDA, re-centre perturbed background on control background)
! L_FGREC : .T. = Re-centre first-guess (in EDA, add total increment from control minimisation NFGREC_MIN)
!                 Re-centre by adding increments from control, cannot use control first-guess directly
!                 because background of vontrol and EDa members differ.
! NFGREC_MIN : Number of control minimisation to take first-guess from. We take the sum
!              of the control analysis increments 0 to NFGREC_MIN.
! L_COMPBGDEP : .T. = Background departures calculated from a separate background trajectory
! LWRICV    : .T. = Write out control vector to file
INTEGER(KIND=JPIM), PARAMETER :: JPNRFT=40
INTEGER(KIND=JPIM),PROTECTED :: NREFTS(0:JPNRFT)
INTEGER(KIND=JPIM),PROTECTED :: NANATS(0:JPNRFT)
INTEGER(KIND=JPIM),PROTECTED :: NGRATS(0:JPNRFT)
CHARACTER (LEN = 80),PROTECTED ::  CFNSIGB
CHARACTER (LEN = 80),PROTECTED ::  CFNSIGA
CHARACTER (LEN = 80),PROTECTED ::  CFNSIGF
CHARACTER (LEN = 80),PROTECTED ::  CFNPCV
INTEGER(KIND=JPIM),PROTECTED :: NITER
INTEGER(KIND=JPIM),PROTECTED :: NITER_MIN
INTEGER(KIND=JPIM),PROTECTED :: NSIMU
INTEGER(KIND=JPIM) :: NINFRA
INTEGER(KIND=JPIM),PROTECTED :: NMIMP
!INTEGER(KIND=JPIM),PROTECTED :: NSIM4D
INTEGER(KIND=JPIM) :: NDIAG
INTEGER(KIND=JPIM),PROTECTED :: NFRREF
INTEGER(KIND=JPIM),PROTECTED :: NFRANA
INTEGER(KIND=JPIM),PROTECTED :: NFRGRA
INTEGER(KIND=JPIM),PROTECTED :: NSIM4DL
INTEGER(KIND=JPIM) :: MSIME
!INTEGER(KIND=JPIM),PROTECTED :: NUPTRA
INTEGER(KIND=JPIM) :: NOBSCURE_NUMBER ! The name explains all
INTEGER(KIND=JPIM),PROTECTED :: MUPTRA
INTEGER(KIND=JPIM) :: NPCVECS ! Changed in PREPPCM
INTEGER(KIND=JPIM) :: NBGVECS ! Changed in BGEVECS
INTEGER(KIND=JPIM) :: NHEVECS ! Changed in XFORMEV
INTEGER(KIND=JPIM) :: NSSBGV_COUNT ! Counter
INTEGER(KIND=JPIM) :: NSSHEV_COUNT ! Counter
INTEGER(KIND=JPIM) :: MCGLVEC !Counter
INTEGER(KIND=JPIM) :: MBGVEC ! Counter
INTEGER(KIND=JPIM) :: NOANEF ! Set in SUANEBUF
INTEGER(KIND=JPIM),PROTECTED :: NFGFCLEN
!INTEGER(KIND=JPIM),PROTECTED :: NITER4D
INTEGER(KIND=JPIM),PROTECTED :: NBGTRUNC
INTEGER(KIND=JPIM),PROTECTED :: NITERGCV
INTEGER(KIND=JPIM),PROTECTED :: NWRIEVEC
INTEGER(KIND=JPIM),PROTECTED :: N_DIAGS_CONVERGENCE
INTEGER(KIND=JPIM),PROTECTED :: N_DIAGS_EIGENVECS
INTEGER(KIND=JPIM),PROTECTED :: NDATE_TIME_WINDOW_END
REAL(KIND=JPRB),PROTECTED :: RDX
REAL(KIND=JPRB),PROTECTED :: ALPHAG
REAL(KIND=JPRB),PROTECTED :: ALPHAV
REAL(KIND=JPRB),PROTECTED :: RXMIN
REAL(KIND=JPRB),PROTECTED :: EVBCGL
REAL(KIND=JPRB),PROTECTED :: RCVGE
REAL(KIND=JPRB),PROTECTED :: R_NORM_REDUCTION_ABORT_LEVEL
REAL(KIND=JPRB) :: GCVJO
REAL(KIND=JPRB),PROTECTED :: ZEPSNEG
REAL(KIND=JPRB),PROTECTED :: R_MAX_CNUM_PC
REAL(KIND=JPRB),PROTECTED :: RTOL_CHECK_GRADIENT
REAL(KIND=JPRB),PROTECTED :: FILTERFACTOR
REAL(KIND=JPRB),PROTECTED :: FILTEREXPO
REAL(KIND=JPRB),PROTECTED :: FILTERRESOL
REAL(KIND=JPRB),PROTECTED :: CTOPBGE
REAL(KIND=JPRB),PROTECTED :: CAMTBGE
REAL(KIND=JPRB),PROTECTED :: RCOEFCO2         ! Cooef. std. CO2 MACC
REAL(KIND=JPRB),PROTECTED :: RCOEFCH4         ! Cooef. std. CH4 MACC
REAL(KIND=JPRB),PROTECTED :: RCOEFCO          ! Cooef. std. CO CHEM
REAL(KIND=JPRB),PROTECTED :: RCOEFNO2         ! Cooef. std. NO2 CHEM
REAL(KIND=JPRB),PROTECTED :: RCOEFGO3         ! Cooef. std. O3 CHEM
REAL(KIND=JPRB),PROTECTED :: DELTA
LOGICAL :: LFCOBS
LOGICAL,PROTECTED :: LFCOBSTEST
LOGICAL,PROTECTED :: LTEST
LOGICAL,PROTECTED :: LREPRO4DVAR
LOGICAL,PROTECTED :: LSLADREP
LOGICAL :: LTRREF
LOGICAL,PROTECTED :: LREFINC
LOGICAL,PROTECTED :: LZOWA
LOGICAL :: LJC
LOGICAL,PROTECTED :: LINITCV
LOGICAL :: LENSCV
LOGICAL,PROTECTED :: LMODERR
LOGICAL :: LVARBC
LOGICAL,PROTECTED :: LPERTMRPAS
LOGICAL :: LTWANA
LOGICAL :: LTWGRA
LOGICAL :: LTWLCZ
LOGICAL,PROTECTED :: LAVCGL
LOGICAL :: LMPCGL
LOGICAL :: LTWINC
LOGICAL,PROTECTED :: LGRASCAL
LOGICAL,PROTECTED :: LSKIPMIN
LOGICAL,PROTECTED :: L_ABS_CONVERGENCE
LOGICAL,PROTECTED :: L_INFO_CONVERGENCE
LOGICAL :: LTWBGV ! Set in BGVECS
LOGICAL :: LTWCGL ! Set in XFORMEV
LOGICAL,PROTECTED :: LWRIBVEC
LOGICAL,PROTECTED :: LWRIBVEC_FULL
LOGICAL,PROTECTED :: LREABVEC
LOGICAL,PROTECTED :: LWRIEVEC
LOGICAL,PROTECTED :: LEVECCNTL
LOGICAL,PROTECTED :: LEVECGRIB
LOGICAL,PROTECTED :: LWRISIGB
LOGICAL,PROTECTED :: LWRISIGA
LOGICAL,PROTECTED :: LWRISIGF
LOGICAL,PROTECTED :: LBGTRUNC
LOGICAL,PROTECTED :: LBGOBS
LOGICAL,PROTECTED :: LBGM
LOGICAL,PROTECTED :: LANOBS
LOGICAL :: L_CHECK_CONVERGENCE
LOGICAL,PROTECTED :: L_CHECK_GRADIENT
LOGICAL :: LJCDFI ! Turned off and on like a jojo
LOGICAL,PROTECTED :: LUSEJCDFI
LOGICAL :: LTOVSCV
LOGICAL,PROTECTED :: LTOVSREP
LOGICAL,PROTECTED :: LGCV
LOGICAL,PROTECTED :: L_FGNOTBG
LOGICAL,PROTECTED :: L_BGREC
LOGICAL,PROTECTED :: L_FGREC
INTEGER(KIND=JPIM),PROTECTED :: NFGREC_MIN
LOGICAL,PROTECTED :: L_COMPBGDEP
LOGICAL,PROTECTED :: LWRICV
LOGICAL :: LGCVJO ! FORECAST_ERROR
LOGICAL,PROTECTED :: LWREINI
LOGICAL,PROTECTED :: LN1CG1
LOGICAL,PROTECTED :: LCONGRAD
LOGICAL :: L3DFGAT ! In TLPROP
LOGICAL,PROTECTED :: LFDBERR
LOGICAL,PROTECTED :: LAMV_REASSIGN_PASSIVE
LOGICAL,PROTECTED :: LAMV_HEIGHT_ADJUST
LOGICAL,PROTECTED :: LREO3_BCOR
LOGICAL,PROTECTED :: LCH4_BCOR
LOGICAL,PROTECTED :: L_INFO_CONTENT
LOGICAL,PROTECTED :: LCLDSINK
LOGICAL,PROTECTED :: L_GUESS_RUNTIME
LOGICAL,PROTECTED :: LPROPTL
LOGICAL,PROTECTED :: LAEOLUSAMD
LOGICAL :: LBACKGE
LOGICAL,PROTECTED :: LBACKGECV
LOGICAL,PROTECTED :: LBACKGERENORM
LOGICAL,PROTECTED :: LUSEWAVRENORM
LOGICAL,PROTECTED :: LCONSTANTZFCE
LOGICAL,PROTECTED :: LWRISB_VPROF
LOGICAL :: LDIAG_LCT ! SET in BGVECS and BGEVECS
LOGICAL,PROTECTED :: LFAMEMBERS
LOGICAL,PROTECTED :: LTOY42
LOGICAL,PROTECTED :: LSUSPQLIM
LOGICAL,PROTECTED :: LSPINT
LOGICAL,PROTECTED :: LJBIMPACT = .FALSE.
LOGICAL,PROTECTED :: LJCIMPACT = .FALSE.
LOGICAL :: LMONITOR_FCDEPAR = .FALSE. ! Changed in shuffle used by Odbtools in ODB - uncomprehensible
LOGICAL,PROTECTED :: LENDA = .FALSE.
LOGICAL,PROTECTED :: LBGPERT = .false.
LOGICAL,PROTECTED :: LCHRESINCR
LOGICAL,PROTECTED :: LINC_TOVSCV
LOGICAL,PROTECTED :: LUSE_EDA_SKT= .false.
LOGICAL,PROTECTED :: LECV
INTEGER(KIND=JPIM),PROTECTED :: NUPTRA_RANGE = 0
INTEGER(KIND=JPIM),PROTECTED :: N_INFO_CONTENT_METHOD
INTEGER(KIND=JPIM),PROTECTED :: N_INFO_CONTENT_SEED
INTEGER(KIND=JPIM),PROTECTED :: N1IMP
INTEGER(KIND=JPIM),PROTECTED :: NPRECO
INTEGER(KIND=JPIM),PROTECTED :: NBFGSB
INTEGER(KIND=JPIM),PROTECTED :: NSELECT
INTEGER(KIND=JPIM),PROTECTED :: NITER_GUESS_RUNTIME
INTEGER(KIND=JPIM),PROTECTED :: NFREQ_GUESS_RUNTIME
INTEGER(KIND=JPIM),PROTECTED :: NMEM_GUESS_RUNTIME
!     ------------------------------------------------------------------
CONTAINS
SUBROUTINE SETUP_VAR(KULOUT)

!**** *SETUP_VAR*   - Routine to initialize variational flags common

!     Purpose.
!     --------
!        Initialize variational control modules: YOMSENS, YOMVRTL, YOMCOSJB, YOMVAR,
!         and also some variables of YOMJQ, YOMMODERR, YOMVCGL, YOMDIAGVAR
!        Reads namelist NAMSENS, NAMVRTL and NAMVAR.

!**   Interface.
!     ----------
!        *CALL* *SETUP_VAR(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!         output modules: YOMSENS, YOMVRTL, YOMCOSJB, YOMVAR
!                         a subset of YOMJQ, YOMMODERR,  YOMVCGL, YOMDIAGVAR

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Authors.
!     -------
!      Philippe Courtier  and Jean-Noel Thepaut *DMN/ECMWF*  90-12-01

!     Modifications.
!     --------------
!      Modified by N. Bormann   : 01-09-21 Initialise lamv_reassign_passive
!      Modified by A. Dethof    : 03-04240 Initialise lreo3_bcor
!      Modified by D. Dee       : 03-10-17 Initialise lvarbc
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y.Tremolet    01-Apr-2004 Setup model error stats
!      G. Radnoti    03-12-2004 : initialize lproptl
!      D. Tan        14-Mar-2005: Initialise laeolusamd, laeolusl2bp
!      A. Benedetti  14-Mar-2005: Initialise laerod
!      G. Desroziers and K. Yessad (sept 2005):
!       - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!       - adapt option LTRAJHR to METEO-FRANCE configurations.
!      Modified by R. Engelen    : Mar 2008 Initialise lch4_bcor
!      G. Desroziers 28-Feb-2008: Possible computation of ensemble bg error variance
!      Y.Tremolet    27-Nov-2008 Setup for long windows
!      G. Desroziers 22-Dec-2008: Enable transf. of ARPEGE file in GRIB format (to be used in femars)
!      H.Varella     15-Nov-2011: Option LWRIBVEC_FULL included
!      M. Rennie     13-Apr-2012: Remove LAEOLUSL2BP
!      K. Yessad (oct 2013): cleaning, re-write in a readable way.
!      LF. Meunier   29-Oct-2013: Remove L_OPENMP_CV
!      K. Yessad (July 2014): use L4DVAR, in order to avoid use of NSTOP.
!      Y. Michel     10-Mar-2015 Local Correlation Tensor
!      V.Chabot      27-Janv-2016 Renormalisation coefficient for wavelet correlation matrix
!      N. Bormann    01-Aug-2016 CVarBC
!      M. Hamrud     May 2018 module routine setup_var, with contents previously in suvar
!      K. Lean       10-Jan-2020 Option LAMV_HEIGHT_ADJUST added
!     -------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMCT0   , ONLY : NCONF, LOBSC1, LBACKG, LFDBOP, LECMWF, LELAM, L4DVAR, L_OOPS
USE YOMMP0   , ONLY : LOPT_SCALAR
USE YOMDIAGVAR,ONLY : DIAG_4DVAR
USE YOMVCGL  , ONLY : NSAVEEV, NSAVEPC
USE YOMMODERR, ONLY : N_COUPLED_WINDOWS
USE YOMJQ    , ONLY : LSTATMERR
USE YOMSENS  , ONLY : LGRVOL   ,NJROPT   ,LBSENS
USE YOMVRTL  , ONLY : L131TL   ,LTLINT   ,LOBSTL   ,LDRYTL   ,LIDMODEL
USE YOMCOSJB , ONLY : LJBZERO, LJPZERO, LJHZERO, LJLZERO, LJTZERO
USE ALGORITHM_STATE_MOD, ONLY : SETUP_ALGORITHM_STATE,SET_NUPTRA,SET_MUPTRA,GET_NSIM4D,SET_NSIM4D,&
 &                              GET_NUPTRA,GET_MUPTRA,GET_ALGOR_TYPE
USE YOMVRTLX, ONLY : LMINI, L801TL
!     -------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT

!     -------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: J
INTEGER(KIND=JPIM) :: NUPTRA
REAL(KIND=JPRB) :: Z
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------------------------

#include "abor1.intfb.h"


!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMVAR:SETUP_VAR',0,ZHOOK_HANDLE)


IF (LHOOK) CALL DR_HOOK('YOMVAR:SETUP_VAR',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_VAR

END MODULE YOMVAR
