! radiation_matrix.F90 - SPARTACUS matrix operations
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
!   2018-10-15  R. Hogan    Added fast_expm_exchange_[23]
!   2020-12-15  P. Ukkonen  Optimized kernels
!
! This module provides the neccessary mathematical functions for the
! SPARTACUS radiation scheme: matrix multiplication, matrix solvers
! and matrix exponentiation, but (a) multiple matrices are operated on
! at once with array access indended to facilitate vectorization, and
! (b) optimization for 2x2 and 3x3 matrices.  There is probably
! considerable scope for further optimization. Note that this module
! is not used by the McICA solver.

module radiation_matrix

  use parkind1, only : jprb

  implicit none
  public

  ! Codes to describe sparseness pattern, where the SHORTWAVE
  ! pattern is of the form:
  ! (x x x)
  ! (x x x)
  ! (0 0 x)
  ! where each element may itself be a square matrix.
  integer, parameter :: IMatrixPatternDense     = 0
  integer, parameter :: IMatrixPatternShortwave = 1

  private :: solve_vec_2, solve_vec_3, solve_mat_2, &
       &     solve_mat_3, lu_factorization, lu_substitution, solve_mat_n, &
       &     diag_mat_right_divide_3

! Allow size of inner dimension (number of g-points) to be known at compile time for some routines
#ifdef NG_SW
    integer, parameter, private :: ng_sw = NG_SW
    integer, parameter, private :: ng_sww = 6*NG_SW
#else
#define ng_sw ng_sw_in
#define ng_sww ng_sw_in
#endif
#ifdef NG_LW
    integer, parameter, private :: ng_lw = NG_LW
#else
#define ng_lw ng_lw_in
#endif

  interface fast_expm_exchange
    module procedure fast_expm_exchange_2, fast_expm_exchange_3, fast_expm_exchange_3_orig
  end interface fast_expm_exchange

contains

  ! --- MATRIX-VECTOR MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A as n m-by-m square matrices (with the n dimension varying
  ! fastest) and b as n m-element vectors, and perform matrix-vector
  ! multiplications on first iend pairs
  function mat_x_vec(n,iend,m,A,b,do_top_left_only_in)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)                   :: n, m, iend
real(jprb), intent(in), dimension(:,:,:) :: A
real(jprb), intent(in), dimension(:,:)   :: b
logical,    intent(in), optional         :: do_top_left_only_in
real(jprb),             dimension(iend,m):: mat_x_vec









end function mat_x_vec

  pure function mat_x_vec_3(ng,A,b)
integer,    intent(in)                    :: ng
real(jprb), intent(in), dimension(ng,3,3) :: A
real(jprb), intent(in), dimension(ng,3)   :: b
real(jprb),             dimension(ng,3):: mat_x_vec_3


end function mat_x_vec_3

  pure function mat_x_vec_3_lw(ng_lw_in,A,b)
integer,    intent(in)                       :: ng_lw_in
real(jprb), intent(in), dimension(ng_lw,3,3) :: A
real(jprb), intent(in), dimension(ng_lw,3)   :: b
real(jprb),             dimension(ng_lw,3):: mat_x_vec_3_lw


end function mat_x_vec_3_lw

  pure function mat_x_vec_3_sw(ng_sw_in,A,b,do_top_left_only_in)
integer,    intent(in)                        :: ng_sw_in
real(jprb), intent(in), dimension(ng_sw,3,3)  :: A
real(jprb), intent(in), dimension(ng_sw,3)    :: b
real(jprb),             dimension(ng_sw,3)    :: mat_x_vec_3_sw
logical,    intent(in), optional              :: do_top_left_only_in




end function mat_x_vec_3_sw

  !---------------------------------------------------------------------
  ! Treat A as an m-by-m square matrix and b as n m-element vectors
  ! (with the n dimension varying fastest), and perform matrix-vector
  ! multiplications on first iend pairs
  function singlemat_x_vec(n,iend,m,A,b)

integer,    intent(in)                    :: n, m, iend
real(jprb), intent(in), dimension(m,m)    :: A
real(jprb), intent(in), dimension(:,:)    :: b
real(jprb),             dimension(iend,m) :: singlemat_x_vec







end function singlemat_x_vec

  pure function singlemat_x_vec_sw(ng_sw_in,A,b)
integer,    intent(in)                      :: ng_sw_in
real(jprb), intent(in), dimension(3,3)      :: A
real(jprb), intent(in), dimension(ng_sw,3)  :: b
real(jprb),             dimension(ng_sw,3)  :: singlemat_x_vec_sw


end function singlemat_x_vec_sw

  pure function singlemat_x_vec_lw(ng_lw_in,A,b)
integer,    intent(in)                      :: ng_lw_in
real(jprb), intent(in), dimension(3,3)      :: A
real(jprb), intent(in), dimension(ng_lw,3)  :: b
real(jprb),             dimension(ng_lw,3)  :: singlemat_x_vec_lw


end function singlemat_x_vec_lw

  ! --- SQUARE MATRIX-MATRIX MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A and B each as n m-by-m square matrices (with the n
  ! dimension varying fastest) and perform matrix multiplications on
  ! all n matrix pairs
  function mat_x_mat(n,iend,m,A,B,i_matrix_pattern)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)                      :: n, m, iend
integer,    intent(in), optional            :: i_matrix_pattern
real(jprb), intent(in), dimension(:,:,:)    :: A, B
real(jprb),             dimension(iend,m,m) :: mat_x_mat










end function mat_x_mat

  pure subroutine mat_x_mat_3_sw(ng_sw_in,A,B,C)
integer,    intent(in)                            :: ng_sw_in
real(jprb), intent(in),     dimension(ng_sw,3,3)  :: A, B
real(jprb), intent(inout),  dimension(ng_sw,3,3)  :: C


end subroutine mat_x_mat_3_sw

  pure subroutine mat_x_mat_3_lw(ng_lw_in,A,B,C)
integer,    intent(in)                            :: ng_lw_in
real(jprb), intent(in),     dimension(ng_lw,3,3)  :: A, B
real(jprb), intent(inout),  dimension(ng_lw,3,3)  :: C


end subroutine mat_x_mat_3_lw

  pure function mat_x_mat_6(ng,A,B)
integer,    intent(in)                    :: ng
real(jprb), intent(in), dimension(ng,6,6) :: A, B
real(jprb),             dimension(ng,6,6) :: mat_x_mat_6


end function mat_x_mat_6

  !---------------------------------------------------------------------
  ! Treat A as an m-by-m matrix and B as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function singlemat_x_mat(n,iend,m,A,B)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)                      :: n, m, iend
real(jprb), intent(in), dimension(m,m)      :: A
real(jprb), intent(in), dimension(:,:,:)    :: B
real(jprb),             dimension(iend,m,m) :: singlemat_x_mat







end function singlemat_x_mat

  pure function singlemat_x_mat_3_sw(ng_sw_in,A,B)
integer,    intent(in)                        :: ng_sw_in
real(jprb), intent(in), dimension(3,3)        :: A
real(jprb), intent(in), dimension(ng_sw,3,3)  :: B
real(jprb),             dimension(ng_sw,3,3)  :: singlemat_x_mat_3_sw



end function singlemat_x_mat_3_sw

  !---------------------------------------------------------------------
  ! Treat B as an m-by-m matrix and A as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function mat_x_singlemat(n,iend,m,A,B)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)                      :: n, m, iend
real(jprb), intent(in), dimension(:,:,:)    :: A
real(jprb), intent(in), dimension(m,m)      :: B
real(jprb),             dimension(iend,m,m) :: mat_x_singlemat







end function mat_x_singlemat

  pure function mat_x_singlemat_3_sw(ng_sw_in, A,B)
integer,    intent(in)                        :: ng_sw_in
real(jprb), intent(in), dimension(ng_sw,3,3)  :: A
real(jprb), intent(in), dimension(3,3)        :: B
real(jprb),             dimension(ng_sw,3,3)  :: mat_x_singlemat_3_sw


end function mat_x_singlemat_3_sw

  !---------------------------------------------------------------------
  ! Compute I-A*B where I is the identity matrix and A & B are n
  ! m-by-m square matrices
  function identity_minus_mat_x_mat(n,iend,m,A,B,i_matrix_pattern)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)                   :: n, m, iend
integer,    intent(in), optional         :: i_matrix_pattern
real(jprb), intent(in), dimension(:,:,:) :: A, B
real(jprb),             dimension(iend,m,m) :: identity_minus_mat_x_mat







end function identity_minus_mat_x_mat

pure subroutine identity_minus_mat_x_mat_3_sw(ng_sw_in,A,B,C)
integer,    intent(in)                            :: ng_sw_in
real(jprb), intent(in),     dimension(ng_sw,3,3)  :: A, B
real(jprb), intent(inout),  dimension(ng_sw,3,3)  :: C




end subroutine identity_minus_mat_x_mat_3_sw

  subroutine identity_minus_mat_x_mat_3_lw(ng_lw_in,A,B,C)
integer,    intent(in)                            :: ng_lw_in
real(jprb), intent(in),     dimension(ng_lw,3,3)  :: A, B
real(jprb), intent(inout),  dimension(ng_lw,3,3)  :: C




end subroutine identity_minus_mat_x_mat_3_lw

  !---------------------------------------------------------------------
  ! Replacement for matmul in the case that the first matrix is sparse
  function sparse_x_dense(sparse, dense)
real(jprb), intent(in) :: sparse(:,:), dense(:,:)
real(jprb) :: sparse_x_dense(size(sparse,1),size(dense,2))
 
 





end function sparse_x_dense


  ! --- REPEATEDLY SQUARE A MATRIX ---

  !---------------------------------------------------------------------
  ! Square m-by-m matrix "A" nrepeat times. A will be corrupted by
  ! this function.
  function repeated_square(m,A,nrepeat,i_matrix_pattern)
integer,    intent(in)           :: m, nrepeat
real(jprb), intent(inout)        :: A(m,m)
integer,    intent(in), optional :: i_matrix_pattern
real(jprb)                       :: repeated_square(m,m)





end function repeated_square

  pure subroutine repeated_square_sw_9(nrepeat,A,B)
integer,    intent(in)            :: nrepeat
real(jprb), intent(inout)         :: A(9,9)
real(jprb), intent(out)           :: B(9,9)


end subroutine repeated_square_sw_9


  ! --- SOLVE LINEAR EQUATIONS ---

  !---------------------------------------------------------------------
  ! Solve Ax=b to obtain x.  Version optimized for 2x2 matrices using
  ! Cramer's method: "A" contains n 2x2 matrices and "b" contains n
  ! 2-element vectors; returns A^-1 b.
  pure subroutine solve_vec_2(n,iend,A,b,x)
integer,    intent(in)  :: n, iend
real(jprb), intent(in)  :: A(:,:,:)
real(jprb), intent(in)  :: b(:,:)
real(jprb), intent(out) :: x(:,:)




end subroutine solve_vec_2


  !---------------------------------------------------------------------
  ! Solve AX=B to obtain X, i.e. the matrix right-hand-side version of
  ! solve_vec_2, with A, X and B all containing n 2x2 matrices;
  ! returns A^-1 B using Cramer's method.
  pure subroutine solve_mat_2(n,iend,A,B,X)
integer,    intent(in)  :: n, iend
real(jprb), intent(in)  :: A(:,:,:)
real(jprb), intent(in)  :: B(:,:,:)
real(jprb), intent(out) :: X(:,:,:)






end subroutine solve_mat_2


  !---------------------------------------------------------------------
  ! Solve Ax=b optimized for 3x3 matrices, using LU
  ! factorization and substitution without pivoting.
  pure subroutine solve_vec_3(n,iend,A,b,x)
integer,    intent(in)  :: n, iend
real(jprb), intent(in)  :: A(:,:,:)
real(jprb), intent(in)  :: b(:,:)
real(jprb), intent(out) :: x(:,:)

























end subroutine solve_vec_3

  pure function solve_vec_3_sw(ng_sw_in,A,b)
integer,    intent(in)  :: ng_sw_in
real(jprb), intent(in)  :: A(ng_sw,3,3)
real(jprb), intent(in)  :: b(ng_sw,3)
real(jprb) :: solve_vec_3_sw(ng_sw,3)




















end function solve_vec_3_sw

  ! Like solve_vec_3 but overwrite b with x
  pure subroutine solve_vec_3_lw_inplace(ng_lw_in,A,b)
integer,    intent(in)  :: ng_lw_in
real(jprb), intent(in)  :: A(ng_lw,3,3)
real(jprb), intent(inout)  :: b(ng_lw,3)





















end subroutine solve_vec_3_lw_inplace

  pure function solve_vec_3_lw(ng_lw_in,A,b)
integer,    intent(in)  :: ng_lw_in
real(jprb), intent(in)  :: A(ng_lw,3,3)
real(jprb), intent(in)  :: b(ng_lw,3)
real(jprb) :: solve_vec_3_lw(ng_lw,3)




















end function solve_vec_3_lw

  pure function solve_vec_3_ng(ng,A,b)
integer,    intent(in)  :: ng
real(jprb), intent(in)  :: A(ng,3,3)
real(jprb), intent(in)  :: b(ng,3)
real(jprb) :: solve_vec_3_ng(ng,3)




















end function solve_vec_3_ng

  !---------------------------------------------------------------------
  ! Solve AX=B optimized for 3x3 matrices, using LU factorization and
  ! substitution with no pivoting.
  pure subroutine solve_mat_3(n,iend,A,B,X)
integer,    intent(in)  :: n, iend
real(jprb), intent(in)  :: A(:,:,:)
real(jprb), intent(in)  :: B(:,:,:)
real(jprb), intent(out) :: X(:,:,:)
















end subroutine solve_mat_3


  !---------------------------------------------------------------------
  ! Return X = B A^-1 = (A^-T B)^T optimized for 3x3 matrices, where B
  ! is a diagonal matrix, using LU factorization and substitution with
  ! no pivoting.
  pure subroutine diag_mat_right_divide_3(n,iend,A,B,X)
integer,    intent(in)  :: n, iend
real(jprb), intent(in)  :: A(iend,3,3)
real(jprb), intent(in)  :: B(iend,3)
real(jprb), intent(out) :: X(n,3,3)















































end subroutine diag_mat_right_divide_3


  !---------------------------------------------------------------------
  ! Treat A as n m-by-m matrices and return the LU factorization of A
  ! compressed into a single matrice (with L below the diagonal and U
  ! on and above the diagonal; the diagonal elements of L are 1). No
  ! pivoting is performed.
  pure subroutine lu_factorization(n, iend, m, A, LU)
integer,    intent(in)  :: n, m, iend
real(jprb), intent(in)  :: A(:,:,:)
real(jprb), intent(out) :: LU(iend,m,m)






end subroutine lu_factorization

  pure subroutine lu_factorization_lw(n, A, LU)
integer,    intent(in)  :: n
real(jprb), intent(in)  :: A(n,6,6)
real(jprb), intent(out) :: LU(n,6,6)






end subroutine lu_factorization_lw

  !---------------------------------------------------------------------
  ! Treat LU as an LU-factorization of an original matrix A, and
  ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
  ! m-element vectors.
  pure subroutine lu_substitution(n,iend,m,LU,b,x)

integer,    intent(in) :: n, m, iend
real(jprb), intent(in) :: LU(iend,m,m)
real(jprb), intent(in) :: b(:,:)
real(jprb), intent(out):: x(iend,m)






end subroutine lu_substitution

! Optimized version, overwrite input vector b with output vector x
  pure subroutine lu_substitution_lw_inplace(n,LU,b) 
integer,    intent(in) :: n
real(jprb), intent(in) :: LU(n,6,6)
real(jprb), intent(inout) :: b(n,6)






































end subroutine lu_substitution_lw_inplace


  !---------------------------------------------------------------------
  ! Return matrix X where AX=B. LU, A, X, B all consist of n m-by-m
  ! matrices.
  pure subroutine solve_mat_n(n,iend,m,A,B,X)
integer,    intent(in) :: n, m, iend
real(jprb), intent(in) :: A(:,:,:)
real(jprb), intent(in) :: B(:,:,:)
real(jprb), intent(out):: X(iend,m,m)




end subroutine solve_mat_n

  ! Solve AX=B, where A, X and B consist of ng m-by-m matrices
  ! To minimize memory use, overwrite B with X, and A is corrupted
  pure subroutine solve_mat_inplace(n,m,A,B)
integer,    intent(in) :: n, m
real(jprb), intent(inout) :: A(n,m,m)
real(jprb), intent(inout) :: B(n,m,m)






end subroutine solve_mat_inplace

  pure subroutine solve_mat_3_inplace(ng,A,B)
integer,    intent(in)  :: ng
real(jprb), intent(in)  :: A(ng,3,3)
real(jprb), intent(inout)  :: B(ng,3,3)






end subroutine solve_mat_3_inplace

  pure subroutine solve_mat_3_sw(ng_sw_in,A,B,X)
integer,    intent(in)  :: ng_sw_in
real(jprb), intent(in)  :: A(ng_sw,3,3)
real(jprb), intent(in)  :: B(ng_sw,3,3)
real(jprb), intent(out) :: X(ng_sw,3,3)











end subroutine solve_mat_3_sw

  pure subroutine solve_mat_3_lw(ng_lw_in,A,B,X)
integer,    intent(in)  :: ng_lw_in
real(jprb), intent(in)  :: A(ng_lw,3,3)
real(jprb), intent(in)  :: B(ng_lw,3,3)
real(jprb), intent(out) :: X(ng_lw,3,3)









 
 





end subroutine solve_mat_3_lw

  !---------------------------------------------------------------------
  ! Solve Ax=b, where A consists of n m-by-m matrices and x and b
  ! consist of n m-element vectors. For m=2 or m=3, this function
  ! calls optimized versions, otherwise it uses general LU
  ! decomposition without pivoting.
  function solve_vec(n,iend,m,A,b)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in) :: n, m, iend
real(jprb), intent(in) :: A(:,:,:)
real(jprb), intent(in) :: b(:,:)
real(jprb)             :: solve_vec(iend,m)





end function solve_vec


  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of n m-by-m matrices. For m=2
  ! or m=3, this function calls optimized versions, otherwise it uses
  ! general LU decomposition without pivoting.
  function solve_mat(n,iend,m,A,B)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)  :: n, m, iend
real(jprb), intent(in)  :: A(:,:,:)
real(jprb), intent(in)  :: B(:,:,:)
real(jprb)              :: solve_mat(iend,m,m)




end function solve_mat

  pure subroutine mat_square_sw(ng_sw_in,nlev_b,A,C)




integer,    intent(in)                                :: ng_sw_in, nlev_b
real(jprb), intent(in),  dimension(ng_sw*nlev_b,9,9)  :: A

real(jprb), intent(out), dimension(ng_sw*nlev_b,9,9)  :: C





end subroutine mat_square_sw

  pure subroutine mat_x_mat_sw(ng_sw_in,nlev_b,A,B,C)
integer,    intent(in)                    :: ng_sw_in, nlev_b
real(jprb), intent(in), dimension(ng_sw*nlev_b,9,9) :: A, B

real(jprb), intent(out),dimension(ng_sw*nlev_b,9,9) :: C










end subroutine mat_x_mat_sw

  pure subroutine mat_x_mat_sw_repeats(ng_sw_in, nlev_b, A, B, C)
integer,    intent(in)                              :: ng_sw_in, nlev_b
real(jprb), intent(in), dimension(ng_sw*nlev_b,9,9) :: A, B

real(jprb), intent(out),dimension(ng_sw*nlev_b,9,9) :: C











 
 
     
end subroutine mat_x_mat_sw_repeats

  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of ng m-by-m matrices
  ! Overwrite B with X. A is corrupted
  pure subroutine solve_mat_sw(ng_sw_in,nlev_b,A,B)
integer,    intent(in)    :: ng_sw_in, nlev_b
real(jprb), intent(inout) :: A(ng_sw*nlev_b,9,9) 

real(jprb), intent(inout) :: B(ng_sw*nlev_b,9,9) 




  
 






















end subroutine solve_mat_sw

  ! --- MATRIX EXPONENTIATION ---

  !---------------------------------------------------------------------
  ! Perform matrix exponential of n m-by-m matrices stored in A (where
  ! index n varies fastest) using the Higham scaling and squaring
  ! method. The result is placed in A. This routine is intended for
  ! speed so is accurate only to single precision.  For simplicity and
  ! to aid vectorization, the Pade approximant of order 7 is used for
  ! all input matrices, perhaps leading to a few too many
  ! multiplications for matrices with a small norm.
  subroutine expm(n,iend,m,A,i_matrix_pattern)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)      :: n, m, iend
real(jprb), intent(inout)   :: A(n,m,m)
integer,    intent(in)      :: i_matrix_pattern







































end subroutine expm

  subroutine expm_lw(n,A)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)      :: n
real(jprb), intent(inout)   :: A(n,6,6)




































end subroutine expm_lw

  !---------------------------------------------------------------------
  ! Like expm, but optimized for the shortwave, which has
  ! a special matrix structure with zeros and repeated elements.
  ! Further assumes nreg = 3  =>  m = 9, and computations are done for all
  ! N elements, where N (the inner dimension) is the number of individual matrices
  ! For performance reasons, N is increased by batching several cloudy layers in the
  ! higher level code: N = nlev_b*ng_sw, where ng_sw is the number of g-points and nlev_b
  ! is the number of batched layers (this is capped at 6 for ecCKD, i.e. nlev_b=1..6)
  subroutine expm_sw(N,ng_sw_in,nlev_b,A)
use yomhook, only : lhook, dr_hook, jphook
integer,    intent(in)      :: N, ng_sw_in, nlev_b
real(jprb), intent(inout)   :: A(ng_sw*nlev_b,9,9)






































end subroutine expm_sw

  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 2x2 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b)
  !   ( a  -b)
  ! and a and b are assumed to be positive or zero.  The solution uses
  ! Putzer's algorithm - see the appendix of Hogan et al. (GMD 2018)
  subroutine fast_expm_exchange_2(n,iend,a,b,R)
use yomhook, only : lhook, dr_hook, jphook
integer,                      intent(in)  :: n, iend
real(jprb), dimension(n),     intent(in)  :: a, b
real(jprb), dimension(n,2,2), intent(out) :: R










end subroutine fast_expm_exchange_2


  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 3x3 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b   0)
  !   ( a -b-c  d)
  !   ( 0   c  -d)
  ! and a-d are assumed to be positive or zero.  The solution uses the
  ! diagonalization method and is a slight generalization of the
  ! solution provided in the appendix of Hogan et al. (GMD 2018),
  ! which assumed c==d.
  subroutine fast_expm_exchange_3_orig(n,iend,a,b,c,d,R)
use yomhook, only : lhook, dr_hook, jphook

integer,                      intent(in)  :: n, iend
real(jprb), dimension(n),     intent(in)  :: a, b, c, d
real(jprb), dimension(n,3,3), intent(out) :: R





































end subroutine fast_expm_exchange_3_orig
  subroutine fast_expm_exchange_3(ng_sw_in, R)
use yomhook, only : lhook, dr_hook, jphook

integer,                      intent(in)  :: ng_sw_in
real(jprb), dimension(ng_sww,3,3), intent(inout) :: R













#ifdef __INTEL_COMPILER



#endif
















end subroutine fast_expm_exchange_3

!  generic :: fast_expm_exchange => fast_expm_exchange_2, fast_expm_exchange_3


end module radiation_matrix
