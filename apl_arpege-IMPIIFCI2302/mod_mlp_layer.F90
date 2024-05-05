module mod_mlp_layer

  ! Defines the layer type and its methods.

  ! Changelog: 
  !   P.Ukkonen, 23.1.2022 : Change the code so that input layer is not included in layers (see mod_mlp)
  !   P.Ukkonen, 14.9.2023: Rename mod_layer to mod_mlp_layer to avoid overlap with OOPS NN modules originally built on same library

  use mod_mlp_activation
  use mo_rte_kind, only: sp
  ! use mod_random, only: randn

  implicit none

  private
  public :: layer_type

  type :: layer_type
    real(sp), allocatable :: b(:) ! biases
    real(sp), allocatable :: w(:,:) ! weights
    real(sp), allocatable :: w_transposed(:,:) ! weights
    !!  !$acc policy<copylayer> copyin(b, w, w_transposed)
    procedure(activation_interface),        pointer, nopass :: activation !=> null()
    procedure(bias_and_activation_interface),pointer, nopass :: bias_and_activation !=> null()

  contains
  
    procedure, public, pass(this) :: set_activation
  end type layer_type

  interface layer_type
    module procedure constructor
  end interface layer_type

contains

  type(layer_type) function constructor(last_size, this_size) result(layer)
    ! Layer class constructor. this_size is the number of neurons in the layer.
    ! next_size is the number of neurons in the next layer, used to allocate
    ! the weights.
    integer, intent(in) :: last_size, this_size
    allocate(layer % w(last_size, this_size))
    allocate(layer % w_transposed(this_size, last_size))
    allocate(layer % b(this_size))
    ! layer % b = randn(this_size)
  end function constructor

  subroutine set_activation(this, activation)
    ! Sets the activation function. Input string must match one of
    ! provided activation functions, otherwise it defaults to sigmoid.
    ! If activation not present, defaults to sigmoid.
    class(layer_type), intent(in out) :: this
    character(len=*), intent(in) :: activation
    select case(trim(activation))
      case('gaussian')
        this % activation => gaussian
        ! this % bias_and_activation => gaussian_b
      case('relu')
        this % activation             => relu
        this % bias_and_activation    => relu_mat_b
      case('sigmoid')
        this % activation             => sigmoid
        this % bias_and_activation    => sigmoid_mat_b
      case('hard_sigmoid')
        this % activation             => hard_sigmoid
        this % bias_and_activation    => hard_sigmoid_mat_b
      case('softsign')
        this % activation             => softsign
        this % bias_and_activation    => softsign_mat_b
      case('tanh')
        this % activation => tanhf
      case('linear')
        this % activation             => linear
        this % bias_and_activation    => linear_mat_b
      case default
        print *, "failed to read activation function, setting to linear"
        this % activation             => linear
        this % bias_and_activation    => linear_mat_b
    end select
    ! if (activation(1:8)=='gaussian') then
    !     this % activation => gaussian
    !     ! this % bias_and_activation => gaussian_b
    ! else if (activation(1:4)=='relu') then
    !     this % activation             => relu
    !     this % bias_and_activation    => relu_mat_b
    ! else if (activation(1:7)=='sigmoid') then
    !     this % activation             => sigmoid
    !     this % bias_and_activation    => sigmoid_mat_b
    ! else if (activation(1:12)=='hard_sigmoid') then
    !     this % activation             => hard_sigmoid
    !     this % bias_and_activation    => hard_sigmoid_mat_b
    ! else if (activation(1:8)=='softsign') then
    !     this % activation             => softsign
    !     this % bias_and_activation    => softsign_mat_b
    ! else if (activation(1:4)=='tanh') then
    !     this % activation => tanhf
    ! else if (activation(1:6)=='linear') then
    !     this % activation             => linear
    !     this % bias_and_activation    => linear_mat_b
    ! else
    !     print *, "failed to read activation function, setting to linear"
    !     this % activation             => linear
    !     this % bias_and_activation    => linear_mat_b
    ! end if
  end subroutine set_activation

end module mod_mlp_layer
