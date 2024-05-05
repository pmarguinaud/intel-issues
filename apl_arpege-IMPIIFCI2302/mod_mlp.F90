module mod_mlp

  ! Changelog:
  !        P. Ukkonen, 2020-2021: Optimized kernels using BLAS library calls for batched inference on CPU/GPU
  !                               Note that the kernels using BLAS/GEMM currently assume "flat" model structure
  !                               where every hidden layer has the same number of neurons
  !        P. Ukkonen, 23.1.2022: Change the code so that input layer is not included in layers,
  !                               remove RRMTPG-specific kernels
  !        P.Ukkonen, 14.9.2023: Rename "network" to "mlp" to avoid overlap with OOPS NN modules originally built on same library

  ! **OLD** NEURAL_FORTRAN view of things: input layer is included, and has weights but not bias or activation
  ! NN with two hidden layers would have allocated 4 layer derived types:
  ! 
  ! lay1      lay2     lay3      lay4
  ! (inp)     (hidden) (hidden)  (outp)   <-- 4 layer objects
  !  7        16       16         224     <-- nodes ("dims") = number of inputs, hidden neurons or outputs
  ! w1(7,16)  w2       w3         w4=0    <-- weights = zero in last layer (not used)
  ! b1=0      b2       b3         b4      <-- biases  = zero in first layer (not used)
  ! A1=lin    A2       A3         A4      <-- activation function = linear in first layer (not used)
  !
  ! **NOW**
  ! Bishop 2006* view of NNs - input layer not counted as a layer, because only the weights and biases matter
  ! The same MLP with two hidden layers now has 3 layers and looks like this:
  ! 
  !  input      hidden      hidden      output
  !             lay1        lay2        lay3  <-- 3 layer objects, each with weights, biases and activation
  !  7          16          16           224  <-- nodes/dims (=number of inputs / hidden neurons / outputs)
  !
  !     w1(7,16)     w2(16,16)    w3(16,224)  <-- weights (one set for each layer, but input layer does not count)
  !     b1(16)       b2(16)       b3(224)     <-- biases
  !     A1           A2           A3          <-- activation function (immediately after adding bias)
  !
  ! *C.M. Bishop (2006): Pattern recognition and Machine Learning
     
     
  use mo_rte_kind, only: sp
  use mod_mlp_layer, only: layer_type

#ifdef USE_OPENACC
  use cublas 
  use openacc
#define sgemm cublassgemm
#endif

  implicit none
 
#ifdef USE_TIMING
  integer, private :: ret, i
#endif

  private

  public :: mlp_type

  type :: mlp_type

    type(layer_type), allocatable :: layers(:)
    integer,          allocatable :: dims(:)
  contains

    procedure, public, pass(this) :: init
    procedure, public, pass(this) :: load
    procedure, public, pass(this) :: output_opt, output_opt_flatmodel       ! Vector input, matrix-vector product
    procedure, public, pass(this) :: output_sgemm_flat, output_sgemm_flat_byrows
    procedure, public, pass(this) :: save
    procedure, public, pass(this) :: set_activation

  end type mlp_type

  interface mlp_type
    module procedure net_constructor
  endinterface mlp_type

contains

  type(mlp_type) function net_constructor(dims, activation) result(net)
    ! network class constructor. Size of input array dims indicates the total
    ! number of layers (hidden + output) plus one (the input layer, whose dim also
    ! also needs to be specified but this pseudo-layer does not contain weights, biases or activation)
    integer, intent(in) :: dims(:)
    character(len=*), intent(in), optional :: activation
    call net % init(dims)
    if (present(activation)) then
      call net % set_activation(activation)
    else
      call net % set_activation('sigmoid')
    end if
    ! call net % sync(1)
  end function net_constructor


    subroutine init(this, dims)
      ! Allocates and initializes the layers with given dimensions dims.
      class(mlp_type), intent(in out) :: this
      integer, intent(in) :: dims(:)
      integer :: n
    allocate(this%dims(size(dims)))
    this % dims = dims

    if (.not. allocated(this % layers)) allocate(this % layers(size(dims)-1))
    do n = 2, size(dims) 
      this % layers(n-1) = layer_type(dims(n-1), dims(n))
    end do
  end subroutine init


  subroutine load(this, filename)
    ! Loads the network from file.
    class(mlp_type), intent(in out) :: this
    character(len=*), intent(in) :: filename
    character(len=20) :: activation_type

    integer :: fileunit, n, num_layers_incl_input, num_layers
    integer, allocatable :: dims(:)
    
    open(newunit=fileunit, file=filename, status='old', action='read')

    read(fileunit, fmt=*) num_layers_incl_input
    allocate(dims(num_layers_incl_input))
    num_layers = num_layers_incl_input - 1

    read(fileunit, fmt=*) dims
    call this % init(dims)

   !$acc enter data copyin(this) 
   !$acc enter data copyin(this % dims)  
   !$acc enter data copyin(this % layers)
    do n = 1, num_layers
      read(fileunit, fmt=*) this % layers(n) % b
      !$acc enter data copyin(this % layers(n) % b) async

      ! print *, "layer", n, "set bias to", this % layers(n) % 
    end do
    !$acc wait
    
    do n = 1, num_layers
      read(fileunit, fmt=*) this % layers(n) % w
      this % layers(n) % w_transposed = transpose(this % layers(n) % w )   
     !$acc enter data copyin(this % layers(n) % w_transposed) async
      ! print *, "layer", n, "set weight to", this % layers(n) % w
   
    end do
    
    do n = 1, num_layers
      read(fileunit, fmt=*) activation_type
      call this % layers(n) % set_activation(activation_type)
      ! print *, "layer", n, "set activation to", activation_type

    end do    

    close(fileunit)
    !$acc wait
  end subroutine load

  pure subroutine output_opt(this, x, output)
    class(mlp_type),    intent(in)  :: this
    real(sp), dimension(:), intent(in)  :: x
    real(sp), dimension(:), intent(out) :: output
    ! Local variables
    real(sp), allocatable   :: a(:)
    integer,  dimension(2)  :: matsize
    integer                 :: n

    associate(layers => this % layers)
      matsize = shape(layers(1) % w_transposed)
      a = matvecmul(layers(1) % w_transposed, x, matsize(1), matsize(1)) + layers(1) % b
      ! sigmoid activation: using an "inout" subroutine to avoid array copy 
      call layers(1) % activation(a)
      ! INTERMEDIATE LAYERS
      do n = 2, size(layers)-1
        matsize = shape(layers(n) % w_transposed)
        a = matvecmul(layers(n) % w_transposed, a, matsize(1), matsize(2)) + layers(n) % b
        call layers(n) % activation(a)
      end do
      ! LAST LAYER (LINEAR ACTIVATION = do nothing, just add biases)
      matsize = shape(layers(n) % w_transposed)
      output = (matvecmul(layers(n) % w_transposed, a, matsize(1), matsize(2)) + layers(n) % b)
      call layers(n) % activation(output)
    end associate
    
  end subroutine

  pure subroutine output_opt_flatmodel(this, x, output)
    ! Use forward propagation to compute the output of the network.
    ! Non-batched computations, so relatively slow but following changes are implemented for speed
    ! 1) Outputs are allocated outside of function, 
    ! 2) use of explicit-shape intermediate array that assumes the number of neurons are the same for all hidden layers,
    ! 3) activation functions are replaced with a subroutine that modifies the arguments (sigmoid), activation from final layer removed (linear activation=redundant 1:1 copy)
    ! 4) matmul replaced by custom function which is faster than matmul for matrix-vector multiplication
    ! 5) weights have been pre-transposed in the load routine.
    ! This procedure is much faster than the original when using gfortran -O3 -march=native or ifort -O3.
    ! For lower optimization levels the custom function (4) may be SLOWER
    class(mlp_type),    intent(in)  :: this
    real(sp), dimension(:), intent(in)  :: x
    real(sp), dimension(:), intent(out) :: output
    ! Local variables
    ! The signal/tensor passing through the network
    real(sp), dimension(size(this % layers(1) % w_transposed,1))  :: a 
    integer :: n, neurons

    neurons = size(this % layers(1) % w_transposed, 1)

    associate(layers => this % layers)
      a = matvecmul(layers(1) % w_transposed, x, neurons, size(x)) + layers(1) % b
      call layers(1) % activation(a)
      ! INTERMEDIATE LAYERS
      do n = 2, size(layers)-1
        a = matvecmul(layers(n) % w_transposed, a, neurons, neurons) + layers(n) % b
        call layers(n) % activation(a)
      end do
      ! LAST LAYER (LINEAR ACTIVATION = do nothing, just add biases)
      output = (matvecmul(layers(n) % w_transposed, a, size(output), neurons) + layers(n) % b)
      call layers(n) % activation(output)
    end associate
  end subroutine

  subroutine output_sgemm_flat(this, nx, ny, nbatch, x, output)
    ! Optimized batched inference function for a multi-output (regression) neural network using BLAS/cuBLAS
    ! Takes a 2D input array where the columns are the input vectors, outer dimension
    ! is the number of samples (nbatch, block size) which could be ncol*nlay
    ! Always in single-precision (sgemm) because more precision is redundant (nets trained in sp)
    ! this version assumes a "flat" network structure with regard to hidden neurons (nneur1=nneur2)    
    !
    !                             Layer Weights (T)     Layer Inputs          Layer Outputs
    ! Input (first hidden) layer    (nneur1 x nx)     * (nx x nbatch )      = (nneur1 x nbatch) 
    ! Other hidden layers           (nneur2 x nneur1) * (nneur1 x nbatch )  = (nneur2 x nbatch) 
    ! output layer                  (ny x nneur2)     * (nneur2 x nbatch )  = (ny x nbatch)  
    ! in GEMM terms:                  A (m x k)       *  B (k x N )         = C (m x N) 
    ! T = weights have been pre-transposed
    !
    class(mlp_type),              intent(in), target  :: this ! a neural network model
    integer, intent(in)                           :: nx, ny, nbatch
    real(sp), dimension(nx, nbatch), intent(in)  :: x            ! Model input
    real(sp), dimension(ny, nbatch), intent(out) :: output       ! Model output
    real(sp), dimension(size(this % layers(1) % w_transposed, 1), nbatch), &
                                          target  :: a1, a2       ! Temporary output/input between layers, of shape (neurons, nbatch)
    real(sp), dimension(:,:), contiguous, pointer :: a, a_next    ! The input to a layer is the output of the previous layer. To avoid memory
                                                                  ! movement, we can use pointers and just switch them around after each layer
    real(sp), dimension(:,:), contiguous, pointer :: wt           ! Weights
    real(sp), dimension(:),   contiguous, pointer :: b            ! BIases
    integer :: n, j, neurons, nlayers, i

    neurons = size(this % layers(1) % w_transposed, 1)
    nlayers = size(this % layers)

    !$acc enter data create(a1, a2)
    associate(layers=>this%layers)    ! so it's easier to read

      ! FIRST LAYER
      wt => layers(1) % w_transposed 
      a  => a1                        
      b  => layers(1) % b            

      !$acc host_data use_device(wt, x, a)
      call sgemm('N','N', neurons, nbatch, nx, 1.0, wt, neurons, x, nx, 0.0, a, neurons)  ! uses GPU version if USE_OPENACC=1
      !$acc end host_data

      call layers(1) % bias_and_activation(a, b)

      ! INTERMEDIATE LAYERS
      a_next => a2

      do n = 2, nlayers-1

        wt => layers(n) % w_transposed
        b  => layers(n) % b

        !$acc host_data use_device(wt, a, a_next)
        call sgemm("N","N", neurons, nbatch, neurons, 1.0, wt, neurons, a, neurons, 0.0, a_next, neurons)
        !$acc end host_data

        call layers(n) % bias_and_activation(a_next, b)

        ! Swap pointers
        if(mod(n,2) .EQ. 0) then
          a       => a2
          a_next  => a1  
        else
          a       => a1
          a_next  => a2
        end if

      end do
      
      ! Output layer
      wt => layers(n) % w_transposed
      b  => layers(n) % b

      !$acc host_data use_device(wt, a, output)
      call sgemm("N","N", ny, nbatch, neurons, 1.0, wt, ny, a, neurons, 0.0, output, ny)
      !$acc end host_data

      call layers(n) % bias_and_activation(output, b)

      end associate
    !$acc exit data detach(a,a_next) delete(a1, a2)

  end subroutine

  subroutine output_sgemm_flat_byrows(this, nx, ny, nbatch, x, output)
    ! Generic inference function using BLAS/cuBLAS for batched prediction (many samples at a time)
    ! In this version, the dimensions of inputs and outputs are reversed so that the individual
    ! NN input/output vectors are strided across the outer dimension
    ! This may be faster for very small networks
    !                               layer input           weights             layer output
    ! input (first hidden) layer    (nbatch x nx)     * (nx x nneur1)       = (nbatch x nneur1)
    ! other hidden layers           (nbatch x nneur1) * (nneur1 x nneur2)   = (nbatch * nneur2)  
    ! output layer                  (nbatch x nneur2) * (nneur2 x ny )      = (nbatch x ny)
    ! in GEMM terms:                  A (m x k)       *  B (k x N )         = C (m x N) 
    ! T = transposed weights
    class(mlp_type),              intent(in), target  :: this ! a neural network model
    integer, intent(in)                           :: nx, ny, nbatch
    real(sp), dimension(nbatch, nx), intent(in)   :: x            ! Model input
    real(sp), dimension(nbatch, ny), intent(out)  :: output       ! Model output
    real(sp), dimension(nbatch, size(this % layers(1) % w_transposed, 1)), &
                                      target      :: a1, a2       ! Temporary output/input between layers
    real(sp), dimension(:,:), contiguous, pointer :: a, a_next    ! The input to a layer is the output of the previous layer. To avoid memory
                                                                  ! movement, we can use pointers and just switch them around after each layer
    real(sp), dimension(:,:), contiguous, pointer :: w            ! Weights
    real(sp), dimension(:),   contiguous, pointer :: b            ! BIases
    integer :: n, j, neurons, nlayers, i

    neurons = size(this % layers(1) % w, 2)
    nlayers = size(this % layers)

    !$acc enter data create(a1, a2)
    associate(layers=>this%layers)    ! so it's easier to read

      ! FIRST LAYER
      w   => layers(1) % w  ! Set the weights to the weights of the first layer
      a   => a1                        
      b   => layers(1) % b            

      !$acc host_data use_device(w, x, a)
      call sgemm('N','N', nbatch, neurons, nx, 1.0, x, nbatch, w, nx, 0.0, a, nbatch)
      !$acc end host_data

      ! INTERMEDIATE LAYERS
      a_next => a2
      
      do n = 2, nlayers-1

        w => layers(n) % w
        b => layers(n) % b

        !$acc host_data use_device(w, a, a_next)
        call sgemm("N","N", nbatch, neurons, neurons, 1.0, a, nbatch, w, neurons, 0.0, a_next, nbatch)
        !$acc end host_data

        ! Swap pointers
        if(mod(n,2) .EQ. 0) then
          a       => a2
          a_next  => a1  
        else
          a       => a1
          a_next  => a2
        end if

      end do

      w => layers(n) % w
      b  => layers(n) % b

      !$acc host_data use_device(w, a, output)
      call sgemm("N","N", nbatch, ny, neurons, 1.0, a, nbatch, w, neurons, 0.0, output, nbatch)
      !$acc end host_data

      call layers(n) % bias_and_activation(output, b)

      end associate
    !$acc exit data detach(a,a_next) delete(a1, a2)

  end subroutine


elemental subroutine reluu(x) 
!$acc routine seq
  !! REctified Linear Unit (RELU) activation subroutine.
  real(sp), intent(inout) :: x
  x = max(0.0_sp, x)
end subroutine reluu

  
elemental subroutine activation_softsign(x)
!$acc routine seq
  real(sp), intent(inout) :: x
  x = x / (abs(x) + 1)
end subroutine



pure function matvecmul(matA,vecB,nrow,ncol)
    implicit none
    integer, intent(in) :: nrow,ncol
    real(sp), intent(in) :: matA(nrow,ncol)
    real(sp), intent(in) :: vecB(ncol)
    integer :: i,j
    real(sp) :: matvecmul(nrow)

    ! each (row,here ncol) element in b (length e.g. 256) is obtained by :
    ! loop through the different elements in vecB (length 50), multiply by the corresponding
    ! column (still 50) element in matA for that particular row (outer loop), add the 50 values together

    matvecmul = 0.0_sp
    do j=1,ncol !   y            X          alpha
        matvecmul = matvecmul + matA(:,j) * vecB(j) !length 256. this is for one column, and then the columns need to be added together ( b = b + ..)
    enddo

end function matvecmul


  subroutine save(this, filename)
    ! Saves the network to a file.
    class(mlp_type), intent(in out) :: this
    character(len=*), intent(in) :: filename
    integer :: fileunit, n
    open(newunit=fileunit, file=filename)
    write(fileunit, fmt=*) size(this % dims)
    write(fileunit, fmt=*) this % dims
    do n = 2, size(this % dims)
      write(fileunit, fmt=*) this % layers(n) % b
    end do
    do n = 1, size(this % dims) - 1
      write(fileunit, fmt=*) this % layers(n) % w
    end do
    close(fileunit)
  end subroutine save

  subroutine set_activation(this, activation)
    ! A thin wrapper around layer % set_activation().
    ! This method can be used to set an activation function
    ! for all layers at once.
    class(mlp_type), intent(in out) :: this
    character(len=*), intent(in) :: activation
    integer :: n
    do n = 1, size(this % layers)
      call this % layers(n) % set_activation(activation)
    end do
  end subroutine set_activation


                                              
  ! end subroutine output_sgemm_tau_sgemmbatched

end module mod_mlp
