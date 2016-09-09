module hci
! Heat-bath Configuration Interaction
! Input two parameters: epsilon1,epsilon2
! Starting with HF, each iteration,
! 1. Find lowest eigenvector in the current set of determinants.
! 2. Add all determinants i for which H_{ij}*c_j>epsilon1 for at least one
!    determinant j in current set.
! Repeat until det list increases by less than 1% on an iteration.
! Once the converged wavefunction is generated, perform second-order perturbation theory to
! improve the energy estimate, using epsilon2 to avoid wasting time on small contributions
! to the sum.
! A Holmes, Started on 7 Feb 2016
use omp_lib
use types, only : rk,ik,ik_vec,i8b
use common_ham, only : nup,ndn,norb,n_core_orb,hamiltonian_type
#ifdef NUM_ORBITALS_GT_127
use overload, only : maskr_vec
#endif

implicit none
save
private
public :: perform_cisdtq,perform_cisd_pt,perform_hci

  type neighbor_table
    integer :: n_hash_functions,n_bins_per_hash
    logical,allocatable :: is_empty(:,:)
    integer,allocatable :: bin_counts(:,:)
    integer,allocatable :: start_index(:,:),end_index(:,:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable :: dets_up(:,:),dets_dn(:,:)
#else
    integer(ik),allocatable :: dets_up(:,:),dets_dn(:,:)
#endif
    real(rk),allocatable :: coeffs(:,:)
  end type neighbor_table

  integer :: RandomHash1(0:255) !,RandomHash2(0:255)
  integer,allocatable :: n_per_bin_singles(:)
  integer,allocatable :: n_per_bin_doubles(:)
 !integer,allocatable :: n_per_bin(:)
  integer,allocatable :: hash_ind_singles(:),hash_info_singles(:)
  integer,allocatable :: hash_ind_doubles(:),hash_info_doubles(:)
 !integer,allocatable :: hash_ind(:),hash_info(:)
  integer,allocatable :: start_index_singles(:),end_index_singles(:)
  integer,allocatable :: start_index_doubles(:),end_index_doubles(:)
 !integer,allocatable :: start_index(:),end_index(:)

 !integer :: n_bins = 1000000007 ! number of bins in hash tables. can change, but this seems reasonable
 !integer :: n_bins = 100000007 ! number of bins in hash tables. can change, but this seems reasonable
  integer :: n_bins = 1000003 ! number of bins in hash tables. can change, but this seems reasonable
  integer :: max_sum_terms

contains

  subroutine perform_cisdtq(min_elem)

    use semistoch, only : hamiltonian

    real(rk),intent(in) :: min_elem

    integer :: iter,max_iters
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable :: old_dets_up(:),old_dets_dn(:),new_dets_up(:),new_dets_dn(:)
    type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik),allocatable :: old_dets_up(:),old_dets_dn(:),new_dets_up(:),new_dets_dn(:)
    integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
    real(rk),allocatable :: old_wts(:),new_wts(:)
    integer :: ndets_old,ndets_new
    real(rk) :: energy,old_energy
    logical :: is_con
    integer :: i
    integer,allocatable :: iorder(:),temp_i_2(:)
   !integer :: cdets = 4000
   !integer :: tdets = 99999999

    stop "CISDTQ not currently up-to-date!"
!    write (6,*) "Performing CISDTQ with min_elem=",min_elem; call flush(6)
!
!    allocate(old_dets_up(tdets))
!    allocate(old_dets_dn(tdets))
!    allocate(old_wts(tdets))
!
!    allocate(new_dets_up(tdets))
!    allocate(new_dets_dn(tdets))
!    allocate(new_wts(tdets))
!
!    ! Start with HF det
!    ndets_old = 1
!#ifdef NUM_ORBITALS_GT_127
!    old_dets_up(1) = maskr_vec(nup)
!    old_dets_dn(1) = maskr_vec(ndn)
!#else
!    old_dets_up(1) = maskr(nup,ik)
!    old_dets_dn(1) = maskr(ndn,ik)
!#endif
!    old_wts(1) = 1._rk
!    call hamiltonian(old_dets_up(1),old_dets_dn(1),old_dets_up(1),old_dets_dn(1),energy,is_con)
!    write (6,*) "After iteration",0,",",1,"dets with energy",energy; call flush(6)
!    old_energy = energy
!
!    allocate(iorder(tdets))
!    allocate(temp_i16_up((tdets+1)/2))
!    allocate(temp_i16_dn((tdets+1)/2))
!    allocate(temp_i_2((tdets+1)/2))
!
!    max_iters = 2
!    do iter=1,max_iters
!      ! First iteration: Full CISD
!      if (iter==1) then
!        call get_next_det_list(cdets,tdets,-1._rk,energy,ndets_old,old_dets_up,old_dets_dn,old_wts,ndets_new,new_dets_up,new_dets_dn,new_wts,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
!      else
!        call get_next_det_list(cdets,tdets,min_elem,energy,ndets_old,old_dets_up,old_dets_dn,old_wts,ndets_new,new_dets_up,new_dets_dn,new_wts,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
!      endif
!
!     !write (6,*) "Before Lanczos:"
!     !do i=1,min(20,ndets_new)
!     !  write (6,*) new_dets_up(i),new_dets_dn(i),new_wts(i)
!     !enddo
!
!      call lanczos_diagonalize(ndets_new,new_dets_up,new_dets_dn,final_wts=old_wts,lowest_energy=energy)
!     !call lanczos_diagonalize(ndets_new,new_dets_up,new_dets_dn,new_wts,old_wts,energy)
!
!      write (6,*) "After iteration",iter,",",ndets_new,"dets with energy",energy; call flush(6)
!
!      ndets_old = ndets_new
!      old_dets_up(1:ndets_old) = new_dets_up(1:ndets_old)
!      old_dets_dn(1:ndets_old) = new_dets_dn(1:ndets_old)
!
!      if (abs(energy-old_energy)<1.e-6_rk)  exit
!
!      old_energy = energy
!
!    enddo
!
!    ! Now, old_dets_up/dn, old_wts, old_energy are the variational wavefunction
!
!    write (6,*) "Final variational wavefunction:"
!    do i=1,min(20,ndets_old)
!      write (6,*) i,old_dets_up(i),old_dets_dn(i),old_wts(i)
!    enddo

  end subroutine perform_cisdtq


  subroutine perform_cisd_pt(pt_eps)

    use semistoch, only : hamiltonian

    real(rk),intent(in) :: pt_eps

    integer :: iter,max_iters
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable :: old_dets_up(:),old_dets_dn(:),new_dets_up(:),new_dets_dn(:)
    type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik),allocatable :: old_dets_up(:),old_dets_dn(:),new_dets_up(:),new_dets_dn(:)
    integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
    real(rk),allocatable :: old_wts(:),new_wts(:)
    integer :: ndets_old,ndets_new
    real(rk) :: energy,old_energy
    logical :: is_con
    integer :: i
    integer,allocatable :: iorder(:),temp_i_2(:)
    real(rk),allocatable :: diag_elems(:)
   !integer :: cdets = 4000
   !integer :: tdets = 99999999

    stop "CISD+2PT not currently up-to-date!"

!    write (6,*) "Performing CISD+2PT with pt_eps=",pt_eps; call flush(6)
!
!    allocate(old_dets_up(tdets))
!    allocate(old_dets_dn(tdets))
!    allocate(old_wts(tdets))
!
!    allocate(new_dets_up(tdets))
!    allocate(new_dets_dn(tdets))
!    allocate(new_wts(tdets))
!
!    ! Start with HF det
!    ndets_old = 1
!#ifdef NUM_ORBITALS_GT_127
!    old_dets_up(1) = maskr_vec(nup)
!    old_dets_dn(1) = maskr_vec(ndn)
!#else
!    old_dets_up(1) = maskr(nup,ik)
!    old_dets_dn(1) = maskr(ndn,ik)
!#endif
!    old_wts(1) = 1._rk
!    call hamiltonian(old_dets_up(1),old_dets_dn(1),old_dets_up(1),old_dets_dn(1),energy,is_con)
!    write (6,*) "After iteration",0,",",1,"dets with energy",energy; call flush(6)
!    old_energy = energy
!
!    allocate(iorder(tdets))
!    allocate(temp_i16_up((tdets+1)/2))
!    allocate(temp_i16_dn((tdets+1)/2))
!    allocate(temp_i_2((tdets+1)/2))
!
!    max_iters = 1
!    do iter=1,max_iters
!      ! First (and only) iteration: Full CISD
!      call get_next_det_list(cdets,tdets,-1._rk,energy,ndets_old,old_dets_up,old_dets_dn,old_wts,ndets_new,new_dets_up,new_dets_dn,new_wts,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
!
!      write (6,*) "Before Lanczos:"
!      do i=1,min(20,ndets_new)
!        write (6,*) new_dets_up(i),new_dets_dn(i),new_wts(i)
!      enddo
!
!      call lanczos_diagonalize(ndets_new,new_dets_up,new_dets_dn,new_wts,old_wts,energy)
!
!      write (6,*) "After iteration",iter,",",ndets_new,"dets with energy",energy; call flush(6)
!
!      ndets_old = ndets_new
!      old_dets_up(1:ndets_old) = new_dets_up(1:ndets_old)
!      old_dets_dn(1:ndets_old) = new_dets_dn(1:ndets_old)
!
!      if (abs(energy-old_energy)<1.e-6_rk)  exit
!
!      old_energy = energy
!
!    enddo
!
!    ! Now, old_dets_up/dn, old_wts, old_energy are the variational wavefunction
!
!    write (6,*) "Final variational wavefunction:"
!    do i=1,min(20,ndets_old)
!      write (6,*) i,old_dets_up(i),old_dets_dn(i),old_wts(i)
!    enddo
!
!    ! Quick hack: compute all diag_elems now
!    allocate(diag_elems(ndets_old))
!    do i=1,ndets_old
!      call hamiltonian(old_dets_up(i),old_dets_dn(i),old_dets_up(i),old_dets_dn(i),diag_elems(i),is_con)
!    enddo
!
!    call second_order_pt(ndets_old,old_dets_up,old_dets_dn,old_wts,diag_elems,energy,pt_eps)
!   !call second_order_pt(ndets_old,old_dets_up,old_dets_dn,old_wts,diag_elems,energy,pt_eps,iorder,temp_i16_up,temp_i16_dn,temp_i_2)

  end subroutine perform_cisd_pt


  subroutine perform_hci(min_elem,pt_eps,hf_up,hf_dn)

    use semistoch, only : hamiltonian
    use common_psi_t, only : trial_wf_iters
    use types, only : i8b

    real(rk),intent(in) :: min_elem,pt_eps

    integer :: iter,max_iters
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),optional,intent(in) :: hf_up,hf_dn
    type(ik_vec),allocatable :: old_dets_up(:),old_dets_dn(:),new_dets_up(:),new_dets_dn(:)
    type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik),optional,intent(in) :: hf_up,hf_dn
    integer(ik),allocatable :: old_dets_up(:),old_dets_dn(:),new_dets_up(:),new_dets_dn(:)
    integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
    real(rk),allocatable :: old_wts(:),new_wts(:)
    integer :: ndets_old,ndets_new
    real(rk) :: energy,old_energy
    logical :: is_con
    integer :: i
    integer,allocatable :: iorder(:),temp_i_2(:)
    real(rk),allocatable :: diag_elems(:)
    integer :: tdets=9999999 ! max number of dets
    integer(i8b) :: n_conn_estimate
    integer :: n_batches
    real(rk) :: rdm(norb,norb)

    write (6,'(''Performing HCI with epsilon_1='',es12.4,'' epsilon_2='',es12.4)') min_elem, pt_eps ; call flush(6)

    allocate(old_dets_up(tdets))
    allocate(old_dets_dn(tdets))
    allocate(old_wts(tdets))

    allocate(new_dets_up(tdets))
    allocate(new_dets_dn(tdets))
    allocate(new_wts(tdets))

    ! Start with HF det
    ndets_old = 1
    if (present(hf_up)) then
      if (hf_up==hf_dn) then
        old_dets_up(1) = hf_up
        old_dets_dn(1) = hf_dn
      else
        ndets_old = 2
        old_dets_up(1) = hf_up
        old_dets_dn(1) = hf_dn
        old_dets_up(2) = hf_dn
        old_dets_dn(2) = hf_up
      endif
    else
#ifdef NUM_ORBITALS_GT_127
      old_dets_up(1) = maskr_vec(nup)
      old_dets_dn(1) = maskr_vec(ndn)
#else
      old_dets_up(1) = maskr(nup,ik)
      old_dets_dn(1) = maskr(ndn,ik)
#endif
    endif
    old_wts(1) = 1._rk
    if (present(hf_up)) then
      if (hf_up.ne.hf_dn)  old_wts(2) = 1._rk
    endif
    call hamiltonian(old_dets_up(1),old_dets_dn(1),old_dets_up(1),old_dets_dn(1),energy,is_con)
    write (6,'(''Iteration   0       '',i4'' dets, energy='',f16.6)') ndets_old,energy; call flush(6)
    old_energy = energy

    allocate(iorder(tdets))
    allocate(temp_i16_up((tdets+1)/2))
    allocate(temp_i16_dn((tdets+1)/2))
    allocate(temp_i_2((tdets+1)/2))

    max_iters = 20
    do iter=1,max_iters
      call get_next_det_list(min_elem,energy,ndets_old,old_dets_up,old_dets_dn,old_wts,ndets_new,new_dets_up,new_dets_dn) !,new_wts,iorder,temp_i16_up,temp_i16_dn,temp_i_2)

     !write (6,*) "Before Lanczos:"
     !do i=1,min(20,ndets_new)
     !  write (6,*) new_dets_up(i),new_dets_dn(i),new_wts(i)
     !enddo

! Warning: tmp
      if (ndets_new<=int(1.01*ndets_old)) then
!     if (ndets_new<=int(1.001*ndets_old)) then
        ndets_new = ndets_old
        exit
      endif
     !if (ndets_old==ndets_new)  exit

      call lanczos_diagonalize(ndets_new,new_dets_up,new_dets_dn,final_wts=old_wts,lowest_energy=energy)
     !call lanczos_diagonalize(ndets_new,new_dets_up,new_dets_dn,new_wts,old_wts,energy)

      write (6,'(''Iteration'',i4,i8,'' dets, energy='',f16.6)') iter, ndets_new, energy; call flush(6)

      ndets_old = ndets_new
      old_dets_up(1:ndets_old) = new_dets_up(1:ndets_old)
      old_dets_dn(1:ndets_old) = new_dets_dn(1:ndets_old)

      if (abs(energy-old_energy)<1.e-6_rk)  exit

      old_energy = energy

    enddo

    write (6,'(/,''Final Iteration'',i4,i8,'' dets, energy='',f16.6)') iter-1, ndets_new, energy; call flush(6)

    ! Now, old_dets_up/dn, old_wts, old_energy are the variational wavefunction

    write (6,'(''Final variational wavefunction (20 most important):'')')
    do i=1,min(20,ndets_old)
      write (6,*) i,old_dets_up(i),old_dets_dn(i),old_wts(i)
    enddo

    call my_second(1,"1rdm")
    ! Quick hack: compute all diag_elems now
    allocate(diag_elems(ndets_old))
    do i=1,ndets_old
      call hamiltonian(old_dets_up(i),old_dets_dn(i),old_dets_up(i),old_dets_dn(i),diag_elems(i),is_con)
    enddo

    n_conn_estimate = estimate_n_connections(ndets_old,old_dets_up,old_dets_dn,old_wts,pt_eps)
    write (6,*) "Estimated number of connections to variational wavefunction=",n_conn_estimate; call flush(6)
    n_batches = max(8,n_conn_estimate/10**8)
    ! n_batches = max(1,n_conn_estimate/10**8)
    write (6,*) "Performing 1-RDM calculation in",n_batches,"batches"; call flush(6)

    call get_1rdm(ndets_old,old_dets_up,old_dets_dn,old_wts,rdm)
    call get_1rdm_with_pt(ndets_old,old_dets_up,old_dets_dn,old_wts,diag_elems,energy,pt_eps,n_batches,rdm)
    call my_second(2,"1rdm")

   !call second_order_pt_ah(ndets_old,old_dets_up,old_dets_dn,old_wts,diag_elems,energy,pt_eps)
    n_conn_estimate = estimate_n_connections(ndets_old,old_dets_up,old_dets_dn,old_wts,pt_eps)
    write (6,*) "Estimated number of connections to variational wavefunction=",n_conn_estimate; call flush(6)
    n_batches = max(8,n_conn_estimate/10**8)
    ! n_batches = max(1,n_conn_estimate/10**8)
    write (6,*) "Performing PT correction in",n_batches,"batches"; call flush(6)
    call second_order_pt_dtm_batches(ndets_old,old_dets_up,old_dets_dn,old_wts,diag_elems,energy,pt_eps,n_batches) !,iorder,temp_i16_up,temp_i16_dn,temp_i_2)

  end subroutine perform_hci


  subroutine get_next_det_list(min_elem,old_energy,ndets_old,old_dets_up,old_dets_dn,old_wts,ndets_new,new_dets_up,new_dets_dn) !,new_wts,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
  ! Gets next determinant list for the next iteration of Heat-bath CI
  ! Warning: Assumes old dets list is sorted in decreasing order by abs wt!
  ! This is because early termination of the loop is possible if a determinant
  ! is reached whose abs wt is less than epsilon/H_doub_max
  ! A Holmes, Feb 2016

    use semistoch, only : find_doubly_excited
   !use generic_sort, only : sort
   !use tools, only : merge_sort2_up_dn

    real(rk),intent(in) :: min_elem,old_energy
    integer,intent(in) :: ndets_old
    integer,intent(out) :: ndets_new
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: old_dets_up(:),old_dets_dn(:)
    type(ik_vec),intent(out) :: new_dets_up(:),new_dets_dn(:)
    type(ik_vec),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
   !type(ik_vec),allocatable :: tmp_dets_up(:),tmp_dets_dn(:)
#else
    integer(ik),intent(in) :: old_dets_up(:),old_dets_dn(:)
    integer(ik),intent(out) :: new_dets_up(:),new_dets_dn(:)
   !integer(ik),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
    integer(ik),allocatable :: tmp_dets_up(:),tmp_dets_dn(:)
#endif
    real(rk),intent(in) :: old_wts(:)
   !real(rk),intent(out) :: new_wts(:)
   !integer,intent(inout) :: iorder(:),temp_i_2(:)
   !real(rk),allocatable :: tmp_wts(:),garbage(:)
    integer :: i,n_ref

    ! Find dets i for which H_ij*c_j > epsilon_1 for at least one j
    call my_second(1,"Generate connections")
    if (min_elem>=0._rk) then
      n_ref = ndets_old
      call find_doubly_excited(n_det=ndets_new,dets_up=tmp_dets_up,dets_dn=tmp_dets_dn,ref_up=old_dets_up(1:n_ref),ref_dn=old_dets_dn(1:n_ref),norb=norb,n_core_orb=n_core_orb,ref_coeffs=old_wts(1:n_ref),ninitiator=n_ref,min_elem=min_elem) !,e_mix_num=tmp_wts,e_mix_den=garbage,pt=old_energy,min_elem=min_elem)
     !call find_doubly_excited(n_det=ndets_new,dets_up=tmp_dets_up,dets_dn=tmp_dets_dn,ref_up=old_dets_up(1:n_ref),ref_dn=old_dets_dn(1:n_ref),norb=norb,n_core_orb=n_core_orb,ref_coeffs=old_wts(1:n_ref),ninitiator=n_ref,e_mix_num=tmp_wts,e_mix_den=garbage,pt=old_energy,min_elem=min_elem)
    else
      stop "Must have epsilon_1 >= 0!"
     !n_ref = min(cdets,ndets_old)
     !call find_doubly_excited(n_det=ndets_new,dets_up=tmp_dets_up,dets_dn=tmp_dets_dn,ref_up=old_dets_up(1:n_ref),ref_dn=old_dets_dn(1:n_ref),norb=norb,n_core_orb=n_core_orb,ref_coeffs=old_wts(1:n_ref),ninitiator=n_ref,e_mix_num=tmp_wts,e_mix_den=garbage,pt=old_energy)
    endif
    call my_second(2,"Generate connections")
    ! tmp_dets sorted by label

   !write (6,*) "After find_doubly_excited:"
   !do i=1,min(20,ndets_new)
   !  write (6,*) i,tmp_dets_up(i),tmp_dets_dn(i),tmp_wts(i)
   !enddo

   !call sort(tmp_wts(1:ndets_new), tmp_dets_up(1:ndets_new), tmp_dets_dn(1:ndets_new))
   !! Now, tmp_dets are sorted by magnitude of coefficients

   !ndets_new = min(tdets,ndets_new)
    new_dets_up(1:ndets_new) = tmp_dets_up(1:ndets_new)
    new_dets_dn(1:ndets_new) = tmp_dets_dn(1:ndets_new)
   !new_wts(1:ndets_new) = tmp_wts(1:ndets_new)

    ! Sort by label, to enable future binary searches for constructing Hamiltonian for Lanczos
   !do i=1,ndets_new
   !  iorder(i) = i
   !enddo
   !call merge_sort2_up_dn(new_dets_up(1:ndets_new),new_dets_dn(1:ndets_new), iorder, ndets_new, temp_i16_up, temp_i16_dn, temp_i_2)
   !new_wts(1:ndets_new) = new_wts(iorder(1:ndets_new)) 
    
    deallocate(tmp_dets_up,tmp_dets_dn)

  end subroutine get_next_det_list


  subroutine lanczos_diagonalize(ndets,dets_up,dets_dn,final_wts,lowest_energy)
  ! A Holmes, 10 Apr 2016

    use chemistry, only : diagonalize_sparse_hamiltonian_chem,time_sym,find_important_connected_dets_chem,find_important_connected_dets_1e_removed,find_important_singly_connected_dets_1e_removed,find_important_doubly_connected_dets_2e_removed,get_new_diag_elem
    use common_ham, only : nelec
    use common_run, only : connected_dets_up,connected_dets_dn,connected_matrix_elements,max_connected_dets,diag_elem_info,connected_diag_elems_info
    use generic_sort, only : sort,sort_by_first_argument

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
#endif
    real(rk),intent(out) :: final_wts(:),lowest_energy
    real(rk) :: tau_out_for_diag_ham
    integer,allocatable :: tmp_ind_singles(:),tmp_info_singles(:)
    integer,allocatable :: tmp_ind_doubles(:),tmp_info_doubles(:)
   !integer,allocatable :: tmp_ind(:),tmp_info(:)
    integer,allocatable :: ref_connection_indices(:)
    integer :: n_allocate_singles, n_allocate_doubles
    integer :: n_hash_singles, n_hash_doubles
    integer :: n_connected_dets
    integer :: i_ref,i,j

    if (hamiltonian_type.ne.'chem')  stop "Only chem implemented so far!"
!
!    ! New method: No need to store entire Hamiltonian; only store O(N_el^2*N_det) partial
!    ! connections. Then, applying Hamiltonian can be done in O(N_nonzero) time. We can
!    ! keep doing Lanczos for now, but may want to do Davidson later.
!
!    ! First, create hash tables
!
!    call init_hash_fn ! Sets up hash function
!
!    n_allocate_singles = ndets*nelec
!    n_allocate_doubles = ndets*(nelec*(nelec-1))/2
!   !n_allocate_doubles = ndets*nelec ! will probably need to be reallocated
!    allocate(hash_ind_singles(n_allocate_singles))
!    allocate(hash_info_singles(n_allocate_singles))
!    allocate(hash_ind_doubles(n_allocate_doubles))
!    allocate(hash_info_doubles(n_allocate_doubles))
!    n_hash_singles = 0
!    n_hash_doubles = 0
!
!    ! Single excitations
!    do i_ref=1,ndets
!!write (6,*) i_ref,dets_up(i_ref), dets_dn(i_ref)
!      call find_important_singly_connected_dets_1e_removed(dets_up(i_ref), dets_dn(i_ref), n_connected_dets, connected_dets_up, connected_dets_dn)
!      ! Now, connected_dets represents all those 1-electron-removed singly excited dets connected to reference
!      ! These are the places in the hash table where the index i_ref needs to be stored.
!      do j=1,n_connected_dets
!        n_hash_singles = n_hash_singles + 1
!        hash_ind_singles(n_hash_singles) = hash_det(connected_dets_up(j),connected_dets_dn(j)) ! index of hash table
!        hash_info_singles(n_hash_singles) = i_ref ! what is stored in hash table
!      enddo ! j
!    enddo ! i_ref
!!write (6,*) "hash_ind_singles="
!!do i=1,n_hash_singles
!!  write (6,*) i,hash_ind_singles(i)
!!enddo
!!call flush(6)
!
!    ! Double excitations
!    do i_ref=1,ndets
!      call find_important_doubly_connected_dets_2e_removed(dets_up(i_ref), dets_dn(i_ref), 0._rk, n_connected_dets, connected_dets_up, connected_dets_dn)
!      ! Now, connected_dets represents all those 2-electron-removed doubly excited dets connected to reference by more than epsilon
!      ! These are the places in the hash table where the index i_ref needs to be stored.
!      if (n_hash_doubles+n_connected_dets>n_allocate_doubles) then ! reallocate
!        write (6,*) "Reallocating hash tables from size",n_allocate_doubles,"to size", 2*n_allocate_doubles; call flush(6)
!        allocate(tmp_ind_doubles(n_allocate_doubles))
!        allocate(tmp_info_doubles(n_allocate_doubles))
!        tmp_ind_doubles(:) = hash_ind_doubles(:)
!        tmp_info_doubles(:) = hash_info_doubles(:)
!        deallocate(hash_ind_doubles,hash_info_doubles)
!        n_allocate_doubles = 2*n_allocate_doubles
!        allocate(hash_ind_doubles(n_allocate_doubles))
!        allocate(hash_info_doubles(n_allocate_doubles))
!        hash_ind_doubles(1:size(tmp_ind_doubles)) = tmp_ind_doubles(:)
!        hash_info_doubles(1:size(tmp_info_doubles)) = tmp_info_doubles(:)
!        deallocate(tmp_ind_doubles,tmp_info_doubles)
!      endif
!      do j=1,n_connected_dets
!        n_hash_doubles = n_hash_doubles + 1
!        hash_ind_doubles(n_hash_doubles) = hash_det(connected_dets_up(j),connected_dets_dn(j)) ! index of hash table
!        hash_info_doubles(n_hash_doubles) = i_ref ! what is stored in hash table
!      enddo ! j
!    enddo ! i_ref
!!write (6,*) "hash_ind_doubles="
!!do i=1,n_hash_doubles
!!  write (6,*) i,hash_ind_doubles(i)
!!enddo
!!call flush(6)
!
!    ! Now, sort the hash tables by hash_ind
!    call sort_by_first_argument(n_hash_singles,hash_ind_singles,hash_info_singles)
!    call sort_by_first_argument(n_hash_doubles,hash_ind_doubles,hash_info_doubles)
!
!    ! Now, store the start_index and end_index for each hash_ind
!    allocate(start_index_singles(n_bins))
!    allocate(n_per_bin_singles(n_bins))
!    start_index_singles(:) = 0
!    start_index_singles(hash_ind_singles(1)) = 1
!    do i=2,n_hash_singles
!      if (hash_ind_singles(i).ne.hash_ind_singles(i-1)) then
!        start_index_singles(hash_ind_singles(i)) = i
!        n_per_bin_singles(hash_ind_singles(i-1)) = start_index_singles(hash_ind_singles(i)) - start_index_singles(hash_ind_singles(i-1))
!      endif
!    enddo
!    n_per_bin_singles(hash_ind_singles(n_hash_singles-1)) = start_index_singles(hash_ind_singles(n_hash_singles)) - start_index_singles(hash_ind_singles(n_hash_singles-1))
!    deallocate(hash_ind_singles)
!
!    allocate(start_index_doubles(n_bins))
!    allocate(n_per_bin_doubles(n_bins))
!    start_index_doubles(:) = 0
!    start_index_doubles(hash_ind_doubles(1)) = 1
!    do i=2,n_hash_doubles
!      if (hash_ind_doubles(i).ne.hash_ind_doubles(i-1)) then
!        start_index_doubles(hash_ind_doubles(i)) = i
!        n_per_bin_doubles(hash_ind_doubles(i-1)) = start_index_doubles(hash_ind_doubles(i)) - start_index_doubles(hash_ind_doubles(i-1))
!      endif
!    enddo
!    n_per_bin_doubles(hash_ind_doubles(n_hash_doubles-1)) = start_index_doubles(hash_ind_doubles(n_hash_doubles)) - start_index_doubles(hash_ind_doubles(n_hash_doubles-1))
!    deallocate(hash_ind_doubles)
!
!    ! At this point, hash_info(start_index(i):end_index(i)), where i=hash_det(det_up,det_dn), are the
!    ! reference dets connected to det_up/dn
!
!    ! Now, do Lanczos without storing Hamiltonian!
!    call matrix_lanczos_partial_connections(ndets,dets_up,dets_dn,final_wts(1:ndets),lowest_energy,hash_info_singles,n_per_bin_singles,start_index_singles,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
!
!    deallocate(hash_info_singles,n_per_bin_singles,start_index_singles,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
!
!    ! Finally, sort dets in order of decreasing abs wt!
!    call sort(final_wts(1:ndets), dets_up(1:ndets), dets_dn(1:ndets)) ! Sort in order of decreasing absolute coef
!
    ! Old method:
    call diagonalize_sparse_hamiltonian_chem(dets_up(1:ndets),dets_dn(1:ndets),ndets,final_wts(1:ndets),lowest_energy,tau_out_for_diag_ham,time_sym,sort_or_not=.true.)

  end subroutine lanczos_diagonalize


  subroutine second_order_pt(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,epsilon)
    ! Compute and print the second order energy lowering
    ! A Holmes, Feb 2016
    use semistoch, only : find_doubly_excited,hamiltonian

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#endif
    real(rk),intent(in) :: wts(:),diag_elems(:),var_energy,epsilon
    real(rk) :: delta_e_2pt
    real(rk),allocatable :: connected_wts(:),new_diag_elems(:),e_mix_den(:)
    integer :: i,ndets_connected

    write(6,*) ; call my_second(1,'Compute 2nd-order energy correction')

    if (epsilon>=0._rk) then
      call find_doubly_excited(n_det=ndets_connected,dets_up=connected_dets_up,dets_dn=connected_dets_dn,ref_up=dets_up(1:ndets),ref_dn=dets_dn(1:ndets),norb=norb,n_core_orb=n_core_orb,ref_coeffs=wts(1:ndets),ninitiator=ndets,e_mix_num=connected_wts,e_mix_den=e_mix_den,min_elem=epsilon,ref_diag_elems=diag_elems(1:ndets),new_diag_elems=new_diag_elems)
    else
      call find_doubly_excited(n_det=ndets_connected,dets_up=connected_dets_up,dets_dn=connected_dets_dn,ref_up=dets_up(1:ndets),ref_dn=dets_dn(1:ndets),norb=norb,n_core_orb=n_core_orb,ref_coeffs=wts(1:ndets),ninitiator=ndets,e_mix_num=connected_wts,e_mix_den=e_mix_den,ref_diag_elems=diag_elems(1:ndets),new_diag_elems=new_diag_elems)
    endif

    write (6,'(''Total number of connected dets='',i12)') ndets_connected

    delta_e_2pt = 0._rk
    do i=1,ndets_connected
      if (e_mix_den(i)==0._rk) then
        delta_e_2pt = delta_e_2pt + connected_wts(i)**2/(var_energy-new_diag_elems(i))
        if (i<=20)  write (6,*) i,connected_dets_up(i),connected_dets_dn(i),connected_wts(i),new_diag_elems(i),connected_wts(i)**2/(var_energy-new_diag_elems(i)),delta_e_2pt
      endif
    enddo

    write (6,'(''Variational energy='',t36,f15.9)') var_energy
    write (6,'(''Second-order PT energy lowering='',t36,f15.9)') delta_e_2pt
    write (6,'(''Total energy='',t36,f15.9)') var_energy+delta_e_2pt
    write (6,'(''ndets, ndets_connected(total), Variational, PT, Total Energies='',i8,i10,9f16.9)') ndets, ndets_connected, var_energy, delta_e_2pt, var_energy+delta_e_2pt
    call flush(6)

    call my_second(2,'Compute 2nd-order energy correction')

  end subroutine second_order_pt


  subroutine second_order_pt_dtm_batches(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,epsilon,n_batches)
    ! Compute and print the second order energy lowering using Cyrus's suggestion of dividing up the
    ! reference into batches and exciting from all pairs of batches.
    ! A Holmes, 20 June 2016
    use semistoch, only : find_doubly_excited,hamiltonian
    use more_tools, only : binary_search
    use tools, only : merge_sort2_up_dn
    use types, only : i8b

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),allocatable :: connected_dets_up_i(:),connected_dets_dn_i(:)
    type(ik_vec),allocatable :: connected_dets_up_j(:),connected_dets_dn_j(:)
    type(ik_vec),allocatable :: sorted_dets_up(:),sorted_dets_dn(:)
    type(ik_vec),allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),allocatable :: connected_dets_up_i(:),connected_dets_dn_i(:)
    integer(ik),allocatable :: connected_dets_up_j(:),connected_dets_dn_j(:)
    integer(ik),allocatable :: sorted_dets_up(:),sorted_dets_dn(:)
    integer(ik),allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
    real(rk),intent(in) :: wts(:),diag_elems(:),var_energy,epsilon
    integer,intent(in) :: n_batches

    real(rk) :: delta_e_2pt, delta_e_2pt_inc
    real(rk),allocatable :: connected_wts_i(:),new_diag_elems_i(:),e_mix_den_i(:)
    real(rk),allocatable :: connected_wts_j(:),new_diag_elems_j(:),e_mix_den_j(:)
    integer :: i,j,ndets_connected,ndets_connected_i,ndets_connected_j
    integer :: ibatch,jbatch
    integer(i8b) :: k
    integer,allocatable :: iorder(:),temp_i_2(:)

    call my_second(1,'Compute 2nd-order energy correction using deterministic batches')

    allocate(sorted_dets_up(ndets))
    allocate(sorted_dets_dn(ndets))

    sorted_dets_up(:) = dets_up(1:ndets)
    sorted_dets_dn(:) = dets_dn(1:ndets)

    ! Sort by label for binary search
    allocate(iorder(ndets))
    allocate(temp_i16_up((ndets+1)/2))
    allocate(temp_i16_dn((ndets+1)/2))
    allocate(temp_i_2((ndets+1)/2))
    do j=1,ndets
      iorder(j)=j
    enddo

    call merge_sort2_up_dn(sorted_dets_up(1:ndets),sorted_dets_dn(1:ndets), iorder, ndets, temp_i16_up, temp_i16_dn, temp_i_2)

    deallocate(iorder)
    deallocate(temp_i16_up)
    deallocate(temp_i16_dn)
    deallocate(temp_i_2)

    delta_e_2pt = 0._rk
    ndets_connected = 0

    write (*, '("master threads: ", I2)') omp_get_num_threads()
        
    call omp_set_num_threads(8)

    !$omp parallel default(shared)
    !$omp single

    do ibatch=1,n_batches

      call find_doubly_excited(n_det=ndets_connected_i,dets_up=connected_dets_up_i,dets_dn=connected_dets_dn_i,ref_up=dets_up(ibatch:ndets:n_batches),ref_dn=dets_dn(ibatch:ndets:n_batches),norb=norb,n_core_orb=n_core_orb,ref_coeffs=wts(ibatch:ndets:n_batches),e_mix_num=connected_wts_i,e_mix_den=e_mix_den_i,min_elem=epsilon,ref_diag_elems=diag_elems(ibatch:ndets:n_batches),new_diag_elems=new_diag_elems_i)
    
      do jbatch=ibatch,n_batches

        call find_doubly_excited(n_det=ndets_connected_j,dets_up=connected_dets_up_j,dets_dn=connected_dets_dn_j,ref_up=dets_up(jbatch:ndets:n_batches),ref_dn=dets_dn(jbatch:ndets:n_batches),norb=norb,n_core_orb=n_core_orb,ref_coeffs=wts(jbatch:ndets:n_batches),e_mix_num=connected_wts_j,e_mix_den=e_mix_den_j,min_elem=epsilon,ref_diag_elems=diag_elems(jbatch:ndets:n_batches),new_diag_elems=new_diag_elems_j)


        !$omp task &
        !$omp& default(none) &
        !$omp& private(j, k, delta_e_2pt_inc) &
        !$omp& firstprivate(n_batches, ndets, ibatch, jbatch, var_energy) &
        !$omp& shared(ndets_connected_i, ndets_connected_j) &
        !$omp& shared(connected_dets_up_i, connected_dets_up_j, connected_dets_dn_i, connected_dets_dn_j) &
        !$omp& shared(sorted_dets_up, sorted_dets_dn, connected_wts_i, connected_wts_j) &
        !$omp& shared(new_diag_elems_i, e_mix_den_i) &
        !$omp& shared(delta_e_2pt)
        write (*, *) omp_get_wtime()
        write (*, '("worker in task: ", I2, "/", I2, "/", I2, " is processing batch: ", I8, I8)'), omp_get_thread_num(), omp_get_num_threads(), omp_get_max_threads(), ibatch, jbatch
        
        j = 1
        delta_e_2pt_inc = 0.0_rk

        write (*, *) "about to start"

        do i=1,ndets_connected_i
          ! Increment j until it is no longer smaller than i
          do while (connected_dets_up_i(i)>connected_dets_up_j(j).or.(connected_dets_up_i(i)==connected_dets_up_j(j).and.connected_dets_dn_i(i)>connected_dets_dn_j(j)))
            if (j==ndets_connected_j)  exit
            j = j+1
          enddo
          write (*, *) "in loop"
          if (connected_dets_up_i(i)==connected_dets_up_j(j).and.connected_dets_dn_i(i)==connected_dets_dn_j(j)) then
            if (n_batches>1) then
              call binary_search(connected_dets_up_i(i),connected_dets_dn_i(i),sorted_dets_up(1:ndets),sorted_dets_dn(1:ndets),k)
              if (k==0) then
                if (ibatch==jbatch) then
                  delta_e_2pt_inc = delta_e_2pt_inc + connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
                else
                  delta_e_2pt_inc = delta_e_2pt_inc + 2*connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
                endif
               !write (6,*) i,connected_dets_up_i(i),connected_dets_dn_i(i),connected_wts_i(i),new_diag_elems_i(i),connected_wts_i(i)**2/(var_energy-new_diag_elems_i(i)),delta_e_2pt
              endif
            else ! n_batches==1
              if (e_mix_den_i(i)==0._rk) then
                if (ibatch==jbatch) then
                  delta_e_2pt_inc = delta_e_2pt_inc + connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
                else
                  delta_e_2pt_inc = delta_e_2pt_inc + 2*connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
                endif
               !write (6,*) i,connected_dets_up_i(i),connected_dets_dn_i(i),connected_wts_i(i),new_diag_elems_i(i),connected_wts_i(i)**2/(var_energy-new_diag_elems_i(i)),delta_e_2pt
              endif
            endif
          endif
        enddo
        write (*, *) "end loop"

        !$omp critical
        delta_e_2pt = delta_e_2pt + delta_e_2pt_inc
        !$omp end critical

        !$omp end task


      enddo

      ndets_connected = ndets_connected + ndets_connected_i

    enddo

    !$omp end single
    !$omp end parallel

    write (6,'(''Variational energy='',t36,f15.9)') var_energy
    write (6,'(''Second-order PT energy lowering='',t36,f15.9)') delta_e_2pt
    write (6,'(''Total energy='',t36,f15.9)') var_energy+delta_e_2pt
    write (6,'(''ndets, ndets_connected, Variational, PT, Total Energies='',i9,i10,f16.9,f15.9,f16.9)') ndets, ndets_connected, var_energy, delta_e_2pt, var_energy+delta_e_2pt
    call flush(6)

    call my_second(2,'Compute 2nd-order energy correction')

  end subroutine second_order_pt_dtm_batches


  subroutine second_order_pt_multi_index_hashing(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,epsilon,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
    ! Compute and print the second order energy lowering
    use semistoch, only : find_doubly_excited,hamiltonian
    use common_run, only : connected_dets_up,connected_dets_dn,connected_matrix_elements,max_connected_dets,diag_elem_info,connected_diag_elems_info
    use chemistry, only : find_important_connected_dets_chem,get_new_diag_elem
    use types, only : i8b

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
    type(ik_vec),allocatable :: potential_connected_dets_up(:),potential_connected_dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
    integer(ik),allocatable :: potential_connected_dets_up(:),potential_connected_dets_dn(:)
#endif
    real(rk),intent(in) :: wts(:),var_energy,epsilon
    integer,intent(inout) :: iorder(:),temp_i_2(:)

    real(rk) :: delta_e_2pt
    integer :: i,j,k
    integer :: n_hash
    real(rk),allocatable :: diag_elems(:)
    type(neighbor_table) :: neighbors,ref_dets
    real(rk) :: H_ii,H_ij
    real(rk) :: old_sum
    integer :: n_connected_dets,n_connections
    real(rk),allocatable :: potential_connected_coeffs(:)
    logical :: is_con
    integer(i8b) :: n_potential, n_actual

    n_potential = 0_i8b
    n_actual = 0_i8b

    call my_second(1,'Compute 2nd-order energy correction')

    n_hash = 5 ! 1 more than the number of orbitals that change occupancy in a double excitation

    call init_hash_fn ! Sets up hash function

    allocate(potential_connected_dets_up(n_hash*ndets))
    allocate(potential_connected_dets_dn(n_hash*ndets))
    allocate(potential_connected_coeffs(n_hash*ndets))

    ! First, create neighbor tables, for putting reference dets into, so that connected dets can quickly check for potential connections
    call create_neighbor_tables(ndets,dets_up,dets_dn,n_hash,n_bins,neighbors) !,orb_group_mat)
    write (6,'(''Done creating neighbor tables'')')

    ! Note: "neighbors" has been allocated but is currently empty!
    ! This is because no reference dets have been visited yet

    ! Also, create a table of just the ref_dets, to enable us to quickly check whether a connected det should be thrown out
    ! (since the PT expression sums over only dets not in reference)
    call create_neighbor_tables(ndets,dets_up,dets_dn,1,n_bins,ref_dets)

    ! This table should be filled up with all references!
    call fill_neighbor_tables(ndets,dets_up,dets_dn,wts,ref_dets)


    ! Next, loop over all references:
    ! For each, compute diagonal elements, tthen generate all connections exceeding epsilon
    ! For each connection, query neighbors to find potential connections
    ! After all connections for a given reference det have been looped over,
    ! add that reference det to neighbors

    delta_e_2pt = 0._rk

    do i=1,ndets

      call find_important_connected_dets_chem(dets_up(i), dets_dn(i), epsilon/abs(wts(i)), n_connected_dets, connected_dets_up, connected_dets_dn,connected_matrix_elements,diag_elems(i),connected_diag_elems_info) ! connected_diag_elems_info is needed for computing diagonal elements in O(N) time, but don't use more than once per connected det!

      if (n_connected_dets==0)  cycle

      do j=2,n_connected_dets

        ! Check by a single hash table whether it is a reference det
        if (is_in_reference(connected_dets_up(j),connected_dets_dn(j),ref_dets))  cycle

        call get_potential_connections(connected_dets_up(j),connected_dets_dn(j),neighbors,n_connections,potential_connected_dets_up,potential_connected_dets_dn,potential_connected_coeffs,iorder,temp_i16_up,temp_i16_dn,temp_i_2) !,orb_group_mat)

        n_potential = n_potential + int(n_connections,i8b)

        old_sum = 0._rk
        if (n_connections>0) then
          if (is_connected(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(1),potential_connected_dets_dn(1))) then
            n_actual = n_actual + 1_i8b
            call hamiltonian(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(1),potential_connected_dets_dn(1),H_ij,is_con)
            if (abs(H_ij*potential_connected_coeffs(1))>epsilon)  old_sum = H_ij*potential_connected_coeffs(1)
          endif
          do k=2,n_connections
            if (potential_connected_dets_up(k)==potential_connected_dets_up(k-1).and.potential_connected_dets_dn(k)==potential_connected_dets_dn(k-1))  cycle
            if (.not.is_connected(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(k),potential_connected_dets_dn(k)))  cycle
            n_actual = n_actual + 1_i8b
            call hamiltonian(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(k),potential_connected_dets_dn(k),H_ij,is_con)
            if (abs(H_ij*potential_connected_coeffs(k))>epsilon)  old_sum = old_sum + H_ij*potential_connected_coeffs(k)
          enddo ! potential connections (in reference det list)
        endif

        ! get H_ii
        if (connected_diag_elems_info(j)%old_diag_elem>1.e50_rk) then ! single excitation; compute the easy way for now
          call hamiltonian(connected_dets_up(j),connected_dets_dn(j),connected_dets_up(j),connected_dets_dn(j),H_ii,is_con)
        else
          call get_new_diag_elem(connected_diag_elems_info(j),connected_dets_up(j),connected_dets_dn(j),H_ii)
        endif

        ! Update contribution from determinant j
        call hamiltonian(connected_dets_up(j),connected_dets_dn(j),dets_up(i),dets_dn(i),H_ij,is_con)
        delta_e_2pt = delta_e_2pt + (-(old_sum**2) + (old_sum + H_ij*wts(i))**2) / (var_energy - H_ii)

      enddo ! n_connected_dets

      ! Add reference i to neighbors
      call add_det_to_neighbor_tables(dets_up(i),dets_dn(i),wts(i),neighbors) !,orb_group_mat)

    enddo ! ndets

    write (6,*) "n_actual=",n_actual,"n_potential=",n_potential,"success rate=",real(n_actual)/real(n_potential)

    write (6,'(''Variational energy='',t36,f15.9)') var_energy
    write (6,'(''Second-order PT energy lowering='',t36,f15.9)') delta_e_2pt
    write (6,'(''Total energy='',t36,f15.9)') var_energy+delta_e_2pt
    call flush(6)

    call my_second(2,'Compute 2nd-order energy correction')

  end subroutine second_order_pt_multi_index_hashing


  subroutine second_order_pt_ah(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,epsilon)
    ! Compute and print the second order energy lowering
    ! For each reference, generate all connections exceeding epsilon.
    ! Then, loop over all connections to see if any are connections to references already passed through; update E_2pt appropriately
    ! When done generating connections to a reference dets, if the number of connections is >0, store all determinants
    ! that can be constructed by removing the most-excited electron that was just excited.
    ! When checking whether a connection is connected to any previous reference dets, loop over all of its electrons,
    ! and look up each resulting determinant to find any potential connections.

    ! A Holmes, 4 Apr 2016

    use semistoch, only : find_doubly_excited,hamiltonian
    use common_run, only : connected_dets_up,connected_dets_dn,connected_matrix_elements,max_connected_dets,diag_elem_info,connected_diag_elems_info
    use chemistry, only : find_important_connected_dets_chem,find_important_connected_dets_1e_removed,find_important_singly_connected_dets_1e_removed,find_important_doubly_connected_dets_2e_removed,get_new_diag_elem
    use types, only : i8b
    use common_ham, only : nelec
    use generic_sort, only : sort

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
    real(rk),intent(in) :: wts(:),var_energy,epsilon

    real(rk) :: delta_e_2pt
    integer :: i,j,k
    real(rk),allocatable :: diag_elems(:)
    type(neighbor_table) :: neighbors,ref_dets
    real(rk) :: H_ii,H_ij
    real(rk) :: old_sum
    integer :: n_connected_dets,n_connections
    real(rk),allocatable :: potential_connected_coeffs(:)
    logical :: is_con
    integer(i8b) :: n_potential, n_actual
    integer :: i_ref
    integer,allocatable :: tmp_ind_singles(:),tmp_info_singles(:)
    integer,allocatable :: tmp_ind_doubles(:),tmp_info_doubles(:)
   !integer,allocatable :: tmp_ind(:),tmp_info(:)
    integer,allocatable :: ref_connection_indices(:)
    integer :: n_allocate_singles, n_allocate_doubles
    integer :: n_hash_singles, n_hash_doubles
   !integer :: n_allocate
   !integer :: n_hash ! number of dets placed in hash table
    integer :: n_sum_terms

    call my_second(1,'Compute 2nd-order energy correction')

    call init_hash_fn ! Sets up hash function
    n_potential = 0_i8b
    n_actual = 0_i8b

    ! First, create hash tables
    n_allocate_singles = ndets*nelec
    n_allocate_doubles = ndets*nelec ! will probably need to be reallocated
    allocate(hash_ind_singles(n_allocate_singles))
    allocate(hash_info_singles(n_allocate_singles))
    allocate(hash_ind_doubles(n_allocate_doubles))
    allocate(hash_info_doubles(n_allocate_doubles))
    n_hash_singles = 0
    n_hash_doubles = 0

    ! Single excitations
    do i_ref=1,ndets
      call find_important_singly_connected_dets_1e_removed(dets_up(i_ref), dets_dn(i_ref), n_connected_dets, connected_dets_up, connected_dets_dn)
      ! Now, connected_dets represents all those 1-electron-removed singly excited dets connected to reference
      ! These are the places in the hash table where the index i_ref needs to be stored.
      do j=1,n_connected_dets
        n_hash_singles = n_hash_singles + 1
        hash_ind_singles(n_hash_singles) = hash_det(connected_dets_up(j),connected_dets_dn(j)) ! index of hash table
       !hash_info(n_hash) = i_ref ! what is stored in hash table
      enddo ! j
    enddo ! i_ref

    ! Double excitations
    do i_ref=1,ndets
      call find_important_doubly_connected_dets_2e_removed(dets_up(i_ref), dets_dn(i_ref), epsilon/abs(wts(i_ref)), n_connected_dets, connected_dets_up, connected_dets_dn)
      ! Now, connected_dets represents all those 2-electron-removed doubly excited dets connected to reference by more than epsilon
      ! These are the places in the hash table where the index i_ref needs to be stored.
      if (n_hash_doubles+n_connected_dets>n_allocate_doubles) then ! reallocate
        write (6,*) "Reallocating hash tables from size",n_allocate_doubles,"to size", 2*n_allocate_doubles; call flush(6)
        allocate(tmp_ind_doubles(n_allocate_doubles))
        allocate(tmp_info_doubles(n_allocate_doubles))
        tmp_ind_doubles(:) = hash_ind_doubles(:)
        tmp_info_doubles(:) = hash_info_doubles(:)
        deallocate(hash_ind_doubles,hash_info_doubles)
        n_allocate_doubles = 2*n_allocate_doubles
        allocate(hash_ind_doubles(n_allocate_doubles))
        allocate(hash_info_doubles(n_allocate_doubles))
        hash_ind_doubles(1:size(tmp_ind_doubles)) = tmp_ind_doubles(:)
        hash_info_doubles(1:size(tmp_info_doubles)) = tmp_info_doubles(:)
        deallocate(tmp_ind_doubles,tmp_info_doubles)
      endif
      do j=1,n_connected_dets
        n_hash_doubles = n_hash_doubles + 1
        hash_ind_doubles(n_hash_doubles) = hash_det(connected_dets_up(j),connected_dets_dn(j)) ! index of hash table
       !hash_info(n_hash) = i_ref ! what is stored in hash table
      enddo ! j
    enddo ! i_ref

    ! Now, sort the hash tables by hash_ind
    call sort(n_hash_singles,hash_ind_singles)
    call sort(n_hash_doubles,hash_ind_doubles)

    ! Now, store the start_index and end_index for each hash_ind
    allocate(start_index_singles(n_bins))
    start_index_singles(:) = 0
    start_index_singles(hash_ind_singles(1)) = 1
    do i=2,n_hash_singles
      if (hash_ind_singles(i).ne.hash_ind_singles(i-1)) then
        start_index_singles(hash_ind_singles(i)) = i
      endif
    enddo
    deallocate(hash_ind_singles)

    allocate(start_index_doubles(n_bins))
    start_index_doubles(:) = 0
    start_index_doubles(hash_ind_doubles(1)) = 1
    do i=2,n_hash_doubles
      if (hash_ind_doubles(i).ne.hash_ind_doubles(i-1)) then
        start_index_doubles(hash_ind_doubles(i)) = i
      endif
    enddo
    deallocate(hash_ind_doubles)

    ! At this point, hash_info(start_index(i):end_index(i)), where i=hash_det(det_up,det_dn), are the
    ! reference dets connected to det_up/dn

    allocate(n_per_bin_singles(n_bins))
    n_per_bin_singles(:) = 0
    hash_info_singles(:) = 0 ! Because we only want to check whether a connection is connected to a reference that has already been looped over!

    allocate(n_per_bin_doubles(n_bins))
    n_per_bin_doubles(:) = 0
    hash_info_doubles(:) = 0 ! Because we only want to check whether a connection is connected to a reference that has already been looped over!

    ! Also, create a table of just the ref_dets, to enable us to quickly check whether a connected det should be thrown out
    ! (since the PT expression sums over only dets not in reference)
    call create_neighbor_tables(ndets,dets_up,dets_dn,1,n_bins,ref_dets)

    ! This table should be filled up with all references!
    call fill_neighbor_tables(ndets,dets_up,dets_dn,wts,ref_dets)

    ! Next, loop over all references:
    ! For each, compute diagonal elements, then generate all connections exceeding epsilon
    ! For each connection, loop over all of the electrons; hash each determinant obtained by removing one electron
    ! This yields the set of reference determinants connected to the query
    ! After all connections for a given reference det have been looped over,
    ! add that reference det to neighbors

    allocate(ref_connection_indices(ndets))

    delta_e_2pt = 0._rk

    max_sum_terms = 0

    do i=1,ndets

      call find_important_connected_dets_chem(dets_up(i), dets_dn(i), epsilon/abs(wts(i)), n_connected_dets, connected_dets_up, connected_dets_dn,connected_matrix_elements,diag_elems(i),connected_diag_elems_info) ! connected_diag_elems_info is needed for computing diagonal elements in O(N) time, but don't use more than once per connected det!

      if (n_connected_dets==0)  cycle

      do j=2,n_connected_dets ! skip 1 because it is the diagonal element

        ! Check by a single hash table whether it is a reference det
        if (is_in_reference(connected_dets_up(j),connected_dets_dn(j),ref_dets))  cycle

        call get_reference_connections_less_memory(ndets,connected_dets_up(j),connected_dets_dn(j),n_connections,ref_connection_indices)
       !call get_reference_connections(connected_dets_up(j),connected_dets_dn(j),n_connections,ref_connection_indices)

        n_potential = n_potential + n_connections

        old_sum = 0._rk
        n_sum_terms = 0

        if (n_connections>0) then
          do k=1,n_connections
            ! Following check has to be in because of hash collisions
            if (.not.is_connected(connected_dets_up(j),connected_dets_dn(j),dets_up(ref_connection_indices(k)),dets_dn(ref_connection_indices(k))))  cycle
            n_actual = n_actual + 1_i8b
            call hamiltonian(connected_dets_up(j),connected_dets_dn(j),dets_up(ref_connection_indices(k)),dets_dn(ref_connection_indices(k)),H_ij,is_con)
            if (abs(H_ij*wts(ref_connection_indices(k)))>epsilon) then
              old_sum = old_sum + H_ij*wts(ref_connection_indices(k))
              n_sum_terms = n_sum_terms + 1
            endif
          enddo ! reference connections (in reference det list)
        endif

        if (n_sum_terms>max_sum_terms)  max_sum_terms = n_sum_terms

        ! get H_ii
        if (connected_diag_elems_info(j)%old_diag_elem>1.e50_rk) then ! single excitation; compute the easy way for now
          call hamiltonian(connected_dets_up(j),connected_dets_dn(j),connected_dets_up(j),connected_dets_dn(j),H_ii,is_con)
        else
          call get_new_diag_elem(connected_diag_elems_info(j),connected_dets_up(j),connected_dets_dn(j),H_ii)
        endif

        ! Update contribution from determinant j
        call hamiltonian(connected_dets_up(j),connected_dets_dn(j),dets_up(i),dets_dn(i),H_ij,is_con)
        delta_e_2pt = delta_e_2pt + (-(old_sum**2) + (old_sum + H_ij*wts(i))**2) / (var_energy - H_ii)

      enddo ! j=2,n_connected_dets

      ! Add reference i to hash table
      call add_reference_to_hash_table(i,dets_up(i),dets_dn(i),wts(i),epsilon)

    enddo ! i=1,ndets

    write (6,*) "Max terms in one sum=",max_sum_terms

    write (6,*) "n_actual=",n_actual,"n_potential=",n_potential,"success rate=",real(n_actual)/real(n_potential)

    write (6,'(''Variational energy='',t36,f15.9)') var_energy
    write (6,'(''Second-order PT energy lowering='',t36,f15.9)') delta_e_2pt
    write (6,'(''Total energy='',t36,f15.9)') var_energy+delta_e_2pt
    call flush(6)

    call my_second(2,'Compute 2nd-order energy correction')

  end subroutine second_order_pt_ah


  subroutine create_neighbor_tables(n_det,dets_up,dets_dn,n_hash,n_bins,neighbors,orb_group_mat)
  ! Create a structure to be used for approximate string matching.
  ! It should be empty for computing 2PT energy, so 
  ! when used for constructing H quickly, call a separate
  ! subroutine to fill it.
  ! orb_group_mat(i,j) is the j'th orbital in the i'th group; if not
  ! present, then just divide up the orbitals into equal-sized groups
  ! by putting every 5th one in a different group
  ! A Holmes, 22 Mar 2016

    integer,intent(in) :: n_det
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
    integer,intent(in) :: n_hash,n_bins ! want n_hash = 5, n_bins ~ n_det
    type(neighbor_table),intent(out) :: neighbors
    integer,optional,intent(in) :: orb_group_mat(:,:)

    integer,allocatable :: bin_counts(:,:)
    integer :: i,j
    integer,allocatable :: hashes(:)

    neighbors%n_hash_functions = n_hash
    neighbors%n_bins_per_hash = n_bins

    allocate(hashes(n_hash))

    allocate(neighbors%is_empty(n_hash,n_bins))
    allocate(neighbors%bin_counts(n_hash,n_bins))
    allocate(neighbors%start_index(n_hash,n_bins))
    allocate(neighbors%end_index(n_hash,n_bins))
    allocate(neighbors%dets_up(n_hash,n_det))
    allocate(neighbors%dets_dn(n_hash,n_det))
    allocate(neighbors%coeffs(n_hash,n_det))

    ! First, count number of dets that hash into each bin
    allocate(bin_counts(n_hash,n_bins))
    bin_counts(:,:) = 0
    if (present(orb_group_mat)) then
      do i=1,n_det
        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes,orb_group_mat)
        do j=1,n_hash
          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
        enddo
      enddo
    else 
      do i=1,n_det
        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes)
        do j=1,n_hash
          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
        enddo
      enddo
    endif
      
    ! Then, "allocate" appropriately, i.e., determine which indices correspond to which bins of hash table
    do j=1,n_hash
      neighbors%is_empty(j,1) = (bin_counts(j,1)==0)
      neighbors%start_index(j,1) = 1
      neighbors%end_index(j,1) = bin_counts(j,1)
      do i=2,n_bins
        neighbors%is_empty(j,i) = (bin_counts(j,i)==0)
        neighbors%start_index(j,i) = neighbors%end_index(j,i-1) + 1
        neighbors%end_index(j,i) = neighbors%end_index(j,i-1) + bin_counts(j,i)
      enddo
    enddo

    neighbors%is_empty(:,:) = .true. ! Because the neighbors table is currently empty! (call fill_neighbor_tables afterwards if needed)
    neighbors%bin_counts(:,:) = 0

  end subroutine create_neighbor_tables


  subroutine fill_neighbor_tables(n_det,dets_up,dets_dn,coeffs,neighbors,orb_group_mat)
  ! Fill neighbor tables with dets_up/dn
  ! Assumes neighbor tables have already been created with create_neighbor_tables
  ! A Holmes, 22 Mar 2016

    integer,intent(in) :: n_det
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
    real(rk),intent(in) :: coeffs(:)
    type(neighbor_table),intent(inout) :: neighbors
    integer,optional,intent(in) :: orb_group_mat(:,:)

    integer,allocatable :: bin_counts(:,:)
    integer :: i,j,ind
    integer :: n_hash
    integer :: n_bins
    integer,allocatable :: hashes(:)

    n_hash = neighbors%n_hash_functions
    n_bins = neighbors%n_bins_per_hash

    allocate(hashes(n_hash))

    allocate(bin_counts(n_hash,n_bins))
    bin_counts(:,:) = 0
    if (present(orb_group_mat)) then
      do i=1,n_det
        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes,orb_group_mat)
        do j=1,n_hash
          ind = neighbors%start_index(j,hashes(j)) + bin_counts(j,hashes(j))
          neighbors%dets_up(j,ind) = dets_up(i)
          neighbors%dets_dn(j,ind) = dets_dn(i)
          neighbors%coeffs(j,ind) = coeffs(i)
          neighbors%is_empty(j,hashes(j)) = .false.
          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
        enddo
      enddo
    else 
      do i=1,n_det
        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes)
        do j=1,n_hash
          ind = neighbors%start_index(j,hashes(j)) + bin_counts(j,hashes(j))
          neighbors%dets_up(j,ind) = dets_up(i)
          neighbors%dets_dn(j,ind) = dets_dn(i)
          neighbors%coeffs(j,ind) = coeffs(i)
          neighbors%is_empty(j,hashes(j)) = .false.
          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
        enddo
      enddo
    endif

    neighbors%bin_counts = bin_counts
      
  end subroutine fill_neighbor_tables


  subroutine add_det_to_neighbor_tables(det_up,det_dn,coeff,neighbors,orb_group_mat)
  ! Add a single det to neighbor tables
  ! Assumes this is one of the dets that was known about when
  ! create_neighbor_tables was called, and that this det has
  ! not yet been placed in the table!
  ! A Holmes, 22 Mar 2016

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
#else
    integer(ik),intent(in) :: det_up,det_dn
#endif
    real(rk),intent(in) :: coeff
    type(neighbor_table),intent(inout) :: neighbors
    integer,optional,intent(in) :: orb_group_mat(:,:)

    integer :: n_hash
    integer :: n_bins
    integer :: ind,j
    integer,allocatable :: hashes(:)

    n_hash = neighbors%n_hash_functions
    n_bins = neighbors%n_bins_per_hash

    allocate(hashes(n_hash))

    if (present(orb_group_mat)) then
      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes,orb_group_mat)
      do j=1,n_hash
        ind = neighbors%start_index(j,hashes(j)) + neighbors%bin_counts(j,hashes(j))
        neighbors%dets_up(j,ind) = det_up
        neighbors%dets_dn(j,ind) = det_dn
        neighbors%coeffs(j,ind) = coeff
        neighbors%is_empty(j,hashes(j)) = .false.
        neighbors%bin_counts(j,hashes(j)) = neighbors%bin_counts(j,hashes(j)) + 1
      enddo
    else 
      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes)
      do j=1,n_hash
        ind = neighbors%start_index(j,hashes(j)) + neighbors%bin_counts(j,hashes(j))
        neighbors%dets_up(j,ind) = det_up
        neighbors%dets_dn(j,ind) = det_dn
        neighbors%coeffs(j,ind) = coeff
        neighbors%is_empty(j,hashes(j)) = .false.
        neighbors%bin_counts(j,hashes(j)) = neighbors%bin_counts(j,hashes(j)) + 1
      enddo
    endif


  end subroutine add_det_to_neighbor_tables


  subroutine get_hashes(det_up,det_dn,n_hash,n_bins,hashes,orb_group_mat)
  ! Get hashes of all the subsets of the orbitals of det_up/dn.
  ! If orb_group_mat not present, then assume the orbitals alternate their groups
  ! A Holmes, 22 Mar 2016

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec),allocatable :: pieces(:)
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik),allocatable :: pieces(:)
#endif
    integer,intent(in) :: n_hash,n_bins
    integer,intent(out) :: hashes(:)
    integer,optional,intent(in) :: orb_group_mat(:,:)

    integer :: i,k,bit

    ! First, divide bits up into different groups
    allocate(pieces(n_hash))
    pieces(:) = 0_ik
    bit = 0
    if (present(orb_group_mat)) then
      stop "orb_group_mat not implemented yet!"
    else
      do i=1,2*norb
        k = mod(i-1,n_hash)+1
        if (k==n_hash)  bit = bit + 1
        if (i<=norb) then
          if (btest(det_up,i-1))  pieces(k) = ibset(pieces(k),bit)
        else
          if (btest(det_dn,i-norb-1))  pieces(k) = ibset(pieces(k),bit)
        endif
      enddo
    endif

    ! Then, apply hash function to each group
    do i=1,n_hash
      call hash(pieces(i),n_bins,hashes(i),RandomHash1)
    enddo

  end subroutine get_hashes


  integer function hash_det(det_up,det_dn)
  ! Hash the determinant det_up/dn, into a hash table with n_bins bins
  ! A Holmes, 5 Apr 2016

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
#else
    integer(ik),intent(in) :: det_up,det_dn
#endif

    call hash_2_dets(det_up,det_dn,n_bins,hash_det,RandomHash1)

  end function hash_det


  subroutine get_potential_connections(det_up,det_dn,neighbors,n_connections,connected_dets_up,connected_dets_dn,connected_coeffs,iorder,temp_i16_up,temp_i16_dn,temp_i_2,orb_group_mat)
  ! Query the neighbor tables (which may not be completely filled in!) with a query
  ! determinant, det_up/dn. Return a set of (potential) connected determinants and
  ! their corresponding coefficients
  ! Returns them sorted by det label, but there may be repeats
  ! A Holmes, 22 Mar 2016

    use tools, only : merge_sort2_up_dn

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec),intent(out) :: connected_dets_up(:),connected_dets_dn(:)
    type(ik_vec),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik),intent(out) :: connected_dets_up(:),connected_dets_dn(:)
    integer(ik),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
    type(neighbor_table),intent(in) :: neighbors
    integer,intent(out) :: n_connections
    real(rk),intent(out) :: connected_coeffs(:)
    integer,intent(inout) :: iorder(:),temp_i_2(:)
    integer,optional,intent(in) :: orb_group_mat(:,:)

    integer :: i,ind
    integer,allocatable :: hashes(:)
    integer :: n_hash
    integer :: n_bins

    n_hash = neighbors%n_hash_functions
    n_bins = neighbors%n_bins_per_hash

    allocate(hashes(n_hash))

    n_connections = 0

    if (present(orb_group_mat)) then
      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes,orb_group_mat)
    else
      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes)
    endif

    do i=1,n_hash
      ind = hashes(i)
      if (neighbors%is_empty(i,ind))  cycle
      connected_dets_up(n_connections+1:n_connections+neighbors%bin_counts(i,ind)) = neighbors%dets_up(i,neighbors%start_index(i,ind):neighbors%start_index(i,ind)+neighbors%bin_counts(i,ind)-1)
      connected_dets_dn(n_connections+1:n_connections+neighbors%bin_counts(i,ind)) = neighbors%dets_dn(i,neighbors%start_index(i,ind):neighbors%start_index(i,ind)+neighbors%bin_counts(i,ind)-1)
      connected_coeffs(n_connections+1:n_connections+neighbors%bin_counts(i,ind)) = neighbors%coeffs(i,neighbors%start_index(i,ind):neighbors%start_index(i,ind)+neighbors%bin_counts(i,ind)-1)
      n_connections = n_connections + neighbors%bin_counts(i,ind)
    enddo

    if (n_connections==0)  return

    ! Finally, sort and merge to remove repeats!
    do i=1,n_connections
      iorder(i) = i
    enddo
    call merge_sort2_up_dn(connected_dets_up(1:n_connections),connected_dets_dn(1:n_connections), iorder, n_connections, temp_i16_up, temp_i16_dn, temp_i_2)
    connected_coeffs(1:n_connections) = connected_coeffs(iorder(1:n_connections)) 
 ! MERGE! (TODO)

  end subroutine get_potential_connections


!  subroutine get_reference_connections(det_up,det_dn,n_connections,ref_connection_indices)
!  ! Return the set of reference determinants connected to det_up/dn
!  ! All determinants will be either single or double excitations (unless there are hash collisions)
!  ! Time complexity is O(N_elec + N_connections*(1+N_hash_collisions))
!  ! A Holmes, 5 Apr 2016
!
!    use chemistry, only : occ_up,occ_dn
!    use tools, only : sort_and_merge
!
!    integer,intent(out) :: n_connections
!    integer,intent(out) :: ref_connection_indices(:)
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec) :: tmp
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!    integer(ik) :: tmp
!#endif
!
!    integer :: i_elec,ind,j
!
!    ! Get occ_up, occ_dn in O(N_elec) time
!    tmp = det_up
!    do i_elec=1,nup
!      occ_up(i_elec) = trailz(tmp) + 1
!      tmp = ibclr(tmp,occ_up(i_elec)-1)
!    enddo
!    tmp = det_dn
!    do i_elec=1,ndn
!      occ_dn(i_elec) = trailz(tmp) + 1
!      tmp = ibclr(tmp,occ_dn(i_elec)-1)
!    enddo
!
!    ! Loop over occupied electron lists; for each, remove that electron and hash to look up connections
!    n_connections = 0
!    do i_elec=1,nup
!      ind = hash_det(ibclr(det_up,occ_up(i_elec)-1),det_dn)
!      if (n_per_bin(ind)==0)  cycle
!      do j=1,n_per_bin(ind)
!        n_connections = n_connections + 1
!        ref_connection_indices(n_connections) = hash_info(start_index(ind)+j-1)
!      enddo ! j
!    enddo ! i_elec
!    do i_elec=1,ndn
!      ind = hash_det(det_up,ibclr(det_dn,occ_dn(i_elec)-1))
!      if (n_per_bin(ind)==0)  cycle
!      do j=1,n_per_bin(ind)
!        n_connections = n_connections + 1
!        ref_connection_indices(n_connections) = hash_info(start_index(ind)+j-1)
!      enddo ! j
!    enddo ! i_elec
!
!    ! Finally, we must sort and merge because there can be up to
!    ! N_elec repeats of each single excitation and up to 
!    ! 2 repeats of each double excitation!
!
!    call sort_and_merge(n_connections,ref_connection_indices)
!
!  end subroutine get_reference_connections


  subroutine get_reference_connections_less_memory(n_ref,det_up,det_dn,n_connections,ref_connection_indices)
  ! Works like get_reference_connections, but uses a factor of N_orb less memory
  ! In worst case, takes a factor of N_elec more time, but since it avoids the problem of having up to
  ! N_elec duplicates that get_reference_connections has, it should not be too much slower in practice
  ! It also avoids the need to sort and merge because there are no duplicates
  ! A Holmes, 6 Apr 2016

    use chemistry, only : occ_up,occ_dn

    integer,intent(in) :: n_ref
    integer,intent(out) :: n_connections
    integer,intent(out) :: ref_connection_indices(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec) :: tmp
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik) :: tmp
#endif

    integer :: i_elec
    integer :: n_single, n_double

    ! Get occ_up, occ_dn in O(N_elec) time
    tmp = det_up
    do i_elec=1,nup
      occ_up(i_elec) = trailz(tmp) + 1
      tmp = ibclr(tmp,occ_up(i_elec)-1)
    enddo
    tmp = det_dn
    do i_elec=1,ndn
      occ_dn(i_elec) = trailz(tmp) + 1
      tmp = ibclr(tmp,occ_dn(i_elec)-1)
    enddo

    call get_singly_connected_ref_dets(det_up,det_dn,n_single,ref_connection_indices,occ_up,occ_dn,hash_info_singles,n_per_bin_singles,start_index_singles)
    call get_doubly_connected_ref_dets(det_up,det_dn,n_double,ref_connection_indices(n_single+1:n_ref),occ_up,occ_dn,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
   !call get_singly_connected_ref_dets(det_up,det_dn,n_single,ref_connection_indices,occ_up,occ_dn)
   !call get_doubly_connected_ref_dets(det_up,det_dn,n_double,ref_connection_indices(n_single+1:n_ref),occ_up,occ_dn)

    n_connections = n_single + n_double

  end subroutine get_reference_connections_less_memory


  subroutine get_singly_connected_ref_dets(det_up,det_dn,n_connections,ref_connection_indices,occ_up,occ_dn,hash_info_singles,n_per_bin_singles,start_index_singles)
  ! Return the set of singly excited reference determinants connected to det_up/dn
  ! All determinants will be single excitations (unless there are hash collisions)
  ! There should be no repeats, so no need to sort and merge
  ! Time complexity is O(N_elec + N_connections*(1+N_hash_collisions))
  ! A Holmes, 6 Apr 2016

    integer,intent(out) :: n_connections
    integer,intent(out) :: ref_connection_indices(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec) :: tmp
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik) :: tmp
#endif
    integer,intent(in) :: occ_up(:),occ_dn(:)
    integer,intent(in) :: hash_info_singles(:),n_per_bin_singles(:),start_index_singles(:)

    integer :: i_elec,ind,j

    ! Loop over occupied electron lists; for each, remove that electron and hash to look up connections
    n_connections = 0
    do i_elec=1,nup
      ind = hash_det(ibclr(det_up,occ_up(i_elec)-1),det_dn)
      if (n_per_bin_singles(ind)==0)  cycle
      do j=1,n_per_bin_singles(ind)
        n_connections = n_connections + 1
        ref_connection_indices(n_connections) = hash_info_singles(start_index_singles(ind)+j-1)
      enddo ! j
    enddo ! i_elec
    do i_elec=1,ndn
      ind = hash_det(det_up,ibclr(det_dn,occ_dn(i_elec)-1))
      if (n_per_bin_singles(ind)==0)  cycle
      do j=1,n_per_bin_singles(ind)
        n_connections = n_connections + 1
        ref_connection_indices(n_connections) = hash_info_singles(start_index_singles(ind)+j-1)
      enddo ! j
    enddo ! i_elec

  end subroutine get_singly_connected_ref_dets


  subroutine get_doubly_connected_ref_dets(det_up,det_dn,n_connections,ref_connection_indices,occ_up,occ_dn,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
  ! Return the set of doubly excited reference determinants connected to det_up/dn
  ! All determinants will be double excitations (unless there are hash collisions)
  ! There should be no repeats, so no need to sort and merge
  ! Time complexity is O(N_elec^2 + N_connections*(1+N_hash_collisions))
  ! Warning: also returns N_elec copies of each singly connected ref det,
  !          which must be filtered out after calling this routine!
  ! A Holmes, 6 Apr 2016

    use chemistry, only : pairs_e1,pairs_e2

    integer,intent(out) :: n_connections
    integer,intent(out) :: ref_connection_indices(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec) :: tmp
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik) :: tmp,new_up,new_dn
#endif
    integer,intent(in) :: occ_up(:),occ_dn(:)
    integer,intent(in) :: hash_info_doubles(:),n_per_bin_doubles(:),start_index_doubles(:)

    integer :: i_elec,ind,j
    integer :: ipair,npairs
    integer :: p,q

    npairs = 0
    ! Same spin, up
    do p=1,nup-1
      do q=p+1,nup
        npairs = npairs + 1
        pairs_e1(npairs) = occ_up(p)
        pairs_e2(npairs) = occ_up(q)
      enddo
    enddo
    ! Same spin, dn
    do p=1,ndn-1
      do q=p+1,ndn
        npairs = npairs + 1
        pairs_e1(npairs) = occ_dn(p)+norb
        pairs_e2(npairs) = occ_dn(q)+norb
      enddo
    enddo
    ! Opposite spin
    do p=1,nup
      do q=1,ndn
        npairs = npairs + 1
        pairs_e1(npairs) = occ_up(p)
        pairs_e2(npairs) = occ_dn(q)+norb
      enddo
    enddo

    ! Loop over occupied electron pairs; for each, remove that electron and hash to look up connections

    n_connections = 0
    do ipair=1,npairs

      p = pairs_e1(ipair)
      q = pairs_e2(ipair)

      ! Commenting out for now, since I think we want to retrieve all partial connections that have been stored:
     !epair = int(combine_2_indices(p,q))
     !if (dtm_hb(epair,1)%absH<min_elem)  exit ! If no double excitations that start from p,q are larger than epsilon, then skip this p,q pair

      ! Generate new determinants by removing p and q
      new_up = det_up
      new_dn = det_dn
      if (p<=norb) then
        new_up = ibclr(new_up,p-1)
      else
        new_dn = ibclr(new_dn,p-norb-1)
      endif
      if (q<=norb) then
        new_up = ibclr(new_up,q-1)
      else
        new_dn = ibclr(new_dn,q-norb-1)
      endif

      ind = hash_det(new_up,new_dn)

      if (n_per_bin_doubles(ind)==0)  cycle

      do j=1,n_per_bin_doubles(ind)
        n_connections = n_connections + 1
        ref_connection_indices(n_connections) = hash_info_doubles(start_index_doubles(ind)+j-1)
      enddo ! j

    enddo ! ipair (p,q)

  end subroutine get_doubly_connected_ref_dets


  logical function is_connected(up1,dn1,up2,dn2)
  ! Returns true if the excitation level between dets 1 and 2 is <=2 and >0
  ! false otherwise
  ! A Holmes, 22 Mar 2016

    use chemistry, only : excitation_level

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: up1,dn1,up2,dn2
#else
    integer(ik),intent(in) :: up1,dn1,up2,dn2
#endif

    integer :: excite_level

    call excitation_level(up1,dn1,up2,dn2,excite_level)
    is_connected = (excite_level>0)

  end function is_connected


  subroutine init_hash_fn
    ! These are random numbers:
    RandomHash1 = (/4735613, 313959, 650067,2918183,2123508,1824383, 368641,4665125, 889838,1653649, 567626, 700470,2072070,3265229,1218533,1765295,1278102,2587196,3192811,4079255,3370036,5073942,1910990,2852420,2908100,1267516,1016462,3757537,4184599, 304348,2123593,4425848,3581496,1791107,3497983,1391199,2198178,2947481,3044171,4776221,2298252,3308124,4255082, 756103,3464666,3571942,3888499,1721390,4148643, 221259,1768488,1669008,1336151, 899703,5072037,2064427,2291053,2385335,3980263,3427668,4590425,3130633,1795017,1298250, 437459,1125424,1810093,5031931,1982167,4491159,1664669,3584140,1229843,1109613,2897092,2456647,2870772,5057078,2580274, 246617,3489605, 743439,2495387,3268699, 945352,2387627,4946247,1554925,1258363,1796338,3509947,3914421,2572562, 907556, 455356,2301574,4622378,2526352,4954820, 850178, 535417, 286038,3498667,3241142,2524137,2937664,2249144,4947862,3411414,4485867,4610744,  17755,4813498,4686397,4337008, 178525,1955778,2093500,3448090,4787430,3377083,  43503,4225437, 761386,5010808,3790117,1625775,5031659,2308473,4029366,1280477,4861902,3954799,2985978,2901621,3061381,1543058, 653181,4504563,2731597,2634156,2414947,1399327, 425706,2280240,2920945,1939056,3003820,2989859,2873210,5001290,3667487,1452012,2119338,4194025,1593590,3594528,3385764,4118903, 308089, 170759,1138565, 903090,4696165,1358916,2896407,1462530,3696410,2300452,3165629,2134285,3437968,1946733,2601552,1185301,3922477,2544259, 260349,2558668,3904972,4254275,3685299,1264351,4273640,2627883,2919797,4986566,1057643,2066436,2962802, 610070,3901649,1381809,2342696,3827048, 138016,1993221,3701099, 658366, 871615,1086936,2323667,1609455,3283476,3220993,4161146,1102928,4279104,1346917,5038082,4026472,1083272,5080232,1091946,4863898,1349339,3686856,2421638, 616475, 147825,  36701,4707242,2266338,3872717, 651830, 404626, 498554,3907313, 855440,2487375,2716318,4383594,1890043,4229224, 506837,1815062,2767149,1773485, 614411,2654510,4654238,5059003,  41437,3902708,1384039,1358399,4635825,3555781,3241332,4683214,1461543,3618113,4767352,5028089,2733360,2020587/)
  end subroutine init_hash_fn


  subroutine hash(det,range,hashx,RandomHash)
  ! hash function borrowed from mpi_routines, which I believe is
  ! a Pearson hash
  ! returns hashx, which is between 1 and range

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det
#else
    integer(ik), intent(in) :: det
#endif
    integer, intent(in) :: range,RandomHash(0:255)
    integer, intent(out) :: hashx
    integer(ik) :: acc
    integer(ik) :: test_nk
    integer :: test_int,i,j
    integer :: val
   !integer, parameter :: n_bits_nk = bit_size(test_nk)
   !integer, parameter :: n_bits = bit_size(test_int)

    acc = 0
#ifdef NUM_ORBITALS_GT_127
    do i=1,num_words    !run over integers defining the bit-string
        do j = 0, norb -1, 8
       !do j = 0, n_bits_nk -1, 8
            val = int(iand(ishft(det%v(i)*368296850_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
            !val is now a number between 0 -> 255
            !Ensure that RandomHash has a mapping for the 0th element too
            !1099511628211_ik = ibset(0_ik,27)
            acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
        enddo
    enddo
#else
    do j = 0, norb -1, 8
   !do j = 0, n_bits_nk -1, 8
        val = int(iand(ishft(det*368296850_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
        !val is now a number between 0 -> 255
        !Ensure that RandomHash has a mapping for the 0th element too
        !1099511628211_ik = ibset(0_ik,27)
        acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
    enddo
#endif
    hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))+1

  end subroutine hash


  subroutine hash_2_dets(det1,det2,range,hashx,RandomHash)
  ! Borrowed from mpi_routines.f90

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det1,det2
#else
    integer(ik), intent(in) :: det1,det2
#endif
    integer, intent(in) :: range,RandomHash(0:255)
    integer, intent(out) :: hashx
    integer(ik) :: acc
    integer(ik) :: test_nk
    integer :: test_int,i,j
    integer :: val
   !integer, parameter :: n_bits_nk = bit_size(test_nk)

    acc = 0
    i = 1
    do j = 0, norb -1, 8
   !do j = 0, n_bits_nk -1, 8
        val = int(iand(ishft(det1*368296850_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
        !val is now a number between 0 -> 255
        !Ensure that RandomHash has a mapping for the 0th element too
        !1099511628211_ik = ibset(0_ik,27)
        acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
    enddo
    i = 2
    do j = 0, norb -1, 8
   !do j = 0, n_bits_nk -1, 8
        val = int(iand(ishft(det2*368296850_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
        !val is now a number between 0 -> 255
        !Ensure that RandomHash has a mapping for the 0th element too
        !1099511628211_ik = ibset(0_ik,27)
        acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
    enddo
    hashx = int(abs(mod(acc,int(range,ik))),kind(test_int)) + 1

  end subroutine hash_2_dets


  logical function is_in_reference(det_up,det_dn,ref_dets)
  ! Check whether det is in ref_dets in O(1+N_hash_collisions) time
  ! A Holmes, 23 Mar 2016

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
#else
    integer(ik), intent(in) :: det_up,det_dn
#endif
    type(neighbor_table),intent(in) :: ref_dets

    integer :: n_bins
    integer :: hashes(1)
    integer :: i,ind

    n_bins = ref_dets%n_bins_per_hash

    call get_hashes(det_up,det_dn,1,n_bins,hashes)

    ind = hashes(1)
    if (ref_dets%is_empty(1,ind)) then
      is_in_reference = .false.
      return
    endif

    do i=1,ref_dets%bin_counts(1,ind) ! Only iterates over more than 1 element if there are hash collisions!
      if (det_up==ref_dets%dets_up(1,ref_dets%start_index(1,ind)+i-1).and.det_dn==ref_dets%dets_dn(1,ref_dets%start_index(1,ind)+i-1)) then
        is_in_reference = .true.
        return
      endif
    enddo

    is_in_reference = .false.

  end function is_in_reference


  subroutine add_reference_to_hash_table(i_ref,det_up,det_dn,wt,epsilon)
  ! Add the reference det_up/dn to the hash table
  ! A Holmes, 5 Apr 2016

    use common_run, only : connected_dets_up, connected_dets_dn
    use chemistry, only : find_important_connected_dets_1e_removed,find_important_singly_connected_dets_1e_removed,find_important_doubly_connected_dets_2e_removed

    integer,intent(in) :: i_ref ! because we only want to store the index of the reference, not the actual configuration
    real(rk),intent(in) :: wt,epsilon
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec) :: tmp
#else
    integer(ik), intent(in) :: det_up,det_dn
    integer(ik) :: tmp
#endif

    integer :: ind,i_elec,j
    integer :: n_connected_dets

      call find_important_singly_connected_dets_1e_removed(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn)
     
      do j=1,n_connected_dets
        ind = hash_det(connected_dets_up(j),connected_dets_dn(j))
        hash_info_singles(start_index_singles(ind)+n_per_bin_singles(ind)) = i_ref
        n_per_bin_singles(ind) = n_per_bin_singles(ind) + 1
      enddo ! j

      call find_important_doubly_connected_dets_2e_removed(det_up, det_dn, epsilon/abs(wt), n_connected_dets, connected_dets_up, connected_dets_dn)
     
      do j=1,n_connected_dets
        ind = hash_det(connected_dets_up(j),connected_dets_dn(j))
        hash_info_doubles(start_index_doubles(ind)+n_per_bin_doubles(ind)) = i_ref
        n_per_bin_doubles(ind) = n_per_bin_doubles(ind) + 1
      enddo ! j

     !call find_important_connected_dets_1e_removed(det_up, det_dn, epsilon/abs(wt), n_connected_dets, connected_dets_up, connected_dets_dn)
     
     !do j=1,n_connected_dets
     !  ind = hash_det(connected_dets_up(j),connected_dets_dn(j))
     !  hash_info(start_index(ind)+n_per_bin(ind)) = i_ref
     !  n_per_bin(ind) = n_per_bin(ind) + 1
     !enddo ! j

  end subroutine add_reference_to_hash_table

!=====================================================================================================================
 subroutine matrix_lanczos_partial_connections(n,dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,hash_info_singles,n_per_bin_singles,start_index_singles,hash_info_doubles,n_per_bin_doubles,start_index_doubles,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
 ! A Holmes, 10 Apr 2016. Similar to Hitesh's Lanczos routine, but for the case where instead of storing the 
 !                        Hamiltonian, we only store the partial connections (i.e., all ways to remove 1 or 2
 !                        electrons) of all the reference determinants
!=====================================================================================================================

 use chemistry, only : occ_up,occ_dn
 use common_ham, only : hamiltonian_type,nelec
 use common_run, only : max_connected_dets
 use semistoch, only : hamiltonian

 implicit none

 ! Dummy
 integer,intent(in)          :: n
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
 type(ik_vec) :: tmp
#else
 integer(ik),intent(in) :: dets_up(:),dets_dn(:)
 integer(ik) :: tmp
#endif
 real(rk),intent(out)        :: lowest_eigenvector(:)
 real(rk),intent(out)        :: lowest_eigenvalue
 integer,intent(in)          :: hash_info_singles(:),n_per_bin_singles(:),start_index_singles(:)
 integer,intent(in)          :: hash_info_doubles(:),n_per_bin_doubles(:),start_index_doubles(:)
 real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue
 real(rk),intent(in),optional :: initial_vector(:)

 ! Local
 integer                    :: i,it,i_elec,j
!real(rk)                   :: rannyu
 real(rk)                   :: energy_shift
 real(rk)                   :: norm,norm_inv
 real(rk),allocatable       :: w(:),v(:,:)
 real(rk),allocatable       :: alphas(:),betas(:)
 integer                    :: iterations
 real(rk)                   :: lowest_eigenvalue_prev
 logical                    :: converged=.false.
 integer                    :: len_work,info
 real(rk),allocatable       :: work(:),eigenvalues(:),tridiag(:,:)
 real(rk)                   :: lanczos_epsilon=1.e-15_rk
 integer,allocatable        :: ref_con_indices(:)
 integer                    :: n_con
 real(rk)                   :: H_ij
 logical                    :: is_con

  allocate(ref_con_indices(max_connected_dets+(nelec*(nelec+1))/2)) ! Because there can be up to nelec(nelec-1)/2 copies of the diagonal element coming from doubles and up to nelec copies coming from singles
  iterations=100          ! User option
  iterations=min(n,iterations)
  allocate (v(n,iterations+1))
  allocate (w(n))
  allocate(alphas(iterations+1))
  allocate(betas(iterations+1))
  w(:)=0._rk

  if (present(initial_vector)) then
    norm = 1._rk/sqrt(dot_product(initial_vector,initial_vector))
    v(:,1) = norm*initial_vector(:)
  else
    v(:,1)=0
    v(1,1)=1
  endif

  energy_shift=0._rk
  betas(1)=0._rk
  allocate (tridiag(iterations,iterations))

  allocate(eigenvalues(iterations))
  len_work = 3*iterations-1
  allocate(work(len_work))

  converged=.false.

  write(6,'(/,''Executing matrix_lanczos_partial_connections in more_tools.f90'')')
  if (n>1) then
      do it=1,iterations
         ! w(:) = H*v(:,it)
         do i=1,n

           ! Get occ_up, occ_dn in O(N_elec) time
           tmp = dets_up(i)
           do i_elec=1,nup
             occ_up(i_elec) = trailz(tmp) + 1
             tmp = ibclr(tmp,occ_up(i_elec)-1)
           enddo
           tmp = dets_dn(i)
           do i_elec=1,ndn
             occ_dn(i_elec) = trailz(tmp) + 1
             tmp = ibclr(tmp,occ_dn(i_elec)-1)
           enddo

           ! Diagonal element
           call hamiltonian(dets_up(i),dets_dn(i),dets_up(i),dets_dn(i),H_ij,is_con)
           w(i) = w(i) + H_ij*v(i,it)

           ! Single excitations
           call get_singly_connected_ref_dets(dets_up(i),dets_dn(i),n_con,ref_con_indices,occ_up,occ_dn,hash_info_singles,n_per_bin_singles,start_index_singles)
           do j=1,n_con
             if (i==ref_con_indices(j))  cycle ! already did diagonal element
             call hamiltonian(dets_up(i),dets_dn(i),dets_up(ref_con_indices(j)),dets_dn(ref_con_indices(j)),H_ij,is_con)
             if (is_con)  w(ref_con_indices(j)) = w(ref_con_indices(j)) + H_ij*v(i,it)
           enddo

           ! Double excitations
           call get_doubly_connected_ref_dets(dets_up(i),dets_dn(i),n_con,ref_con_indices,occ_up,occ_dn,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
           do j=1,n_con
             if (i==ref_con_indices(j))  cycle ! already did diagonal element
             ! Now filter out single excitations:
             if (dets_up(i)==dets_up(ref_con_indices(j))) then
               if (popcnt(iand(dets_dn(i),not(dets_dn(ref_con_indices(j)))))==1)  cycle
             elseif (dets_dn(i)==dets_dn(ref_con_indices(j))) then
               if (popcnt(iand(dets_up(i),not(dets_up(ref_con_indices(j)))))==1)  cycle
             endif
             call hamiltonian(dets_up(i),dets_dn(i),dets_up(ref_con_indices(j)),dets_dn(ref_con_indices(j)),H_ij,is_con)
             if (is_con)  w(ref_con_indices(j)) = w(ref_con_indices(j)) + H_ij*v(i,it)
           enddo

         enddo ! i

         if (it .gt. 1) w(:)=w(:)-betas(it)*v(:,it-1)
         alphas(it)=dot_product(w,v(:,it))
         w(:)=w(:)-alphas(it)*v(:,it)
         norm=dot_product(w,w)
         if (norm<(1.e-12_rk))  converged=.true.
         betas(it+1)=norm**(0.5_rk)
         norm_inv=1._rk/betas(it+1)
         v(:,it+1)=w(:)*norm_inv
         w(:)=v(:,it+1)
         do i=1,it        ! Reorthogonalization
             norm=dot_product(v(:,it+1),v(:,i))
             call flush(6)
             w(:)=w(:)-norm*v(:,i)
         enddo
         v(:,it+1)=w(:)
         w(:)=0._rk
         norm=dot_product(v(:,it+1),v(:,it+1))
         norm_inv=1._rk/(norm**(0.5_rk))
         v(:,it+1)=v(:,it+1)*norm_inv
         tridiag(:,:)=0._rk

         eigenvalues(:)=0._rk

         do i=1,it
             tridiag(i,i)=alphas(i)
             if (i<it) then
                 tridiag(i,i+1)=betas(i+1)
                 tridiag(i+1,i)=betas(i+1)
             endif
         enddo

         !diagonalize with lapack routine
         len_work = 3*it-1

         call dsyev('V', 'U', it, tridiag(1:it,1:it), it, eigenvalues, work, len_work, info)

         !deallocate(tridiag) !deallocate(work)
         lowest_eigenvalue=eigenvalues(1)
         if (present(highest_eigenvalue))  highest_eigenvalue=eigenvalues(it)
         if (present(second_lowest_eigenvalue).and.it>1)  second_lowest_eigenvalue=eigenvalues(2)
         !call print_real_matrix(size(eigenvalues),1,eigenvalues)
         if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<lanczos_epsilon) then
             converged=.true.
             exit
         else
             lowest_eigenvalue_prev=lowest_eigenvalue
             write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue
             call flush(6)
         endif
         if (converged)  exit

      enddo

      it=min(it,iterations)
      write(6,'(''matrix_lanczos_partial_connections: n, Lowest eigenvalue ='',i10, f16.10)') n, lowest_eigenvalue
      call flush(6)

      v(:,1)=matmul(v(:,1:it),tridiag(1:it,1))

      if (allocated(eigenvalues)) deallocate(eigenvalues)
      if (allocated(tridiag))     deallocate(tridiag)
      if (allocated(work))        deallocate(work)

  else
      write (6,'(''Diagonalization attempted with n=1'')')
      call hamiltonian(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),lowest_eigenvalue,is_con)
  endif

  lowest_eigenvector(1:n)=v(1:n,1)

  end subroutine matrix_lanczos_partial_connections


  subroutine get_1rdm(ndet,dets_up_in,dets_dn_in,coeffs_in,rdm)
  ! Gets the 1-RDM of the variational wavefunction
  ! A Holmes, 20 Jun 2016

    use tools, only : merge_sort2_up_dn
    use more_tools, only : get_occ_orbs,binary_search

    integer,intent(in) :: ndet
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up_in(:),dets_dn_in(:)
    type(ik_vec),allocatable :: dets_up(:),dets_dn(:)
    type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#else
    integer(ik),intent(in) :: dets_up_in(:),dets_dn_in(:)
    integer(ik),allocatable :: dets_up(:),dets_dn(:)
    integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#endif
    real(rk),intent(in) :: coeffs_in(:)
    real(rk),intent(out) :: rdm(norb,norb)

    real(rk),allocatable :: coeffs(:)
    integer,allocatable :: iorder(:),temp_i_2(:)
    integer :: idet,i_elec,p,r
    integer :: j,k
    integer,allocatable :: occ_up(:),occ_dn(:)

    allocate(dets_up(ndet))
    allocate(dets_dn(ndet))
    allocate(coeffs(ndet))
    dets_up(1:ndet) = dets_up_in(1:ndet)
    dets_dn(1:ndet) = dets_dn_in(1:ndet)
    coeffs(1:ndet) = coeffs_in(1:ndet)

    ! Sort by label for binary search
    allocate(iorder(ndet))
    allocate(temp_i16_up((ndet+1)/2))
    allocate(temp_i16_dn((ndet+1)/2))
    allocate(temp_i_2((ndet+1)/2))
    do j=1,ndet
      iorder(j)=j
    enddo

    call merge_sort2_up_dn(dets_up(1:ndet),dets_dn(1:ndet), iorder, ndet, temp_i16_up, temp_i16_dn, temp_i_2)
    coeffs(1:ndet) = coeffs(iorder(1:ndet))

    deallocate(iorder)
    deallocate(temp_i16_up)
    deallocate(temp_i16_dn)
    deallocate(temp_i_2)

    allocate(occ_up(nup))
    allocate(occ_dn(ndn))

    write (6,'(/,''Computing 1-RDM'')')

    do idet=1,ndet
      call get_occ_orbs(dets_up(idet),dets_dn(idet),occ_up,occ_dn)
      do i_elec=1,nup
        p = occ_up(i_elec)
        do r=1,norb
          if (p.ne.r.and.btest(dets_up(idet),r-1))  cycle
          call binary_search(ibset(ibclr(dets_up(idet),p-1),r-1),dets_dn(idet),dets_up,dets_dn,k)
          if (k>0) then
            rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k)
          endif
        enddo ! r
      enddo ! i_elec
      do i_elec=1,ndn
        p = occ_dn(i_elec)
        do r=1,norb
          if (p.ne.r.and.btest(dets_dn(idet),r-1))  cycle
          call binary_search(dets_up(idet),ibset(ibclr(dets_dn(idet),p-1),r-1),dets_up,dets_dn,k)
          if (k>0) then
            rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k)
          endif
        enddo ! r
      enddo ! i_elec
    enddo ! idet

    write (6,'(/,''1-RDM:'')')
    do p=1,norb
      write(6,'(1000es16.8)') rdm(p,:)
    enddo
    write (6,*)

   !write (6,*) "Writing 1-RDM to file rdm_file"
   !open(8, file='rdm_file',status='new')

   !do p=1,norb
   !  write(8,'(1000es13.5)') rdm(p,:)
   !enddo ! p

   !close(8)

   !write (6,*) "Done computing and printing 1-RDM"

  end subroutine get_1rdm


  subroutine get_1rdm_with_pt(ndet,dets_up_in,dets_dn_in,coeffs_in,diag_elems,var_energy,epsilon,n_batches,rdm)
  ! Gets the 1-RDM to lowest nonzero order in PT, i.e., if psi = psi0+psi1+..., then
  ! <psi|rho|psi> ~ <psi0|rho|psi0> + 2 <psi0|rho|psi1>
  ! Just like in the HCI energy correction, only uses terms in the numerator for which |H_{ij} c_j| > epsilon
  ! A Holmes, 27 Jun 2016

    use tools, only : merge_sort2_up_dn
    use more_tools, only : get_occ_orbs,binary_search
    use semistoch, only : find_doubly_excited,hamiltonian
   !use common_run, only : connected_dets_up,connected_dets_dn

    integer,intent(in) :: ndet
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up_in(:),dets_dn_in(:)
    type(ik_vec),allocatable :: dets_up(:),dets_dn(:)
    type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up_in(:),dets_dn_in(:)
    integer(ik),allocatable :: dets_up(:),dets_dn(:)
    integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#endif
    real(rk),intent(in) :: coeffs_in(:),diag_elems(:),var_energy,epsilon
    integer,intent(in) :: n_batches
    real(rk),intent(out) :: rdm(:,:)

    real(rk),allocatable :: coeffs(:)
    integer,allocatable :: iorder(:),temp_i_2(:)
    integer :: idet,i_elec,p,r
    integer :: j,k,k2
    integer,allocatable :: occ_up(:),occ_dn(:)
    real(rk) :: norm ! because 1st-order corrected wavefunction is not normalized
    integer :: ibatch,i
    real(rk),allocatable :: e_mix_den(:)
    real(rk),allocatable :: connected_wts(:),new_diag_elems(:)
    integer :: ndets_connected

    allocate(dets_up(ndet))
    allocate(dets_dn(ndet))
    allocate(coeffs(ndet))
    dets_up(1:ndet) = dets_up_in(1:ndet)
    dets_dn(1:ndet) = dets_dn_in(1:ndet)
    coeffs(1:ndet) = coeffs_in(1:ndet)

    ! Sort by label for binary search
    allocate(iorder(ndet))
    allocate(temp_i16_up((ndet+1)/2))
    allocate(temp_i16_dn((ndet+1)/2))
    allocate(temp_i_2((ndet+1)/2))
    do j=1,ndet
      iorder(j)=j
    enddo

    call merge_sort2_up_dn(dets_up(1:ndet),dets_dn(1:ndet), iorder, ndet, temp_i16_up, temp_i16_dn, temp_i_2)
    coeffs(1:ndet) = coeffs(iorder(1:ndet))

    deallocate(iorder)
    deallocate(temp_i16_up)
    deallocate(temp_i16_dn)
    deallocate(temp_i_2)

    allocate(occ_up(nup))
    allocate(occ_dn(ndn))

    write (6,'(/,''Computing 1-RDM to lowest nonzero order in PT'')')


    rdm(:,:) = 0._rk

    do ibatch=1,n_batches

      call find_doubly_excited(n_det=ndets_connected,dets_up=connected_dets_up,dets_dn=connected_dets_dn,ref_up=dets_up(ibatch:ndet:n_batches),ref_dn=dets_dn(ibatch:ndet:n_batches),norb=norb,n_core_orb=n_core_orb,ref_coeffs=coeffs(ibatch:ndet:n_batches),e_mix_num=connected_wts,e_mix_den=e_mix_den,min_elem=epsilon,ref_diag_elems=diag_elems(ibatch:ndet:n_batches),new_diag_elems=new_diag_elems)

      ! Now loop over single excitations from psi0, binary searching each one

      do idet=1,ndet
        call get_occ_orbs(dets_up(idet),dets_dn(idet),occ_up,occ_dn)
        do i_elec=1,nup
          p = occ_up(i_elec)
          do r=1,norb
            if (p.ne.r.and.btest(dets_up(idet),r-1))  cycle
            call binary_search(ibset(ibclr(dets_up(idet),p-1),r-1),dets_dn(idet),connected_dets_up(1:ndets_connected),connected_dets_dn(1:ndets_connected),k)
            if (k>0) then
              call binary_search(ibset(ibclr(dets_up(idet),p-1),r-1),dets_dn(idet),dets_up(1:ndet),dets_dn(1:ndet),k2)
              if (k2>0) then ! connected to variational det
                rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k2)
              else ! connected to connected det
                ! Have to do both because <psi0|rho_{pr}|psi1> != <psi1|rho_{pr}|psi0>
                rdm(p,r) = rdm(p,r) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
                rdm(r,p) = rdm(r,p) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
              endif
            endif
          enddo ! r
        enddo ! i_elec
        do i_elec=1,ndn
          p = occ_dn(i_elec)
          do r=1,norb
            if (p.ne.r.and.btest(dets_dn(idet),r-1))  cycle
            call binary_search(dets_up(idet),ibset(ibclr(dets_dn(idet),p-1),r-1),connected_dets_up(1:ndets_connected),connected_dets_dn(1:ndets_connected),k)
            if (k>0) then
              call binary_search(dets_up(idet),ibset(ibclr(dets_dn(idet),p-1),r-1),dets_up(1:ndet),dets_dn(1:ndet),k2)
              if (k2>0) then ! connected to variational det
                rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k2)
              else ! connected to connected det
                ! Have to do both because <psi0|rho_{pr}|psi1> != <psi1|rho_{pr}|psi0>
                rdm(p,r) = rdm(p,r) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
                rdm(r,p) = rdm(r,p) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
              endif
            endif
          enddo ! r
        enddo ! i_elec

      enddo ! idet

    enddo

    write (6,'(/,''PT-corrected 1-RDM:'')')
    do p=1,norb
      write(6,'(1000es13.5)') rdm(p,:)
    enddo
    write (6,*)

   !write (6,*) "Writing 1-RDM to file rdm_file"
   !open(8, file='rdm_file',status='new')

   !do p=1,norb
   !  write(8,'(1000es13.5)') rdm(p,:)
   !enddo ! p

   !close(8)

   !write (6,*) "Done computing and printing 1-RDM"

  end subroutine get_1rdm_with_pt


  integer(i8b) function estimate_n_connections(ndets,dets_up,dets_dn,coeffs,epsilon2)
  ! Estimate the number of connections that exceed eps2 for the reference wavefunction
  ! A Holmes, 24 June 2016

    use chemistry, only : find_important_connected_dets_chem
    use common_run, only : connected_dets_up,connected_dets_dn

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
    real(rk),intent(in) :: coeffs(:),epsilon2

    integer :: j
    integer :: n_connected_dets

      estimate_n_connections = 0_i8b
      do j=1,min(100,ndets)
        call find_important_connected_dets_chem(dets_up((j*ndets)/min(100,ndets)), dets_dn((j*ndets)/min(100,ndets)), epsilon2/abs(coeffs((j*ndets)/min(100,ndets))), n_connected_dets, connected_dets_up, connected_dets_dn)
        estimate_n_connections = estimate_n_connections + int(n_connected_dets,i8b)
      enddo
      estimate_n_connections = int(real(estimate_n_connections)*real(ndets)/real(min(100,ndets)),i8b)

  end function estimate_n_connections


end module hci
