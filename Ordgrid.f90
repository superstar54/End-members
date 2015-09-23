!====================================================================
!  Formula to calculate the number and permutations of unique endmembers 
!  with q*m sublattices and  n constituents in each.
!====================================================================
! Xing Wang,  xingwang@csu.edu.cn, Central South University, China
! 2015-1-21, FCC
! 2015-1-25, BCC and general
!====================================================================
! Warning: calculate more than 8 sublattices and 8 constituents, 
! the size of array 2000 shoud be increased.
!====================================================================
!
!1 All the sublattices are identical (FCC structure belongs to this case)
!
!number of unique endmembers
!recurrence formula for a(m,n)
! n>=m
!a(m,n) = C(n,1)*a(m-1,1)+...+C(n,k)*a(m-k,k)+...C(n,m)*a(0,m)
! n<m
!a(m,n) = C(n,1)*a(m-1,1)+...+C(n,k)*a(m-k,k)+...C(n,m)*a(0,m)
!
!Permutations  of unique endmembers are also calculated in recurrence formula and 
!coorrespond to the calculation of number.
!e.g. for permutations p(m,n), the C(n,1)*a(m-1,1) term 
!in the first sublattice is one of the n constituents, the rest of the m-1 sublattices
!is the permutations of p(m-1,1)
!
!2 Two kinds of sublattices in symmetric position (BCC structure belongs to this case)
!
!Suppose we have 2*m sublattices (2 kinds of sublattices in symmetric position, each 
!kind have m sublattices) and n constituents in each sublattices, and the number of unique endmembers 
!for this structure is b(2*m,n). The formula of general term of b(2*m,n) can be written into a 
!recurrence formula as follows:
! b(2*m,n) = a(2,a(m,n))
!
!Suppose p(m:m,n) are the permutations of structure with 2m sublattices 
!(2 kinds of sublattices in symmetric position, each kind have m sublattices) 
!and n constituents in each, the strategy to obtain p(m:m,n) is as follows:
!p(m:m,n) are the permutations of 2 sublattices, in each sublattices is the 
!permutations of unique endmembers in m identical sublattices and n constituents in each.

!====================================================================

!.....................................................................
program odgrid
  implicit none
  integer :: i, j
  integer :: m, n, q, qm, k, nsub, nsub1, amn
  integer, allocatable :: set(:), set1(:), subset(:,:), subset1(:,:)

  character :: yn
  character*80 :: filename

  write(*, *) '        Calculate the unique endmembers'
  write(*,*)
  write(*, *) 'For example'
  write(*,*)  '-------------------------------------------------------------'
  write(*, *) '   Kinds of sublattices    number of sublattices in each kind'
  write(*, *) 'FCC        1                         4'
  write(*, *) 'BCC        2                         2'
  write(*,*)  '-------------------------------------------------------------'
  write(*, *)
  write(*,'(A,2X, $)') 'How many kinds of sublatitces: '
  read(*,*) q
  write(*,'(A,2X, $)') 'How many sublatitces in each kinds: '
  read(*,*) m
  write(*,'(A,2X, $)') 'Number of constituents in each sublattices: '
  read(*,*) n
  !m = 4
  !n = 7
  allocate(set1(n), subset1(m,2000))
  do i = 1,n,1
    set1(i) = i
  end do
  !call calcamn(amn, m,n)

  ! Calculate the unique endmembers of same kind of sublattices
  call calcperm(subset1, nsub1, set1, m,n,m)

  allocate(set(nsub1), subset(q,2000))

  do i = 1,nsub1,1
    set(i) = i
  end do
  ! Calculate the unique endmembers of different kinds of sublattices
  call calcperm(subset, nsub, set, q,nsub1,q)

  
10 format('There are total ', i4 '  unique endmembers')
20 format('Permutation ', i4 ':   ', $)

  do i = 1, nsub
    write(*, 20) i
    do j = 1, q
      do k = 1, m
        write(*, '(i3,$)') subset1(k,subset(j,i))
      end do
    end do
    write(*, *)
  end do
  write(*,*)
  write(*, 10) nsub

  write(*, *)
  write(*,'(A,2X, $)') 'Save the results?    Y/N   '
  read(*,*) yn
  if(yn.eq.'Y'.or.yn.eq.'y') then
    write(*,'(A,2X, $)') 'filename '
    read(*,*) filename
    open(11, file=filename, form='formatted', access='sequential')
      write(11, 10) nsub
      do i = 1, nsub
        write(11, 20) i
        do j = 1, q
          do k = 1, m
            write(11, '(i3,$)') subset1(k,subset(j,i))
          end do
        end do
        write(11, *)
      end do
    close(11)
  end if
  

end program odgrid


  recursive subroutine calcamn(amn, m,n)
  ! Calculate the number of unique endmembers
  !recurrence formula
  ! n>=m
  !a(m,n) = C(n,1)*a(m-1,1)+...+C(n,k)*a(m-k,k)+...C(n,m)*a(0,m)
  ! n<m
  !a(m,n) = C(n,1)*a(m-1,1)+...+C(n,k)*a(m-k,k)+...C(n,m)*a(0,m)
    integer :: i, m, n
    integer :: amn, amn1, cmn

    amn = 0
    ! the first term in the 2-dimensional matirx
    if(m.eq.0) then  ! the first row
      amn = 1
      return
    elseif(m.eq.1) then ! the second row
      amn = n
      return
    elseif(n.eq.1) then ! the first column
      amn = 1
      return
    end if
    !
    if(n.lt.m) then
      do i = 1, n, 1
        call calcamn(amn1, m-i, i)
        call calccmn(cmn, n,i)
        amn = amn + cmn*amn1  !recurrence formula
      end do
    else
      do i = 1, m, 1
        call calcamn(amn1, m-i, i)
        call calccmn(cmn, n,i)
        amn = amn + cmn*amn1
      end do
    endif
  end subroutine calcamn

  recursive subroutine calcperm(subset, nsub, set, m,n,m0)
  ! Calculate the number of unique endmembers and the permutation
  ! p(m,n)
    integer :: i, m, n
    integer :: set(n), subset(m0,2000), subtemp1(m0,2000), subtemp2(m0,2000)
    integer :: amn, cmn
    integer :: nsub, nsub1, nsub2

    nsub = 0
    ! the first term in the 2-dimensional matirx
    if(m.eq.0) then
      nsub = 0
      return
    elseif(m.eq.1) then
      do i = 1, n, 1
        nsub = nsub + 1
        subset(1,nsub) = set(i)
      end do
      return
    elseif(n.eq.1) then
      nsub = nsub + 1
      do i = 1, m, 1
        subset(i,nsub) = set(n)
      end do
      return
    end if
    !
    if(n.le.m) then
      do i = 1, n, 1

        call combine(subtemp1, nsub1, set, n, i, m0)
        call calccmn(cmn, n,i)
        call calcamn(amn, m-i, i)
        do j = 1, cmn, 1

          call calcperm(subtemp2, nsub2, subtemp1(1:i,j), m-i,i,m0)
          ! nsub2=0 means m=0, the a(0,n)=1, however the number of permutation is zero
          if(nsub2.eq.0) then  
            do k = 1, amn, 1
              nsub = nsub + 1
              subset(1:i, nsub) = subtemp1(1:i,j)
            end do
          else
            do k = 1, amn, 1
              nsub = nsub + 1
              ! the 1~i sublattices
              subset(1:i, nsub) = subtemp1(1:i,j) 
              ! the rest of the m-i sublattices is the same with p(m-i,k)
              subset((i+1):m, nsub) = subtemp2(1:(m-i), k)
            end do
          end if
        end do
      end do
    else
      do i = 1, m, 1
        call combine(subtemp1, nsub1, set, n, i, m0)
        call calccmn(cmn, n,i)
        call calcamn(amn, m-i, i)
        do j = 1, cmn, 1
          call calcperm(subtemp2, nsub2, subtemp1(1:i,j), m-i,i,m0)
          if(nsub2.eq.0) then
            do k = 1, amn, 1
              nsub = nsub + 1
              subset(1:i, nsub) = subtemp1(1:i,j)
            end do
          else
            do k = 1, amn, 1
              nsub = nsub + 1
              ! the 1~i sublattices
              subset(1:i, nsub) = subtemp1(1:i,j)
              ! the rest of the m-i sublattices is the same with p(m-i,k)
              subset((i+1):m, nsub) = subtemp2(1:(m-i), k)
            end do
          end if
        end do
      end do
    endif
  end subroutine calcperm

  subroutine combine(subset, nsub, set, n,k, m)
  ! select k element from total n element
  ! calculate all the combinations
    integer :: i, j, count
    integer :: n, k
    integer :: set(n), vec(n), subset(m, 2000)
	logical :: flag

    ! build the 0-1 vector
    nsub = 0
    do i = 1,n,1
      if(i.le.k) then
        vec(i) = 1
      else
        vec(i) = 0
      end if
    end do
    !begin scan
    flag = .TRUE.
    do while(flag)
      ! get choosen
      nsub = nsub + 1
      j = 0
      do i = 1,n,1
        if(vec(i).eq.1) then
          j = j + 1
          subset(j, nsub) = set(i)
        end if
      end do
      !call sprint(subset(nsub,:),k)
      flag = .FALSE.
      do i = 1,(n-1),1
        if(vec(i).eq.1.and.vec(i+1).eq.0) then
          vec(i) = 0
          vec(i+1) =1
          ! move all 1 to left most side.
          count = 0
          do j = 1,i,1
            if(vec(j).eq.1) then
              count = count + 1
            end if
          end do
          if(count.le.i) then
            do j = 1, count, 1
              vec(j) = 1
            end do
            do j = (count + 1), i, 1
              vec(j) = 0
            end do
          end if
          flag = .TRUE.
          exit
        end if
      end do
    end do
  end subroutine combine

  subroutine calccmn(nsub, n,k)
  ! select k element from total n element
  ! calculate the number of combinations
    integer :: i, j, count
    integer :: cmn, n, k
    integer :: vec(n)
	logical :: flag

    ! build the 0-1 vector
    nsub = 0
    do i = 1,n,1
      if(i.le.k) then
        vec(i) = 1
      else
        vec(i) = 0
      end if
    end do
    !begin scan
    flag = .TRUE.
    do while(flag)
      nsub = nsub + 1
      flag = .FALSE.
      do i = 1,(n-1),1
        if(vec(i).eq.1.and.vec(i+1).eq.0) then
          vec(i) = 0
          vec(i+1) =1
          ! move all 1 to left most side.
          count = 0
          do j = 1,i,1
            if(vec(j).eq.1) then
              count = count + 1
            end if
          end do
          if(count.le.i) then
            do j = 1, count, 1
              vec(j) = 1
            end do
            do j = (count + 1), i, 1
              vec(j) = 0
            end do
          end if
          flag = .TRUE.
          exit
        end if
      end do
    end do
  end subroutine calccmn
