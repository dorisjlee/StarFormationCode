program ReadPractice
implicit none
real :: sum,mm1,xx1,xx2,xx3,vv1,vv2,vv3,ll1,ll2,ll3
integer::io,x
sum = 0
open(10,file="ic_sink",form='formatted')
DO
   read(10,*,iostat=io)mm1,xx1,xx2,xx3,vv1,vv2,vv3,ll1,ll2,ll3
   write(*,*)"after read"
!   READ(*,*,IOSTAT=io)  x
   IF (io > 0) THEN
      WRITE(*,*) 'Check input.  Something was wrong'
      EXIT
   ELSE IF (io < 0) THEN
      WRITE(*,*)  'The total is ', sum
      EXIT
   ELSE
      sum = sum + x
   END IF
END DO
      write(*,*)"after do"
       write(*,*)"IOstat:",io
     write(*,*)"Info: ", mm1,xx1,xx2,xx3,vv1,vv2,vv3,ll1,ll2,ll3
close(12)
end program ReadPractice


!program ReadPractice
!implicit none
!real, dimension(32,32) :: array1, array2
!open(12, file="density.txt")
!! read in values
!read(12,*) array1
!array1 = transpose(array1)
!call printMatrix(array1,32,32)
!close(12)
!end program ReadPractice

!subroutine printMatrix(array, n, m)
!implicit none
!real, intent(in) :: array(n,m)
!integer, intent(in) :: n,m
!integer :: i
!do i = 1,n
!print*, array(i,:)
!end do
!end subroutine printMatrix
