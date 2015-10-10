program ReadPractice
implicit none
real, dimension(32,32) :: array1, array2
open(12, file="density.txt")
! read in values
read(12,*) array1
array1 = transpose(array1)
call printMatrix(array1,32,32)
close(12)
end program ReadPractice

subroutine printMatrix(array, n, m)
implicit none
real, intent(in) :: array(n,m)
integer, intent(in) :: n,m
integer :: i
do i = 1,n
print*, array(i,:)
end do
end subroutine printMatrix
