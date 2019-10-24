program mymain

    USE procedures
    INTEGER pos, no_lines, iostat_value, nn, t_count, skip
    REAL b
    REAL, dimension(:), allocatable :: values
    CHARACTER(30) :: buffer
    no_lines = 0

! **************************************************************************************
!     Added new variable skip. Stores value when skip=1,2; and Skips value when skip=3,4
! **************************************************************************************
    skip = 0


    open (25, file= './Patient results A Interval/patient68.csv')

    ! Count the number of lines in the file, by incrementing 'no_lines'
    ! everytime a line is read into 'buffer'
    print *, "-------------------------------"
    open(10, file='./Patient SIgnals/a0068.txt')
    do
       read(10, "(A)", IOSTAT=iostat_value) buffer
       if (iostat_value > 0)  then
          print *, "Something wrong"
      else if (iostat_value < 0) then
          print *, "Reached end of file"
          EXIT
      else
          no_lines = no_lines + 1
      end if
    end do
   close(10)
   no_lines = no_lines - 2
   ! **************************************************************************************
   !     Dividing no_lines by 2 since we are discarding half the data
   ! **************************************************************************************
   no_lines = no_lines /  2
   call lengthOfData(no_lines, nn)

   ! Allocate the size of the array = 'no_lines' minus the 2 header lines
   allocate(values(no_lines))

   ! Open file again, and skip the first 2 lines
   ! Set 'no_lines' = 0 and use as an index for inserting into the 'values' array
   t_count = 1
   open(15, file='./Patient Signals/a0068.txt')
   ! **************************************************************************************
   !     Skipping the first 3 now instead of 2, since we want to start from .001
   ! **************************************************************************************
   do i=1,3
       read(15,*)
   end do
   do
      read(15, "(A)", IOSTAT=iostat_value) buffer
      if (iostat_value > 0)  then
         print *, "Something wrong"
     else if (iostat_value < 0) then
         print *, "Reached end of file"
         EXIT
     else
         ! **************************************************************************************
         !   Check skip. IF < 2: store in values then increment skip; ELSE: increment no_skip
         !      when no_skip >=2 then reset skip and no_skip
         ! **************************************************************************************
         if (skip < 2) then
             ! get the position of the comma in 'buffer'
             pos = index(buffer, ",")
             ! read integer into 'b' from the index after the comma
             read(buffer(pos+1:), *) b
             ! save value into the 'values' array in index 't_count'
             values(t_count) = b
             t_count = t_count + 1
             skip = skip + 1
             !write (20,*) b
         else
             skip = skip + 1;
             if (skip >= 4) then
                 skip = 0
             end if
         end if

     end if

     if (t_count > no_lines) then
         EXIT
      end if

    end do
  close(15)

  print *, no_lines, " values"
  call calculate_rescaled_range(values,no_lines,nn)

end program mymain

module procedures
    implicit none

contains

    subroutine lengthOfData (no_lines, nn)
      INTEGER, INTENT(INOUT) :: no_lines
      REAL tempo
      INTEGER, INTENT(OUT) :: nn
      tempo = REAL(no_lines)
      nn = floor(log10(tempo)/log10(2.0))
      no_lines = 2**nn
      nn = nn
    end subroutine lengthOfData

    ! @values: the array containing the readings
    ! @data_size: the size of the values array
    subroutine calculate_rescaled_range(values,data_size,n)
        INTEGER, INTENT(IN) :: data_size
        REAL, dimension(:), INTENT(IN) :: values
        INTEGER, INTENT(IN) :: n
        REAL, dimension(1:n) :: RSave, logRSave
        REAL, dimension(n,data_size) :: mean_values, rescaled_ranges
        REAL temp_sum, mean, temp_variable, rescaled_range
        INTEGER i, j, k, no_of_regions, loop_count
        INTEGER size_of_set ! The size of the data in each region
        INTEGER lower_bound, upper_bound
        REAL,dimension(1:n) :: sizeOfSet,logsizeOfSet !Array for the size of set

        ! Iterating from 1 to n.
        ! Iteration 1 does 1 region
        ! 2 breaks it into 2 regions, 3 does 4 regions, 4 does 8 regions...
        ! 'no_of_regions' = 2 ^ (Iteration -1)
        do i=1,n
            print *, "---  REGION SIZE: ", i, "   ---"
            no_of_regions = 2**(i-1)
            temp_variable = data_size / no_of_regions
            size_of_set = FLOOR(temp_variable)
            sizeOfSet(i) = size_of_set !To store each each value for size of set in the Array
            lower_bound = 1
            upper_bound = size_of_set
            loop_count = 1

            print *, "No of regions: ", no_of_regions
            print *, "Size of set: ", size_of_set
            print *, "Lower bound: ", lower_bound, ": ", values(lower_bound)
            print *, "Upper bound: ", upper_bound, ": ", values(upper_bound)


            ! Loop for the regions
            do while (loop_count <= no_of_regions) !TODO: Changed from <= i


                ! Set 'actual_size_of_set' because of uneven divisions
                !actual_size_of_set = upper_bound - lower_bound + 1

                ! Iterate from 'lower_bound' to 'upper_bound' to find the mean
                ! of the values in the region
                temp_sum = 0
                do j=lower_bound, upper_bound
                    temp_sum = temp_sum + values(j)
                enddo
                mean = temp_sum / size_of_set

                ! Store the mean multidimensional array 'mean_values'
                mean_values(i,loop_count) = mean
                ! Find the standard deviation of the region

                call calculate_deviation(values,mean,lower_bound,upper_bound,rescaled_range)

                rescaled_ranges(i,loop_count) = rescaled_range
                lower_bound = lower_bound + size_of_set
                upper_bound = upper_bound + size_of_set
                loop_count = loop_count + 1
            enddo


        end do

        ! Calculate the averages for the values in rescaled_ranges
        do j=1,n
            temp_sum = 0
            do k=1,loop_count+1
              if (k <= 2**(j-1)) then
                temp_sum = temp_sum + rescaled_ranges(j,k)
              else
                EXIT
              endif
            enddo
            k = k -1
            RSave(j) = temp_sum / k
        enddo

        write(*,501)
        write (25,505)
        write (*,*) "----------------------------------------"
        do i = 1,n
          write (*,500) sizeOfSet(i),  RSave(i)
          !write (25,500) sizeOfSet(i),  RSave(i)
          write (25,'(1F7.1,A,1F9.4)') sizeOfSet(i),',',RSave(i)
        end do
        write (25,*)

        call calculate_hurst_exponent(sizeOfSet, RSave,logsizeOfSet, logRSave, n)

        500 Format (5X, 1F7.1, 15x, 1F9.4)
        505 Format ("Size of Set",",","Average Rescaled Range")
        501 Format (2X,"Size of Set", 3X, "Average Rescaled Range")

    end subroutine calculate_rescaled_range

    ! @values: the array containing the readings
    ! @mean: the mean of the region to calculate the deviation for
    ! @lower_bound: the starting index of the region in 'values'
    ! @upper_bound: the last index in the region in 'values'
    ! @rescaled_range: value to be calculated and returned by the subroutine
    subroutine calculate_deviation(values, mean, lower_bound, upper_bound, rescaled_range)
        INTEGER, INTENT(IN) :: lower_bound, upper_bound
        REAL, INTENT(IN) :: mean
        REAL, dimension(:), INTENT(IN) :: values
        REAL, INTENT(OUT) :: rescaled_range
        INTEGER i, counter, size
        REAL running_total, difference, max_value, min_value, range, variance
        REAL stdDev, runnig_sum_sq, temp_square
        REAL, dimension(1:upper_bound-lower_bound+1) :: deviation

        running_total = 0
        max_value = 0
        min_value = 0
        runnig_sum_sq = 0
        temp_square = 0
        counter = lower_bound
        size = upper_bound-lower_bound+1

        ! Iterate through the region and find the deviation
        do i=1,size
            difference = values(counter) - mean
            running_total = running_total + difference
            deviation(i) = difference
            if ( running_total > max_value ) then
                max_value = running_total
            end if
            if ( running_total < min_value ) then
                min_value = running_total
            end if
            counter = counter + 1
        enddo

        range = max_value - min_value
        do i = 1,size
            temp_square = deviation(i) ** 2
            runnig_sum_sq = runnig_sum_sq + temp_square
        end do
        variance = runnig_sum_sq / (size)
        stdDev = sqrt (variance)
        if (stdDev == 0) then
          rescaled_range = 1
        else
        rescaled_range = range / stdDev
      end if
    end subroutine calculate_deviation

    subroutine calculate_hurst_exponent (x,y,logx,logy,n)
      INTEGER, INTENT(IN) :: n
      Real, dimension(n) :: x,y
      Real, dimension (n) :: logx, logy
      integer :: i,ipt
      real :: sum_logx=0, sum_logy=0, sqsumlogx=0, sum_logx_logy=0
      real :: logx_ave, logy_ave, m, c
      do 26 i = 1,n
      logx = log10(x)
      logy = log10(y)
    !  write (*,*)
    !  write(*,500) logx(i), logy(i)
    !  write(*,*)
      !500 format (1F7.4, 2x, 1F6.4)

        sum_logx  = sum_logx  + logx(i)
        sum_logy  = sum_logy  + logy(i)
        sqsumlogx = sqsumlogx + logx(i)*logx(i)
        sum_logx_logy = sum_logx_logy + logx(i)*logy(i)
     26 continue
     ipt = i-1
      logx_ave = sum_logx / float(ipt)
      logy_ave = sum_logy / float(ipt)

      m = (sum_logx_logy - sum_logx*logy_ave) / (sqsumlogx - sum_logx*logx_ave)
      c = logy_ave - m*logx_ave

    !  write(*,*)' intercept =      ', c
      write(*,*)' hurst exponent = ', m
      write(25,*)
      write(25,*)' Hurst exponent = ', m
   end subroutine calculate_hurst_exponent

end module
