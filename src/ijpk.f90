! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Packs two integers i and j in the lower triangular format. Assumes
! i ge. j

    integer*4 function ijpk(i,j)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    integer::i,j

    ijpk=(i*(i-1))/2+j

    return
    end function

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

