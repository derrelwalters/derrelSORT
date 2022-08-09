! Derrel Walters © 2019
! These sort routines were written by Derrel Walters ~ 2019-01-23


SUBROUTINE iSORT(arrA)
  ! This implementation of derrelSORT is for integers,
  ! but the same principles apply for other datatypes.
  ! This routine executes linearly with respect to the
  ! number of elements in the passed array.
  !
  ! ~ Derrel Walters
  IMPLICIT NONE

!  INTEGER(KIND=8),INTENT(IN) :: nA
!  INTEGER(KIND=4),INTENT(IN) :: nA
  INTEGER,DIMENSION(:),INTENT(INOUT) :: arrA

  INTEGER(KIND=4) :: nA 
  INTEGER,DIMENSION(:),ALLOCATABLE :: arrB 
!  INTEGER(KIND=8) :: lowIDX, highIDX, midIDX
  INTEGER(KIND=4) :: lowIDX, highIDX, midIDX
  INTEGER ::  iStat
!  INTEGER(KIND=8) :: i, j, A, B, C, thisHigh, mergeSize, nLoops
  INTEGER(KIND=4) :: i, j, A, B, C, thisHigh, mergeSize, nLoops
  INTEGER,DIMENSION(:),ALLOCATABLE :: iterMark
  LOGICAL,DIMENSION(:),ALLOCATABLE :: moreToGo
  
  nA = SIZE(arrA)
  ALLOCATE(arrB(nA), STAT=iStat)
  arrB = arrA
  mergeSize = 2
  lowIDX = 1 - mergeSize
  highIDX = 0

  nLoops = INT(LOG(REAL(nA))/LOG(2.0))
  ALLOCATE(iterMark(nLoops), moreToGo(nLoops), STAT=iStat)
  moreToGo = .FALSE.
  iterMark = 0

  DO i = 1, nLoops
    iterMark(i) = FLOOR(REAL(nA)/2**i) 
    IF (MOD(nA, 2**i) > 0) THEN
      moreToGo(i) = .TRUE.
      iterMark(i) = iterMark(i) + 1
    END IF
  END DO

  DO i = 1, nLoops
      DO j = 1, iterMark(i)
        A = 0
        B = 1
        C = 0
        lowIDX = lowIDX + mergeSize
        highIDX = highIDX + mergeSize
        midIDX = (lowIDX + highIDX + 1) / 2
        thisHigh = highIDX
        IF (j == iterMark(i).AND.moreToGo(i)) THEN 
          lowIDX = lowIDX - mergeSize
          highIDX = highIDX - mergeSize
          midIDX = (lowIDX + highIDX + 1) / 2
          A = midIDX - lowIDX
          B = 2
          C = nA - 2*highIDX + midIDX - 1
          thisHigh = nA
        END IF
!! The traditional merge can also be used (see subroutine for comment). !!
!                                                                        !
!        CALL imerge(arrA(lowIDX:midIDX-1+A), B*(midIDX-lowIDX),   &     !
!                    arrA(midIDX+A:thisHigh), highIDX-midIDX+1+C, &      !
!                    arrB(lowIDX:thisHigh), thisHigh-lowIDX+1)           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL imerge2(arrA(lowIDX:midIDX-1+A), B*(midIDX-lowIDX),   &
                    arrA(midIDX+A:thisHigh), highIDX-midIDX+1+C,   &
                    arrB(lowIDX:thisHigh), thisHigh-lowIDX+1)
        arrA(lowIDX:thisHigh) = arrB(lowIDX:thisHigh)
      END DO
      mergeSize = 2*mergeSize
      lowIDX = 1 - mergeSize
      highIDX = 0
  END DO

END SUBROUTINE iSORT

SUBROUTINE imerge(arrA, nA, arrB, nB, arrC, nC)
  ! This merge is a traditional merge that places
  ! the lowest element first. The form that the
  ! time complexity takes, O(n), is not affected
  ! by the merge routine - yet this routine
  ! does not run as fast as the merge used in
  ! imerge2.  
  !
  ! ~Derrel Walters
  IMPLICIT NONE

!  INTEGER(KIND=8),INTENT(IN) :: nA, nB , nC
  INTEGER(KIND=4),INTENT(IN) :: nA, nB , nC

  INTEGER,DIMENSION(nA),INTENT(IN) :: arrA
  INTEGER,DIMENSION(nB),INTENT(IN) :: arrB
  INTEGER,DIMENSION(nC),INTENT(INOUT) :: arrC

!  INTEGER(KIND=8) :: i, j, k
  INTEGER(KIND=4) :: i, j, k

  arrC = 0
  i = 1
  j = 1
  k = 1

  DO
    IF (i > nA .OR. j > NB) EXIT
    IF (arrB(j) < arrA(i)) THEN
      arrC(k) = arrB(j)
      j = j + 1
    ELSE
      arrC(k) = arrA(i)
      i = i + 1
    END IF
    k = k + 1
  END DO

  IF (i <= nA) THEN
    DO
      IF (i > nA) EXIT
        arrC(k) = arrA(i)
        i = i + 1
        k = k + 1
    END DO
  ELSEIF (j <= nB) THEN
    DO
      IF (j > nB) EXIT
        arrC(k) = arrB(j)
        j = j + 1
        k = k + 1
    END DO
  END IF

END SUBROUTINE imerge

SUBROUTINE imerge2(arrA, nA, arrB, nB, arrC, nC)
  ! This merge is a faster merge.  Array A arrives
  ! just to the left of Array B, and Array C is
  ! filled from both ends simultaneously - while
  ! still preserving the stability of the sort.
  ! The derrelSORT routine is so fast, that
  ! the merge does not affect the O(n) time 
  ! complexity of the sort in practice 
  ! (perhaps, making its execution more linear
  ! at small numbers of elements). 
  !
  ! ~ Derrel Walters
  IMPLICIT NONE

!  INTEGER(KIND=8),INTENT(IN) :: nA, nB , nC
  INTEGER(KIND=4),INTENT(IN) :: nA, nB , nC

  INTEGER,DIMENSION(nA),INTENT(IN) :: arrA
  INTEGER,DIMENSION(nB),INTENT(IN) :: arrB
  INTEGER,DIMENSION(nC),INTENT(INOUT) :: arrC

!  INTEGER(KIND=8) :: i, j, k, x, y, z
  INTEGER(KIND=4) :: i, j, k, x, y, z

  arrC = 0
  i = 1
  j = 1
  k = 1
  x = nA
  y = nB
  z = nC

  DO
    IF (i > x .OR. j > y) EXIT
    IF (arrB(j) < arrA(i)) THEN
      arrC(k) = arrB(j)
      j = j + 1
    ELSE
      arrC(k) = arrA(i)
      i = i + 1
    END IF
    IF (arrA(x) > arrB(y)) THEN
      arrC(z) = arrA(x)
      x = x - 1
    ELSE
      arrC(z) = arrB(y)
      y = y - 1
    END IF
    k = k + 1
    z = z - 1
  END DO

  IF (i <= x) THEN
    DO
      IF (i > x) EXIT
        arrC(k) = arrA(i)
        i = i + 1
        k = k + 1
    END DO
  ELSEIF (j <= y) THEN
    DO
      IF (j > y) EXIT
        arrC(k) = arrB(j)
        j = j + 1
        k = k + 1
    END DO
  END IF
END SUBROUTINE imerge2
