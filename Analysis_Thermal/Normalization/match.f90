
function match( d1, d2 )
  implicit none
  double precision, intent(in) :: d1, d2
  !character(len=1), intent(out) :: match
  double precision, parameter :: tolerance=1.0d-6
  character(len=1) :: match

  if( 0.0d0 <= d1 .and. 0.0d0 <= d2 )then
    if( d2*( 1.0d0 -tolerance ) <= d1 .and. d1 <= d2*( 1.0d0 +tolerance ) )then
      match='y'
    else
      !write(*,*) 'd1=', d1
      !write(*,*) 'd2=', d2
      !write(*,*) 'd2*( 1.0d0 -tolerance )=', d2*( 1.0d0 -tolerance ) 
      !write(*,*) 'd2*( 1.0d0 +tolerance )=', d2*( 1.0d0 +tolerance ) 
      !stop
      match='n'
    end if
  else if( d1 < 0.0d0 .and. d2 < 0.0d0 )then
    if( d1 <= d2*( 1.0d0 -tolerance ) .and. d2*( 1.0d0 +tolerance ) <= d1 )then
      match='y'
    else
      match='n'
    end if
  else if( d1*d2 < 0.0d0 )then
    match='n'
  else if( d1*d2 == 0.0d0 )then
    if( d1==0.0d0 .and. d2==0.0d0 )then
       match='y'
    else
       match='n'
    end if
  else
    write(*,*)'d1=', d1, 'd2=', d2
    write(*,*) d2*( 1.0d0 -tolerance ), '<=', d1, '<=', d2*( 1.0d0 +tolerance ), '<-- OK ?'
    !match='n'
    write(*,*)' stop at match.f90'
    stop
  end if

  !write(*,*)'d1=', d1, 'd2=', d2
  !write(*,*) d2*( 1.0d0 -tolerance ), '<=', d1, '<=', d2*( 1.0d0 +tolerance ), '<-- OK ?'

end function

