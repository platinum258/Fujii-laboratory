Module Objective_function_subprogram
        implicit none
contains
        subroutine function_inf(fn, fname, min_val, max_val)
                double precision, intent(out) :: min_val, max_val
                integer, intent(in) :: fn
                character(100), intent(out) :: fname

                if(fn == 1) then
                        min_val = -3.0d0
                        max_val = 4.0d0
                        write(fname, '(a)') 'McCorminck_function'
                else if(fn == 2) then
                        min_val = -512.0d0
                        max_val = 0.0d0
                        write(fname, '(a)') 'Eggholder_function'
                else if(fn == 3) then
                        min_val = 2.0d0
                        max_val = 4.0d0
                        write(fname, '(a)')'Easom_function'
                else if(fn == 4) then
                        min_val = -2.0d0
                        max_val = 2.0d0
                        write(fname, '(a)')'Three_hunp_camel_function'
                else if(fn == 5) then
                        min_val = -20.0d0
                        max_val = 20.0d0
                        write(fname, '(a)')'Five-well_potential_function'
                else if (fn == 6) then
                        min_val = -400.0d0
                        max_val = 100.0d0
                        write(fname,'(a)')'Schwefel_function'
                else if (fn == 7) then
                        min_val = -100.0d0
                        max_val = 100.0d0
                        write(fname,'(a)')'Sphere_function'
                else if (fn == 8) then
                        min_val = -5.12d0
                        max_val = 5.12d0
                        write(fname,'(a)')'Rastrign_function'
                end if

                write(*, *) trim(fname)


        end subroutine function_inf

        function Objective_function( xin, fn)
                double precision, intent(in) :: xin(:)
                double precision :: Objective_function
                integer, intent(in) :: fn
                integer :: i, j, n
                double precision :: pi, x, y

                pi = 2.0d0 * acos(0.0d0)
                n = size(xin)

                Objective_function = 0.0d0
                
                x = xin(1)
                y = xin(2) 

                if(fn == 1) then
                        ! ==================================
                        !        McCorminck_function
                        ! ==================================
                        Objective_function = sin(x + y) &
                                + (x - y) ** 2.0d0 &
                                - 1.5d0 * x &
                                + 2.5d0 * y &
                                + 1.0d0
                else if(fn == 2) then
                        ! ==================================
                        !        Eggholder_function
                        ! ==================================
                        Objective_function = - (y+47) &
                                * sin(sqrt(abs(y + (x/2.0d0) + 47))) &
                                - x * sin(sqrt(abs(x-(y+47))))
                else if(fn == 3) then
                        ! ==================================
                        !          Easom_function
                        ! ==================================
                        Objective_function = - (cos(x) * cos(y) &
                                * exp(-((x-pi)**2.0d0 + (y-pi)**2.0d0)))
                else if(fn == 4) then
                        ! ==================================
                        !     Three_hunp_camel_function
                        ! ==================================
                        Objective_function = (2.0d0 * (x ** 2.0d0)) &
                                - (1.05d0 * (x ** 4.0d0)) &
                                + ((x ** 6) / 6.0d0) &
                                + x * y &
                                + (y ** 2.0d0)
                else if(fn == 5) then
                        ! ==================================
                        !   Five-well_potential_function
                        ! ==================================
                        Objective_function = (1.0d0 &
                                           - (1.0d0 / (1.0d0 + 0.05d0 &
                                           * ((x ** 2.0d0) + (y - 10.0d0) ** 2.0d0))) &
                                           - (1.0d0 / (1.0d0 + 0.05d0 &
                                           * (((x - 10.0d0) ** 2.0d0) + y ** 2.0d0))) &
                                           - (1.5d0 / (1.0d0 + 0.03d0 &
                                           * (((x + 10.0d0) ** 2.0d0) + y ** 2.0d0))) &
                                           - (2.0d0 / (1.0d0 + 0.05d0 &
                                           * ((x - 5.0d0) ** 2.0d0 + (y + 10.0d0) ** 2.0d0))) &
                                           - (1.0d0 / (1.0d0 + 0.1d0 &
                                           * ((x + 5.0d0) ** 2.0d0 + (y + 10.0d0) ** 2.0d0)))) &
                                           * (1.0d0 + 0.0001 * (x ** 2.0d0 + y ** 2.0d0) ** 1.2d0)
                else if(fn == 6) then
                        ! ==================================
                        !         Schwefel_function
                        ! ==================================
                        do i = 1, n
                                Objective_function = Objective_function - (xin(i) * sin(sqrt(abs(xin(i)))))
                        end do
                else if(fn == 7) then
                        ! ==================================
                        !         Sphere_function
                        ! ==================================
                        do i = 1, n
                                Objective_function = Objective_function + (xin(i)-50.0d0) ** 2.0d0
                        end do
                else if(fn == 8) then
                        ! ==================================
                        !        Rastrigin_function
                        ! ==================================
                        Objective_function = 10.0d0 * dble(n)
                        do i = 1, n
                                Objective_function = Objective_function &
                                        + ( ((xin(i)-2) ** 2.0d0) - (10.0d0 * cos( 2.0d0 * pi * (xin(i)-2))))
                        end do
                end if

                return

        end function Objective_function
end module Objective_function_subprogram
