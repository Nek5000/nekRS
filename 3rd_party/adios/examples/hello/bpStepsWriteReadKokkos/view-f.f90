module view_f_mod use, intrinsic ::iso_c_binding use,
    intrinsic ::iso_fortran_env

            use ::flcl_mod

                implicit none

        public

        interface subroutine f_init_view(x, val) &
        &bind(c, name = 'c_init_view') use,
    intrinsic ::iso_c_binding use ::flcl_mod type(c_ptr), intent(in)::x integer(c_int),
    intent(in)::val end subroutine f_init_view

        subroutine f_print_view(x) &
        &bind(c, name = 'c_print_view') use,
    intrinsic ::iso_c_binding use ::flcl_mod type(c_ptr),
    intent(in)::x end subroutine f_print_view end interface

    contains

    subroutine init_view(x, val) use,
    intrinsic ::iso_c_binding use ::flcl_mod implicit none type(view_i32_2d_t),
    intent(inout)::x integer(c_int),
    intent(in)::val

    call f_init_view(x % ptr(), val)

        end subroutine init_view

    subroutine print_view(x) use,
    intrinsic ::iso_c_binding use ::flcl_mod implicit none type(view_i32_2d_t),
    intent(in)::x

    call f_print_view(x % ptr())

        end subroutine print_view

    end module view_f_mod
