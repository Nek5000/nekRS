#include "flcl-cxx.hpp"
#include <Kokkos_Core.hpp>

using view_type = flcl::view_i32_2d_t;

extern "C" {

void c_init_view(view_type **v_x, int *val)
{
    using flcl::view_from_ndarray;

    view_type x = **v_x;
    int d_val = *val;
    Kokkos::parallel_for("initBuffer", x.extent(1), KOKKOS_LAMBDA(const size_t idy) {
        for (size_t idx = 0; idx < x.extent(0); idx++)
            x(idx, idy) = d_val + idx;
    });
    Kokkos::fence();

    return;
}

void c_print_view(view_type **v_x)
{
    using flcl::view_from_ndarray;

    view_type x = **v_x;

    Kokkos::parallel_for("printBuffer", 1, KOKKOS_LAMBDA(const size_t id) {
        for (size_t idx = 0; idx < x.extent(0); idx++)
        {
            for (size_t idy = 0; idy < x.extent(1); idy++)
                Kokkos::printf("%d ", x(idx, idy));
            Kokkos::printf("\n");
        }
    });
    Kokkos::fence();

    return;
}
}
