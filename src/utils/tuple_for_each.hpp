#ifndef TUPLE_TOOLS_HPP_
#define TUPLE_TOOLS_HPP_

#include <tuple>
#include <utility>

namespace impl {
template <int N> struct tag_type {
  static constexpr int value = N;
};
} // namespace impl

// adapted from
// https://stackoverflow.com/questions/1198260/how-can-you-iterate-over-the-elements-of-an-stdtuple
template <std::size_t I = 0, typename FuncT, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
tuple_for_each(std::tuple<Tp...> &, FuncT) // Unused arguments are given no names.
{
}

template <std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if <
    I<sizeof...(Tp), void>::type tuple_for_each(std::tuple<Tp...> &t, FuncT f)
{
  const auto T = impl::tag_type<I>{};
  f(T);
  tuple_for_each<I + 1, FuncT, Tp...>(t, f);
}

namespace impl {
// adapted from https://stackoverflow.com/questions/67477403/making-an-index-sequence-tuple
template <typename T, typename I> struct n_tuple_helper {
};

template <typename T, std::size_t... I> struct n_tuple_helper<T, std::index_sequence<I...>> {
  using type = std::tuple<std::enable_if_t<I || true, T>...>;
};
} // namespace impl

template <typename T, std::size_t N>
using n_tuple = typename impl::n_tuple_helper<T, std::make_index_sequence<N>>::type;

#endif