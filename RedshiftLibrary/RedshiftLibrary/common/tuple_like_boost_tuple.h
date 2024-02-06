// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef TUPLE_LIKE_BOOST_TUPLE_TUPLE_LIKE_H
#define TUPLE_LIKE_BOOST_TUPLE_TUPLE_LIKE_H

#include <boost/version.hpp>

#if BOOST_VERSION < 107400

#include <boost/tuple/tuple.hpp>
#include <cstddef>
#include <tuple>

/**
 * This helper is needed to make boost::tuple std::tuple-like. It is derived
 * from this: https://github.com/boostorg/tuple/pull/8
 *
 * This makes it possible to use boost::tuple with C++17 structured bindings
 * with recent versions of boost (like 1.65 which is the version Ubuntu 18.04
 * chips with).
 */
namespace boost_helper {
template <class T>
struct tuples_length
    : boost::integral_constant<
          int, 1 + tuples_length<typename T::tail_type>::value> {};

template <>
struct tuples_length<boost::tuple<>> : boost::integral_constant<int, 0> {};

template <>
struct tuples_length<boost::tuple<> const> : boost::integral_constant<int, 0> {
};

template <>
struct tuples_length<boost::tuples::null_type>
    : boost::integral_constant<int, 0> {};

template <>
struct tuples_length<boost::tuples::null_type const>
    : boost::integral_constant<int, 0> {};
} // namespace boost_helper

namespace std {
// std::tuple_size
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class T8, class T9, class T10>
struct tuple_size<boost::tuples::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>>
    : boost_helper::tuples_length<
          boost::tuples::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>> {};

template <class H, class T>
struct tuple_size<boost::tuples::cons<H, T>>
    : boost_helper::tuples_length<boost::tuples::cons<H, T>> {};

template <>
struct tuple_size<boost::tuples::null_type>
    : boost_helper::tuples_length<boost::tuples::null_type> {};

// std::tuple_element
template <std::size_t I, class T1, class T2, class T3, class T4, class T5,
          class T6, class T7, class T8, class T9, class T10>
struct tuple_element<
    I, boost::tuples::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>>
    : boost::tuples::element<
          I, boost::tuples::tuple<T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>> {};

template <std::size_t I, class H, class T>
struct tuple_element<I, boost::tuples::cons<H, T>>
    : boost::tuples::element<I, boost::tuples::cons<H, T>> {};
} // namespace std

#endif // boost version

#endif // TUPLE_LIKE_BOOST_TUPLE_H
