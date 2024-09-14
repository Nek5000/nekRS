/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosYAML.h basic YAML parsing functionality for ADIOS config file schema
 *
 *  Created on: Dec 19, 2019
 *      Author: Axel Huebl <axelhuebl@lbl.gov>
 */

#ifndef ADIOS2_HELPER_ADIOSJSONCOMPLEX_H_
#define ADIOS2_HELPER_ADIOSJSONCOMPLEX_H_

#include <complex>
#include <nlohmann_json.hpp>

// JSON std::complex handling
namespace std
{
template <class T>
inline void to_json(nlohmann::json &j, const std::complex<T> &p)
{
    j = nlohmann::json{p.real(), p.imag()};
}

template <class T>
inline void from_json(const nlohmann::json &j, std::complex<T> &p)
{
    p.real(j.at(0));
    p.imag(j.at(1));
}
} // end namespace std

#endif /* ADIOS2_HELPER_ADIOSJSONCOMPLEX_H_ */
