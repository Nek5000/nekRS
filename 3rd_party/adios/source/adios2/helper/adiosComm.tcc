/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosComm.tcc : specialization of template functions defined in
 * adiosComm.h
 */

#ifndef ADIOS2_HELPER_ADIOSCOMM_TCC_
#define ADIOS2_HELPER_ADIOSCOMM_TCC_

#include "adiosComm.h"

namespace adios2
{
namespace helper
{

// BroadcastValue full specializations forward-declared in 'adiosComm.inl'.
template <>
std::string Comm::BroadcastValue(const std::string &input, const int rankSource) const
{
    const size_t inputSize = input.size();
    const size_t length = this->BroadcastValue(inputSize, rankSource);
    std::string output;

    if (rankSource == this->Rank())
    {
        output = input;
    }
    else
    {
        output.resize(length);
    }

    this->Bcast(const_cast<char *>(output.data()), length, rankSource);

    return output;
}

// Comm::Impl::GetDatatype full specializations forward-declared in
// 'adiosComm.inl'.
template <>
CommImpl::Datatype CommImpl::GetDatatype<signed char>()
{
    return CommImpl::Datatype::SignedChar;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<char>()
{
    return CommImpl::Datatype::Char;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<short>()
{
    return CommImpl::Datatype::Short;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<int>()
{
    return CommImpl::Datatype::Int;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<long>()
{
    return CommImpl::Datatype::Long;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned char>()
{
    return CommImpl::Datatype::UnsignedChar;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned short>()
{
    return CommImpl::Datatype::UnsignedShort;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned int>()
{
    return CommImpl::Datatype::UnsignedInt;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned long>()
{
    return CommImpl::Datatype::UnsignedLong;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<unsigned long long>()
{
    return CommImpl::Datatype::UnsignedLongLong;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<long long>()
{
    return CommImpl::Datatype::LongLong;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<double>()
{
    return CommImpl::Datatype::Double;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<long double>()
{
    return CommImpl::Datatype::LongDouble;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<int, int>>()
{
    return CommImpl::Datatype::Int_Int;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<float, int>>()
{
    return CommImpl::Datatype::Float_Int;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<double, int>>()
{
    return CommImpl::Datatype::Double_Int;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<long double, int>>()
{
    return CommImpl::Datatype::LongDouble_Int;
}

template <>
CommImpl::Datatype CommImpl::GetDatatype<std::pair<short, int>>()
{
    return CommImpl::Datatype::Short_Int;
}

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSCOMM_TCC_ */
