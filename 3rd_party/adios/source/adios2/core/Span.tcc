/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Span.tcc
 *
 *  Created on: Apr 17, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#ifndef ADIOS2_CORE_SPAN_TCC_
#define ADIOS2_CORE_SPAN_TCC_

#include "Span.h"

#include "adios2/core/Engine.h"

namespace adios2
{
namespace core
{

template <class T>
Span<T>::Span(Engine &engine, const size_t size) : m_Engine(engine), m_Size(size)
{
}

template <class T>
size_t Span<T>::Size() const noexcept
{
    return m_Size;
}

template <class T>
T *Span<T>::Data() const noexcept
{
    return m_Engine.BufferData<T>(m_BufferIdx, m_PayloadPosition);
}

template <class T>
T &Span<T>::At(const size_t position)
{
    if (position > m_Size)
    {
        helper::Throw<std::invalid_argument>("Core", "Span", "At",
                                             "position " + std::to_string(position) +
                                                 " is out of bounds for span of size " +
                                                 std::to_string(m_Size));
    }

    return (*this)[position];
}

template <class T>
T &Span<T>::operator[](const size_t position)
{
    T &data = *m_Engine.BufferData<T>(m_BufferIdx, m_PayloadPosition + position * sizeof(T));
    return data;
}

} // end namespace core
} // end namespace adios2

#endif // ADIOS2_CORE_SPAN_TCC_
