/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ExampleReadPlugin.tcc
 *
 *  Created on: Jul 5, 2021
 *      Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#ifndef EXAMPLEREADPLUGIN_TCC_
#define EXAMPLEREADPLUGIN_TCC_

#include "ExampleReadPlugin.h"

namespace adios2
{
namespace plugin
{

template <typename T>
void ExampleReadPlugin::AddVariable(const std::string &name, Dims shape, Dims start, Dims count)
{
    core::Variable<T> *v = m_IO.InquireVariable<T>(name);
    if (!v)
    {
        m_IO.DefineVariable<T>(name, shape, start, count);
    }
}

template <class T>
inline void ExampleReadPlugin::ReadVariable(core::Variable<T> &variable, T *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        auto delimPos = line.find(",");
        auto name = line.substr(0, delimPos);
        size_t step = std::stoul(line.substr(delimPos + 1));
        if (name == variable.m_Name && m_CurrentStep == step)
        {
            std::string vals;
            std::getline(m_DataFile, vals);
            if (vals.find(",") == vals.npos)
            {
                values[0] = helper::StringTo<T>(vals, "");
            }
            else
            {
                std::istringstream ss(vals);
                std::string val;
                int i = 0;
                while (std::getline(ss, val, ','))
                {
                    values[i] = helper::StringTo<T>(val, "");
                    i++;
                }
            }
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<std::string> &variable,
                                            std::string *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        if (line == variable.m_Name)
        {
            std::getline(m_DataFile, values[0]);
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<char> &variable, char *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        if (line == variable.m_Name)
        {
            for (size_t i = 0; i < variable.SelectionSize(); ++i)
            {
                m_DataFile.get(&values[i], 1, ',');
            }
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<unsigned char> &variable,
                                            unsigned char *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        if (line == variable.m_Name)
        {
            for (size_t i = 0; i < variable.SelectionSize(); ++i)
            {
                char val;
                m_DataFile.get(&val, 1, ',');
                values[i] = static_cast<unsigned char>(val);
            }
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<signed char> &variable,
                                            signed char *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        if (line == variable.m_Name)
        {
            for (size_t i = 0; i < variable.SelectionSize(); ++i)
            {
                char val;
                m_DataFile.get(&val, 1, ',');
                values[i] = static_cast<signed char>(val);
            }
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<short> &variable, short *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        if (line == variable.m_Name)
        {
            for (size_t i = 0; i < variable.SelectionSize(); ++i)
            {
                std::string val;
                std::getline(m_DataFile, val, ',');
                values[i] = static_cast<short>(std::stoi(val));
            }
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<unsigned short> &variable,
                                            unsigned short *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        if (line == variable.m_Name)
        {
            for (size_t i = 0; i < variable.SelectionSize(); ++i)
            {
                std::string val;
                std::getline(m_DataFile, val, ',');
                values[i] = static_cast<unsigned short>(std::stoi(val));
            }
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<long double> &variable,
                                            long double *values)
{
    while (m_DataFile.good())
    {
        std::string line;
        std::getline(m_DataFile, line);
        if (line == variable.m_Name)
        {
            for (size_t i = 0; i < variable.SelectionSize(); ++i)
            {
                std::string val;
                std::getline(m_DataFile, val, ',');
                try
                {
                    values[i] = std::strtold(val.c_str(), nullptr);
                }
                catch (...)
                {
                    std::throw_with_nested(
                        std::invalid_argument("ERROR: could not cast " + val + " to long double "));
                }
            }
            break;
        }
    }
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<std::complex<float>> &variable,
                                            std::complex<float> *values)
{
    throw std::invalid_argument("ERROR: std::complex<float> not supported in this engine");
}

template <>
inline void ExampleReadPlugin::ReadVariable(core::Variable<std::complex<double>> &variable,
                                            std::complex<double> *values)
{
    throw std::invalid_argument("ERROR: std::complex<double> not supported in this engine");
}
} // end namespace plugin
} // end namespace adios2
#endif /* EXAMPLEREADPLUGIN_TCC_ */
