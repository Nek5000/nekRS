/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * heatTransfer.cpp
 *
 * Recreates heat_transfer.f90 (Fortran) in C++
 *
 * Created on: Feb 2017
 *     Author: Norbert Podhorszki
 *
 */

#include <cstring>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdexcept>
#include <string>

#include "HeatTransfer.h"

HeatTransfer::HeatTransfer(const Settings &settings)
: m_s(settings), m_T1Buf(new double[(m_s.ndx + 2) * (m_s.ndy + 2)]),
  m_T2Buf(new double[(m_s.ndx + 2) * (m_s.ndy + 2)]), m_T1(new double *[m_s.ndx + 2]),
  m_T2(new double *[m_s.ndx + 2])
{
    m_T1[0] = m_T1Buf.get();
    m_T2[0] = m_T2Buf.get();
    for (size_t i = 1; i < m_s.ndx + 2; ++i)
    {
        m_T1[i] = m_T1[0] + i * (m_s.ndy + 2);
        m_T2[i] = m_T2[0] + i * (m_s.ndy + 2);
    }
    m_TCurrent = m_T1.get();
    m_TNext = m_T2.get();
}

void HeatTransfer::init(bool init_with_rank)
{
    if (init_with_rank)
    {
        std::fill_n(m_T1Buf.get(), (m_s.ndx + 2) * (m_s.ndy + 2), static_cast<double>(m_s.rank));
    }
    else
    {
        const double hx = 2.0 * 4.0 * atan(1.0) / m_s.ndx;
        const double hy = 2.0 * 4.0 * atan(1.0) / m_s.ndy;

        double x, y;
        for (unsigned int i = 0; i < m_s.ndx + 2; i++)
        {
            x = 0.0 + hx * (i - 1);
            for (unsigned int j = 0; j < m_s.ndy + 2; j++)
            {
                y = 0.0 + hy * (j - 1);
                m_T1[i][j] = cos(8 * x) + cos(6 * x) - cos(4 * x) + cos(2 * x) - cos(x) +
                             sin(8 * y) - sin(6 * y) + sin(4 * y) - sin(2 * y) + sin(y);
            }
        }
    }
    m_TCurrent = m_T1.get();
    m_TNext = m_T2.get();
}

void HeatTransfer::printT(std::string message, MPI_Comm comm) const
{
    int rank, size;
    int tag = 1;
    int token;
    MPI_Status status;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (rank > 0)
    {
        MPI_Recv(&token, 1, MPI_INT, rank - 1, tag, comm, &status);
    }

    std::cout << "Rank " << rank << " " << message << std::endl;
    for (unsigned int i = 0; i < m_s.ndx + 2; i++)
    {
        std::cout << "  T[" << i << "][] = ";
        for (unsigned int j = 0; j < m_s.ndy + 2; j++)
        {
            std::cout << std::setw(6) << m_TCurrent[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::flush << std::endl;

    if (rank < size - 1)
    {
        MPI_Send(&token, 1, MPI_INT, rank + 1, tag, comm);
    }
}

void HeatTransfer::iterate()
{
    for (unsigned int i = 1; i <= m_s.ndx; ++i)
    {
        for (unsigned int j = 1; j <= m_s.ndy; ++j)
        {
            m_TNext[i][j] = omega / 4 *
                                (m_TCurrent[i - 1][j] + m_TCurrent[i + 1][j] +
                                 m_TCurrent[i][j - 1] + m_TCurrent[i][j + 1]) +
                            (1.0 - omega) * m_TCurrent[i][j];
        }
    }
    std::swap(m_TCurrent, m_TNext);
}

void HeatTransfer::heatEdges()
{
    // Heat the whole global edges
    if (m_s.posx == 0)
        std::fill_n(m_TCurrent[0], m_s.ndy + 2, edgetemp);

    if (m_s.posx == m_s.npx - 1)
        std::fill_n(m_TCurrent[m_s.ndx + 1], m_s.ndy + 2, edgetemp);

    if (m_s.posy == 0)
        for (unsigned int i = 0; i < m_s.ndx + 2; ++i)
            m_TCurrent[i][0] = edgetemp;

    if (m_s.posy == m_s.npy - 1)
        for (unsigned int i = 0; i < m_s.ndx + 2; ++i)
            m_TCurrent[i][m_s.ndy + 1] = edgetemp;
}

void HeatTransfer::exchange(MPI_Comm comm)
{
    // Build a custom MPI type for the column vector to allow strided access
    MPI_Datatype tColumnVector;
    MPI_Type_vector(m_s.ndx + 2, 1, m_s.ndy + 2, MPI_DOUBLE, &tColumnVector);
    MPI_Type_commit(&tColumnVector);

    // Exchange ghost cells, in the order left-right-up-down

    // send to left + receive from right
    int tag = 1;
    MPI_Status status;
    if (m_s.rank_left >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send left to rank "
        //          << m_s.rank_left << std::endl;
        MPI_Send(m_TCurrent[0] + 1, 1, tColumnVector, m_s.rank_left, tag, comm);
    }
    if (m_s.rank_right >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from right from rank "
        //          << m_s.rank_right << std::endl;
        MPI_Recv(m_TCurrent[0] + (m_s.ndy + 1), 1, tColumnVector, m_s.rank_right, tag, comm,
                 &status);
    }

    // send to right + receive from left
    tag = 2;
    if (m_s.rank_right >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send right to rank "
        //          << m_s.rank_right << std::endl;
        MPI_Send(m_TCurrent[0] + m_s.ndy, 1, tColumnVector, m_s.rank_right, tag, comm);
    }
    if (m_s.rank_left >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from left from rank "
        //          << m_s.rank_left << std::endl;
        MPI_Recv(m_TCurrent[0], 1, tColumnVector, m_s.rank_left, tag, comm, &status);
    }

    // Cleanup the custom column vector type
    MPI_Type_free(&tColumnVector);

    // send down + receive from above
    tag = 3;
    if (m_s.rank_down >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send down to rank "
        //          << m_s.rank_down << std::endl;
        MPI_Send(m_TCurrent[m_s.ndx], m_s.ndy + 2, MPI_DOUBLE, m_s.rank_down, tag, comm);
    }
    if (m_s.rank_up >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from above from rank "
        //          << m_s.rank_up << std::endl;
        MPI_Recv(m_TCurrent[0], m_s.ndy + 2, MPI_DOUBLE, m_s.rank_up, tag, comm, &status);
    }

    // send up + receive from below
    tag = 4;
    if (m_s.rank_up >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send up to rank " <<
        // m_s.rank_up
        //          << std::endl;
        MPI_Send(m_TCurrent[1], m_s.ndy + 2, MPI_DOUBLE, m_s.rank_up, tag, comm);
    }
    if (m_s.rank_down >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from below from rank "
        //          << m_s.rank_down << std::endl;
        MPI_Recv(m_TCurrent[m_s.ndx + 1], m_s.ndy + 2, MPI_DOUBLE, m_s.rank_down, tag, comm,
                 &status);
    }
}

/* Copies the internal ndx*ndy section of the ndx+2 * ndy+2 local array
 * into a separate contiguous vector and returns it.
 * @return A vector with ndx*ndy elements
 */
std::vector<double> HeatTransfer::data_noghost() const
{
    std::vector<double> d(m_s.ndx * m_s.ndy);
    for (unsigned int i = 1; i <= m_s.ndx; ++i)
    {
        std::memcpy(&d[(i - 1) * m_s.ndy], m_TCurrent[i] + 1, m_s.ndy * sizeof(double));
    }
    return d;
}
