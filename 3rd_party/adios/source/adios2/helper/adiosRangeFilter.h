/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosRangeFilter.h a simple range filter
 *
 *  Created on: Feb 4, 2021
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#ifndef ADIOS2_HELPER_RANGEFILTER_H_
#define ADIOS2_HELPER_RANGEFILTER_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <string>
#include <vector>
/// \endcond

namespace adios2
{
namespace helper
{

constexpr size_t MaxAllowedNumber = 9999999ULL;

class RangeFilter
{
public:
    RangeFilter(){};
    ~RangeFilter() = default;

    /** Create a rangefilter from a list of range definitions.
     * @param selection string of space-separated list of range definitions in
     * the form of "start:end:step". Indexing starts from 0. If 'end' is 'n' or
     * 'N', then it is an unlimited range expression. Range definitions are
     * adding up.
     * Examples:
     * "0 6 3 2" selected 4 steps indexed 0,2,3 and 6
     * "1:5" selects 5 consecutive steps, skipping step 0, and starting from 1
     * "2:n" selects all steps from step 2
     * "0:n:2" selects every other steps from the beginning (0,2,4,6...)
     * "0:n:3  10:n:5" selects every third step from the beginning and
     * additionally every fifth steps from step 10.
     * RangeFilters that are not initialized will return true for everything
     * @return Will throw an invalid_argument exception if the selection
     * does not conform to the specification
     */
    void ParseSelection(const std::string &selection);

    /**
     * Checks if an element is present in the selection
     * @param n non-negative integer
     * @return true if 'element' was in the selection
     */
    bool IsSelected(size_t n);

private:
    // true: element was in selection, false: was not in selection
    std::vector<bool> m_Selection;
    // unlimited selection rules (pair of Start and Steps)
    std::vector<std::pair<size_t, size_t>> m_UnlimitedRules;
    size_t ToSizeT(const std::string &input);
}; // class

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_RANGEFILTER_H_ */
