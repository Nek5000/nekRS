#!/bin/sh
# Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
# HYPRE Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

TNAME=`basename $0 .sh`
RTOL=$1
ATOL=$2

#=============================================================================
# compare with baseline case
#=============================================================================

echo  "# ${TNAME}.out.default"               > ${TNAME}.out
tail -13 ${TNAME}.out.default     | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.no_orthchk"           >> ${TNAME}.out
tail -13 ${TNAME}.out.no_orthchk  | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.pcgitr.0"             >> ${TNAME}.out
tail -13 ${TNAME}.out.pcgitr.0    | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.pcgitr.1"             >> ${TNAME}.out
tail -13 ${TNAME}.out.pcgitr.1    | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.pcgitr.2"             >> ${TNAME}.out
tail -13 ${TNAME}.out.pcgitr.2    | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.pcgtol.01"            >> ${TNAME}.out
tail -13 ${TNAME}.out.pcgtol.01   | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.pcgtol.05"            >> ${TNAME}.out
tail -13 ${TNAME}.out.pcgtol.05   | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.seed"                 >> ${TNAME}.out
tail -13 ${TNAME}.out.seed        | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.seed.repeat"          >> ${TNAME}.out
tail -13 ${TNAME}.out.seed.repeat | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.solver.none"          >> ${TNAME}.out
tail -13 ${TNAME}.out.solver.none | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.verb.1"               >> ${TNAME}.out
tail -13 ${TNAME}.out.verb.1      | head -3 >> ${TNAME}.out

echo  "# ${TNAME}.out.gen.1"                >> ${TNAME}.out
tail -14 ${TNAME}.out.gen.1       | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.gen.2"                >> ${TNAME}.out
tail -14 ${TNAME}.out.gen.2       | head -3 >> ${TNAME}.out
echo  "# ${TNAME}.out.orthchk"              >> ${TNAME}.out
tail -14 ${TNAME}.out.orthchk     | head -3 >> ${TNAME}.out

echo  "# ${TNAME}.out.itr.100"              >> ${TNAME}.out
tail -15 ${TNAME}.out.itr.100     | head -5 >> ${TNAME}.out
echo  "# ${TNAME}.out.itr.2"                >> ${TNAME}.out
tail -15 ${TNAME}.out.itr.2       | head -5 >> ${TNAME}.out
echo  "# ${TNAME}.out.vrand.2"              >> ${TNAME}.out
tail -15 ${TNAME}.out.vrand.2     | head -5 >> ${TNAME}.out

echo  "# ${TNAME}.out.verb.0"               >> ${TNAME}.out
tail -40 ${TNAME}.out.verb.0      | head -2 >> ${TNAME}.out

echo  "# ${TNAME}.out.verb.2"               >> ${TNAME}.out
tail -11 ${TNAME}.out.verb.2      | head -3 >> ${TNAME}.out

# Make sure that the output files are reasonable
CHECK_LINE="Complexity"
OUT_COUNT=`grep "$CHECK_LINE" ${TNAME}.out | wc -l`
SAVED_COUNT=`grep "$CHECK_LINE" ${TNAME}.saved | wc -l`
if [ "$OUT_COUNT" != "$SAVED_COUNT" ]; then
   echo "Incorrect number of \"$CHECK_LINE\" lines in ${TNAME}.out" >&2
fi

if [ -z $HYPRE_NO_SAVED ]; then
   #diff -U3 -bI"time" ${TNAME}.saved ${TNAME}.out >&2
   (../runcheck.sh ${TNAME}.out ${TNAME}.saved $RTOL $ATOL) >&2
fi

#=============================================================================
# remove temporary files
#=============================================================================

# rm -f ${TNAME}.testdata*
