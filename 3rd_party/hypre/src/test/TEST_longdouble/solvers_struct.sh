#!/bin/sh
# Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
# HYPRE Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

TNAME=`basename $0 .sh`
CONVTOL=$1

# Set default check tolerance
if [ x$CONVTOL = "x" ];
then
    CONVTOL=0.0
fi
#echo "tol = $CONVTOL"
#=============================================================================
# compare with baseline case
#=============================================================================

FILES="\
 ${TNAME}.out.0\
 ${TNAME}.out.1\
 ${TNAME}.out.2\
 ${TNAME}.out.3\
 ${TNAME}.out.4\
"

for i in $FILES
do
  echo "# Output file: $i"
  tail -3 $i
done > ${TNAME}.out

FILES="\
 ${TNAME}.out.10.lobpcg\
 ${TNAME}.out.11.lobpcg\
 ${TNAME}.out.17.lobpcg\
 ${TNAME}.out.18.lobpcg\
 ${TNAME}.out.19.lobpcg\
"

for i in $FILES
do
  echo "# Output file: $i"
  tail -3 $i
  echo "# Output file: $i.1"
  tail -13 $i.1 | head -3
  echo "# Output file: $i.5"
  tail -21 $i.5 | head -11
done >> ${TNAME}.out

# Make sure that the output files are reasonable
CHECK_LINE="Iterations"
OUT_COUNT=`grep "$CHECK_LINE" ${TNAME}.out | wc -l`
SAVED_COUNT=`grep "$CHECK_LINE" ${TNAME}.saved | wc -l`
if [ "$OUT_COUNT" != "$SAVED_COUNT" ]; then
   echo "Incorrect number of \"$CHECK_LINE\" lines in ${TNAME}.out" >&2
fi

if [ -z $HYPRE_NO_SAVED ]; then
   (../runcheck.sh ${TNAME}.out ${TNAME}.saved $CONVTOL) >&2
fi

#=============================================================================
# remove temporary files
#=============================================================================

# rm -f ${TNAME}.testdata*
