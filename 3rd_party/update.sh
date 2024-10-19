#!/bin/bash
set -a

if [ "$1" == "nek" ]; then

git rm -rf 3rd_party/nek5000 3rd_party/nek5000_parRSB
rm -rf 3rd_party/nek5000 3rd_party/nek5000_parRSB
git commit -m 'remove nek'
git subtree add --prefix 3rd_party/nek5000 https://github.com/Nek5000/nek5000.git master --squash
git subtree add --prefix 3rd_party/nek5000_parRSB https://github.com/Nek5000/parRSB.git master --squash
rm -rf 3rd_party/nek5000/tools 3rd_party/nek5000/run 3rd_party/nek5000/examples 3rd_party/nek5000/short_tests
git reset HEAD~3 --soft
git add -u
git commit -m 'import latest nek'

elif [ "$1" == "gslib" ]; then

git rm -rf 3rd_party/gslib
rm -rf 3rd_party/gslib
git commit -m 'remove gslib'
git subtree add --prefix 3rd_party/gslib https://github.com/Nek5000/gslib.git master --squash
git reset HEAD~2 --soft
git add -u
git commit -m 'import latest gslib'

elif [ "$1" == "adios" ]; then

git rm -rf 3rd_party/adios
rm -rf 3rd_party/adios
git commit -m 'remove adios'
git subtree add --prefix 3rd_party/adios https://github.com/ornladios/ADIOS2.git v2.10.1 --squash
git reset HEAD~2 --soft
git add -u
git commit -m 'import latest adios'

elif [ "$1" == "hypre" ]; then

git rm -rf 3rd_party/hypre
git commit -m 'remove hypre'
rm -rf 3rd_party/hypre
git subtree add --prefix 3rd_party/hypre https://github.com/hypre-space/hypre.git v2.32.0 --squash
rm -rf 3rd_party/hypre/src/examples 3rd_party/hypre/src/docs 3rd_party/hypre/src/test
git reset HEAD~2 --soft
git add -u
git commit -m 'import latest hypre'

elif [ "$1" == "AMGX" ]; then

git rm -rf 3rd_party/AMGX
git commit -m 'remove AMGX'
git subtree add --prefix 3rd_party/AMGX https://github.com/NVIDIA/AMGX.git main --squash
git reset HEAD~2 --soft
git add -u
git commit -m 'import latest AMGX'

elif [ "$1" == "occa" ]; then

git rm -rf 3rd_party/occa
git commit -m 'remove occa'
git subtree add --prefix 3rd_party/occa https://github.com/libocca/occa.git development --squash
git reset HEAD~2 --soft
git add -u
git commit -m 'import latest occa'

else

echo "$1 not found"

fi
