###  CORI
ADIOS_DIR=/global/homes/p/pnorbert/adios/master/login
#ADIOS_LIBS=-Wl,-rpath=${ADIOS_DIR}/lib64 -L${ADIOS_DIR}/lib64 -ladios2 
#ADIOS_INC=-I${ADIOS_DIR}/include
MEXLIBS="LDFLAGS=${ADIOS_LIBS}"

###   TITAN
#ADIOS_DIR=/ccs/proj/e2e/pnorbert/ADIOS/ADIOS2/build.titan.gnu-nompi/install
#MEXLIBS=${ADIOS_LIBS}

###  Kathleen's PC
#ADIOS_DIR=/opt/adios/2.2-lean
#ADIOS_LIBS=-Wl,-rpath=${ADIOS_DIR}/lib64 -L${ADIOS_DIR}/lib -ladios2
#ADIOS_INC=-I${ADIOS_DIR}/include
#MEXLIBS="LDFLAGS=${ADIOS_LIBS}"

### General - this should work everywhere
ADIOS_INC=-I${ADIOS_DIR}/include
ADIOS_LIBS=`${ADIOS_DIR}/bin/adios2-config --c-libs`

### MacOS - example using homebrew installed ADIOS2 and Xcode 15 clang
###      1) Install homebrew (https://brew.sh/) and Xcode (App Store)
###      2) brew install adios2
###      OR
###      2) Compile Adios2 from scratch and update ADIOS_DIR below to match install directory
#ADIOS_DIR=/opt/homebrew/opt/adios2
#ADIOS_INC=-I${ADIOS_DIR}/include
#ADIOS_LIBS=-Wl,-rpath,${ADIOS_DIR}/lib -shared -L${ADIOS_DIR}/lib -ladios2_c -ladios2_core
#MEXLIBS="LDFLAGS=${ADIOS_LIBS}"


MEXOPTS=-largeArrayDims -DDEBUG CFLAGS="-g -std=c99 -fPIC -O0"
default:
	@echo "ADIOS reader methods for MATLAB"
	@echo "Have 'mex' compiler in the path and run"
	@echo "> make matlab"


matlab: adiosopenc adiosclosec adiosreadc

INC:=$(shell echo ${ADIOSREAD_SEQ_INC})
LIB:=$(shell echo ${ADIOSREAD_SEQ_CLIB})

adiosopenc:  
	mex adiosopenc.c  ${MEXOPTS} -I. ${ADIOS_INC} ${MEXLIBS}

adiosreadc: 
	mex adiosreadc.c  ${MEXOPTS} -I. ${ADIOS_INC} ${MEXLIBS}

adiosclosec:
	mex adiosclosec.c ${MEXOPTS} -I. ${ADIOS_INC} ${MEXLIBS}

#adiosdbgc:
#       mex adiosdbgc.c ${MEXOPTS} -I. ${ADIOS_INC} ${ADIOS_LIBS}


clean:
	rm -f *.mexa64 core* *.o 

install:
	@echo "Install *.m and *.mexa64 into MATLABPATH"
	@echo "Do not have 'adios' anywhere in the path where you install the files"
	@echo "otherwise 'help adios' will not work in Matlab"
