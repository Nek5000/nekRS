FROM ornladios/adios2:ci-fedora-sanitizers-base

# Set up some arguments
ONBUILD ARG INSTALL_PREFIX
ONBUILD ARG TOOLCHAIN_FILE
ONBUILD ARG CFLAGS
ONBUILD ARG CXXFLAGS
ONBUILD ARG LLVM_USE_SANITIZER

# Build and install libc++
# Have to build in debug mode due to a bug in LLVM:
# https://bugs.llvm.org/show_bug.cgi?id=44586
ONBUILD WORKDIR /root/llvm
ONBUILD RUN git clone --branch llvmorg-9.0.0 --depth 1 \
        https://github.com/llvm/llvm-project.git source && \
    mkdir -p build && \
    cd build && \
    /opt/cmake/bin/cmake \
        -GNinja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
        -DCMAKE_C_COMPILER=/usr/bin/clang \
        -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
        -DLLVM_ENABLE_PROJECTS="libcxx;libcxxabi" \
        -DLLVM_USE_SANITIZER=${LLVM_USE_SANITIZER} \
        ../source/llvm && \
    ninja install-cxxabi install-cxx

# Copy in the toolchain
ONBUILD COPY \
    ${TOOLCHAIN_FILE} \
    ${INSTALL_PREFIX}/toolchain.cmake

# Build and install zlib
ONBUILD WORKDIR /root/zlib
ONBUILD RUN git clone --branch v1.2.11 --depth 1 \
        https://github.com/madler/zlib.git source && \
    mkdir -p build && \
    cd build && \
    /opt/cmake/bin/cmake \
        -GNinja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_TOOLCHAIN_FILE=${INSTALL_PREFIX}/toolchain.cmake \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
        ../source && \
    ninja install

# Build and install bzip2
ONBUILD WORKDIR /root/bzip2
ONBUILD RUN git clone --branch bzip2-1.0.8 \
        https://sourceware.org/git/bzip2.git source && \
    cd source && \
    sed -e "s_^CC=.*\$_CC=/usr/bin/clang_" \
        -e "s_^CFLAGS=.*\$_CFLAGS=-fpic -fPIC -Wall -Winline -O2 ${CFLAGS} \$(BIGFILES)_" \
        -i Makefile-libbz2_so && \
    make -f Makefile-libbz2_so && \
    sed -e "s_^CC=.*\$_CC=/usr/bin/clang_" \
        -e "s_^PREFIX=.*\$_PREFIX=${INSTALL_PREFIX}_" \
        -e "s_^CFLAGS=.*\$_CFLAGS=-Wall -Winline -O2 ${CFLAGS} \$(BIGFILES)_" \
        -i Makefile && \
    make install && \
    install libbz2.so.1.0.8 ${INSTALL_PREFIX}/lib && \
    ln -s -T libbz2.so.1.0.8 ${INSTALL_PREFIX}/lib/libbz2.so.1.0 && \
    ln -s -T libbz2.so.1.0 ${INSTALL_PREFIX}/lib/libbz2.so

# Build and install zeromq
ONBUILD WORKDIR /root/zeromq
ONBUILD RUN git clone --branch v4.3.2 \
        https://github.com/zeromq/libzmq.git source && \
    mkdir -p build && \
    cd build && \
    /opt/cmake/bin/cmake \
        -GNinja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_TOOLCHAIN_FILE=${INSTALL_PREFIX}/toolchain.cmake \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
        ../source && \
    ninja install

# Build and install libpng
ONBUILD WORKDIR /root/libpng
ONBUILD RUN git clone --branch v1.6.9 \
        https://git.code.sf.net/p/libpng/code.git source && \
    mkdir -p build && \
    cd build && \
    /opt/cmake/bin/cmake \
        -GNinja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_TOOLCHAIN_FILE=${INSTALL_PREFIX}/toolchain.cmake \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
        ../source && \
    ninja install

# Build and install hdf5
ONBUILD WORKDIR /root/hdf5
ONBUILD RUN git clone --branch hdf5-1_10_6 \
        https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git source && \
    mkdir -p build && \
    cd build && \
    /opt/cmake/bin/cmake \
        -GNinja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_TOOLCHAIN_FILE=${INSTALL_PREFIX}/toolchain.cmake \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
        ../source && \
    ninja install

# Build and install blosc
ONBUILD WORKDIR /root/blosc
ONBUILD RUN git clone --branch v1.17.1 \
        https://github.com/Blosc/c-blosc.git source && \
    mkdir -p build && \
    cd build && \
    /opt/cmake/bin/cmake \
        -GNinja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_TOOLCHAIN_FILE=${INSTALL_PREFIX}/toolchain.cmake \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
        ../source && \
    ninja install

# Build and install libfabric
ONBUILD WORKDIR /root/libfabric
ONBUILD RUN mkdir -p source build && \
    curl -L \
        https://github.com/ofiwg/libfabric/releases/download/v1.9.0/libfabric-1.9.0.tar.bz2 | \
    tar -C source -xj --strip-components=1 && \
    cd build && \
    ../source/configure \
        --prefix=${INSTALL_PREFIX} \
        CC="/usr/bin/clang ${CFLAGS} -L${INSTALL_PREFIX}/lib -Wl,-rpath,${INSTALL_PREFIX}/lib -Wno-unused-command-line-argument" \
        CXX="/usr/bin/clang++ ${CXXFLAGS} -L${INSTALL_PREFIX}/lib -Wl,-rpath,${INSTALL_PREFIX}/lib -Wno-unused-command-line-argument -nostdinc++ -isystem ${INSTALL_PREFIX}/include/c++/v1 -stdlib=libc++" && \
    make -j4 install

# Build and install libffi
ONBUILD WORKDIR /root/libffi
ONBUILD RUN mkdir -p source build && \
    curl -L \
        https://github.com/libffi/libffi/releases/download/v3.3/libffi-3.3.tar.gz | \
    tar -C source -xz --strip-components=1 && \
    cd build && \
    ../source/configure \
        --prefix=${INSTALL_PREFIX} \
        CC="/usr/bin/clang ${CFLAGS} -L${INSTALL_PREFIX}/lib -Wl,-rpath,${INSTALL_PREFIX}/lib -Wno-unused-command-line-argument" \
        CXX="/usr/bin/clang++ ${CXXFLAGS} -L${INSTALL_PREFIX}/lib -Wl,-rpath,${INSTALL_PREFIX}/lib -Wno-unused-command-line-argument -nostdinc++ -isystem ${INSTALL_PREFIX}/include/c++/v1 -stdlib=libc++" && \
    make -j4 install
