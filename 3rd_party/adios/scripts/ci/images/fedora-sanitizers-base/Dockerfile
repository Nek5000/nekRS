FROM fedora:latest

# Install core dev packages
RUN dnf -y install gcc gcc-c++ gcc-gfortran git make curl file patch \
        zlib-devel bzip2 bzip2-libs bzip2-devel blosc-devel vim valgrind \
        zeromq-devel libpng-devel libfabric-devel libffi-devel \
        python3-devel python3-numpy \
        clang llvm-devel libcxxabi-devel libcxx-devel

# Install the most recent CMake nightly binary
WORKDIR /opt/cmake
RUN curl -L https://cmake.org/files/dev/$(curl https://cmake.org/files/dev/ | sed -n '/Linux-x86_64.tar.gz/s/.*>\(cmake[^<]*\)<.*/\1/p' | sort | tail -1) | tar --strip-components=1 -xzv
ENV PATH=/opt/cmake/bin:${PATH}

# Installing Kitware Ninja
RUN mkdir -p /usr/local/bin && \
    curl -L https://github.com/Kitware/ninja/releases/download/v1.9.0.g99df1.kitware.dyndep-1.jobserver-1/ninja-1.9.0.g99df1.kitware.dyndep-1.jobserver-1_x86_64-linux-gnu.tar.gz | \
        tar -C /usr/local/bin --strip-components=1 -xzv

# Misc cleanup
RUN dnf clean all && \
    rm -rfv /tmp/* /var/cache/dnf
