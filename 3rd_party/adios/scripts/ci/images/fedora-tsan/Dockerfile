ARG INSTALL_PREFIX=/opt/tsan
ARG TOOLCHAIN_FILE=toolchain.cmake
ARG CFLAGS="-fsanitize=thread"
ARG CXXFLAGS="-fsanitize=thread"
ARG LLVM_USE_SANITIZER=Thread

FROM ornladios/adios2:ci-fedora-full-stack-onbuild AS tmp-stage

FROM ornladios/adios2:ci-fedora-sanitizers-base
COPY --from=tmp-stage /opt/tsan/ /opt/tsan/
