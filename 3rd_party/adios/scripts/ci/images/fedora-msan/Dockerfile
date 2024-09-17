ARG INSTALL_PREFIX=/opt/msan
ARG TOOLCHAIN_FILE=toolchain.cmake
ARG CFLAGS="-fsanitize=memory -fsanitize-memory-track-origins"
ARG CXXFLAGS="-fsanitize=memory -fsanitize-memory-track-origins"
ARG LLVM_USE_SANITIZER=MemoryWithOrigins

FROM ornladios/adios2:ci-fedora-full-stack-onbuild AS tmp-stage

FROM ornladios/adios2:ci-fedora-sanitizers-base
COPY --from=tmp-stage /opt/msan/ /opt/msan/
