ARG baseos
FROM ornladios/adios2:spack-dependencies-${baseos}

# Install from CI source
ARG ci_source_dir=ci-source
COPY ${ci_source_dir} /opt/adios2/source
RUN spack dev-build \
        -j$(grep -c '^processor' /proc/cpuinfo) \
        -d /opt/adios2/source \
        --skip-patch \
        adios2@develop && \
    spack clean -a

# Build the environment
RUN spack env create --without-view adios2 && \
    spack -e adios2 add $(spack find --format "/{hash}") && \
    spack -e adios2 install -v && \
    rm -rf /root/.spack && \
    spack env activate adios2 && \
    spack env deactivate && \
    spack -e adios2 env loads

# Initialize the first-use-of-environment stuff for the adios user
WORKDIR /home/adios
USER adios
RUN spack env activate adios2 && \
    spack env deactivate && \
    echo "source /opt/spack/share/spack/setup-env.sh" >> ~/.bash_profile && \
    echo "module use ${SPACK_ROOT}/share/spack/modules/linux-ubuntu18.04-haswell" >> ~/.bash_profile && \
    echo "source ${SPACK_ROOT}/var/spack/environments/adios2/loads" >> ~/.bash_profile

ENTRYPOINT []
CMD ["bash", "--login"]
