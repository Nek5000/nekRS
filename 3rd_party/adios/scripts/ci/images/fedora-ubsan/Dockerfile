FROM ornladios/adios2:ci-fedora-sanitizers-base

# Install core dev packages
RUN dnf -y install libubsan python3-devel python3-numpy hdf5-devel

# Misc cleanup
RUN dnf clean all && \
    rm -rfv /tmp/* /var/cache/dnf
