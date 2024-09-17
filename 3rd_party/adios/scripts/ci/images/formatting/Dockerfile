FROM fedora:38

RUN dnf update -y && \
    dnf install -y \
    clang-tools-extra \
    curl \
    git \
    python3-flake8 \
    python3-pylint \
    python3-black \
    ShellCheck \
    && \
    dnf clean all

# Check that the target programs has been installed
RUN command -v clang-format > /dev/null && \
    command -v flake8 > /dev/null && \
    command -v shellcheck > /dev/null
