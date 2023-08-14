FROM gipert/bat:1.0.0

LABEL maintainer.name="Luigi Pertoldi"
LABEL maintainer.email="gipert@pm.me"

ARG HMIXFIT_VERSION="1.0.3"

USER root
WORKDIR /root

# install hmixfih
RUN mkdir -p src /opt/hmixfit && \
    wget -q -O- "https://github.com/gipert/hmixfit/archive/refs/tags/v${HMIXFIT_VERSION}.tar.gz" \
        | tar --strip-components 1 -xz -C src && \
    cd src && \
    make PREFIX=/opt/hmixfit install && \
    cd .. && rm -r src

ENV PATH="/opt/hmixfit/bin:$PATH"

ENTRYPOINT ["/opt/hmixfit/bin/hmixfit"]
CMD ["--help"]
