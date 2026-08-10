FROM ubuntu:24.04@sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea

LABEL org.opencontainers.image.title="genome-population-simulator F1 distribution test"
LABEL org.opencontainers.image.version="0.1.0"
LABEL org.opencontainers.image.base.name="ubuntu:24.04@sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea"

COPY bin/genome-population-simulator /usr/local/bin/genome-population-simulator
COPY schemas /opt/genome-population-simulator/schemas
COPY licenses /opt/genome-population-simulator/licenses
COPY examples /opt/genome-population-simulator/examples

ENTRYPOINT ["/usr/local/bin/genome-population-simulator"]
CMD ["--version"]
