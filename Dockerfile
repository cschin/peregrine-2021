#ref: https://alexbrand.dev/post/how-to-package-rust-applications-into-minimal-docker-containers/
FROM rust:1.58.1 AS build
WORKDIR /usr/src

# since we use the mimalloc crate, it does not work well with musl.
#RUN rustup target add x86_64-unknown-linux-musl

# Create a dummy project and build the app's dependencies.
# If the Cargo.toml or Cargo.lock files have not changed,
# we can use the docker build cache and skip these (typically slow) steps.
RUN USER=root cargo new peregrine-r 
WORKDIR /usr/src/peregrine-r
#COPY Cargo.toml Cargo.lock build.rs ./
COPY Cargo.toml build.rs ./
COPY .git ./.git
RUN apt-get update && apt-get install -y build-essential cmake 
RUN cargo build --release

# Copy the source and build the application.
COPY src ./src
#RUN cargo install --target x86_64-unknown-linux-musl --path .
RUN cargo install --path .

# Copy the statically-linked binary into a scratch container.
#FROM scratch
#COPY --from=build /usr/local/cargo/bin/pg_* ./
#USER 1000
CMD ["./pg_asm"]

FROM ubuntu:20.04
COPY --from=build /usr/local/cargo/bin/pg_* /usr/local/bin/
RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y parallel pigz awscli
RUN apt-get install -y samtools exonerate jq wget
COPY run_script_from_s3.sh /usr/local/bin
RUN chmod u+x /usr/local/bin/run_script_from_s3.sh
CMD ["pg_asm"]
