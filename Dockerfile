FROM python:3.6

RUN apt-get update \
    && apt-get install -y \
    libopenblas-dev=0.3.13+ds-3 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy libraries
COPY --from=quay.io/biocontainers/htslib:1.14--h9093b5e_0 /usr/local/lib/libhts.so.3 /usr/local/lib/libtinfow.so.6 /usr/local/lib/
COPY --from=quay.io/biocontainers/bcftools:1.14--h88f3f91_0 /usr/local/lib/libgsl.so.25 /usr/local/lib/libcblas.so.3 /usr/local/lib/

# Copy binaries
COPY --from=quay.io/biocontainers/htslib:1.14--h9093b5e_0 /usr/local/bin/bgzip /usr/local/bin/htsfile /usr/local/bin/tabix /usr/local/bin/
COPY --from=quay.io/biocontainers/spoa:4.0.7--h9a82719_1 /usr/local/bin/spoa /usr/local/bin/
COPY --from=quay.io/biocontainers/racon:1.4.20--h9a82719_1 /usr/local/bin/racon /usr/local/bin/
COPY --from=quay.io/biocontainers/minimap2:2.23--h5bf99c6_0 /usr/local/bin/minimap2 /usr/local/bin/
COPY --from=quay.io/biocontainers/samtools:1.14--hb421002_0 /usr/local/bin/samtools /usr/local/bin/
COPY --from=quay.io/biocontainers/bcftools:1.14--h88f3f91_0 /usr/local/bin/bcftools /usr/local/bin/

RUN pip install medaka==1.5.0 NGSpeciesID

ENTRYPOINT [ "NGSpeciesID" ]
