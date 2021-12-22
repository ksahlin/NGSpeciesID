FROM quay.io/biocontainers/htslib:1.14--h9093b5e_0 AS htslib
FROM quay.io/biocontainers/spoa:4.0.7--h9a82719_1 AS spoa
FROM quay.io/biocontainers/racon:1.4.20--h9a82719_1 AS racon
FROM quay.io/biocontainers/minimap2:2.23--h5bf99c6_0 AS minimap2
FROM quay.io/biocontainers/samtools:1.14--hb421002_0 AS samtools
FROM quay.io/biocontainers/bcftools:1.14--h88f3f91_0 AS bcftools

FROM python:3.6

RUN apt-get update \
    && apt-get install -y \
    libopenblas-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip install medaka NGSpeciesID

# Copy libraries
COPY --from=htslib /usr/local/lib/libhts.so.3 /usr/local/lib/libtinfow.so.6 /usr/local/lib/
COPY --from=bcftools /usr/local/lib/libgsl.so.25 /usr/local/lib/libcblas.so.3 /usr/local/lib/

# Copy binaries
COPY --from=htslib /usr/local/bin/bgzip /usr/local/bin/htsfile /usr/local/bin/tabix /usr/local/bin/
COPY --from=spoa /usr/local/bin/spoa /usr/local/bin/
COPY --from=racon /usr/local/bin/racon /usr/local/bin/
COPY --from=minimap2 /usr/local/bin/minimap2 /usr/local/bin/
COPY --from=samtools /usr/local/bin/samtools /usr/local/bin/
COPY --from=bcftools /usr/local/bin/bcftools /usr/local/bin/

ENTRYPOINT [ "NGSpeciesID" ]