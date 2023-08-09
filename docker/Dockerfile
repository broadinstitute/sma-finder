FROM python:3.9-slim-bullseye

RUN apt-get purge ca-certificates-java

# install core tools and libraries
RUN apt update && apt-get install -y --no-install-recommends \
    autoconf \
    build-essential \
    bzip2 \
    ca-certificates \
    curl \
    dpkg-dev \
    g++ \
    gcc \
    git \
    make \
    openjdk-17-jre \
    pkg-config \
	python3-pip \
    libbz2-dev \
    libc6-dev \
    libcurl4-openssl-dev \
    libffi-dev \
    libgdbm-dev \
    liblzma-dev \
    libncursesw5-dev \
    libsqlite3-dev \
    libssl-dev \
    openssl \
    unzip \
    wget \
    zlib1g-dev

# install htslib
ENV SAMTOOLS_VERSION="1.15.1"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# install gcloud
ENV GCLOUD_SDK_VERISON="441.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

# install GATK in order to use PrintReads to localize a subset of the CRAM file
ENV GATK_VERSION="4.4.0.0"
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
	&& unzip gatk-${GATK_VERSION}.zip \
	&& rm gatk-${GATK_VERSION}.zip \
	&& ln -s /gatk-${GATK_VERSION}/gatk /usr/local/bin/gatk

# install sma_finder.py
COPY requirements.txt /
RUN python3 -m pip install --upgrade -r requirements.txt

COPY sma_finder.py /

ENTRYPOINT /bin/bash

