# Use the Ubuntu 20.04 as the base image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# Install core dependencies
RUN apt-get update && \
    apt-get install -y \
        wget \
        curl \
        unzip \
        openjdk-8-jdk \
        python3-pip \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        bzip2 && \
    rm -rf /var/lib/apt/lists/*

# Set Java 1.8 as def.
RUN update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java

# Install GATK 4.2.6.1
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip && \
    unzip gatk-4.2.6.1.zip -d /opt && \
    rm gatk-4.2.6.1.zip && \
    ln -s /opt/gatk-4.2.6.1/gatk /usr/local/bin/gatk

# Install BWA 0.7.17
RUN wget -q https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download -O bwa-0.7.17.tar.bz2 && \
    tar -xjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    mv bwa /usr/local/bin && \
    cd .. && \
    rm -rf bwa-0.7.17 bwa-0.7.17.tar.bz2

# Install FastQC 0.11.9
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip -d /opt && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.11.9.zip

# Update the package list
RUN apt-get update && apt-get install -y \
    build-essential \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    wget \
    nano \
    bzip2

# Install Samtools 1.15.1
RUN wget -q https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
    tar -xjf samtools-1.15.1.tar.bz2 && \
    cd samtools-1.15.1 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.15.1 samtools-1.15.1.tar.bz2

# Install MultiQC 1.13
RUN pip3 install multiqc==1.13

# Create a symlink for python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Download and install SRA Toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    mv sratoolkit.3.1.1-ubuntu64/bin/* /usr/local/bin/ && \ 
    rm -rf sratoolkit.3.1.1-ubuntu64 sratoolkit.current-ubuntu64.tar.gz

# Optional: Set the PATH for the SRA Toolkit
ENV PATH="/usr/local/bin:${PATH}"

# Update the package lists
RUN apt-get update && \
    apt-get install -y \
    openjdk-11-jre \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Download Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -O /opt/Trimmomatic.zip && \
    unzip /opt/Trimmomatic.zip -d /opt/ && \
    rm /opt/Trimmomatic.zip

# Set environment variables for Trimmomatic
ENV TRIMMOMATIC_PATH="/opt/Trimmomatic-0.39"
ENV PATH="${PATH}:${TRIMMOMATIC_PATH}"

# Create wrapper for Trimmomatic
RUN echo '#!/bin/bash\njava -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /usr/local/bin/trimmomatic && \
    chmod +x /usr/local/bin/trimmomatic

# Install R and dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('ggplot2', 'dplyr'), repos='http://cran.rstudio.com/')"

# Set up PATH
ENV PATH="/usr/local/bin:$PATH"

WORKDIR "/app/pipeline_folder"

# Set default command to open a shell
CMD ["/bin/bash"]