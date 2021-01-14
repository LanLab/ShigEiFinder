FROM ncbi/blast:2.10.0

RUN apt-get update -y && apt-get install -y \
	wget \
	build-essential\
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev\
    gcc-multilib\
    bwa\
    python3.7

# #Install samtools
RUN wget -P / "https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2"
RUN tar -xvf samtools-1.10.tar.bz2
RUN cd samtools-1.10 && ./configure --prefix=/usr/local/
RUN cd samtools-1.10 && make 
RUN cd samtools-1.10 && make install

# # Install BWA 
# RUN apt-get install --yes git
# WORKDIR /tmp
# RUN git clone https://github.com/lh3/bwa.git
# WORKDIR /tmp/bwa
# RUN git checkout v0.7.15
# RUN make
# RUN cp -p bwa /usr/local/bin && rm -rf /tmp/bwa
# RUN apt-get clean

# Set working dir
WORKDIR /app

