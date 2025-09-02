FROM ubuntu:24.04

RUN apt update && apt upgrade -y
RUN apt install -y build-essential git libcairo2-dev libpango1.0-dev wget

WORKDIR /opt/

RUN wget http://genometools.org/pub/genometools-1.6.2.tar.gz
RUN tar genometools-1.6.2.tar.gz

WORKDIR /opt/genometools-1.6.2/
RUN sed -i.bak 's/-Werror//g' Makefile
RUN make -j8
RUN make install

RUN wget https://github.com/BrendelGroup/AEGeAn/archive/refs/tags/v0.16.0.tar.gz
RUN tar xzvf v0.16.0.tar.gz