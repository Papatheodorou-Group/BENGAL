# syntax=docker/dockerfile:1


FROM python:3.7.12-slim-buster


RUN apt-get update && apt-get -y install gcc && apt-get -y install g++ \
&& apt-get -y install wget && apt-get -y install autoconf && apt-get -y install automake && apt-get -y install libxml2-dev && \
wget https://github.com/Kitware/CMake/releases/download/v3.22.3/cmake-3.22.3-linux-x86_64.sh && \
cp cmake-3.22.3-linux-x86_64.sh /opt/ && chmod +x /opt/cmake-3.22.3-linux-x86_64.sh && \
cd /opt/ &&yes y | bash /opt/cmake-3.22.3-linux-x86_64.sh && ln -s /opt/cmake-3.22.3-linux-x86_64/bin/* /usr/local/bin && \
apt-get install -y git

RUN pip3 install click matplotlib pandas anndata sklearn

RUN git clone https://github.com/YY-SONG0718/sccaf && cd sccaf && pip3 install .

RUN pip3 install scanpy==1.9.1

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
