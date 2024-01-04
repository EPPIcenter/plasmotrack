FROM centos:8
LABEL authors="Maxwell Murphy"

ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
EXPOSE 22

RUN cd /etc/yum.repos.d/
RUN sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-*
RUN sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-*
# Install dependencies
RUN dnf update -y && \
    dnf --enablerepo=powertools install -y \
    gcc-toolset-11-gcc \
    gcc-toolset-11-gcc-c++ \
    git \
    make \
    ninja-build \
    gdb \
    python3 \
    python3-pip \
    python3-devel \
    wget \
    openssh \
    openssh-server \
    openssh-clients \
    passwd \
    rsync \
    && dnf clean all

RUN curl -L "https://github.com/Kitware/CMake/releases/download/v3.27.9/cmake-3.27.9-linux-x86_64.sh" -o cmake.sh && \
    chmod +x cmake.sh && \
    ./cmake.sh --skip-license --prefix=/usr/local && \
    rm cmake.sh

RUN echo "source /opt/rh/gcc-toolset-11/enable" >> /etc/bashrc
SHELL ["/bin/bash", "--login", "-c"]

RUN pip3 install --upgrade pip && pip3 install conan

ADD start.sh /start.sh

RUN chmod +x /start.sh /tini

ENTRYPOINT ["/tini", "-g", "--"]

CMD /start.sh