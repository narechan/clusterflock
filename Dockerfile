## Dockerfile
## Clusterflock
## version 0.1
## Clusterflock is a clustering algorithm that uses the AI of flocking organisms

# the base image
FROM ubuntu:14.04
USER root

# install program dependencies
RUN apt-get update && \
    apt-get install -y \
    g++ \
    make \
    git \
    cpanminus \
    default-jre \
    gifsicle \
    gnuplot \
    imagemagick

# install all CPAN dependencies
RUN cpanm \
    List::MoreUtils \
    List::Util \
    Getopt::Long \
    Chart::Gnuplot \
    Parallel::ForkManager \
    Statistics::Descriptive

# pull clusterflock code
RUN git clone https://github.com/narechan/clusterflock.git /home/clusterflock

# update container environment
ENV PATH /home/clusterflock:$PATH
WORKDIR /home/clusterflock
