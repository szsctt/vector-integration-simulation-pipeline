# Dockerfile for viral integration simulation pipeline

FROM mambaorg/micromamba:2.3.0-ubuntu24.04

USER root

#  $ docker build . -t szsctt/simvi:latest -t szsctt/simvi:1
#  $ docker run --rm -it szsctt/simvi:latest /bin/bash
#  $ docker push szsctt/simvi:latest
#  $ docker push szsctt/simvi:1


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

# install conda stuff
ADD scripts/consolidate_envs.py /opt/simvi/scripts/
ADD envs /opt/simvi/envs/
RUN micromamba install -n base -c conda-forge pip pyyaml=6 python=3.11 -y &&\
	/opt/conda/bin/python3 /opt/simvi/scripts/consolidate_envs.py /opt/simvi/envs/*yml /opt/simvi/envs/simvi.yml &&\
	micromamba env update -n base -f /opt/simvi/envs/simvi.yml -y &&\
	micromamba clean --all -y 	

# include simvi scripts, etc
ADD scripts /opt/simvi/scripts/
ADD snakemake_rules /opt/simvi/snakemake_rules
ADD Snakefile /opt/simvi/Snakefile
ADD run_singularity.sh /opt/simvi/run_singularity.sh

# add test files
ADD tests /opt/simvi/tests
RUN mkdir -p /opt/simvi/tests/out .snakemake

WORKDIR /opt/simvi

CMD tests/snakemake/runMe.sh
