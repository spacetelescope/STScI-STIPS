FROM continuumio/miniconda3
MAINTAINER Space Telescope Science Institute <help@stsci.edu>

#########
# Setup #
#########

ENV HOME /root
WORKDIR $HOME

# Clone STScI-STIPS
RUN git clone https://github.com/spacetelescope/STScI-STIPS.git
ENV STIPSDIR $HOME/STScI-STIPS

##########################
# Basic apt Requirements #
##########################

WORKDIR $STIPSDIR

RUN conda env create -f environment.yml
RUN echo "conda activate stips" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

##################
# Reference Data #
##################

WORKDIR $HOME

RUN echo -e "import stips\nstips.DownloadReferenceData()" >> activate_stips.py
RUN python activate_stips.py
