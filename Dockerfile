FROM continuumio/miniconda3
MAINTAINER Space Telescope Science Institute <help@stsci.edu>

# RUN "sh" "-c" "echo nameserver 8.8.8.8 >> /etc/resolv.conf"
# RUN "sh" "-c" "echo forcing APT update"


#########
# Setup #
#########

ENV HOME /root
WORKDIR $HOME

# Place to store data and source
RUN mkdir -p /opt

# Clone STScI-STIPS
RUN git clone https://github.com/spacetelescope/STScI-STIPS.git
ENV STIPSDIR $HOME/STScI-STIPS

##########################
# Basic apt Requirements #
##########################

WORKDIR $HOME

RUN apt-get update
RUN apt-get install -y python3.7
RUN apt-get install -y python-pip
RUN apt-get install -y wget

# matplotlib requirement
RUN apt-get install -y pkg-config
RUN apt-get install -y libfreetype6-dev
RUN apt-get install -y rsync
RUN apt-get install -y git

# SciPy
RUN apt-get install -y libblas-dev
RUN apt-get install -y liblapack-dev
RUN apt-get install -y gfortran

# Update pip
RUN pip install -U pip
RUN pip install -U setuptools


##################
# Reference Data #
##################

# Install reference data under /opt

WORKDIR /opt

# Extract STIPS reference data
RUN wget -qO- https://stsci.box.com/shared/static/4nebx2ndxr7c77lgocfbvxo7c2hyd3in.tgz | tar xvz
ENV stips_data /opt/stips_data

# Extract PySynphot reference data
RUN wget -qO- https://ssb.stsci.edu/trds/tarfiles/synphot1.tar.gz | tar xvz
RUN wget -qO- https://ssb.stsci.edu/trds/tarfiles/synphot2.tar.gz | tar xvz
RUN wget -qO- https://ssb.stsci.edu/trds/tarfiles/synphot5.tar.gz | tar xvz
RUN wget -qO- https://ssb.stsci.edu/trds/tarfiles/synphot6.tar.gz | tar xvz
ENV PYSYN_CDBS /opt/grp/hst/cdbs

# Extract Pandeia reference data
RUN wget -qO- https://stsci.box.com/shared/static/7voehzi5krrpml5wgyg8bo954ew7arh2.gz | tar -xvz
ENV pandeia_refdata /opt/pandeia_data-1.5.2rc2_roman

# Extract WebbPSF reference data
# (note: version number env vars are declared close to where they are used
# to prevent unnecessary rebuilds of the Docker image)
ENV WEBBPSF_DATA_VERSION 0.9.0
RUN wget -qO- https://stsci.box.com/shared/static/qcptcokkbx7fgi3c00w2732yezkxzb99.gz | tar xvz
ENV WEBBPSF_PATH /opt/webbpsf-data

# Prepare environment variables
ENV PYHOME /opt/conda
ENV PYTHON_VERSION 3.7
ENV PATH $HOME/bin:$PATH
ENV LD_LIBRARY_PATH $HOME/lib:$LD_LIBRARY_PATH


########################################
# Setup Conda and Install Requirements #
########################################

WORKDIR $HOME

#RUN conda env update --file $STIPSDIR/environment.yml
RUN conda env create -f $STIPSDIR/environment.yml
ENV PATH /opt/conda/envs/stips/bin:$PATH
RUN echo "conda activate stips" >> /root/.bashrc
ENV CONDA_DEFAULT_ENV stips

# Prepare environment variables
ENV WEBBPSF_SKIP_CHECK 1


##############################
# Install Other Requirements #
##############################

WORKDIR $HOME

# Montage
WORKDIR /opt
RUN git clone https://github.com/Caltech-IPAC/Montage.git
# ADD initdistdata.c ./Montage/lib/src/two_plane_v1.1/initdistdata.c
WORKDIR /opt/Montage
RUN make
RUN pip install montage-wrapper
ENV PATH /opt/Montage/bin:$PATH
WORKDIR $HOME


#################
# Install STIPS #
#################

WORKDIR $STIPSDIR
RUN python setup.py develop

WORKDIR $HOME
