# Dockerfile for building third-party applications for the
# BioSimSpace notebook server.

FROM siremol/sire-devel:latest

USER root

RUN yum -y install patch vim

WORKDIR $HOME

USER $FN_USER

ARG par_url
ENV PAR_URL=$par_url

# Add Sire bin directory to PATH.
ENV PATH=$HOME/plumed/bin:$HOME/sire.app/bin:$PATH
ENV PLUMED_KERNEL=$HOME/plumed/lib/libplumedKernel.so
ENV LD_LIBARY_PATH=$HOME/plumed/lib:$HOME/sire.app/lib:$LD_LIBARY_PATH

# Install the additional Conda packages that are required for the build process.
RUN $HOME/sire.app/bin/conda install -y -c conda-forge binutils c-compiler compilers \
    cxx-compiler fortran-compiler gsl openblas openmpi-mpicc

# Download and build PLUMED.
RUN git clone https://github.com/plumed/plumed2.git && \
    cd plumed2 && \
    ./configure --prefix=$HOME/plumed && \
    make CC=mpicc LDFLAGS="-L$HOME/sire.app/lib -lgsl -lopenblas -lz" -j 4 && \
    make install && \
    cd $HOME && \
    tar -cjf  plumed.tar.bz2 plumed

# Download and build GROMACS, patched with PLUMED.
# Note tha CMake fails to link against various libraries, hence the use of CMAKE_EXE_LINKER_FLAGS.
RUN wget http://ftp.gromacs.org/pub/gromacs/gromacs-2019.3.tar.gz -O gromacs-src.tar.gz && \
    mkdir gromacs-build && \
    tar -xzf gromacs-src.tar.gz -C gromacs-build --strip 1 && \
    cd gromacs-build && \
    export LD_LIBRARY_PATH=$LD_LIBARY_PATH && \
    echo 2 | $HOME/plumed/bin/plumed patch -p && \
    cmake . -DGMX_BUILD_OWN_FFTW=ON \
            -DGMX_SIMD=AVX_512 \
            -DCMAKE_INSTALL_PREFIX=$HOME/gromacs \
            -DCMAKE_EXE_LINKER_FLAGS="-L$HOME/plumed/lib -L$HOME/sire.app/lib /home/sireuser/sire.app/lib/libopenblas.so.0 -lgsl -lmpi_cxx -lplumed -lplumedKernel -lz" && \
    make -j 4 && \
    make install && \
    cd $HOME && \
    tar -cjf gromacs.tar.bz2 gromacs

# Download AmberTools19 from Anaconda cloud.
RUN mkdir amber && \
    wget https://anaconda.org/AmberMD/ambertools/19.0/download/linux-64/ambertools-19.0-0.tar.bz2 -O amber-conda.tar.bz2 && \
    tar -xjf amber-conda.tar.bz2 -C amber && \
    tar -cjf amber.tar.bz2 amber

ENTRYPOINT ["bash"]
