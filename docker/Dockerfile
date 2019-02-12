# Use an official Python runtime as a parent image
FROM gadgetron/ubuntu_1804_cuda92_cudnn7_base:latest

# set arguments
ARG GADGETRON_URL=https://github.com/gadgetron/gadgetron
ARG ISMRMRD_URL=https://github.com/ismrmrd/ismrmrd.git
ARG ISMRMRD_PYTHON_API_URL=https://github.com/ismrmrd/ismrmrd-python.git
ARG ISMRMRD_PYTHON_TOOLS_URL=https://github.com/ismrmrd/ismrmrd-python-tools.git
ARG BART_URL=https://github.com/mrirecon/bart.git
ARG SIEMENS_TO_ISMRMRD_URL=https://github.com/ismrmrd/siemens_to_ismrmrd.git
ARG PHILIPS_TO_ISMRMRD_URL=https://github.com/ismrmrd/philips_to_ismrmrd.git
ARG ITK_URL=http://itk.org/ITK.git
ARG ITK_BRANCH=v4.13.1
ARG ELASTIX_URL=https://github.com/SuperElastix/elastix.git
ARG CSMOCO_URL=https://github.com/thomaskuestner/CS_MoCo_LAB.git

# Copy the current directory contents into the container at /app
ADD . /tmp

# install newest version of cmake
RUN	apt-get -y purge cmake \
	&& wget https://cmake.org/files/v3.13/cmake-3.13.0.tar.gz \
	&& tar -xzvf cmake-*.tar.gz \
	&& cd cmake-3*/ \
	&& ./bootstrap \
	&& make -j$(expr $(nproc) + 1)\
	&& make install \
	&& cd .. \
	&& rm -rf cmake-3* \
	# update system
	&& apt-get -y update \
	&& apt-get -y upgrade \
	# install NumPy dependency
	&& apt-get -y install python-numpy \
	# install other software
	&& apt-get -y install gcc-6 g++-6 gfortran \
	# and cleanup
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# install ISMRMRD
RUN	cd /opt/code \
	&& git clone ${ISMRMRD_URL} \
	&& cd ismrmrd \
	&& mkdir build \
	&& cd build \
	&& cmake ../ \
	&& make -j$(expr $(nproc) + 1) \
	&& make install \
	&& cd ../.. \
	&& rm -rf ismrmrd/build

# Fix compiler error with CUDA for now
RUN	sed -i "s;#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;//#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;g" /usr/local/cuda/include/host_config.h

# install BART
RUN	cd /opt/code \
	&& git clone ${BART_URL} \
	&& cd bart \
	&& mkdir build \
	&& cd build \
	&& cmake -DBART_FPIC=ON -DBART_ENABLE_MEM_CFL=ON -DBART_REDEFINE_PRINTF_FOR_TRACE=ON -DBART_LOG_BACKEND=ON -DBART_LOG_GADGETRON_BACKEND=ON ../ \
	&& make -j$(expr $(nproc) + 1) \
	&& make install \
	&& cd ../.. \
	&& rm -rf bart/build

# install GADGETRON
RUN	cd /opt/code \
	&& git clone ${GADGETRON_URL} \
	&& cd gadgetron \
	&& mkdir build \
	&& cd build \
	&& cmake -DCUDA_NVCC_FLAGS="-ccbin gcc-6" -DBUILD_WITH_PYTHON3=ON ../ \
	&& make -j$(expr $(nproc) + 1) \
	&& make install \
	&& /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` \
	&& cp ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml.example ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml \
	&& cp /opt/code/gadgetron/docker/start_supervisor /opt/ \
	&& cp /opt/code/gadgetron/docker/supervisord.conf /opt/ \
	&& cd ../.. \
	&& rm -rf gadgetron/build

# HASH for ISMRMRD
RUN	cd /opt/code/ismrmrd \
	&& /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrd.sha1 --value `git rev-parse HEAD`

# install ISMRMRD PYTHON API
RUN	cd /opt/code \
	&& git clone ${ISMRMRD_PYTHON_API_URL} \
	&& cd ismrmrd-python \
	&& python3 setup.py install \
	&& /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpython.sha1 --value `git rev-parse HEAD`

# install ISMRMRD PYTHON TOOLS
RUN	cd /opt/code \
	&& git clone ${ISMRMRD_PYTHON_TOOLS_URL} \
	&& cd ismrmrd-python-tools \
	&& python3 setup.py install \
	&& /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpythontools.sha1 --value `git rev-parse HEAD`

# install SIEMENS_TO_ISMRMRD
RUN	cd /opt/code \
	&& git clone ${SIEMENS_TO_ISMRMRD_URL} \
	&& cd siemens_to_ismrmrd \
	&& mkdir build \
	&& cd build \
	&& cmake ../ \
	&& make -j$(expr $(nproc) + 1) \
	&& make install \
	&& /opt/code/gadgetron/docker/manifest --key .io.gadgetron.siemens_to_ismrmrd.sha1 --value `git rev-parse HEAD` \
	&& cd ../.. \
	&& rm -rf siemens_to_ismrmrd/build

# install PHILIPS_TO_ISMRMRD
RUN	cd /opt/code \
	&& git clone ${PHILIPS_TO_ISMRMRD_URL} \
	&& cd philips_to_ismrmrd \
	&& mkdir build \
	&& cd build \
	&& cmake ../ \
	&& make -j $(nproc) \
	&& make install \
	&& /opt/code/gadgetron/docker/manifest --key .io.gadgetron.philips_to_ismrmrd.sha1 --value `git rev-parse HEAD` \
	&& cd ../.. \
	&& rm -rf philips_to_ismrmrd/build

# install ITK
RUN	cd /opt/code \
	&& git clone -b ${ITK_BRANCH} ${ITK_URL} \
	&& mkdir ITK/build \
	&& cd ITK/build \
	&& cmake ../ \
	&& make -j$(expr $(nproc) + 1)\
	&& make install \
	&& cd ../.. \
	&& rm -rf ITK/build

# install elastix (twice)
RUN	cd /opt/code \
	&& git clone ${ELASTIX_URL} \
	&& mkdir elastix/build \
	&& cd elastix/build \
	&& cmake ../ \
	&& make -j$(expr $(nproc) + 1)\
	&& make install \
	&& cd ../.. \
	&& patch elastix/CMakeLists.txt < /tmp/elastix.patch \
	&& rm -rf elastix/build \
	&& mkdir elastix/build \
	&& cd elastix/build \
	&& cmake ../ \
	&& make -j$(expr $(nproc) + 1)\
	&& make install \
	&& cd ../.. \
	&& mkdir /usr/local/lib/cmake/Elastix \
	&& ln -s /opt/code/elastix/build/ElastixConfig.cmake /usr/local/lib/cmake/Elastix/ElastixConfig.cmake

# install CS_MoCo_LAB
RUN	cd /opt/code \
	&& git clone ${CSMOCO_URL} \
	&& cd CS_MoCo_LAB/reconstruction/gadgetron/CS_LAB_Gadget \
	&& mkdir build \
	&& cd build \
	&& cmake ../ \
	&& make -j$(expr $(nproc) + 1)\
	&& make install \
	&& cd .. \
	&& rm -rf build

# Set the shell entry point to the container
WORKDIR /opt/data

# set CMD
CMD ["/opt/start_supervisor"]
