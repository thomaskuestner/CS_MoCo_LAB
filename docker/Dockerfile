# Use an official Python runtime as a parent image
FROM gadgetron/ubuntu_1604_cuda75:latest

# Set the working directory
WORKDIR /opt/code

# Copy the current directory contents into the container at /app
ADD . /tmp

# Fetch, compile & install any needed packages
RUN	git clone https://github.com/mrirecon/bart.git \
	&& cd bart \
	&& make \
	&& make install \
	&& cd ..

RUN	git clone -b v4.13.0 http://itk.org/ITK.git \
	&& mkdir ITK/build \
	&& cd ITK/build \
	&& cmake ../ \
	&& make \
	&& make install \
	&& cd ../..

RUN	git clone https://github.com/SuperElastix/elastix.git \
	&& cp -R /tmp/elastix . \
	&& rm -rf /tmp/elastix \
	&& mkdir elastix/build \
	&& cd elastix/build \
	&& cmake ../ \
	&& make \
	&& make install \
	&& cd ../.. \
	&& mkdir /usr/local/lib/cmake/Elastix \
	&& ln -s /opt/code/elastix/build/ElastixConfig.cmake /usr/local/lib/cmake/Elastix/ElastixConfig.cmake \
	&& ln -s libelastix-4.8.a /usr/local/lib/libelastix.a \
	&& ln -s libtransformix-4.8.a /usr/local/lib/libtransformix.a

RUN	git clone https://github.com/thomaskuestner/CS_MoCo_LAB.git \
    && cd CS_MoCo_LAB/reconstruction/gadgetron/CS_LAB_Gadget \
	&& mkdir build \
	&& cd build \
	&& cmake ../ \
	&& make \
	&& make install \
	&& cd ../../../../..

# Set the shell entry point to the container
WORKDIR /opt/data
