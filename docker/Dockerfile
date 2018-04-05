# Use an official Python runtime as a parent image
FROM gadgetron/ubuntu_1604_cuda75:latest

# Set the working directory
WORKDIR /opt/code

# Copy the current directory contents into the container at /app
ADD . /tmp

# Fetch, compile & install any needed packages
RUN	git clone https://github.com/mrirecon/bart.git \
	&& cd bart \
	&& make -j$(expr $(nproc) + 1)\
	&& make install \
	&& cd ..

RUN	git clone -b v4.13.0 http://itk.org/ITK.git \
	&& mkdir ITK/build \
	&& cd ITK/build \
	&& cmake ../ \
	&& make -j$(expr $(nproc) + 1)\
	&& make install \
	&& cd ../.. \
	&& rm -rf ITK/build

RUN	git clone https://github.com/SuperElastix/elastix.git \
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

RUN	git clone https://github.com/thomaskuestner/CS_MoCo_LAB.git \
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
