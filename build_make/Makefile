
# fdemoullin
# generate this makefile instead of using cmake to run on GPGPU-Sim

PROJECT_NAME = CWBVH_make
BUILD_DIR = ./
SRC_DIR = ../src/

# OptiX
OPTIX_INCLUDE_DIR = ~/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/include
OPTIX_LIB = /home/francois/Documents/RayTracing/NVIDIA-OptiX-SDK-5.1.0-linux64/lib64/ # make sure this is in LD_LIBRARY_PATH

# Embree3
EMBREE_INCLUDE_DIR = /usr/local/include
EMBREE_LIB = /usr/local/lib/

# CUDA
# plain cuda, use this option for running on hardware not on gpgpu-sim
# CUDAPATH = ~/cuda/include
NVCC = nvcc
NVCCFLAGS = -c -g -std=c++11
CUDAARCH = -gencode arch=compute_75,code=compute_75
LFLAGS = -lcuda --cudart shared
GPUDEBUGFLAG = -G # Only use it when debugging the kernel in cuda gdb, don't set this when running with GPGPU-Sim!

all: build clean

build: gpu cpu
	$(NVCC) $(LFLAGS) $(CUDAARCH) -o $(BUILD_DIR)/$(PROJECT_NAME) *.o -L $(OPTIX_LIB) -L $(EMBREE_LIB) -lembree3 -ltbb -loptix_prime -loptix

gpu:
	$(NVCC) --keep $(NVCCFLAGS) $(CUDAARCH) $(SRC_DIR)*.cu

cpu:			
	$(NVCC) $(NVCCFLAGS) $(CUDAARCH) $(SRC_DIR)*.cpp -I $(OPTIX_INCLUDE_DIR) -I $(EMBREE_INCLUDE_DIR)

clean:
	rm *.o
