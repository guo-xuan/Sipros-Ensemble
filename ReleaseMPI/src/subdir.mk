################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/directoryStructure.cpp \
../src/isotopologue.cpp \
../src/mpimain.cpp \
../src/ms2scan.cpp \
../src/ms2scanvector.cpp \
../src/peptide.cpp \
../src/proNovoConfig.cpp \
../src/proteindatabase.cpp \
../src/ptm.cpp \
../src/tokenvector.cpp 

OBJS += \
./src/directoryStructure.o \
./src/isotopologue.o \
./src/mpimain.o \
./src/ms2scan.o \
./src/ms2scanvector.o \
./src/peptide.o \
./src/proNovoConfig.o \
./src/proteindatabase.o \
./src/ptm.o \
./src/tokenvector.o 

CPP_DEPS += \
./src/directoryStructure.d \
./src/isotopologue.d \
./src/mpimain.d \
./src/ms2scan.d \
./src/ms2scanvector.d \
./src/peptide.d \
./src/proNovoConfig.d \
./src/proteindatabase.d \
./src/ptm.d \
./src/tokenvector.d 

override CXXFLAGS = -Wextra -static -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -I$(MSTOOLKIT)/include


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: C++ Compiler'
	-$(MCC) $(MOPTS) $(CXXFLAGS) -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


