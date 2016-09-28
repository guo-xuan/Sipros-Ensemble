################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Scores/CometSearch.cpp \
../src/Scores/MVH.cpp 

OBJS += \
./src/Scores/CometSearch.o \
./src/Scores/MVH.o 

CPP_DEPS += \
./src/Scores/CometSearch.d \
./src/Scores/MVH.d 


# Each subdirectory must supply rules for building sources it contributes
src/Scores/%.o: ../src/Scores/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	/home/xgo/local/bin/mpiCC -fopenmp -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


