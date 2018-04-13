################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Scores/CometSearchMod.cpp \
../src/Scores/MVH.cpp 

OBJS += \
./src/Scores/CometSearchMod.o \
./src/Scores/MVH.o 

CPP_DEPS += \
./src/Scores/CometSearchMod.d \
./src/Scores/MVH.d 

# Each subdirectory must supply rules for building sources it contributes
# g++ -std=c++11 -fopenmp -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
src/Scores/%.o: ../src/Scores/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: C++ Compiler'
	-$(CC) $(OPTS) $(CXXFLAGS) -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


