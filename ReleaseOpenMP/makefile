################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

MSTOOLKIT = ../src/MSToolkit
COMETSEARCH = ../src/CometSearch
override CXXFLAGS = -Wall -Wextra -static -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -I$(MSTOOLKIT)/include 


# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include src/Scores/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs
-include ../compiler.mk

# Add inputs and outputs from these tool invocations to the build variables 


OBJS += SiprosReader.o


# All Target
all: Sipros_OpenMP

# Tool invocations
# g++ -std=c++11 -fopenmp -g3 -o "Sipros_OpenMP" $(OBJS) $(USER_OBJS) $(LIBS)
Sipros_OpenMP: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker'
	cd ./src/MSToolkit ; make lite
	-$(CC) $(OPTS) $(CXXFLAGS) -L./src/MSToolkit -o "Sipros_OpenMP" $(OBJS) $(USER_OBJS) $(LIBS)
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS)
	@echo 'Finished building target: $@'
	@echo ' '
	
SiprosReader.o: ../src/SiprosReader.cpp $(DEPS)
	${CXX} ${CXXFLAGS} ../src/SiprosReader.cpp -c

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) Sipros_OpenMP
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
