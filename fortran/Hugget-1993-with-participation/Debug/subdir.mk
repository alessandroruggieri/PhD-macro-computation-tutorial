################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F95_SRCS += \
../globals.f95 \
../main.f95 \
../mathutils.f95 \
../sub_equilibrium.f95 \
../sub_model.f95 \
../sub_print.f95 \
../sub_read.f95 \
../sub_sim.f95 

OBJS += \
./globals.o \
./main.o \
./mathutils.o \
./sub_equilibrium.o \
./sub_model.o \
./sub_print.o \
./sub_read.o \
./sub_sim.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f95
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

globals.o: ../globals.f95

main.o: ../main.f95 globals.o

mathutils.o: ../mathutils.f95

sub_equilibrium.o: ../sub_equilibrium.f95 globals.o mathutils.o

sub_model.o: ../sub_model.f95 globals.o mathutils.o

sub_print.o: ../sub_print.f95 globals.o

sub_read.o: ../sub_read.f95 globals.o

sub_sim.o: ../sub_sim.f95 globals.o mathutils.o


