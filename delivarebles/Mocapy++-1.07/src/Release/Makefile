INCL = -I$(BOOSTPATH)

OBJS =	abstractinfengine.o dbn.o discretedensities.o discreteess.o discretepriors.o emengine.o essbase.o gaussiandensities.o gaussianess.o gibbsrandom.o infenginemcmc.o mcmc.o multigauss.o node.o parentmap.o random_data.o utils.o vonmises.o vonmisesdensities.o vonmisesess.o \
poissondensities.o poissoness.o MultinomialESS.o MultinomialDensity.o MultinomialDensities.o LogFactorial.o MultinomialESS.o moreabstractinfengine.o sampleinfengine.o forwardbacktracker.o infenginehmm.o

LIBS = 

TARGET = libMocapy.a

CPP_DEPS += \
./abstractinfengine.d \
./childnode.d \
./discretedensities.d \
./discreteess.d \
./dbn.d \
./emengine.d \
./essbase.d \
./gaussiandensities.d \
./gaussianess.d \
./gibbsrandom.d \
./infenginemcmc.d \
./mcmc.d \
./multigauss.d \
./node.d \
./parentmap.d \
./random_data.d \
./utils.d \
./vonmises.d \
./vonmisesdensities.d \
./vonmisesess.d \
./poissondensities.d \
./poissoness.d \
./LogFactorial.d \
./MultinomialDensities.d \
./MultinomialDensity.d \
./MultinomialESS.d \
./moreabstractinfengine.d \
./sampleinfengine.d \
./forwardbacktracker.d \
./infenginehmm.d

-include $(CPP_DEPS)

%.o:	../%.cpp 
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DNDEBUG $(INCL) -O3 -Wall -c -fmessage-length=0 -MMD -MP  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


$(TARGET):	$(OBJS)
	ar -r $(TARGET) $(OBJS) $(LIBS)
	mv $(TARGET) ../../libs

all:	$(TARGET)

clean:
	rm -f *.o *.d *~ $(TARGET)
