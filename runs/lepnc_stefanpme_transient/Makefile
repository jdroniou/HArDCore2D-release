MAIN = compute_rates

# variables
OUTPUT = $(shell echo $(MAIN) | sed  s/"\(.\)"/\\U\\1/g)
COMP = gfortran
OPTS = -O5 -fbackslash -fintrinsic-modules-path /usr/local/lib/mods/
OBJS_f90 =
OBJS =
LIBS = -lhsl

all: $(OUTPUT)

###
## EXECUTION
#.PHONY: run
run: $(OUTPUT)
	./$(OUTPUT)

##
# Compute rates
compute_rates: COMPUTE_RATES
COMPUTE_RATES: compute_rates.f90
	$(COMP) $(OPTS) $< -o COMPUTE_RATES

compute_rates_run: COMPUTE_RATES
	@ ./COMPUTE_RATES

##
# NETTOYAGE
.PHONY : clean
clean: 
	@echo "Cleanup";
	rm -f *.o *.mod;
  
