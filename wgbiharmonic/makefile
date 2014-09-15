FCOMP    = ifort 
OPTS     = -c -O2 
LINKOPTS = -O2 -o 
#ALLOC = -mcmodel=medium -shared-intel
OBJS = main.o setup.o \
         assemblestiff.o exactu.o gauss.o \
	 gb.o gradient.o inverse.o \
	 unique.o mesh.o quad_rule.o \
	 rhs.o write2tec.o

BOUND_GEN:$(OBJS)
	$(FCOMP) $(LINKOPTS) ./rundir/run $(OBJS)

setup.o:setup.f90
	$(FCOMP) $(OPTS) setup.f90

main.o: main.f90 setup.o 
	$(FCOMP) $(OPTS) main.f90

assemblestiff.o:assemblestiff.f90
	$(FCOMP) $(OPTS) assemblestiff.f90

exactu.o:exactu.f90
	$(FCOMP) $(OPTS) exactu.f90

gauss.o:gauss.f90
	$(FCOMP) $(OPTS) gauss.f90

gb.o:gb.f90
	$(FCOMP) $(OPTS) gb.f90

gradient.o:gradient.f90
	$(FCOMP) $(OPTS) gradient.f90

inverse.o:inverse.f90
	$(FCOMP) $(OPTS) inverse.f90

unique.o:unique.f90
	$(FCOMP) $(OPTS) unique.f90

mesh.o:mesh.f90
	$(FCOMP) $(OPTS) mesh.f90

quad_rule.o:quad_rule.f90
	$(FCOMP) $(OPTS) quad_rule.f90

rhs.o:rhs.f90
	$(FCOMP) $(OPTS) rhs.f90

write2tec.o:write2tec.f90
	$(FCOMP) $(OPTS) write2tec.f90
clean: 
	rm *.o *.mod
