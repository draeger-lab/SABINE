machine = $(shell uname -m)
installpath = ./bin/$(machine)/.

compileoption = -c

ifeq ($(DEBUG),true)
	compileoption += -g -O0
else
	compileoption += -O3
endif

ifneq ($(machine),x86_64)
	compiler = g++
	linkoption = -o
else
	compiler = /usr/bin/g++
	linkoption = -o
endif

compile = $(compiler) $(compileoption)
link = $(compiler) $(linkoption)

progs = cstat sstat scluster pfmic bsanno cwords costat pfmqual

sources = $(shell ls *.cpp)
objs = $(sources:.cpp=.o)
deps = .deps

#implicit compile rule
%.o : %.cpp
	$(compile) $< -o $@

bs_annotator.o: bs_annotator.cpp bs_annotator.h pfm.h pfm_helper.h \
  convolution.h stringutils.h sequences.h exceptions.h
bsanno.o: bsanno.cpp bsanno.h pfm.h pfm_helper.h convolution.h \
  stringutils.h bs_annotator.h sequences.h exceptions.h
clustermatrix.o: clustermatrix.cpp clustermatrix.h pfm.h pfm_helper.h \
  convolution.h stringutils.h simstat.h countpars.h mytimer.h \
  similaritymatrix.h sge.h exceptions.h
convolution.o: convolution.cpp convolution.h
countpars.o: countpars.cpp countpars.h pfm.h pfm_helper.h convolution.h \
  stringutils.h
countstat.o: countstat.cpp countstat.h pfm.h pfm_helper.h convolution.h \
  stringutils.h countpars.h mytimer.h
cstat.o: cstat.cpp cstat.h pfm.h pfm_helper.h convolution.h stringutils.h \
  countstat.h countpars.h mytimer.h
exceptions.o: exceptions.cpp exceptions.h
helper.o: helper.cpp helper.h
mytimer.o: mytimer.cpp mytimer.h
pfm.o: pfm.cpp pfm.h pfm_helper.h convolution.h stringutils.h
pfm_helper.o: pfm_helper.cpp pfm_helper.h
pfmic.o: pfmic.cpp pfm.h pfm_helper.h convolution.h stringutils.h
scluster.o: scluster.cpp scluster.h pfm_helper.h clustermatrix.h pfm.h \
  convolution.h stringutils.h simstat.h countpars.h mytimer.h \
  similaritymatrix.h sge.h exceptions.h
sequences.o: sequences.cpp sequences.h exceptions.h stringutils.h
sge.o: sge.cpp sge.h stringutils.h mytimer.h exceptions.h
similaritymatrix.o: similaritymatrix.cpp similaritymatrix.h pfm.h \
  pfm_helper.h convolution.h stringutils.h simstat.h countpars.h \
  mytimer.h sge.h exceptions.h
simstat.o: simstat.cpp simstat.h pfm.h pfm_helper.h convolution.h \
  stringutils.h countpars.h mytimer.h
sstat.o: sstat.cpp sstat.h similaritymatrix.h pfm.h pfm_helper.h \
  convolution.h stringutils.h simstat.h countpars.h mytimer.h sge.h \
  exceptions.h
stringutils.o: stringutils.cpp stringutils.h

#some global make rules
all:
	make $(progs)

install:
	make all
	mkdir -p $(installpath)
	cp $(progs) $(installpath)

#link options for different programs
pfmqual: pfmqual.o pfm.o pfm_helper.o convolution.o stringutils.o pfmloader.o pfm_quality.o sequences.o simstat.o countpars.o convolution.o mytimer.o countparshetero.o
	$(link) $@ $?

bsanno: bsanno.o bs_annotator.o pfm.o convolution.o sequences.o exceptions.o pfm_helper.o stringutils.o pfmloader.o countpars.o countstat.o mytimer.o sequences.o
	$(link) $@ $?

pfmic: pfmic.o pfm.o pfm_helper.o convolution.o stringutils.o pfmloader.o sequences.o
	$(link) $@ $?

cwords: cwords.o pfm.o pfm_helper.o convolution.o stringutils.o pfmloader.o sequences.o
	$(link) $@ $?

sstat: sstat.o pfm.o pfm_helper.o simstat.o countpars.o convolution.o mytimer.o similaritymatrix.o stringutils.o sge.o exceptions.o pfmloader.o countparshetero.o sequences.o
	$(link) $@ $?

costat: costat.o pfm.o pfm_helper.o coocstat.o countstat.o countpars.o countparshetero.o convolution.o mytimer.o stringutils.o pfmloader.o sequences.o
	$(link) $@ $?

cstat: cstat.o pfm.o pfm_helper.o countstat.o countpars.o convolution.o mytimer.o stringutils.o pfmloader.o sequences.o
	$(link) $@ $?

scluster: scluster.o pfm.o pfm_helper.o simstat.o countpars.o convolution.o mytimer.o clustermatrix.o sge.o stringutils.o similaritymatrix.o exceptions.o stringutils.o pfmloader.o countparshetero.o sequences.o
	$(link) $@ $?

clean:
	rm -f *.o
	rm -f $(deps)

mrproper:
	make clean
	rm -f $(progs)
