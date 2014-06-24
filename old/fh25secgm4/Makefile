#Makefile made for the system of convertion of POTLIB code
#Written by Aaron Tagliaboschi <aaron.tagliaboschi@gmail.com>

SYSTEM=fh2#Put your system name here
FULLNAME=fh25secgm4#Put full name of library here)
FC=gfortran#Put fortran copiler command here

origintest: $(FULLNAME)_original.f origintest.f
	$(FC) $(FULLNAME)_original.f origintest.f -o origintest
	./origintest >origintest.out

potlibtest: $(FULLNAME).f
	$(FC) utility.f $(FULLNAME).f tester.f -o potlibtest
	./potlibtest <$(SYSTEM).dat >potlibtest.out

clean:
	mv junk* native* obsolete
