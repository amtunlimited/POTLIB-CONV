#Makefile made for the system of convertion of POTLIB code
#Written by Aaron Tagliaboschi <aaron.tagliaboschi@gmail.com>

SYSTEM=nafh#Put your system name here
FULLNAME=nafh#Put full name of library here)
FC=gfortran#Put fortran copiler command here

origintest: $(FULLNAME)_original.f origintest.f
	$(FC) $(FULLNAME)_original.f origintest.f -o origintest
	./origintest >origintest.out

potlibtest: $(FULLNAME).f
	$(FC) utility.f $(FULLNAME).f tester.f -o potlibtest
	./potlibtest <$(SYSTEM).dat >potlibtest.out

clean:
	mv junk* native* obsolete
