SYSTEM=clh2
FULLNAME=clh2g3ta4
FC=gfortran#Put fortran copiler command here

nativetest:
	$(FC) $(FULLNAME)_native.f nativetest.f -o nativetest
	./nativetest > nativetest.out
	
potlib:
	$(FC) utility.f $(FULLNAME).f tester.f -o potlib
	./potlib < $(SYSTEM).dat > $(FULLNAME).out
