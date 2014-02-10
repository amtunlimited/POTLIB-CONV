ifort ../../Version2/Utilities/tester.f junk6.f ../../Version2/Utilities/utility.f -o junk6.x
junk6.x <tester.dat >junk6.out
diff junk6.out junk5.out
