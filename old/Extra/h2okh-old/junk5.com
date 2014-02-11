ifort tester.f junk5.f ../../Section_1/Utilities/utility.f -o junk5.x
junk5.x <tester.dat >junk5.out
cat junk5.out

