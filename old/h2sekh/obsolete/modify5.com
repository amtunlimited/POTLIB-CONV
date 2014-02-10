gfortran ../../Version2/Utilities/tester.f junk5.f ../../Version2/Utilities/utility.f -o junk5.x
./junk5.x <tester.dat >junk5.out
cp junk5.out tester.std
cat junk5.out
