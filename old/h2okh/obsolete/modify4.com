gfortran nativetest4.f junk4.f ../../Version2/Utilities/utility.f -o junk4.x
./junk4.x <nativetest.datstd >junk4.out
diff junk4.out nativetest.outstd
