gfortran nativetest.f junk2.f ../../Version2/Utilities/utility.f -o junk2.x
./junk2.x <nativetest.datstd >junk2.out
diff junk2.out nativetest.outstd
