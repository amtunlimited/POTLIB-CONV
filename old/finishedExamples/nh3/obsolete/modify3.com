ifort nativetest3.f junk3.f ../../Version2/Utilities/utility.f -o junk3.x
junk3.x <nativetest.datstd >junk3.out
diff junk3.out nativetest.outstd
