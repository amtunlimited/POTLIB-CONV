modify1.x < native.f > junk1.f
ifort nativetest.f junk1.f ../../Version2/Utilities/utility.f -o junk1.x
junk1.x <nativetest.datstd >junk1.out
diff junk1.out nativetest.outstd
