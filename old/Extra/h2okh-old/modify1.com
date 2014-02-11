modify1.x < h2okh.f > junk1.f
ifort h2okhtest.f junk1.f ../../Section_1/Utilities/utility.f -o junk1.x
junk1.x <h2okhtest.datstd >junk1.out
diff junk1.out h2okhtest.outstd
