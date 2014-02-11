ifort h2okhtest3.f junk3.f ../../Section_1/Utilities/utility.f -o junk3.x
junk3.x <h2okhtest.datstd >junk3.out
diff junk3.out h2okhtest.outstd
