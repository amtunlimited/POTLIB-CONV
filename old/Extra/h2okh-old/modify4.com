ifort h2okhtest4.f junk4.f ../../Section_1/Utilities/utility.f -o junk4.x
junk4.x <h2okhtest.datstd >junk4.out
diff junk4.out h2okhtest.outstd
