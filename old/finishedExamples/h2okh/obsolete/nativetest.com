ifort nativetest.f native.f -o nativetest.x
nativetest.x <nativetest.datstd >nativetest.outstd
