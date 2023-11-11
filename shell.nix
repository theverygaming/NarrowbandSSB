with import <nixpkgs> { };
let
  gccForLibs = stdenv.cc.cc;
in
stdenv.mkDerivation {
  name = "NarrowbandSSB";
  buildInputs = [
    cmake
    pkg-config
    volk
    fftwFloat
    libsndfile
  ];
}
