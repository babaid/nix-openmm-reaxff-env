{ pkgs ? import <nixpkgs> {config = {allowUnfree=true; } ;} }:

let
  python = pkgs.python3;

  openmm-reaxff = import ./openmm.nix {
  inherit (pkgs) stdenv lib fetchFromGitHub cmake gfortran fftwSinglePrec doxygen swig;
  enablePython = true;
  python3Packages = pkgs.python3Packages;
  enableOpencl = true;
  opencl-headers = pkgs.opencl-headers;
  ocl-icd = pkgs.ocl-icd;
  config = pkgs.config;
  enableCuda = true;  # or true if you want to enable it
  cudaPackages = pkgs.cudaPackages or {};
  addDriverRunpath = pkgs.addDriverRunpath or (_: "");
  fetchgit = pkgs.fetchgit;
  };



  openmmforcefields = pkgs.stdenv.mkDerivation {
    name = "openmmforcefields";
    src = pkgs.fetchgit {
      url = "https://github.com/openmm/openmmforcefields";
      rev = "HEAD";
      sha256 = "sha256-halp6QCf6C4hByMQ8vvxCVyg936O/1QQs6laF4yryqA=";
    };

    nativeBuildInputs = with pkgs.python3Packages; [ pip setuptools wheel build ];

    installPhase = ''
      ${pkgs.python3.interpreter} -m pip install . --prefix=$out --no-deps --no-build-isolation
    '';
  };

  openff-units = pkgs.stdenv.mkDerivation {
    name = "openff-units";
    src = pkgs.fetchgit {
      url = "https://github.com/openforcefield/openff-units";
      rev = "HEAD";
      sha256 = "sha256-CJynoRHmLwpwx4bniiSq/cBMvjLOBvhrE8SBpQlpbmg=";  # replace after nix-build
    };

    nativeBuildInputs = with pkgs.python3Packages; [ pip setuptools wheel build ];

    installPhase = ''
      ${pkgs.python3.interpreter} -m pip install . --prefix=$out --no-deps --no-build-isolation
    '';
  };

  openff-utilities = pkgs.stdenv.mkDerivation {
    name = "openff-utilities";
    src = pkgs.fetchgit {
      url = "https://github.com/openforcefield/openff-utilities";
      rev = "HEAD";
      sha256 = "sha256-EBWz5G3T0DYxIwm2vE38fNDq4g7pd2KFseqisgUg+Lk=";  # replace after nix-build
    };

    nativeBuildInputs = with pkgs.python3Packages; [ pip setuptools wheel build ];

    installPhase = ''
      ${pkgs.python3.interpreter} -m pip install . --prefix=$out --no-deps --no-build-isolation
    '';
  };

  openff-toolkit = pkgs.stdenv.mkDerivation {
    name = "openff-toolkit";
    src = pkgs.fetchgit {
      url = "https://github.com/openforcefield/openff-toolkit";
      rev = "HEAD";
      sha256 = "sha256-22OfbGqhlFHfcuYWdo/cSsIcePSdFAyYxxvUOe8P0G4=";
    };

    nativeBuildInputs = with pkgs.python3Packages; [ pip setuptools wheel build ];
    buildInputs = [ openff-units openff-utilities ];

    installPhase = ''
      ${pkgs.python3.interpreter} -m pip install . --prefix=$out --no-deps --no-build-isolation
    '';
  };

    pymbar = pkgs.python3Packages.buildPythonPackage rec {
    pname = "pymbar";
    version = "4.0.1";

    src = pkgs.python3Packages.fetchPypi {
      inherit pname version;
      sha256 = "sha256-bpqyiAUHKqKJmsSxTTQdGA55X25B5q6hL3Yq2ob/zpM=";  # Update if needed
    };

    propagatedBuildInputs = with pkgs.python3Packages; [
      numpy
      scipy
      matplotlib
      networkx
      pandas
      versioneer
      pkgs.cudatoolkit
      pkgs.libglvnd
    ];

    doCheck = false;
  };


in {
  inherit openmm-reaxff openmmforcefields openff-units openff-utilities openff-toolkit pymbar;
}

