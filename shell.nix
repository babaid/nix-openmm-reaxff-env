{ pkgs ? import <nixpkgs> {config = {allowUnfree=true;}; } }:

let
  python = pkgs.python3;
  pyPkgs = python.pkgs;

  packages = import ./default.nix { inherit pkgs; };

  pythonEnv = python.withPackages (ps: with ps; [
    ipython
    jupyterlab
    matplotlib
    numpy
    rdkit
    pint
    cachetools
    networkx
    pandas
    scipy
    openpyxl
    lxml
    alive-progress
    pycuda
  ]);

  pythonSitePackages = builtins.concatStringsSep ":" [
    "${packages.openmm-reaxff}/lib/${python.libPrefix}/site-packages"
    "${packages.openmmforcefields}/lib/${python.libPrefix}/site-packages"
    "${packages.openff-toolkit}/lib/${python.libPrefix}/site-packages"
    "${packages.openff-units}/lib/${python.libPrefix}/site-packages"
    "${packages.openff-utilities}/lib/${python.libPrefix}/site-packages"
  ];
in

pkgs.mkShell {
  name = "openmm-env";

  buildInputs = [
    pythonEnv
    packages.openmm-reaxff
    packages.openmmforcefields
    packages.openff-toolkit
    packages.openff-units
    packages.openff-utilities
    packages.pymbar
    pkgs.git
    pkgs.openbabel
    pkgs.cudatoolkit
    pkgs.libglvnd
  ];

  shellHook = ''
    export PYTHONPATH=${pythonSitePackages}:$PYTHONPATH
    export LD_LIBRARY_PATH=${pkgs.libGL}/lib:${pkgs.libGLU}/lib:${pkgs.freeglut}/lib:${pkgs.xorg.libX11}/lib:${pkgs.stdenv.cc.cc.lib}/lib:${pkgs.cudatoolkit}/lib:${pkgs.cudatoolkit.lib}/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$(nixGLNvidia printenv LD_LIBRARY_PATH):$LD_LIBRARY_PATH
  '';
}

