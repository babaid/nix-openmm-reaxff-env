#!/bin/bash

nix-channel --add https://github.com/nix-community/nixGL/archive/main.tar.gz nixgl && nix-channel --update
nix-env -iA nixgl.auto.nixGLDefault  

exname=$(ls ~/.nix-profile/bin/ | grep nixGLNvidia)
echo "To be able to use the GPU in the Nix environment, use following command when entering it: ${exname} nix-shell"
