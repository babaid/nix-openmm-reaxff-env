This repo sets up an isolated and reproducible environment for the reaxFF plugin for OpenMM implemented by me.



In order to use nix you need to execute following command:

´´´bash
sh <(curl --proto '=https' --tlsv1.2 -L https://nixos.org/nix/install) --no-daemon
´´´

This will install the nix package manager. The next step is to decide whether you want to use an NVidia GPU for computations. For GPU support run the following command: 

´´´bash
bash setup-nvidia.sh
´´´

This will provide an a command one has to use in order to use a GPU inside the environment.

Finally type:

´´´bash
nix-build
´´´

and 

´´´bash
nix-shell #or nixGLNvidia-..... nix-shell
´´´


