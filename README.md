# modular_social_networks


A C++ implementation of the model of games on graphs proposed in [Santos et al., 2006](https://doi.org/10.1098/rspb.2005.3272).
It was used to conduct the computational experiments in:

Wechsler, D., and J. Bascompte. 2019. Thresholds in the resilience of modular
social networks to invasion by defectors. *Journal of Theoretical Biology* 460:56–63 [10.1016/j.jtbi.2018.10.018](https://doi.org/10.1016/j.jtbi.2018.10.018).


## Install and run

To compile and install the program, run the following  commands 
in a terminal (assuming the current directory is `modular_social_networks`):

`$ cmake .`

`$ make`

`$ make install`

A new directory called `bin` is created. It contains the executable 
`games_on_graphs` which can be started using the following command
(from within the `bin` directory):

`$ ./games_on_graphs -c PATH_TO_CONFIG_FILE`

While `PATH_TO_CONFIG_FILE` is the path to a configuration file containing
the desired parameters (e.g. the network). An example file is provided
in the root folder of the git project. Specific files can be created
using the python script `create_config.py`.












