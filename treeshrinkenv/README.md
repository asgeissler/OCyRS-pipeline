See ../renv/README.md for more details.

This folder contains the info for the `treeshrink` tool.
The image was created as described below for platform locked 
conda enviornment (on Friday 26th 2022).

        conda lock -p linux-64 -f treeshrinkenv.yml
        sudo singularity build treeshrinkenv.sif Singularity
