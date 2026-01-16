#!/bin/bash
# Activation script for demux environment
# Sets up library paths - does NOT install packages
# For WarpDemuX installation, run: pixi run -e demux install-warpdemux

# Ensures conda libraries take precedence over system libraries for PyTorch compatibility
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
