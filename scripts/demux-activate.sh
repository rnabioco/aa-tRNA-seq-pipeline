#!/bin/bash
# Activation script for demux environment
# Ensures conda libraries take precedence over system libraries for PyTorch compatibility
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
