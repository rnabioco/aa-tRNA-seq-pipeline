#!/bin/bash
# Activation script for demux environment
# Ensures conda libraries take precedence over system libraries for PyTorch compatibility
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# Auto-install WarpDemuX if not present
if ! python -c "import warpdemux" 2>/dev/null; then
    echo "Installing WarpDemuX..."
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_DIR="$SCRIPT_DIR/../resources/tools/WarpDemuX"
    if [ ! -d "$REPO_DIR" ]; then
        git clone --recursive https://github.com/KleistLab/WarpDemuX.git "$REPO_DIR"
    fi
    uv pip install -e "$REPO_DIR"
fi
