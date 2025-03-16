# Determine OS-specific and architecture-specific Dorado URL
import platform

system = platform.system().lower()
machine = platform.machine().lower()

# Detect architecture
if machine in ["arm64", "aarch64"]:
    arch = "arm64"
elif machine in ["x86_64", "amd64", "x64"]:
    arch = "x64"
else:
    raise ValueError(f"Unsupported architecture: {machine}")

# Construct OS-specific suffix and file extension
if system == "linux":
    os_suffix = f"linux-{arch}"
    file_ext = "tar.gz"
elif system == "darwin":
    os_suffix = f"osx-{arch}"
    file_ext = "zip"
else:
    raise ValueError(f"Unsupported operating system: {system}")

DORADO_URL = f"https://cdn.oxfordnanoportal.com/software/analysis/dorado-{DORADO_VERSION}-{os_suffix}.{file_ext}"


rule setup_dorado:
    output:
        dorado_bin=f"{DORADO_DIR}/bin/dorado",
    params:
        dorado_url=DORADO_URL,
        dorado_dir=DORADO_DIR,
        file_ext=file_ext,
    shell:
        """
        # Create directory structure if it doesn't exist
        mkdir -p {params.dorado_dir}

        # Download Dorado tarball
        curl -L -o dorado.{params.file_ext} {params.dorado_url}

        # Extract based on file type
        if [ "{params.file_ext}" = "tar.gz" ]; then
            tar -xzf dorado.{params.file_ext} -C {params.dorado_dir} --strip-components=1
        elif [ "{params.file_ext}" = "zip" ]; then
            unzip -o dorado.{params.file_ext} -d {params.dorado_dir}_temp
            # Find the extracted directory and move its contents
            mv {params.dorado_dir}_temp/*/* {params.dorado_dir}/
            rm -rf {params.dorado_dir}_temp
        fi

        # Remove the tarball
        rm dorado.{params.file_ext}

        # Make the binary executable (just to be sure)
        chmod +x {output.dorado_bin}
        """


rule dorado_model:
    """
    Download dorado base-calling model
    """
    output:
        directory(os.path.join("resources", "models", config["dorado_model"])),
    log:
        os.path.join(outdir, "logs", "dorado_model.log"),
    params:
        model=config["dorado_model"],
        model_dir=os.path.join("resources", "models"),
    shell:
        """
        mkdir -p {params.model_dir}

        # Run Dorado download
        dorado download --model {params.model} --models-directory {params.model_dir} > {log} 2>&1

        # Create a marker file if needed
        if [ ! -d "{output}" ]; then
            echo "Error: Model directory not created: {output}" >> {log}
            exit 1
        fi

        # List the contents of the model directory to help with debugging
        echo "Contents of {output}:" >> {log}
        ls -la {output} >> {log} 2>&1

        # Touch a marker file in case the download doesn't create any files
        # This ensures the output directory is not empty
        touch {output}/.downloaded
        """


rule setup_modkit:
    """
    Install modkit from pre-built binary or source depending on the system
    """
    output:
        modkit_bin=f"{MODKIT_DIR}/bin/modkit",
    params:
        modkit_dir=MODKIT_DIR,
        modkit_version=MODKIT_VERSION,
        modkit_repository="https://github.com/nanoporetech/modkit",
        linux_binary_url=f"https://github.com/nanoporetech/modkit/releases/download/v{MODKIT_VERSION}/modkit_v{MODKIT_VERSION}_u16_x86_64.tar.gz",
    log:
        os.path.join(outdir, "logs", "setup_modkit.log"),
    shell:
        """
        if ! command -v rustc &> /dev/null; then
            echo "Rust not found, installing..." >> {log} 2>&1
            curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y >> {log} 2>&1
            source "$HOME/.cargo/env"
        fi

        export PATH="$HOME/.cargo/bin:$PATH"

        export CARGO_NET_GIT_FETCH_WITH_CLI=true
        cargo install --git {params.modkit_repository} \\
                     --tag v{params.modkit_version} \\
                     --root {params.modkit_dir} \\
                     --jobs 4 >> {log} 2>&1

        """
