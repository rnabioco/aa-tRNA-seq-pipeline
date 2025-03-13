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
    download dorado base-calling model
    """
    output:
        # Create a dummy flag file to track completion
        touch(os.path.join("resources/models", f"{config['dorado_model']}.done")),
    log:
        os.path.join(outdir, "logs", "dorado"),
    params:
        model=config["dorado_model"],
        model_dir=os.path.join("resources/models"),
    shell:
        """
        dorado download --model {params.model} --models-directory {params.model_dir} > {log} 2>&1
        """
