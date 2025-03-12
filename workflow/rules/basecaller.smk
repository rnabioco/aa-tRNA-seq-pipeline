DORADO_URL = f"https://cdn.oxfordnanoportal.com/software/analysis/dorado-{DORADO_VERSION}-linux-x64.tar.gz"

rule setup_dorado:
    output:
        dorado_bin = f"{DORADO_DIR}/bin/dorado"
    params:
        dorado_url = DORADO_URL,
        dorado_dir = DORADO_DIR
    shell:
        """
        # Create directory structure if it doesn't exist
        mkdir -p {params.dorado_dir}
        
        # Download Dorado tarball
        curl -L -o dorado.tar.gz {params.dorado_url}
       
        # Unpack the tarball to the specified directory
        tar -xzf dorado.tar.gz -C {params.dorado_dir} --strip-components=1
        
        # Remove the tarball
        rm dorado.tar.gz
        
        # Make the binary executable (just to be sure)
        chmod +x {output.dorado_bin}
        """

rule dorado_model:
  """
  download dorado base-calling model
  """
    output:
        os.path.join("resources/models")
    log:
        os.path.join(outdir, "logs", "dorado")
    params:
        model: config["dorado_model"]    
    shell:
        """
    dorado download --model {params.model} --models-directory {output}
    """
