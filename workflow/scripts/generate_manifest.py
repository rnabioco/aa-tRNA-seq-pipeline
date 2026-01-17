"""Generate pipeline manifest file for reproducibility tracking."""

import datetime
import getpass
import json
import os
import re
import socket
import subprocess

from git import Repo


def get_pipeline_version(pipeline_dir):
    """
    Get comprehensive pipeline version information from git.

    Returns dict with commit, tag, branch, and dirty status.
    """
    repo = Repo(pipeline_dir)
    commit = repo.head.commit

    # Find tags pointing to current commit
    tags = [tag.name for tag in repo.tags if tag.commit == commit]

    # Get branch name (None if detached HEAD)
    try:
        branch = repo.active_branch.name
    except TypeError:
        branch = None

    return {
        "git_commit": str(commit),
        "git_tag": tags[0] if tags else None,
        "git_branch": branch,
        "git_dirty": repo.is_dirty(),
    }


def parse_pixi_lock_version(lock_path, package_name):
    """
    Parse version from pixi.lock for a conda package.

    Extracts version from conda URL patterns like:
    - conda: https://conda.anaconda.org/bioconda/linux-64/samtools-1.23-h96c455f_0.conda

    Returns version string or None if not found.
    """
    if not os.path.exists(lock_path):
        return None

    # Pattern to match package-version-build in conda URLs
    pattern = rf"/{re.escape(package_name)}-([0-9][0-9a-zA-Z\.\-]*)-[^/]+\.(?:conda|tar\.bz2)"

    with open(lock_path) as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                return match.group(1)
    return None


def get_command_version(cmd, version_flag="--version"):
    """
    Get version by running a command with version flag.

    Returns version string or None if command fails.
    """
    try:
        result = subprocess.run(
            [cmd, version_flag],
            capture_output=True,
            text=True,
            timeout=30,
        )
        output = result.stdout.strip() or result.stderr.strip()
        # Extract version from common patterns
        # dorado outputs: "dorado 1.3.1"
        # bwa outputs: version info in stderr
        if output:
            # Try to find version pattern
            version_match = re.search(r"(\d+\.\d+(?:\.\d+)?)", output)
            if version_match:
                return version_match.group(1)
        return output.split("\n")[0] if output else None
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return None


def get_python_package_version(package_name):
    """Get version of an installed Python package."""
    try:
        result = subprocess.run(
            ["python", "-c", f"import {package_name}; print({package_name}.__version__)"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass
    return None


def get_tool_versions(pipeline_dir, config):
    """
    Collect versions of all tools used in the pipeline.

    Uses pixi.lock for conda tools and runtime commands for others.
    """
    lock_path = os.path.join(pipeline_dir, "pixi.lock")

    versions = {}

    # From pixi.lock (conda packages)
    conda_tools = {
        "samtools": "samtools",
        "bwa": "bwa",
        "modkit": "ont-modkit",
        "snakemake": "snakemake",
        "bedtools": "bedtools",
    }

    for tool_name, package_name in conda_tools.items():
        version = parse_pixi_lock_version(lock_path, package_name)
        if version:
            versions[tool_name] = version

    # Dorado - from config or runtime
    dorado_version = config.get("dorado_version")
    if dorado_version:
        versions["dorado"] = dorado_version
    else:
        # Fallback to runtime check
        dorado_bin = os.path.join(
            pipeline_dir,
            "resources",
            "tools",
            "dorado",
            config.get("dorado_version", ""),
            "bin",
            "dorado",
        )
        if os.path.exists(dorado_bin):
            versions["dorado"] = get_command_version(dorado_bin)

    # Remora - Python package
    remora_version = get_python_package_version("remora")
    if remora_version:
        versions["remora"] = remora_version

    return versions


def extract_config_params(config):
    """
    Extract key configuration parameters for the manifest.

    Filters to only include reproducibility-relevant parameters.
    """
    params = {}

    # Reference files
    if "fasta" in config:
        params["fasta"] = config["fasta"]

    # Model paths and names
    if "base_calling_model" in config:
        params["base_calling_model"] = config["base_calling_model"]
    if "dorado_model" in config:
        params["dorado_model"] = config["dorado_model"]
    if "dorado_version" in config:
        params["dorado_version"] = config["dorado_version"]
    if "remora_cca_classifier" in config:
        params["remora_cca_classifier"] = config["remora_cca_classifier"]

    # Adapter sequences
    if "adapters" in config:
        params["adapters"] = config["adapters"]

    # Modkit settings
    if "modkit" in config:
        params["modkit"] = config["modkit"]

    # Output directory
    if "output_directory" in config:
        params["output_directory"] = config["output_directory"]

    # Command options
    if "opts" in config:
        params["opts"] = config["opts"]

    # WarpDemuX settings
    if "warpdemux" in config:
        params["warpdemux"] = config["warpdemux"]

    return params


def extract_sample_info(samples):
    """Extract sample information for manifest."""
    sample_names = list(samples.keys())
    input_paths = {}

    for sample, info in samples.items():
        # Convert set to list for JSON serialization
        paths = info.get("path", set())
        if isinstance(paths, set):
            paths = list(paths)
        input_paths[sample] = paths[0] if len(paths) == 1 else paths

    return {
        "count": len(sample_names),
        "names": sample_names,
        "input_paths": input_paths,
    }


def generate_manifest(
    config,
    samples,
    status,
    start_time,
    output_dir,
    pipeline_dir,
):
    """
    Generate and write the pipeline manifest JSON file.

    Args:
        config: Snakemake config dictionary
        samples: Parsed samples dictionary
        status: Pipeline status ("success" or "failed")
        start_time: ISO format timestamp when pipeline started
        output_dir: Output directory path
        pipeline_dir: Pipeline repository root path
    """
    end_time = datetime.datetime.now(datetime.timezone.utc).isoformat()

    manifest = {
        "manifest_version": "1.0",
        "pipeline": {
            "name": "aa-tRNA-seq-pipeline",
            **get_pipeline_version(pipeline_dir),
        },
        "execution": {
            "timestamp_start": start_time,
            "timestamp_end": end_time,
            "status": status,
            "hostname": socket.gethostname(),
            "user": getpass.getuser(),
            "working_directory": os.getcwd(),
        },
        "config": extract_config_params(config),
        "samples": extract_sample_info(samples),
        "tools": get_tool_versions(pipeline_dir, config),
    }

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    manifest_path = os.path.join(output_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Pipeline manifest written to: {manifest_path}")
