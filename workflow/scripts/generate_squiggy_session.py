"""Generate squiggy session JSON file for loading samples in Positron."""

import argparse
import datetime
import hashlib
import json
import os


def calculate_md5(filepath, chunk_size=8192):
    """Calculate MD5 hash of a file."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def get_file_info(filepath):
    """Get file metadata including MD5 checksum."""
    stat = os.stat(filepath)
    return {
        "md5": calculate_md5(filepath),
        "size": stat.st_size,
        "lastModified": datetime.datetime.fromtimestamp(
            stat.st_mtime, tz=datetime.timezone.utc
        ).isoformat(),
    }


def make_relative_path(filepath, base_dir):
    """Convert absolute path to relative path from base directory."""
    abs_path = os.path.abspath(filepath)
    abs_base = os.path.abspath(base_dir)
    return os.path.relpath(abs_path, abs_base)


def generate_squiggy_session(
    sample_names,
    output_dir,
    fasta_path,
    pod5_paths=None,
    bam_paths=None,
    session_name=None,
    compute_checksums=True,
):
    """
    Generate squiggy session JSON for loading pipeline outputs.

    Args:
        sample_names: List of sample names
        output_dir: Pipeline output directory (paths will be relative to this)
        fasta_path: Path to reference FASTA file
        pod5_paths: List of pod5 file paths (must match order of sample_names)
        bam_paths: List of BAM file paths (must match order of sample_names)
        session_name: Optional session name (defaults to directory name)
        compute_checksums: Whether to compute MD5 checksums for files

    Returns:
        dict: Session data structure
    """
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()

    if session_name is None:
        session_name = f"aa-tRNA-seq: {os.path.basename(os.path.abspath(output_dir))}"

    samples = {}
    file_checksums = {}

    # Make fasta path relative to output directory
    fasta_rel = make_relative_path(fasta_path, output_dir)

    # Compute fasta checksum once (shared across samples)
    if compute_checksums and os.path.exists(fasta_path):
        file_checksums[fasta_rel] = get_file_info(fasta_path)

    for i, sample in enumerate(sample_names):
        # Use explicit paths if provided, otherwise build default paths
        if pod5_paths and i < len(pod5_paths):
            pod5_path = pod5_paths[i]
        else:
            pod5_path = os.path.join(output_dir, "pod5", sample, f"{sample}.pod5")

        if bam_paths and i < len(bam_paths):
            bam_path = bam_paths[i]
        else:
            bam_path = os.path.join(output_dir, "bam", "final", sample, f"{sample}.bam")

        # Convert to relative paths
        pod5_rel = make_relative_path(pod5_path, output_dir)
        bam_rel = make_relative_path(bam_path, output_dir)

        samples[sample] = {
            "pod5Paths": [pod5_rel],
            "bamPath": bam_rel,
            "fastaPath": fasta_rel,
        }

        # Compute checksums
        if compute_checksums:
            if os.path.exists(pod5_path):
                file_checksums[pod5_rel] = get_file_info(pod5_path)
            if os.path.exists(bam_path):
                file_checksums[bam_rel] = get_file_info(bam_path)

    session = {
        "version": "1.0.0",
        "timestamp": timestamp,
        "sessionName": session_name,
        "samples": samples,
        "plotOptions": {
            "mode": "EVENTALIGN",
            "normalization": "ZNORM",
            "showDwellTime": False,
            "showBaseAnnotations": True,
            "scaleDwellTime": False,
            "downsample": 5,
            "showSignalPoints": False,
        },
        "ui": {
            "expandedSamples": list(sample_names),
            "selectedSamplesForComparison": [],
        },
    }

    if compute_checksums and file_checksums:
        session["fileChecksums"] = file_checksums

    return session


def main():
    parser = argparse.ArgumentParser(
        description="Generate squiggy session JSON for pipeline outputs"
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Sample names to include in session",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Pipeline output directory (paths will be relative to this)",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to reference FASTA file",
    )
    parser.add_argument(
        "--pod5s",
        nargs="+",
        default=None,
        help="Pod5 file paths (must match order of --samples)",
    )
    parser.add_argument(
        "--bams",
        nargs="+",
        default=None,
        help="BAM file paths (must match order of --samples)",
    )
    parser.add_argument(
        "--session-name",
        default=None,
        help="Optional session name",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output JSON file path",
    )
    parser.add_argument(
        "--no-checksums",
        action="store_true",
        help="Skip computing MD5 checksums (faster but no integrity verification)",
    )

    args = parser.parse_args()

    session = generate_squiggy_session(
        sample_names=args.samples,
        output_dir=args.output_dir,
        fasta_path=args.fasta,
        pod5_paths=args.pod5s,
        bam_paths=args.bams,
        session_name=args.session_name,
        compute_checksums=not args.no_checksums,
    )

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    with open(args.output, "w") as f:
        json.dump(session, f, indent=2)

    print(f"Squiggy session written to: {args.output}")


if __name__ == "__main__":
    main()
