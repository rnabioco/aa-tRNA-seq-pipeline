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


def read_samples_file(path):
    """Sample names, one per line (blank lines and `#` comments ignored)."""
    with open(path) as f:
        lines = (line.strip() for line in f)
        return [line for line in lines if line and not line.startswith("#")]


def read_pod5_map(path):
    """{sample: [pod5 path, ...]} from a two-column `sample<TAB>path` file.

    File order is kept, so a sample's POD5s list in the order they were written.
    """
    pod5_paths = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            sample, sep, pod5 = line.partition("\t")
            if not sep or not sample or not pod5:
                raise ValueError(f"{path}:{n}: expected SAMPLE<TAB>PATH, got {line!r}")
            pod5_paths.setdefault(sample, []).append(pod5)
    return pod5_paths


def generate_squiggy_session(
    sample_names,
    output_dir,
    fasta_path,
    session_name=None,
    compute_checksums=True,
    pod5_paths=None,
):
    """
    Generate squiggy session JSON for loading pipeline outputs.

    Args:
        sample_names: List of sample names
        output_dir: Pipeline output directory
        fasta_path: Path to reference FASTA file
        session_name: Optional session name (defaults to directory name)
        compute_checksums: Whether to compute MD5 checksums for files
        pod5_paths: {sample: [pod5 file, ...]}. The pipeline no longer writes a
            POD5 per sample -- the signal stays in the raw run -- so the caller
            says where each sample's files are. A sample not in the mapping falls
            back to the pre-v0.7 layout, pod5/{sample}/{sample}.pod5.

    Returns:
        dict: Session data structure
    """
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()

    if session_name is None:
        session_name = f"aa-tRNA-seq: {os.path.basename(os.path.abspath(output_dir))}"

    samples = {}
    file_checksums = {}

    fasta_abs = os.path.abspath(fasta_path)

    # Compute fasta checksum once (shared across samples)
    if compute_checksums and os.path.exists(fasta_path):
        file_checksums[fasta_abs] = get_file_info(fasta_path)

    for sample in sample_names:
        # Build absolute paths to pipeline outputs
        sample_pod5s = (pod5_paths or {}).get(sample) or [
            os.path.join(output_dir, "pod5", sample, f"{sample}.pod5")
        ]
        sample_pod5s = [os.path.abspath(p) for p in sample_pod5s]
        bam_path = os.path.abspath(
            os.path.join(output_dir, "bam", "final", sample, f"{sample}.bam")
        )

        samples[sample] = {
            "pod5Paths": sample_pod5s,
            "bamPath": bam_path,
            "fastaPath": fasta_abs,
        }

        # Compute checksums
        if compute_checksums:
            for pod5_path in sample_pod5s:
                if os.path.exists(pod5_path):
                    file_checksums[pod5_path] = get_file_info(pod5_path)
            if os.path.exists(bam_path):
                file_checksums[bam_path] = get_file_info(bam_path)

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
            "expandedSamples": sample_names,
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
        default=[],
        help=(
            "Sample names to include in session. Fine for a handful; a large "
            "run overflows the OS argument limit, so use --samples-file there."
        ),
    )
    parser.add_argument(
        "--samples-file",
        help="File of sample names, one per line (added to any --samples)",
    )
    parser.add_argument(
        "--pod5-map",
        help=(
            "Two-column TSV, SAMPLE<TAB>PATH, one row per POD5 file (added to "
            "any --pod5). Use this instead of --pod5 on a large run."
        ),
    )
    parser.add_argument(
        "--pod5",
        action="append",
        default=[],
        metavar="SAMPLE=PATH",
        help=(
            "A POD5 file belonging to SAMPLE. Repeatable; a sample with several "
            "files gets several. Samples given none fall back to "
            "pod5/{sample}/{sample}.pod5 under --output-dir."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Pipeline output directory",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to reference FASTA file",
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

    sample_names = list(args.samples)
    if args.samples_file:
        sample_names += read_samples_file(args.samples_file)
    if not sample_names:
        parser.error("no samples: give --samples and/or --samples-file")

    pod5_paths = {}
    if args.pod5_map:
        try:
            pod5_paths = read_pod5_map(args.pod5_map)
        except ValueError as e:
            parser.error(str(e))
    for item in args.pod5:
        sample, sep, path = item.partition("=")
        if not sep or not sample or not path:
            parser.error(f"--pod5 expects SAMPLE=PATH, got {item!r}")
        pod5_paths.setdefault(sample, []).append(path)

    session = generate_squiggy_session(
        sample_names=sample_names,
        output_dir=args.output_dir,
        fasta_path=args.fasta,
        session_name=args.session_name,
        compute_checksums=not args.no_checksums,
        pod5_paths=pod5_paths,
    )

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    with open(args.output, "w") as f:
        json.dump(session, f, indent=2)

    print(f"Squiggy session written to: {args.output}")


if __name__ == "__main__":
    main()
