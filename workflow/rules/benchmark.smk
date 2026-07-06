"""
Migration result-parity benchmarks (dorado 2.0.1 / escpod / leech).

These rules are NOT part of `rule all` (none of their outputs are in
pipeline_outputs()), so they never run with a normal pipeline invocation. Run
them by naming the target explicitly, e.g.:

    pixi run snakemake benchmark_fingerprint --configfile config/config-test.yml --cores 4
    pixi run snakemake benchmark_pod5_equivalence_all --configfile config/config-test.yml --cores 4
    pixi run snakemake benchmark_classifier_equivalence_all --profile cluster/slurm

See benchmark/README.md for the fingerprint/compare workflow across two refs.
All outputs land under {outdir}/benchmark/.
"""

BENCH_DIR = os.path.join(outdir, "benchmark")
BENCH_SRC = os.path.join(PIPELINE_DIR, "benchmark")


def _bench_cfg(key, default=None):
    """Read benchmark.<key> from config, honoring a flat benchmark_<key> override
    (handy for `--config benchmark_baseline=...` on the CLI)."""
    flat = config.get(f"benchmark_{key}")
    if flat is not None:
        return flat
    return config.get("benchmark", {}).get(key, default)


rule benchmark_fingerprint:
    """
    Snapshot this run's result-bearing outputs into a comparable fingerprint.
    Depends on the full pipeline, then summarizes it. Run once per ref/config;
    diff two fingerprints with `benchmark_compare` (or benchmark/compare.py).
    """
    input:
        pipeline_outputs(),
    output:
        fp=os.path.join(BENCH_DIR, "fingerprint", "fingerprint.json"),
    log:
        os.path.join(outdir, "logs", "benchmark", "fingerprint"),
    params:
        src=BENCH_SRC,
        outdir=outdir,
        fp_dir=os.path.join(BENCH_DIR, "fingerprint"),
        label=_bench_cfg("label", "run"),
        # precompute the whole flag so an empty ref contributes nothing
        gitref_arg=(
            f"--git-ref {_bench_cfg('git_ref')}" if _bench_cfg("git_ref") else ""
        ),
    shell:
        """
        python {params.src}/fingerprint.py \
            {params.outdir} {params.fp_dir} \
            --label {params.label} {params.gitref_arg} \
            >{log} 2>&1
        """


rule benchmark_compare:
    """
    Diff two fingerprints under a tolerance profile (strict | aggregate).
    Provide the two fingerprint dirs via config:
        benchmark:
          baseline: /path/to/old/fingerprint
          candidate: /path/to/new/fingerprint
          profile: aggregate
    or on the CLI: --config benchmark_baseline=... benchmark_candidate=... benchmark_profile=aggregate

    """
    output:
        report=os.path.join(BENCH_DIR, "compare", "report.txt"),
        json=os.path.join(BENCH_DIR, "compare", "report.json"),
    log:
        os.path.join(outdir, "logs", "benchmark", "compare"),
    params:
        src=BENCH_SRC,
        baseline=lambda w: _bench_cfg("baseline")
        or sys.exit(
            "benchmark_compare: set benchmark.baseline (fingerprint dir) in config"
        ),
        candidate=lambda w: _bench_cfg("candidate")
        or sys.exit(
            "benchmark_compare: set benchmark.candidate (fingerprint dir) in config"
        ),
        profile=_bench_cfg("profile", "aggregate"),
    shell:
        """
        python {params.src}/compare.py \
            {params.baseline} {params.candidate} \
            --profile {params.profile} \
            --json {output.json} \
            --report-only \
            >{output.report} 2>{log}
        """


rule benchmark_pod5_equivalence:
    """
    escpod vs pod5 merge losslessness for one sample: the pipeline's escpod-merged
    POD5 (merge_pods output) is compared read-for-read against a fresh `pod5 merge`
    of the same raw inputs. Uses the benchmark pixi env (ships the pod5 CLI + the
    pod5 python oracle).
    """
    input:
        escpod_pod5=os.path.join(outdir, "pod5", "{sample}", "{sample}.pod5"),
        raw=get_raw_inputs,
    output:
        report=os.path.join(BENCH_DIR, "pod5_equivalence", "{sample}.txt"),
    log:
        os.path.join(outdir, "logs", "benchmark", "pod5_equivalence", "{sample}"),
    params:
        src=BENCH_SRC,
        ref_pod5=os.path.join(BENCH_DIR, "pod5_equivalence", "{sample}.pod5ref.pod5"),
    shell:
        # `set -e` fails the rule on a real tool error (e.g. pod5 merge crash);
        # the trailing `|| true` on the comparator keeps the report (which states
        # PASS/FAIL) even when the two POD5s diverge.
        """
        pixi run -e benchmark bash -c '
            set -e
            rm -f {params.ref_pod5}
            pod5 merge {input.raw} -o {params.ref_pod5}
            python {params.src}/equivalence/pod5_equivalence.py compare \
                --a {input.escpod_pod5} --b {params.ref_pod5} || true
        ' >{output.report} 2>{log}
        rm -f {params.ref_pod5}
        """


rule benchmark_pod5_equivalence_all:
    """Run pod5 equivalence for every sample."""
    input:
        expand(
            os.path.join(BENCH_DIR, "pod5_equivalence", "{sample}.txt"),
            sample=samples.keys(),
        ),


rule benchmark_classifier_equivalence:
    """
    leech vs remora charging equivalence for one sample. Runs both engines on the
    same tagged BAM + POD5 (same cca_classifier.pt) and compares per-read scores.
    GPU rule (leech): give it GPU resources like classify_charging.
    """
    input:
        pod5=get_classification_pod5,
        bam=os.path.join(outdir, "bam", "tagged", "{sample}", "{sample}.tagged.bam"),
    output:
        report=os.path.join(BENCH_DIR, "classifier_equivalence", "{sample}.txt"),
    log:
        os.path.join(outdir, "logs", "benchmark", "classifier_equivalence", "{sample}"),
    threads: 4
    params:
        src=BENCH_SRC,
        model=config["remora_cca_classifier"],
        workdir=os.path.join(BENCH_DIR, "classifier_equivalence", "{sample}.work"),
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            export CUDA_VISIBLE_DEVICES
        fi

        # `|| true` keeps the report (RESULT: PASS/FAIL) when the engines diverge;
        # a hard error from leech/remora still surfaces in {log}.
        python {params.src}/equivalence/classifier_equivalence.py run \
            --pod5 {input.pod5} \
            --bam {input.bam} \
            --model {params.model} \
            --workdir {params.workdir} \
            >{output.report} 2>{log} || true

        rm -rf {params.workdir}
        """


rule benchmark_classifier_equivalence_all:
    """Run classifier equivalence for every sample."""
    input:
        expand(
            os.path.join(BENCH_DIR, "classifier_equivalence", "{sample}.txt"),
            sample=samples.keys(),
        ),


rule benchmark_isolation_all:
    """Both isolation checks (escpod-vs-pod5 + leech-vs-remora) for all samples."""
    input:
        rules.benchmark_pod5_equivalence_all.input,
        rules.benchmark_classifier_equivalence_all.input,
