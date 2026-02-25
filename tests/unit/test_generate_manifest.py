"""Unit tests for generate_manifest.py."""

import pytest

from generate_manifest import (
    extract_config_params,
    extract_sample_info,
    parse_pixi_lock_version,
)


class TestParsePixiLockVersion:
    def test_finds_package(self, temp_dir):
        lock = temp_dir / "pixi.lock"
        lock.write_text(
            "- conda: https://conda.anaconda.org/bioconda/linux-64/samtools-1.23-h96c455f_0.conda\n"
            "- conda: https://conda.anaconda.org/bioconda/linux-64/bwa-0.7.18-he4a0461_1.conda\n"
        )
        assert parse_pixi_lock_version(str(lock), "samtools") == "1.23"
        assert parse_pixi_lock_version(str(lock), "bwa") == "0.7.18"

    def test_missing_package(self, temp_dir):
        lock = temp_dir / "pixi.lock"
        lock.write_text(
            "- conda: https://conda.anaconda.org/bioconda/linux-64/samtools-1.23-h96c455f_0.conda\n"
        )
        assert parse_pixi_lock_version(str(lock), "nonexistent") is None

    def test_missing_file(self):
        assert parse_pixi_lock_version("/nonexistent/pixi.lock", "samtools") is None


class TestExtractConfigParams:
    def test_extracts_known_keys(self):
        config = {
            "fasta": "/path/to/ref.fa",
            "base_calling_model": "sup",
            "dorado_model": "dna_r10.4.1",
            "opts": {"bwa": "-x ont2d"},
            "unknown_key": "should_be_ignored",
        }
        result = extract_config_params(config)
        assert result["fasta"] == "/path/to/ref.fa"
        assert result["base_calling_model"] == "sup"
        assert "opts" in result
        assert "unknown_key" not in result

    def test_ignores_unknown_keys(self):
        config = {"random": "value", "other": 123}
        result = extract_config_params(config)
        assert len(result) == 0


class TestExtractSampleInfo:
    def test_count_and_names(self):
        samples = {
            "sample1": {"path": {"/data/run1"}},
            "sample2": {"path": {"/data/run2"}},
        }
        result = extract_sample_info(samples)
        assert result["count"] == 2
        assert set(result["names"]) == {"sample1", "sample2"}

    def test_set_to_list_conversion(self):
        samples = {"s1": {"path": {"/data/run1"}}}
        result = extract_sample_info(samples)
        # Single path should be extracted (not a list)
        assert isinstance(result["input_paths"]["s1"], str)

    def test_multiple_paths(self):
        samples = {"s1": {"path": ["/data/run1", "/data/run2"]}}
        result = extract_sample_info(samples)
        assert isinstance(result["input_paths"]["s1"], list)
