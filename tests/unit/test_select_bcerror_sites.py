"""Unit tests for select_bcerror_sites.py."""

import pandas as pd
import pytest

from select_bcerror_sites import select_sites


def write_bcerror(path, rows):
    """Write a minimal bcerror table with the columns the selector reads."""
    pd.DataFrame(
        rows,
        columns=["Reference", "Position", "Spanning_Reads", "BCErrorFreq"],
    ).to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def one_sample(temp_dir):
    return write_bcerror(
        temp_dir / "a.bcerror.tsv",
        [
            ("ref1", 1, 100, 0.02),  # below the error threshold
            ("ref1", 2, 100, 0.30),
            ("ref1", 3, 5, 0.90),  # below the coverage threshold
            ("ref2", 4, 100, 0.15),
        ],
    )


class TestSelectSites:
    def test_applies_both_thresholds(self, one_sample):
        sites = select_sites([one_sample], min_error=0.1, min_cov=20)

        assert list(zip(sites["ref"], sites["pos"])) == [("ref1", 2), ("ref2", 4)]

    def test_output_columns(self, one_sample):
        sites = select_sites([one_sample], min_error=0.1, min_cov=20)

        assert list(sites.columns) == [
            "ref",
            "pos",
            "n_samples",
            "max_error",
            "mean_error",
        ]

    def test_unions_across_samples(self, temp_dir, one_sample):
        other = write_bcerror(
            temp_dir / "b.bcerror.tsv",
            [
                ("ref1", 2, 100, 0.40),
                ("ref1", 9, 100, 0.50),  # only this sample clears it
            ],
        )
        sites = select_sites([one_sample, other], min_error=0.1, min_cov=20)

        assert list(zip(sites["ref"], sites["pos"])) == [
            ("ref1", 2),
            ("ref1", 9),
            ("ref2", 4),
        ]

    def test_summarizes_error_across_samples(self, temp_dir, one_sample):
        other = write_bcerror(
            temp_dir / "b.bcerror.tsv",
            [("ref1", 2, 100, 0.40)],
        )
        sites = select_sites([one_sample, other], min_error=0.1, min_cov=20)
        at_two = sites[(sites["ref"] == "ref1") & (sites["pos"] == 2)].iloc[0]

        assert at_two["n_samples"] == 2
        assert at_two["max_error"] == pytest.approx(0.40)
        assert at_two["mean_error"] == pytest.approx(0.35)

    def test_min_samples_requires_reproducibility(self, temp_dir, one_sample):
        other = write_bcerror(
            temp_dir / "b.bcerror.tsv",
            [("ref1", 2, 100, 0.40)],
        )
        sites = select_sites(
            [one_sample, other], min_error=0.1, min_cov=20, min_samples=2
        )

        assert list(zip(sites["ref"], sites["pos"])) == [("ref1", 2)]

    def test_sorted_by_reference_and_position(self, temp_dir):
        path = write_bcerror(
            temp_dir / "a.bcerror.tsv",
            [
                ("ref2", 9, 100, 0.5),
                ("ref1", 7, 100, 0.5),
                ("ref1", 3, 100, 0.5),
            ],
        )
        sites = select_sites([path], min_error=0.1, min_cov=20)

        assert list(zip(sites["ref"], sites["pos"])) == [
            ("ref1", 3),
            ("ref1", 7),
            ("ref2", 9),
        ]

    def test_no_passing_sites_returns_empty_frame(self, one_sample):
        sites = select_sites([one_sample], min_error=0.99, min_cov=20)

        assert len(sites) == 0
        assert list(sites.columns) == [
            "ref",
            "pos",
            "n_samples",
            "max_error",
            "mean_error",
        ]
