"""Tests for benchmark aggregation helpers."""

from __future__ import annotations

from pathlib import Path

from seqnado.cli.benchmark_helpers import (
    compute_assay_output_sizes,
    discover_alignment_processing_logs,
    discover_benchmark_files,
    discover_snakemake_logs,
    format_bytes,
    format_seconds,
    format_compact_number,
    load_benchmark_table,
    parse_alignment_processing_logs,
    parse_snakemake_log_timeline,
    parse_snakemake_logs_timeline,
    summarize_benchmarks,
    write_html_report,
)
from seqnado.cli.commands.benchmark import (
    _infer_output_root,
    _infer_run_root,
    _resolve_benchmark_dir,
    _resolve_report_outputs,
)


def _write_benchmark(
    path: Path,
    runtime: float,
    cpu_time: float,
    max_rss: float,
    io_in: float = 0,
    io_out: float = 0,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "s\th:m:s\tmax_rss\tcpu_time\tio_in\tio_out\tmean_load\n"
        f"{runtime}\t0:00:00\t{max_rss}\t{cpu_time}\t{io_in}\t{io_out}\t0.50\n",
        encoding="utf-8",
    )


def test_discover_benchmark_files(tmp_path: Path) -> None:
    _write_benchmark(tmp_path / ".benchmark" / "align" / "sample_a.tsv", 10, 20, 100, 5, 7)
    _write_benchmark(tmp_path / ".benchmark" / "peaks" / "sample_b.tsv", 30, 40, 200, 11, 13)
    (tmp_path / "counts.tsv").write_text("sample\tcount\nA\t1\n", encoding="utf-8")

    files = discover_benchmark_files(tmp_path)

    assert [f.name for f in files] == ["sample_a.tsv", "sample_b.tsv"]


def test_load_benchmark_table_adds_metadata_and_sorts(tmp_path: Path) -> None:
    _write_benchmark(tmp_path / ".benchmark" / "align" / "sample_a.tsv", 10, 20, 100, 5, 7)
    _write_benchmark(tmp_path / ".benchmark" / "peaks" / "sample_b.tsv", 30, 40, 200, 11, 13)

    table = load_benchmark_table(tmp_path)

    assert list(table["benchmark_file"]) == ["peaks/sample_b.tsv", "align/sample_a.tsv"]
    assert list(table["group"]) == ["peaks", "align"]
    assert list(table["job"]) == ["sample_b", "sample_a"]


def test_load_benchmark_table_prefers_extended_rule_name(tmp_path: Path) -> None:
    benchmark_file = tmp_path / ".benchmark" / "qc" / "sample_a.tsv"
    benchmark_file.parent.mkdir(parents=True, exist_ok=True)
    benchmark_file.write_text(
        "\t".join(
            [
                "s",
                "max_rss",
                "cpu_time",
                "io_in",
                "io_out",
                "jobid",
                "rule_name",
                "threads",
                "cpu_usage",
                "input_size_mb",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "10",
                "100",
                "20",
                "5",
                "7",
                "1",
                "atac_fastqc_raw_paired",
                "4",
                "774.03",
                "0.18",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    table = load_benchmark_table(tmp_path)

    assert list(table["assay"]) == ["atac"]
    assert list(table["display_rule"]) == ["fastqc_raw_paired"]
    assert list(table["jobid"]) == [1]
    assert list(table["cpu_usage"]) == [774.03]
    assert list(table["input_size_mb"]) == [0.18]


def test_summarize_benchmarks(tmp_path: Path) -> None:
    _write_benchmark(tmp_path / ".benchmark" / "align" / "sample_a.tsv", 10, 20, 100, 5, 7)
    _write_benchmark(tmp_path / ".benchmark" / "peaks" / "sample_b.tsv", 30, 40, 200, 11, 13)

    summary = summarize_benchmarks(load_benchmark_table(tmp_path))

    assert summary.jobs == 2
    assert summary.groups == 2
    assert summary.total_runtime_seconds == 40
    assert summary.total_cpu_time_seconds == 60
    assert summary.peak_rss_mb == 200
    assert summary.total_io_in == 16
    assert summary.total_io_out == 20


def test_write_html_report(tmp_path: Path) -> None:
    _write_benchmark(tmp_path / ".benchmark" / "align" / "sample_a.tsv", 10, 20, 100, 121254.87, 1740695.35)
    _write_benchmark(tmp_path / ".benchmark" / "peaks" / "sample_b.tsv", 30, 40, 200, 11, 13)
    alignment_log = tmp_path / "seqnado_output" / "chip" / "logs" / "alignment_post_process" / "chip-rx_MLL.log"
    alignment_log.parent.mkdir(parents=True, exist_ok=True)
    alignment_log.write_text(
        "Step\tSample\tReads Before\tReads After\n"
        "Aligned (align_paired)\tchip-rx_MLL\t0\t29346\n"
        "Sort\tchip-rx_MLL\t29346\t29346\n",
        encoding="utf-8",
    )
    table = load_benchmark_table(tmp_path)
    output_html = tmp_path / "benchmark_report.html"

    write_html_report(
        table,
        benchmark_dir=tmp_path,
        output_file=output_html,
        top_n=5,
        read_counts_df=parse_alignment_processing_logs(tmp_path / "seqnado_output"),
    )

    html = output_html.read_text(encoding="utf-8")
    assert "SeqNado Benchmark Report" in html
    assert "Contents" in html
    assert "report-toc" in html
    assert "href='#run-timeline'" in html
    assert "<details class='report-section' id='run-timeline'>" in html
    assert "Run Timeline" in html
    assert "Runtime By Rule Box Plot" in html
    assert "Max RSS By Rule Box Plot" in html
    assert "Read Counts By Rule Box Plot" in html
    assert "href='#read-counts'" in html
    assert "<details class='report-section' id='read-counts'>" in html
    assert "href='#io'" in html
    assert "<details class='report-section' id='io'>" in html
    assert "I/O By Rule Box Plot" in html
    assert "input/output volume" in html
    assert "Output Size By Assay" in html
    assert "Samples Per Assay" in html
    assert "RNA-seq" in html
    assert "Report Provenance" in html
    assert "generated this report with SeqNado version" in html
    assert "Total IO In" in html
    assert "Total IO Out" in html
    assert "Export CSV" in html
    assert "All Assays" in html
    assert "gantt-legend" in html
    assert "gantt-bar" in html
    assert "Timeline of logged jobs parsed from Snakemake run logs." in html
    assert "maximum resident set size (peak memory use) in megabytes" in html
    assert "input/output read volume" in html
    assert "input/output write volume" in html
    assert "data:text/csv" in html
    assert "121.3k" in html
    assert "1.74M" in html
    assert "align" in html
    assert "All Benchmarks" not in html


def test_parse_alignment_processing_logs(tmp_path: Path) -> None:
    log_file = tmp_path / "seqnado_output" / "chip" / "logs" / "alignment_post_process" / "chip-rx_MLL.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text(
        "\n".join(
            [
                "Step\tSample\tReads Before\tReads After",
                "Aligned (align_paired)\tchip-rx_MLL\t0\t29346",
                "Sort\tchip-rx_MLL\t29346\t29346",
                "Blacklist\tchip-rx_MLL\t29346\t29346",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    logs = discover_alignment_processing_logs(tmp_path / "seqnado_output")
    assert logs == [log_file]

    parsed = parse_alignment_processing_logs(tmp_path / "seqnado_output")
    assert list(parsed["assay"]) == ["chip", "chip", "chip"]
    assert list(parsed["sample"]) == ["chip-rx_MLL", "chip-rx_MLL", "chip-rx_MLL"]
    assert list(parsed["display_rule"]) == [
        "Aligned (align_paired)",
        "Sort",
        "Blacklist",
    ]
    assert list(parsed["read_count"]) == [29346, 29346, 29346]


def test_write_html_report_includes_extended_benchmark_sections(tmp_path: Path) -> None:
    benchmark_file = tmp_path / ".benchmark" / "qc" / "sample_a.tsv"
    benchmark_file.parent.mkdir(parents=True, exist_ok=True)
    benchmark_file.write_text(
        "\t".join(
            [
                "s",
                "h:m:s",
                "max_rss",
                "max_vms",
                "max_uss",
                "max_pss",
                "io_in",
                "io_out",
                "mean_load",
                "cpu_time",
                "jobid",
                "rule_name",
                "wildcards",
                "params",
                "threads",
                "cpu_usage",
                "resources",
                "input_size_mb",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "6.13",
                "0:00:06",
                "327.00",
                "6225.16",
                "303.58",
                "306.89",
                "0.00",
                "4.63",
                "126.24",
                "9.02",
                "1",
                "atac_fastqc_raw_paired",
                "{'sample': 'atac'}",
                "{'extra': '--quiet'}",
                "4",
                "774.03",
                "{'_cores': 4}",
                "0.18",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    table = load_benchmark_table(tmp_path)
    output_html = tmp_path / "benchmark_report.html"

    write_html_report(table, benchmark_dir=tmp_path, output_file=output_html, top_n=5)

    html = output_html.read_text(encoding="utf-8")
    assert "Max USS" in html
    assert "Max PSS" in html
    assert "CPU Usage" in html
    assert "Input Size" in html
    assert "maximum unique set size in megabytes" in html
    assert "maximum proportional set size in megabytes" in html
    assert "CPU usage percentage" in html
    assert "input size in megabytes" in html


def test_format_seconds() -> None:
    assert format_seconds(59) == "0:00:59"
    assert format_seconds(3661) == "1:01:01"


def test_format_compact_number() -> None:
    assert format_compact_number(999) == "999"
    assert format_compact_number(121254.87) == "121.3k"
    assert format_compact_number(1740695.35) == "1.74M"


def test_format_bytes() -> None:
    assert format_bytes(512) == "512 B"
    assert format_bytes(2048) == "2.00 KiB"
    assert format_bytes(5 * 1024**3) == "5.00 GiB"


def test_compute_assay_output_sizes(tmp_path: Path) -> None:
    chip_dir = tmp_path / "seqnado_output" / "chip"
    rna_dir = tmp_path / "seqnado_output" / "rna"
    (chip_dir / ".benchmark").mkdir(parents=True)
    (rna_dir / ".benchmark").mkdir(parents=True)
    (chip_dir / "results.bigWig").write_bytes(b"a" * 1024)
    (rna_dir / "counts.tsv").write_bytes(b"a" * 2048)
    (rna_dir / ".benchmark" / "ignored.tsv").write_bytes(b"a" * 9999)

    sizes = compute_assay_output_sizes(tmp_path / "seqnado_output")

    assert list(sizes["assay"]) == ["rna", "chip"]
    assert list(sizes["output_size_bytes"]) == [2048, 1024]


def test_discover_snakemake_logs(tmp_path: Path) -> None:
    log_dir = tmp_path / ".snakemake" / "log"
    log_dir.mkdir(parents=True)
    (log_dir / "2026-01-01T000000.000000.snakemake.log").write_text("", encoding="utf-8")
    (log_dir / "2026-01-02T000000.000000.snakemake.log").write_text("", encoding="utf-8")

    logs = discover_snakemake_logs(tmp_path)

    assert [p.name for p in logs] == [
        "2026-01-01T000000.000000.snakemake.log",
        "2026-01-02T000000.000000.snakemake.log",
    ]


def test_parse_snakemake_log_timeline(tmp_path: Path) -> None:
    log_file = tmp_path / "run.snakemake.log"
    log_file.write_text(
        "\n".join(
            [
                "[Fri Apr 17 11:06:27 2026]",
                "Job 209: Trimming reads for sample rna-spikein-treated-rep2 using Trim Galore",
                "[Fri Apr 17 11:06:30 2026]",
                "Job 140: Running FastQC on raw FASTQ files for sample rna",
                "[Fri Apr 17 11:06:37 2026]",
                "Finished jobid: 209 (Rule: rna_trimgalore_paired)",
                "[Fri Apr 17 11:06:39 2026]",
                "Finished jobid: 140 (Rule: rna_fastqc_raw_paired)",
                "[Fri Apr 17 11:06:40 2026]",
                "Job 999: all",
                "[Fri Apr 17 11:06:41 2026]",
                "Finished jobid: 999 (Rule: all)",
            ]
        ),
        encoding="utf-8",
    )

    timeline = parse_snakemake_log_timeline(log_file)

    assert list(timeline["jobid"]) == [209, 140]
    assert list(timeline["rule"]) == ["rna_trimgalore_paired", "rna_fastqc_raw_paired"]
    assert list(timeline["display_rule"]) == ["trimgalore_paired", "fastqc_raw_paired"]
    assert list(timeline["assay"]) == ["rna", "rna"]
    assert list(timeline["entity"]) == ["rna-spikein-treated-rep2", "rna"]
    assert timeline.loc[0, "label"].startswith("Trimming reads for sample")
    assert timeline.loc[0, "endtime"] > timeline.loc[0, "starttime"]


def test_parse_snakemake_logs_timeline_prefers_most_recent_duplicate(tmp_path: Path) -> None:
    log_file_1 = tmp_path / "2026-04-17T110000.000000.snakemake.log"
    log_file_2 = tmp_path / "2026-04-17T120000.000000.snakemake.log"

    log_file_1.write_text(
        "\n".join(
            [
                "[Fri Apr 17 11:00:00 2026]",
                "Job 1: Running FastQC on raw FASTQ files for sample rna",
                "[Fri Apr 17 11:00:05 2026]",
                "Finished jobid: 1 (Rule: rna_fastqc_raw_paired)",
            ]
        ),
        encoding="utf-8",
    )
    log_file_2.write_text(
        "\n".join(
            [
                "[Fri Apr 17 12:00:00 2026]",
                "Job 2: Running FastQC on raw FASTQ files for sample rna",
                "[Fri Apr 17 12:00:08 2026]",
                "Finished jobid: 2 (Rule: rna_fastqc_raw_paired)",
                "[Fri Apr 17 12:00:10 2026]",
                "Job 3: Trimming reads for sample chip-rx_MLL using Trim Galore",
                "[Fri Apr 17 12:00:20 2026]",
                "Finished jobid: 3 (Rule: chip_trimgalore_paired)",
            ]
        ),
        encoding="utf-8",
    )

    timeline = parse_snakemake_logs_timeline([log_file_1, log_file_2])

    assert list(timeline["rule"]) == ["rna_fastqc_raw_paired", "chip_trimgalore_paired"]
    assert list(timeline["entity"]) == ["rna", "chip-rx_MLL"]
    assert list(timeline["jobid"]) == [2, 3]


def test_resolve_benchmark_dir_falls_back_to_seqnado_output_root_for_multiomics(tmp_path: Path, monkeypatch) -> None:
    (tmp_path / "seqnado_output" / "chip" / ".benchmark").mkdir(parents=True)
    (tmp_path / "seqnado_output" / "rna" / ".benchmark").mkdir(parents=True)
    monkeypatch.chdir(tmp_path)

    resolved = _resolve_benchmark_dir(Path(".benchmark"))

    assert resolved == Path("seqnado_output")


def test_infer_output_root() -> None:
    assert _infer_output_root(Path("seqnado_output")) == Path("seqnado_output")
    assert _infer_output_root(Path("seqnado_output/chip/.benchmark")) == Path("seqnado_output")
    assert _infer_output_root(Path(".benchmark")) is None


def test_infer_run_root() -> None:
    assert _infer_run_root(Path("seqnado_output")) == Path(".")
    assert _infer_run_root(Path("seqnado_output/chip/.benchmark")) == Path(".")
    assert _infer_run_root(Path("/tmp/run/seqnado_output/chip/.benchmark")) == Path("/tmp/run")
    assert _infer_run_root(Path(".benchmark")) == Path(".")


def test_resolve_report_outputs_defaults_to_seqnado_output() -> None:
    html_path, tsv_path = _resolve_report_outputs(
        Path("seqnado_output"),
        Path("benchmark_report.html"),
        Path("benchmark_summary.tsv"),
    )

    assert html_path == Path("seqnado_output/seqnado_benchmark.html")
    assert tsv_path == Path("seqnado_output/seqnado_benchmark.tsv")


def test_resolve_report_outputs_preserves_explicit_paths() -> None:
    html_path, tsv_path = _resolve_report_outputs(
        Path("seqnado_output"),
        Path("custom/report.html"),
        Path("custom/report.tsv"),
    )

    assert html_path == Path("custom/report.html")
    assert tsv_path == Path("custom/report.tsv")
