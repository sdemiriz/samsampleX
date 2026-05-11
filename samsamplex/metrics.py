"""Quality metrics for comparing depth distributions."""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from typing import Callable, Sequence, TextIO

import numpy as np

from .depth import DepthArray

# Output header line (TSV)
STATS_HEADER = (
    "chrom\tstart\tend\tmean_depth_temp\tdepth_diff\trel_depth_diff\tmean_depth_res\n"
)


def depth_from_path(path: str, contig: str, start: int, end: int) -> DepthArray:
    """Load a DepthArray from *path*, auto-detecting format by extension.

    .bed files are read via bed_read_depths; everything else is treated as BAM.
    """
    if path.endswith(".bed"):
        from .bed import bed_read_depths

        return bed_read_depths(path, contig, start, end)
    from .depth import depth_from_bam

    return depth_from_bam(path, contig, start, end)


@dataclass
class WindowMetricsRow:
    chrom: str
    start: int
    end: int
    mean_depth_temp: float
    depth_diff: float
    rel_depth_diff: float
    mean_depth_res: float


def _mean_absolute_error(a: np.ndarray, b: np.ndarray) -> float:
    """Mean absolute per-base depth difference between aligned arrays a and b."""
    if len(a) == 0 and len(b) == 0:
        return 0.0
    assert len(a) > 0 and len(b) > 0, "Non-empty arrays are required for stats calculations"
    return float(np.mean(np.abs(a.astype(np.float64) - b.astype(np.float64))))


def _mean_signed_depth_difference(a: np.ndarray, b: np.ndarray) -> float:
    """Mean (b − a) per base; positive when b is deeper on average."""
    if len(a) == 0 and len(b) == 0:
        return 0.0
    assert len(a) > 0 and len(b) > 0, "Non-empty arrays are required for stats calculations"
    return float(np.mean(b.astype(np.float64) - a.astype(np.float64)))


def metrics_window_rows(
    depth_temp: DepthArray,
    depth_res: DepthArray,
    window: int,
) -> list[WindowMetricsRow]:
    """
    Split the shared region into non-overlapping windows and compute metrics per window.
    """
    if depth_temp.length != depth_res.length:
        raise ValueError(
            f"Depth arrays have different lengths ({depth_temp.length} vs {depth_res.length})"
        )
    if depth_temp.contig != depth_res.contig or depth_temp.start != depth_res.start:
        raise ValueError("Depth arrays must share contig and start")
    if window < 1:
        raise ValueError("window must be >= 1")

    t = depth_temp.depths.astype(np.float64)
    r = depth_res.depths.astype(np.float64)
    chrom = depth_temp.contig
    base = depth_temp.start
    n = depth_temp.length
    rows: list[WindowMetricsRow] = []

    for ws in range(0, n, window):
        we = min(ws + window, n)
        seg_t = t[ws:we]
        seg_r = r[ws:we]
        mean_temp = float(np.mean(seg_t))
        mean_res = float(np.mean(seg_r))
        depth_diff = float(np.mean(seg_r - seg_t))
        if mean_temp > 0.0:
            rel_depth_diff = depth_diff / mean_temp
        else:
            rel_depth_diff = float("nan")

        rows.append(
            WindowMetricsRow(
                chrom=chrom,
                start=base + ws,
                end=base + we,
                mean_depth_temp=mean_temp,
                depth_diff=depth_diff,
                rel_depth_diff=rel_depth_diff,
                mean_depth_res=mean_res,
            )
        )
    return rows


def top_k_window_indices_by_abs(
    rows: Sequence[WindowMetricsRow],
    value_getter: Callable[[WindowMetricsRow], float],
    *,
    k: int = 10,
) -> list[int]:
    """
    Indices of up to *k* windows with largest absolute value; ties broken by lower genomic start. Windows with mean == 0 values are excluded.
    """
    n = len(rows)
    if n == 0:
        return []

    vals = np.array([value_getter(rows[i]) for i in range(n)], dtype=np.float64)
    starts = np.array([rows[i].start for i in range(n)], dtype=np.int64)
    finite = np.isfinite(vals)
    if not np.any(finite):
        return []

    idx = np.nonzero(finite)[0]
    v = vals[finite]
    s = starts[finite]
    av = np.abs(v)
    order = np.lexsort((s, -av))
    idx_sorted = idx[order]
    take = min(k, len(idx_sorted))
    return idx_sorted[:take].tolist()


def any_zero_depth_window(rows: Sequence[WindowMetricsRow]) -> bool:
    """True if any window has mean template or mean result depth exactly zero."""
    for row in rows:
        if row.mean_depth_temp == 0.0 or row.mean_depth_res == 0.0:
            return True
    return False


def write_stats_tsv(fp: TextIO, rows: Sequence[WindowMetricsRow]) -> None:
    """Write header and one row per window."""
    fp.write(STATS_HEADER)
    for row in rows:
        fp.write(
            f"{row.chrom}\t{row.start}\t{row.end}\t"
            f"{row.mean_depth_temp:.6g}\t{row.depth_diff:.6g}\t{row.rel_depth_diff:.6g}\t"
            f"{row.mean_depth_res:.6g}\n"
        )


def _fmt_cmdline(argv: Sequence[str]) -> str:
    return " ".join(f"{a}" for a in argv)


def _fmt_stderr_metric(val: float, *, signed: bool) -> str:
    """Format a metric for stderr tail lines; signed metrics use explicit +/-."""
    if not math.isfinite(val):
        return "nan" if isinstance(val, float) and math.isnan(val) else str(val)
    if signed:
        return f"{val:+.6g}"
    return f"{val:.6g}"


def print_stats_stderr(
    argv: Sequence[str],
    rows: Sequence[WindowMetricsRow],
    *,
    warn_zero_depth: bool,
    stderr_rows: int = 10,
) -> None:
    """stderr: command line, optional WARNING, top windows by absolute value of a metric."""
    print(_fmt_cmdline(argv), file=sys.stderr)
    if warn_zero_depth:
        print(
            "WARNING: Zero mean template or result depth caused a NaN in at least one window.",
            file=sys.stderr,
        )

    n_w = len(rows)
    if n_w == 0:
        return

    sections: list[tuple[str, Callable[[WindowMetricsRow], float]]] = [
        ("depth_diff", lambda r: r.depth_diff),
        ("rel_depth_diff", lambda r: r.rel_depth_diff),
    ]

    for name, getter in sections:
        idxs = top_k_window_indices_by_abs(rows, getter, k=stderr_rows)
        print("", file=sys.stderr)
        print(
            f"Top {len(idxs)} windows by absolute {name} ({len(idxs)}/{n_w} windows)",
            file=sys.stderr,
        )
        print("=" * 80, file=sys.stderr)
        for i in idxs:
            row = rows[i]
            v = getter(row)
            print(
                f"{row.chrom}\t{row.start}\t{row.end}\t{_fmt_stderr_metric(v, signed=True)}",
                file=sys.stderr,
            )
