"""
summarize_errors.py  –  Aggregate errors from all Results_N_*_Gam_*_Omega_*.csv files.

Usage:
    python summarize_errors.py              # searches current directory
    python summarize_errors.py /path/to/dir

Outputs:
    error_summary.txt   – formatted, human-readable table
    error_summary.csv   – flat CSV for import into Excel / plotting tools
"""

import sys
import glob
from pathlib import Path
from process_results import parse_file   # reuse the parser we already have


# ── collect all CSV files ─────────────────────────────────────────────────────

def collect_files(search_dir: str) -> list[str]:
    pattern = str(Path(search_dir) / "Results_N_*_Gam_*_Omega_*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No matching CSV files found in '{search_dir}'")
    return files


def _to_pct(val):
    try:
        return float(val) * 100
    except (ValueError, TypeError):
        return val


# ── flatten parsed data into rows ─────────────────────────────────────────────

def extract_rows(data: dict) -> list[dict]:
    """One row per (run × case)."""
    rows = []
    for run_idx, run in enumerate(data["runs"], 1):
        sv = run["starting_values"]
        for case in run["cases"]:
            errs = case["errors"]
            rows.append({
                "N":            data["N"],
                "gamma":        data["gamma"],
                "omega":        data["omega"],
                "run":          run_idx,
                "dt":           run["dt"],
                "tau_mag":      sv.get("tau mag.", "?"),
                "t_target":     sv.get("t target", "?"),
                "fit_sw":       case["fit_switch"],
                "proc_sw":      case["proc_data_switch"],
                "err_omega":    _to_pct(errs.get("err omega",    "N/A")),
                "err_gamma":    _to_pct(errs.get("err gamma",    "N/A")),
                "err_integral": _to_pct(errs.get("err integral", "N/A")),
            })
    return rows


# ── text report ───────────────────────────────────────────────────────────────

# Column widths for the text table
_CW = {
    "N": 7, "gamma": 7, "omega": 7, "run": 4, "dt": 14,
    "tau_mag": 8, "t_target": 12,
    "fit_sw": 6, "proc_sw": 7,
    "err_omega": 13, "err_gamma": 13, "err_integral": 13,
}
_COLS = list(_CW.keys())
_HEADERS = {
    "N": "N", "gamma": "gamma", "omega": "omega", "run": "run",
    "dt": "dt", "tau_mag": "tau_mag", "t_target": "t_target",
    "fit_sw": "fit_sw", "proc_sw": "proc_sw",
    "err_omega": "err_omega (%)", "err_gamma": "err_gamma (%)", "err_integral": "err_integ (%)",
}

def _fmt(key, val) -> str:
    w = _CW[key]
    if key in ("err_omega", "err_gamma", "err_integral"):
        try:
            return f"{float(val):{w}.6g}"
        except (ValueError, TypeError):
            return f"{str(val):>{w}}"
    elif key == "dt":
        try:
            return f"{float(val):{w}.8g}"
        except (ValueError, TypeError):
            return f"{str(val):>{w}}"
    else:
        return f"{str(val):>{w}}"

def _divider(char="-") -> str:
    return char * (sum(_CW.values()) + 2 * len(_CW) + 1)

def _header_row() -> str:
    cells = [f"{_HEADERS[c]:>{_CW[c]}}" for c in _COLS]
    return "| " + "  ".join(cells) + " |"

def _data_row(row: dict) -> str:
    cells = [_fmt(c, row[c]) for c in _COLS]
    return "| " + "  ".join(cells) + " |"


def write_text_report(all_rows: list[dict], outpath: Path, n_files: int) -> None:
    # Group rows by (N, gamma, omega) for sectioned output
    from itertools import groupby

    groups = {}
    for r in all_rows:
        key = (r["N"], r["gamma"], r["omega"])
        groups.setdefault(key, []).append(r)

    with open(outpath, "w") as f:
        f.write("=" * len(_divider("=")) + "\n")
        f.write(f"  CONSOLIDATED ERROR SUMMARY  —  {n_files} file(s)\n")
        f.write(f"  Errors: err_omega | err_gamma | err_integral (torque)\n")
        f.write("=" * len(_divider("=")) + "\n\n")

        for (N, gamma, omega), rows in groups.items():
            f.write(f"  N = {N}   gamma = {gamma}   omega = {omega}\n")
            f.write(_divider("─") + "\n")
            f.write(_header_row() + "\n")
            f.write(_divider("─") + "\n")

            prev_run = None
            for row in rows:
                if prev_run is not None and row["run"] != prev_run:
                    f.write(_divider("·") + "\n")   # light separator between dt runs
                f.write(_data_row(row) + "\n")
                prev_run = row["run"]

            f.write(_divider("─") + "\n\n")

        f.write(f"  Total rows: {len(all_rows)}\n")

    print(f"  Text report → {outpath}")


# ── CSV export ────────────────────────────────────────────────────────────────

def write_csv(all_rows: list[dict], outpath: Path) -> None:
    import csv
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_COLS)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"  CSV export  → {outpath}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    search_dir = sys.argv[1] if len(sys.argv) > 1 else "."

    files = collect_files(search_dir)
    print(f"\nFound {len(files)} file(s):")
    for f in files:
        print(f"  {f}")
    print()

    all_rows = []
    for fp in files:
        data = parse_file(fp)
        all_rows.extend(extract_rows(data))

    out_dir = Path(search_dir)
    write_text_report(all_rows, out_dir / "error_summary.txt", len(files))
    write_csv(all_rows,         out_dir / "error_summary.csv")

    print(f"\nDone — {len(all_rows)} rows written.")


if __name__ == "__main__":
    main()
