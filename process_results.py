"""
process_results.py  –  Parse & plot simulation results for a given N.

Usage:
    python process_results.py <N>
    python process_results.py 2048

Outputs (written to  output_N<N>/):
    errors_N<N>.txt          – tabular error summary (omega, gamma, torque integral)
    phi_run_<k>_dt_<dt>.png  – phi_gen vs exact phi  (one per dt run)
    tau_run_<k>_dt_<dt>.png  – torque vs exact tau   (one per dt run)
"""

import sys
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ═══════════════════════════════════════════════════════════════════════════════
#  PARSER
# ═══════════════════════════════════════════════════════════════════════════════

def _floats(line: str) -> list[float]:
    return [float(x) for x in line.strip().split(",")[1:]]


def parse_file(filepath: str) -> dict:
    """
    Returns a dict with:
        N            – int
        gamma, omega – str
        exact_t      – list[float]
        exact_phi    – list[float]
        exact_tau    – list[float]
        runs         – list of dicts:
            dt             – float
            starting_values – dict
            t              – list[float]
            cases          – list of dicts:
                fit_switch, proc_data_switch – int
                generated_phi, torque        – list[float]
                errors                       – dict[str, float]
    """
    fp = Path(filepath)

    # Extract N, gamma, omega from filename  e.g. Results_N_2048_Gam_0.345_Omega_7.650.csv
    stem = fp.stem  # Results_N_2048_Gam_0.345_Omega_7.650
    parts = stem.split("_")
    try:
        N     = int(parts[parts.index("N") + 1])
        gamma = parts[parts.index("Gam") + 1]
        omega = parts[parts.index("Omega") + 1]
    except (ValueError, IndexError):
        N, gamma, omega = "?", "?", "?"

    with open(fp) as f:
        lines = f.readlines()

    result = {
        "N": N, "gamma": gamma, "omega": omega,
        "exact_t": [], "exact_phi": [], "exact_tau": [],
        "runs": [],
    }

    i = 0
    # ── global header ──────────────────────────────────────────────────────────
    while i < len(lines):
        raw = lines[i].strip()
        if raw.startswith("exact t,"):
            result["exact_t"] = _floats(raw)
        elif raw.startswith("exact phi,"):
            result["exact_phi"] = _floats(raw)
        elif raw.startswith("exact tau,"):
            result["exact_tau"] = _floats(raw)
        elif raw == "STARTING VALUES":
            break
        i += 1

    # ── run sections ───────────────────────────────────────────────────────────
    while i < len(lines):
        raw = lines[i].strip()

        if raw != "STARTING VALUES":
            i += 1
            continue

        sv_keys = lines[i + 1].strip().split(",")
        sv_vals = lines[i + 2].strip().split(",")
        sv = dict(zip(sv_keys, sv_vals))

        run: dict = {
            "dt": float(sv.get("dt", 0)),
            "starting_values": sv,
            "t": [],
            "cases": [],
        }
        i += 3

        while i < len(lines):
            raw = lines[i].strip()
            if raw == "STARTING VALUES":
                break

            if raw.startswith("t,") and raw.split(",")[0] == "t":
                run["t"] = _floats(raw)

            elif raw.startswith("fit switch,"):
                p = raw.split(",")
                case: dict = {
                    "fit_switch":       int(p[1]),
                    "proc_data_switch": int(p[3]),
                    "generated_phi":    [],
                    "torque":           [],
                    "errors":           {},
                }

                j = i + 1
                while j < len(lines):
                    jraw = lines[j].strip()
                    if jraw.startswith("torque,"):
                        case["torque"] = _floats(jraw)
                    elif jraw.startswith("generated phi,"):
                        case["generated_phi"] = _floats(jraw)
                    elif jraw.startswith("err omega,"):
                        err_keys = jraw.split(",")
                        err_vals = lines[j + 1].strip().split(",")
                        case["errors"] = dict(zip(err_keys, err_vals))
                        j += 2
                        break
                    j += 1

                run["cases"].append(case)
                i = j
                continue

            i += 1

        result["runs"].append(run)

    return result


# ═══════════════════════════════════════════════════════════════════════════════
#  TEXT REPORT
# ═══════════════════════════════════════════════════════════════════════════════

def write_error_report(data: dict, outpath: Path) -> None:
    N, gamma, omega = data["N"], data["gamma"], data["omega"]
    runs = data["runs"]

    col_w = 13

    def fmt(v):
        try:
            return f"{float(v) * 100:{col_w}.6f}"
        except (ValueError, TypeError):
            return f"{str(v):>{col_w}}"

    with open(outpath, "w") as f:
        f.write("=" * 76 + "\n")
        f.write(f"  ERROR SUMMARY  —  N={N}   gamma={gamma}   omega={omega}\n")
        f.write(f"  Columns: err_omega (%) | err_gamma (%) | err_integral (%) torque\n")
        f.write("=" * 76 + "\n\n")

        for k, run in enumerate(runs, 1):
            sv = run["starting_values"]
            f.write(f"Run {k}\n")
            f.write("-" * 76 + "\n")
            f.write(
                f"  dt = {run['dt']:.8g}   |   "
                f"tau_mag = {sv.get('tau mag.','?')}   |   "
                f"tau_start = {sv.get('tau start time','?')}   |   "
                f"tau_end = {sv.get('tau end time','?')}   |   "
                f"t_target = {sv.get('t target','?')}\n\n"
            )
            f.write(
                f"  {'fit_sw':>6}  {'proc_sw':>7}  "
                f"{'err_omega (%)':>{col_w}}  {'err_gamma (%)':>{col_w}}  {'err_integ (%)':>{col_w}}\n"
            )
            f.write(f"  {'-'*6}  {'-'*7}  {'-'*col_w}  {'-'*col_w}  {'-'*col_w}\n")

            for case in run["cases"]:
                errs = case["errors"]
                f.write(
                    f"  {case['fit_switch']:>6}  {case['proc_data_switch']:>7}  "
                    f"{fmt(errs.get('err omega','N/A'))}  "
                    f"{fmt(errs.get('err gamma','N/A'))}  "
                    f"{fmt(errs.get('err integral','N/A'))}\n"
                )
            f.write("\n")

        total_cases = sum(len(r["cases"]) for r in runs)
        f.write("=" * 76 + "\n")
        f.write(f"  Total runs: {len(runs)}   |   Total cases: {total_cases}\n")

    print(f"  Report  → {outpath}")


# ═══════════════════════════════════════════════════════════════════════════════
#  PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

# Presentation palette — ref line is thick black; cases use distinct colours
REF_STYLE = dict(color="black",    lw=2.4, ls="-",           zorder=1, label="Exact (reference)")
CASE_STYLES = [
    dict(color="#2166ac", lw=1.8, ls="--",            zorder=2),   # fit=0 proc=0  blue
    dict(color="#d6604d", lw=1.8, ls="-.",            zorder=2),   # fit=0 proc=1  red-orange
    dict(color="#1a9850", lw=1.8, ls=":",             zorder=2),   # fit=1 proc=0  green
    dict(color="#9970ab", lw=2.0, ls=(0,(4,1,1,1)),  zorder=2),   # fit=1 proc=1  purple
]

def _apply_style(ax, xlabel, ylabel, title, note):
    ax.set_xlabel(xlabel, fontsize=13, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=13, labelpad=8)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.text(
        0.01, 0.01, note,
        transform=ax.transAxes,
        fontsize=9, color="#555555",
        va="bottom", ha="left",
    )
    ax.legend(
        loc="best", fontsize=10.5,
        framealpha=0.92, edgecolor="#cccccc",
        handlelength=2.8, handletextpad=0.6,
    )
    ax.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.5)
    ax.tick_params(labelsize=11)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="#aaaaaa")


def make_plots(data: dict, out_dir: Path) -> None:
    N, gamma, omega = data["N"], data["gamma"], data["omega"]
    exact_t   = data["exact_t"]
    exact_phi = data["exact_phi"]
    exact_tau = data["exact_tau"]

    for k, run in enumerate(data["runs"], 1):
        sv   = run["starting_values"]
        dt   = run["dt"]
        t    = run["t"]
        note = (f"N={N}   dt={dt:.6g}   γ={gamma}   ω={omega}   "
                f"τ_mag={sv.get('tau mag.','?')}   "
                f"t_target={sv.get('t target','?')}")

        # ── phi plot ─────────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(13, 6))
        fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.12)

        ax.plot(exact_t, exact_phi, **REF_STYLE)

        for idx, case in enumerate(run["cases"]):
            gp = case["generated_phi"]
            if not gp or not t:
                continue
            fs, ps = case["fit_switch"], case["proc_data_switch"]
            ax.plot(t, gp,
                    label=f"phi_gen   fit={fs}  proc={ps}",
                    **CASE_STYLES[idx % len(CASE_STYLES)])

        _apply_style(
            ax,
            xlabel="Time  (s)",
            ylabel="φ  (rad)",
            title=f"Run {k}  —  Generated φ vs Exact φ",
            note=note,
        )

        tag  = f"N{N}_Gam{gamma}_Omega{omega}"
        path = out_dir / f"phi_{tag}_run{k:02d}_dt{dt:.6g}.png"
        fig.savefig(path, dpi=600, bbox_inches="tight")
        plt.close(fig)
        print(f"  Plot    → {path}")

        # ── torque plot ───────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(13, 6))
        fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.12)

        ax.plot(exact_t, exact_tau, **REF_STYLE)

        for idx, case in enumerate(run["cases"]):
            tor = case["torque"]
            if not tor or not t:
                continue
            fs, ps = case["fit_switch"], case["proc_data_switch"]
            ax.plot(t, tor,
                    label=f"torque   fit={fs}  proc={ps}",
                    **CASE_STYLES[idx % len(CASE_STYLES)])

        _apply_style(
            ax,
            xlabel="Time  (s)",
            ylabel="Torque  (N·m)",
            title=f"Run {k}  —  Fitted Torque vs Exact τ",
            note=note,
        )

        path = out_dir / f"tau_{tag}_run{k:02d}_dt{dt:.6g}.png"
        fig.savefig(path, dpi=600, bbox_inches="tight")
        plt.close(fig)
        print(f"  Plot    → {path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

def find_files(N: int, gamma: str = "*", omega: str = "*", search_dir: str = ".") -> list[str]:
    gam_pat   = gamma if gamma != "*" else "*"
    omega_pat = omega if omega != "*" else "*"
    pattern   = str(Path(search_dir) / f"Results_N_{N}_Gam_{gam_pat}_Omega_{omega_pat}.csv")
    matches   = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No CSV files found matching: {pattern}")
    return matches


def _out_tag(data: dict) -> str:
    return f"N{data['N']}_Gam{data['gamma']}_Omega{data['omega']}"


def process_one(filepath: str) -> None:
    print(f"\nProcessing: {filepath}")
    data    = parse_file(filepath)
    tag     = _out_tag(data)
    out_dir = Path(f"output_{tag}")
    out_dir.mkdir(exist_ok=True)
    print(f"Output dir: {out_dir}/\n")

    write_error_report(data, out_dir / f"errors_{tag}.txt")
    make_plots(data, out_dir)

    n_plots = len(data["runs"]) * 2
    print(f"Done — 1 report + {n_plots} plots written to '{out_dir}/'")


def main():
    if len(sys.argv) < 2:
        print("Usage: python process_results.py <N> [gamma] [omega] [search_dir]")
        print("       python process_results.py 2048")
        print("       python process_results.py 2048 0.345")
        print("       python process_results.py 2048 0.345 7.650")
        sys.exit(1)

    N          = int(sys.argv[1])
    gamma      = sys.argv[2] if len(sys.argv) > 2 else "*"
    omega      = sys.argv[3] if len(sys.argv) > 3 else "*"
    search_dir = sys.argv[4] if len(sys.argv) > 4 else "."

    files = find_files(N, gamma, omega, search_dir)
    print(f"Found {len(files)} file(s) matching N={N}"
          + (f"  gamma={gamma}" if gamma != "*" else "")
          + (f"  omega={omega}" if omega != "*" else ""))

    for fp in files:
        process_one(fp)

    print(f"\nAll done — processed {len(files)} file(s).")


if __name__ == "__main__":
    main()
