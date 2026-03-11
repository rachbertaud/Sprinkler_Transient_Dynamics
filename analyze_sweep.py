"""
analyze_sweep.py - Reads sweep_results.txt and reports the best case
                   (lowest error) for each error metric.
"""

import re
import os

results_path = os.path.join(os.path.dirname(__file__), 'results_sweep.txt')

with open(results_path, 'r') as f:
    text = f.read()

num = r'([\d.eE+\-]+)'

case_pattern = re.compile(
    r'fit_switch=(\d),\s+proc_data_switch=(\d).*?'
    r'Error of analytical case with noise:\s+'     + num + r'.*?'
    r'Error of analytical case with no noise:\s+'  + num + r'.*?'
    r'Error of torque int:\s+'                     + num + r'.*?'
    r'Error of signal with noise:\s+'              + num + r'.*?'
    r'Error of signal with no noise:\s+'           + num + r'.*?'
    r'Error of forward case with noise:\s+'        + num + r'.*?'
    r'Error of forward case with no noise:\s+'     + num,
    re.DOTALL
)

cases = []
for m in case_pattern.finditer(text):
    cases.append({
        'fit_switch':           int(m.group(1)),
        'proc_data_switch':     int(m.group(2)),
        'err_an_noise':         float(m.group(3)),
        'err_an_no_noise':      float(m.group(4)),
        'err_torque_int':       float(m.group(5)),
        'err_sig_noise':        float(m.group(6)),
        'err_sig_no_noise':     float(m.group(7)),
        'err_fwd_noise':        float(m.group(8)),
        'err_fwd_no_noise':     float(m.group(9)),
    })

if not cases:
    print("No complete cases found in sweep_results.txt.")
    raise SystemExit

# ---- Summary table ----
cols = [
    ('err_an_noise',     'an/noise'),
    ('err_an_no_noise',  'an/clean'),
    ('err_torque_int',   'torq_int'),
    ('err_sig_noise',    'sig/noise'),
    ('err_sig_no_noise', 'sig/clean'),
    ('err_fwd_noise',    'fwd/noise'),
    ('err_fwd_no_noise', 'fwd/clean'),
]

lines = []

header = f"{'fit':>4} {'proc':>4} | " + " ".join(f"{lbl:>10}" for _, lbl in cols)
lines.append(header)
lines.append("-" * len(header))
for c in cases:
    row = f"{c['fit_switch']:>4} {c['proc_data_switch']:>4} | "
    row += " ".join(f"{c[key]:>10.6f}" for key, _ in cols)
    lines.append(row)

lines.append("")

# ---- Best case per metric ----
metrics = [
    ('err_an_noise',     'Analytical  (with noise)'),
    ('err_an_no_noise',  'Analytical  (clean)     '),
    ('err_torque_int',   'Torque integral         '),
    ('err_sig_noise',    'Signal      (with noise)'),
    ('err_sig_no_noise', 'Signal      (clean)     '),
    ('err_fwd_noise',    'Forward     (with noise)'),
    ('err_fwd_no_noise', 'Forward     (clean)     '),
]

lines.append("Best case per metric:")
lines.append("-" * 55)
for key, label in metrics:
    best = min(cases, key=lambda c, k=key: c[k])
    lines.append(f"  {label}  fit={best['fit_switch']}  proc={best['proc_data_switch']}  ->  {best[key]:.6f}")

# Print to console and write to file
output = "\n".join(lines)
print(output)

out_path = os.path.join(os.path.dirname(__file__), 'analysis_results.txt')
with open(out_path, 'w') as f:
    f.write(output + "\n")

print(f"\nResults written to: {out_path}")
