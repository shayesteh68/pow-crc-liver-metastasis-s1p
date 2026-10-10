from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
FIG_DIR = PROJECT_ROOT / "figures"
PLOTS_DIR = RESULTS_DIR / "plots"

for d in (DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR):
    d.mkdir(parents=True, exist_ok=True)
