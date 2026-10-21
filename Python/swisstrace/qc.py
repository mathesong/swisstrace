"""Two-panel raw/corrected QC visualization."""
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure


def swisstrace_qc(res, xlim=None, col_raw="0.4", col_corr="firebrick"):
    """Return a matplotlib Figure (save with figure.savefig or display in a notebook)."""
    raw, tac = res["raw"], res["tac"]
    relative = raw.time - res["t0_seconds"]
    onset = res["lead"] if res["t0_detected"] else 0
    if xlim is None:
        xlim = (relative.min(), tac.time.max())
    fig = Figure(figsize=(1000 / 110, 750 / 110), dpi=110, layout="constrained")
    top, bottom = fig.subplots(2, 1, sharex=True)
    top.plot(relative, raw.coincidence, ".", markersize=3, color=col_raw)
    if np.isfinite(onset):
        top.axvspan(xlim[0], onset, color="steelblue", alpha=.12)
        top.axvline(onset, color="darkgreen", linestyle=":", label="onset (t0 + lead)")
    top.axhline(res["background"], color="steelblue", linestyle="--",
                label=f'background = {res["background"]:.1f} (n = {res["n_background"]})')
    top.axvline(0, color="black", label="time 0")
    top.set(ylabel="coincidences (counts/sec)",
            title=f'{Path(res["file"]).name} | {res["date"]} | {res["isotope"]}, calib = {res["calibration_factor"]:g}')
    top.legend(fontsize="small")
    bottom.plot(tac.time, tac.activity, ".-", markersize=3, color=col_corr)
    bottom.axhline(0, color="0.7")
    bottom.axvline(0, color="black")
    bottom.set(xlim=xlim, xlabel="time relative to time 0 (s)",
               ylabel="activity (kBq/cc)", title="Corrected & calibrated")
    return fig
