"""Figures from the saved Monte Carlo runs (run examples/make_mc.py first).

1. assets/failure_map.png  pulse-SNR loss (dB) vs ink density for each ink colour and LED wavelength
2. assets/spo2_bias.png     apparent SpO2 of an oximeter calibrated on ink-free skin, for coloured inks
3. assets/separation.png    which source-detector spacing keeps the most SNR under black ink
4. assets/ppg.gif           a green-LED PPG waveform losing its pulse as ink density rises
Prints the numbers used in the README.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from inkppg import MCResult, Skin, apparent_spo2, measure, optics, ratio_of_ratios, spo2_calibration

MC = {wl: MCResult.load(f"data/mc_{int(wl)}.npz") for wl in optics.WAVELENGTHS}
RHO = 2.0
DENS = np.linspace(0, 4, 21)
COL = {"black": "#222222", "red": "#c0392b", "green": "#1e8449", "blue": "#1f4e9c"}
LED = {optics.GREEN: "green 530 nm", optics.RED: "red 660 nm", optics.IR: "infrared 940 nm"}
Path("assets").mkdir(exist_ok=True)
skin = Skin()

# ---- 1. SNR loss map
fig, ax = plt.subplots(1, 3, figsize=(11, 3.2), sharey=True)
print("SNR loss at ink density 2/mm, rho 2 mm:")
for a, wl in zip(ax, optics.WAVELENGTHS):
    base = measure(MC[wl], skin, wl, RHO).snr_db
    for ink in optics.INKS:
        loss = [measure(MC[wl], skin.with_ink(ink, d), wl, RHO).snr_db - base for d in DENS]
        a.plot(DENS, loss, color=COL[ink], label=f"{ink} ink", lw=2)
        print(f"  {LED[wl]:16s} {ink:6s} {np.interp(2.0, DENS, loss):6.1f} dB")
    a.axhline(-20, ls="--", color="grey", lw=1)
    a.set_title(LED[wl], fontsize=10)
    a.set_xlabel("ink density (peak absorption, 1/mm)", fontsize=8)
    a.tick_params(labelsize=8)
ax[0].set_ylabel("pulse SNR change (dB)", fontsize=9)
ax[0].text(0.1, -19, "20 dB loss (illustrative dropout margin)", fontsize=7, color="grey", va="bottom")
ax[2].legend(fontsize=8)
fig.suptitle("How much pulse signal survives a tattoo, by ink colour and LED wavelength", fontsize=11)
fig.tight_layout()
fig.savefig("assets/failure_map.png", dpi=120)

# ---- 2. SpO2 bias: an oximeter calibrated on ink-free skin, reading through ink
calib = spo2_calibration(MC[optics.RED], MC[optics.IR], skin, 3.0)
fig, axs = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True)
for ax, true in zip(axs, (0.90, 0.98)):
    sk = Skin(spo2=true)
    print(f"Apparent SpO2 at true {true * 100:.0f} %, rho 3 mm:")
    for ink in optics.INKS:
        app = [100 * apparent_spo2(ratio_of_ratios(MC[optics.RED], MC[optics.IR], sk.with_ink(ink, d), 3.0), calib)
               for d in DENS]
        ax.plot(DENS, app, color=COL[ink], lw=2, label=f"{ink} ink")
        print(f"  {ink:6s} density 0.5: {np.interp(0.5, DENS, app):5.1f} %   1: {np.interp(1.0, DENS, app):5.1f} %"
              f"   2: {np.interp(2.0, DENS, app):5.1f} %")
    ax.axhline(true * 100, ls=":", color="grey")
    ax.set_title(f"true SpO$_2$ {true * 100:.0f} %", fontsize=10)
    ax.set_xlabel("ink density (peak absorption, 1/mm)", fontsize=9)
    ax.set_ylim(85, 100.5)
axs[0].set_ylabel("SpO$_2$ the oximeter reports (%)", fontsize=9)
axs[1].legend(fontsize=8)
fig.suptitle("Coloured inks that absorb red light make the oximeter read HIGH", fontsize=11)
fig.tight_layout()
fig.savefig("assets/spo2_bias.png", dpi=120)

# ---- 3. separation
fig, ax = plt.subplots(figsize=(5.2, 3.6))
seps = np.arange(1.0, 5.01, 0.5)
for wl in optics.WAVELENGTHS:
    clean = [measure(MC[wl], skin, wl, s).snr_db for s in seps]
    inked = [measure(MC[wl], skin.with_ink("black", 2.0), wl, s).snr_db for s in seps]
    ax.plot(seps, clean, color=COL["green"] if wl == optics.GREEN else ("#c0392b" if wl == optics.RED else "#7d3c98"),
            lw=2, label=f"{LED[wl]}, no ink")
    ax.plot(seps, inked, ls="--", color=ax.lines[-1].get_color(), lw=2, label=f"{LED[wl]}, black ink")
ax.set_xlabel("source-detector separation (mm)", fontsize=9)
ax.set_ylabel("pulse SNR (dB, arbitrary photon budget)", fontsize=9)
ax.set_title("Sensor spacing vs tattoo", fontsize=10)
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig("assets/separation.png", dpi=120)

# ---- 4. waveform GIF (green LED, black ink density ramp). Absolute SNR depends on LED power and
# perfusion, which are assumptions; the animation uses an illustrative low-power case with 20 dB pulse SNR
# on ink-free skin and applies the model's SNR loss for each ink density.
rng = np.random.default_rng(0)
t = np.linspace(0, 4, 400)
beat = np.maximum(0, np.sin(2 * np.pi * 1.2 * t)) ** 1.5 + 0.35 * np.maximum(0, np.sin(2 * np.pi * 1.2 * t - 2.2)) ** 2
beat = (beat - beat.mean()) / np.abs(beat - beat.mean()).max()
dens = np.concatenate([np.linspace(0, 4, 50), np.full(10, 4.0)])
BASE_DB = 20.0
fig, ax = plt.subplots(figsize=(5.6, 2.6))
(line,) = ax.plot([], [], color="#1e8449", lw=1.2)
ax.set_xlim(0, 4)
ax.set_ylim(-14, 14)
ax.set_xlabel("time (s)", fontsize=8)
ax.set_yticks([])
title = ax.set_title("", fontsize=9)
base = measure(MC[optics.GREEN], skin, optics.GREEN, RHO)


def frame(i):
    p = measure(MC[optics.GREEN], skin.with_ink("black", dens[i]), optics.GREEN, RHO)
    snr_db = BASE_DB + p.snr_db - base.snr_db
    y = 10 ** (snr_db / 20) * beat + rng.normal(size=t.size)        # noise has unit standard deviation
    line.set_data(t, y)
    title.set_text(f"green-LED PPG, black ink {dens[i]:.1f}/mm: pulse SNR {snr_db:4.1f} dB (illustrative)")
    return line, title


fig.tight_layout()
FuncAnimation(fig, frame, frames=len(dens), blit=False).save("assets/ppg.gif", writer=PillowWriter(fps=8), dpi=90)
print("wrote assets/failure_map.png, spo2_bias.png, separation.png, ppg.gif")
