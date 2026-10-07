"""Run the white Monte Carlo for each LED wavelength and save data/mc_<nm>.npz."""
import sys
import time
from pathlib import Path

from tattooptics import mc, optics

wls = [float(sys.argv[1])] if len(sys.argv) > 1 else list(optics.WAVELENGTHS)
n = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
Path("data").mkdir(exist_ok=True)
for wl in wls:
    t = time.time()
    r = mc.run(optics.scattering(wl), n_photons=n, seed=int(wl))
    r.save(f"data/mc_{int(wl)}.npz")
    print(f"{wl:.0f} nm: mu_s {r.mus:.1f}/mm, {len(r.rho)} of {n} photons exit within 8 mm ({time.time() - t:.0f} s)")
