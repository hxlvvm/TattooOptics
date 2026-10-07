# How inkppg works, in plain language

## What a smartwatch actually measures

The green or red light on the back of a smartwatch is a **photoplethysmography (PPG)** sensor. An LED
shines into the skin. Most of the light scatters around in the tissue, and a little of it comes back out
to a photodiode a few millimetres away. With every heartbeat the small blood vessels swell slightly. Blood
absorbs light, so slightly less light comes back.

That tiny rhythmic dip, often less than 1 % of the signal, is the pulse:

- **DC**: the steady amount of light that returns.
- **AC**: the part that pulses with the heartbeat.
- **AC/DC** (the *perfusion index*): how strong the pulse is relative to the light level.

Pulse oximeters compare the AC/DC of **red** light (660 nm) with **infrared** light (940 nm). Oxygenated
and deoxygenated blood absorb these two colours differently. The ratio of the two pulse strengths, the
*ratio of ratios*, is mapped to an oxygen saturation (SpO₂) through a calibration curve.

## Why tattoos are a problem

Tattoo pigment sits in the dermis, about 0.5–1 mm deep. That is right in the path the light takes
between the LED and the detector, and the pigment is a strong absorber.

- **Black ink** (carbon) absorbs every colour. Less light comes back, the shot noise of the photodiode
  stays, and the pulse can drown in noise. This is the heart-rate *dropout* that a 2025 study observed in
  about a third of participants at rest.
- **Coloured inks** absorb some colours much more than others. A red tattoo looks red because it absorbs
  green light, which is the colour most heart-rate sensors use. A blue or green tattoo absorbs red light.
  That changes the red channel without changing the infrared one, so it can bias the oxygen reading itself,
  not only add noise.

## The model

### Light transport: Monte Carlo

The model follows millions of simulated photons through the skin:

- Each photon travels a random distance and then scatters, changing direction according to the
  Henyey–Greenstein law that describes skin (strongly forward-scattering).
- Photons that reach the surface either escape or are reflected back inside, following Fresnel's equations
  for skin (refractive index 1.4) meeting air.
- For every photon that escapes, the model records **how far it went from the LED** (the source–detector
  distance) and **how much distance it travelled at each depth**, in 50 µm slices.

### The trick that makes it fast: white Monte Carlo

Absorption is **not** applied while the photons move. Instead, absorption is applied afterwards: a photon
that travelled path length `L_k` in a slice with absorption `μa_k` survives with probability
`exp(−Σ μa_k L_k)`. That is Beer–Lambert's law along its exact path.

One transport run per wavelength therefore answers every question afterwards: any skin tone, any blood
oxygen, any ink colour, density or depth. Each question is just a re-weighting of the stored photons and
takes milliseconds.

### From photons to a PPG

- **DC** is the summed weight of the photons that exit in a ring at the detector distance.
- **AC/DC**: a heartbeat raises blood volume by a small fraction. To first order, the dip in the signal is
  that change in absorption times the average distance the detected light travelled through blood. Both
  numbers come straight out of the weighted photons.
- **Noise**: a photodiode counting `N` photons has shot noise √N. So the pulse **signal-to-noise ratio**
  is `AC/DC × √N`.
- **Oxygen bias**: run the ratio of ratios on ink-free skin over a range of true SpO₂ values; that is the
  calibration curve. Then compute the ratio on tattooed skin and read the curve backwards. The difference
  is the error an oximeter would report.

## How it is checked

- The photon transport is compared with **diffusion theory** (Farrell 1992; Kienle & Patterson 1997), an
  independent analytic solution. They agree within 15 % at 1–4 mm, which is the expected accuracy of
  diffusion theory near the source.
- Unit checks cover Fresnel reflection, the scattering law, photon bookkeeping and the direction of every
  ink effect.

## What this is not

- **Inks are hypothetical.** Real tattoo pigments vary by brand and their optical properties are poorly
  documented, so the inks here are parametric what-if absorbers. Treat the results as maps of *how much
  could go wrong*, not as predictions for a specific device.
- The skin is a flat layered slab with the same scattering at every depth. There are no hair follicles,
  no motion and no ambient light.
- The photon budget that sets the absolute SNR is an assumption. Compare configurations with each other.
