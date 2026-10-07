"""tattooptics: how tattoo ink in the skin affects wrist photoplethysmography (PPG) and pulse oximetry."""
from . import mc, optics
from .mc import MCResult, run
from .ppg import PPG, Skin, absorption_profile, apparent_spo2, measure, ratio_of_ratios, spo2_calibration

__all__ = ["mc", "optics", "MCResult", "run", "Skin", "PPG", "absorption_profile", "measure",
           "ratio_of_ratios", "spo2_calibration", "apparent_spo2"]
__version__ = "0.1.0"
