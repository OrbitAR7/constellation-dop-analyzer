# LEO-PNT DOP Analysis

A lightweight Python tool for analyzing **Dilution of Precision (DOP)** performance using real **TLE data** from different satellite constellations such as **GPS (MEO)** and **Starlink (LEO)**.

## Overview

This tool computes and visualizes DOP metrics to evaluate how satellite geometry affects positioning accuracy across latitudes and mask angles. It’s designed for quick constellation geometry testing and comparison between LEO- and MEO-based PNT systems.

## Features

* Parse and propagate **TLE data** using `sgp4`
* Compute **HDOP**, **VDOP**, **PDOP**, and **GDOP**
* Analyze performance at multiple **elevation mask angles**
* Generate **global DOP contour maps**
* Compare **different constellations** (e.g., GPS vs Starlink)

## Example Data

The `/data` folder includes:

* `GPS_TLEs.txt` — Current GPS constellation (MEO, ~20,200 km, ~55° inclination)
* `STARLINK_DTC_TLES.txt` — Starlink DTC satellites (LEO, ~550 km, ~53° inclination)

## Visualization

You can generate contour plots showing how DOP varies across the globe for different constellations and mask angles.
<img width="3821" height="1951" alt="example_gdop_mask_30deg" src="https://github.com/user-attachments/assets/8c8821b7-08f9-4650-8414-95e3ceb3e0a0" />

## Notes

* DOP results depend heavily on constellation geometry and observer location.
* The provided TLEs are for demonstration and testing; users can replace them with updated or custom TLE datasets.
