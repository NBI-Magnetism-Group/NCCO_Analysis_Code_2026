# Data Analysis for - Investigation of the effect of reductive annealing on spin fluctuations in superconducting electron-doped Nd$_{1.85}$Ce$_{0.15}$CuO$_{4-\delta}$

## Overview
This repository contains the scripts, raw data, and figures used for the data analysis presented in a scientific article. The datasets are derived from experiments conducted using the IN20 instrument at the Institut Laue-Langevin (ILL), the TAIPAN instrument at ANSTO, and the SQUID magnetometer measurements from Utrecht University.

---

## Repository Structure

The repository contains raw data organized by instrument and measurement type:

- **`Data_IN20/`** – Raw IN20 neutron scattering experiment data from Institut Laue-Langevin (ILL).  
- **`Data_TAIPAN/`** – Raw TAIPAN neutron scattering experiment data from ANSTO.  
- **`SQUID_Utrecht_data/`** – SQUID magnetometry raw data from Utrecht University.  
- **`(NCCO Laue/)`** – X-ray Laue images (not included by default; available upon request).  
- **`(Student_analysis/)`** – Supplemental student analysis files (not included by default). :contentReference[oaicite:3]{index=3}

---

## Figures

Figures generated from the analysis scripts are organized by output format and instrument:

- **`Figures_png/`**
  - `SQUID/`, `IN20/`, `TAIPAN/`
- **`Figures_eps/`**
  - `SQUID/`, `IN20/`, `TAIPAN/`
- **`Figures_svg/`**
  - `SQUID/`, `IN20/`, `TAIPAN/`  

Each subfolder contains the corresponding plots used for visualizing data trends and results. :contentReference[oaicite:4]{index=4}

---

## Scripts

Analysis code is provided primarily as **Jupyter Notebooks** and supporting Python modules:

- **Jupyter Notebooks**
  - `DataAnalysis_IN20.ipynb` – Analysis workflow for IN20 instrument data.  
  - `DataAnalysis_TAIPAN.ipynb` – Analysis workflow for TAIPAN data.  
  - `SQUID_analysis.ipynb` – Analysis workflow for SQUID magnetometer data.  
- **Python Modules**
  - `TASDataObjectIN20.py`  
  - `TASDataObjectIN20SQ.py`  
  - `TASDataObjectTAIPANSQ.py`  
  - `TASDataObjectTAIPAN_CHI.py`  

These modules define `TASDataObject` classes and helper functions used within notebooks for data handling, fitting, and plotting. :contentReference[oaicite:5]{index=5}

---

## How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/KristineKrighaar/NCCO_Article_Code
   ```
2. Create the Anaconda environment:
   ```bash
   conda env create -f NCCO_Env.yml
   conda activate NCCO_Env
   ```
3. Open and run the analysis scripts using Jupyter Notebook:
   ```bash
   jupyter notebook DataAnalysis_IN20.ipynb
   ```
   or
   ```bash
   jupyter notebook DataAnalysis_TAIPAN.ipynb
   ```
   or
   ```bash
   jupyter notebook SQUID_analysis.ipynb
   ```

---

## Dependencies
The project relies on the following Python packages (managed via the Anaconda environment in `NCCO_Env.yml` which can be provided upon request):
- `numpy`
- `matplotlib`
- `pandas`
- `scipy`
- `iminuit`
- `math`
- `sys`
- `scipy.optimize`
- `glob`
- `matplotlib.collections.LineCollection`
- `matplotlib.colors`

Kristine's Golden Standard (KGS)

The repository uses KGS (Kristine's Golden Standard), a personal library that incorporates robust data-fitting techniques. It leverages the iminuit library to perform:

- Chi-squared minimization
- Customizable regression analysis

Matplotlib configurations included:
- Line width: `mpl.rcParams['lines.linewidth'] = 0.3`
- Error bar cap size: `mpl.rcParams['errorbar.capsize'] = 2`
- Marker size: `mpl.rcParams['lines.markersize'] = 15`
- Font size: `mpl.rcParams['font.size'] = 14`
- Font family: `mpl.rcParams['font.serif'] = 'Computer Modern Roman'`
- Ticks: `mpl.rcParams['ytick.right'] = True`, `mpl.rcParams['xtick.top'] = True`
- Tick direction: `mpl.rcParams['xtick.direction']='in'`, `mpl.rcParams['ytick.direction']='in'`

---

## Contributing
This repository contains the finalized code for a published scientific article. In order to maintain documentation, contributions are not welcome.

---

## License
This project is licensed under the [MIT License](LICENSE). Feel free to use the materials as needed.

---

## Acknowledgments
This work incorporates data from the following sources:

IN20 instrument at Institut Laue-Langevin (ILL): under proposal number TEST-3346 (doi:10.5291/ILL-DATA.TEST-3346)
TAIPAN instrument at ANSTO under proposal P13914
SQUID measurements from Utrecht University

