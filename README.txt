# Data Analysis for - Investigation of the effect of reductive annealing on spin fluctuations in superconducting electron-doped Nd$_{1.85}$Ce$_{0.15}$CuO$_{4-\delta}$

## Overview
This repository contains the scripts, raw data, and figures used for the data analysis presented in a scientific article. The datasets are derived from experiments conducted using the IN20 instrument at the Institut Laue-Langevin (ILL), the TAIPAN instrument at ANSTO, and the SQUID magnetometer measurements from Utrecht University.

---

## Repository Structure

### **Data Folders**
- **`Data_IN20/`**: Contains raw data files from the IN20 instrument at ILL.
- **`Data_TAIPAN/`**: Contains raw data files from the TAIPAN instrument at ANSTO.
- **`SQUID_Utrecht_data/`**: Contains raw SQUID data from measurements conducted at Utrecht University.

### **Figure Folders**
The repository includes generated figures in multiple formats, organized by instrument:
- **`Figures_png/`**: Figures in PNG format.
  - `SQUID/`: Figures from SQUID measurements.
  - `IN20/`: Figures from IN20 experiments.
  - `TAIPAN/`: Figures from TAIPAN experiments.
- **`Figures_eps/`**: Figures in EPS format.
  - `SQUID/`: Figures from SQUID measurements.
  - `IN20/`: Figures from IN20 experiments.
  - `TAIPAN/`: Figures from TAIPAN experiments.
- **`Figures_svg/`**: Figures in SVG format.
  - `SQUID/`: Figures from SQUID measurements.
  - `IN20/`: Figures from IN20 experiments.
  - `TAIPAN/`: Figures from TAIPAN experiments.

### **Scripts**
- **`DataAnalysis_IN20.ipynb`**: Jupyter Notebook for analyzing data from the IN20 experiment.
- **`DataAnalysis_TAIPAN.ipynb`**: Jupyter Notebook for analyzing data from the TAIPAN experiment.
- **`SQUID_analysis.ipynb`**: Jupyter Notebook for analyzing the SQUID measurement data.
- **`TASDataObjectIN20.py`**: Contains a dataclass that aids in the IN20 analysis notebooks, with functions for data handling and processing.
- **`TASDataObjectIN20SQ.py`**: Similar to `TASDataObjectIN20.py` but performs fitting and analysis in units of S(q,omega) instead of χ''.
- **`TASDataObjectTAIPANSQ.py`**: Modified for TAIPAN data structure, providing fitting and analysis in units of S(q,omega).
- **`NCCO_Env.yml`**: Anaconda environment file containing all dependencies for this project.

---

## How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/your-repo-name.git
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
The project relies on the following Python packages (managed via the Anaconda environment in `NCCO_Env.yml`):
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
This repository contains the finalized code for a published scientific article. Contributions are not required.

---

## License
This project is licensed under the [MIT License](LICENSE). Feel free to use and adapt the materials as needed.

---

## Acknowledgments
This work incorporates data from the following sources:
- IN20 instrument at Institut Laue-Langevin (ILL)
- TAIPAN instrument at ANSTO
- SQUID measurements from Utrecht University

