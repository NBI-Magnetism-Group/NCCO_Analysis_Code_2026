import numpy as np
import pandas as pd
from KGS import *

class Dataset:
    def __init__(self, qh, qk, cnt, cnt_err, mn, en, tt, **kwargs):
        self.QH = qh
        self.QK = qk
        self.CNT = cnt
        self.CNT_err = cnt_err
        self.MN = mn
        self.EN = en
        self.TT = tt
        self.metadata = kwargs
        self.I = None
        self.I_err = None
        self.SQ = None
        self.SQ_err = None
        self.Chi = None
        self.Chi_err = None
        self.fq = None
        self.res_vol = None
        self.amp = None
        self.amp_err = None
        self.sig_area = None
        self.sig_area_err = None
        self.best_fit_obj = None
        self.fit_type = None

    def __repr__(self):
        return (f"EN={self.EN}, TT={self.TT}, "
                f"data_length={len(self.QH)}, metadata={self.metadata})")
    
    def calc_I(self):
        """
        Normalized raw counts with respect to the monitor
        """
        self.I = self.CNT/self.MN
        self.I_err = self.CNT_err/self.MN
    
    def calc_SQ(self):
        """
        Calculating the correlation functions

        RETURNS: in units of [abs. unt.]
        """
        self.SQ = (13.77*self.I)/(2**2*self.fq**2*1*self.res_vol)
        self.SQ_err = (13.77*self.I_err)/(2**2*self.fq**2*1*self.res_vol)

    def calc_Chi(self):
        """
        Calculating the correlation functions

        RETURNS: in units of $mu_B^2 eV^-1$
        """
        self.Chi = 10**(3)*np.pi/2*(1-np.exp(-self.EN/(0.08617*self.TT)))*(13.77*self.I)/(self.fq**2*1*self.res_vol)
        self.Chi_err = 10**(3)*np.pi/2*(1-np.exp(-self.EN/(0.08617*self.TT)))*(13.77*self.I_err)/(self.fq**2*1*self.res_vol)

    def p3_scan_amplitude_chi(self):
        """
        Extracts the signal amplitude from a 3-point inelastic neutron scattering measurement.

        Parameters:
        I (array-like): Array of intensities [I(E1), I(E2), I(E3)]
        I_error (array-like): Array of errors in intensities [error(E1), error(E2), error(E3)]

        Returns:
        I_signal (float): The extracted signal amplitude.
        I_signal_error (float): The error in the signal amplitude.
        """

        # Calculate the background as the average of the off-peak points (I(E1) and I(E3))
        I_bg = (self.Chi[0] + self.Chi[2]) / 2

        # Extract the signal amplitude by subtracting the background from the central point (I(E2))
        self.amp = self.Chi[1] - I_bg

        # Calculate the error in the signal amplitude
        self.amp_err = np.sqrt(self.Chi_err[1]**2 + ((self.Chi_err[0]**2 + self.Chi_err[2]**2) / 4))

    def FindBestFit(self, model1, model2, initial_guess1, initial_guess2, fixed_params1=None, limits1=None, fixed_params2=None, limits2=None):
        
        m_gauss = fit(self.QK, self.Chi, self.Chi_err, model1, initial_guess1, fixed_params=fixed_params1, limits=limits1)

        m_const = fit(self.QK, self.Chi, self.Chi_err, model2, initial_guess2, fixed_params=fixed_params2, limits=limits2)

        red_chi2_gauss = m_gauss.fval/(len(self.QK)-4)
        #print('gauss chi = ', red_chi2_gauss)
        red_chi2_const = m_const.fval/(len(self.QK)-1)
        #print('const chi = ', red_chi2_const)

        # Compare reduced chi-squared and return the better fit
        if red_chi2_gauss < red_chi2_const:
             self.fit_type = "gauss"
             self.best_fit_obj = m_gauss
        else:
            self.fit_type = "const"
            self.best_fit_obj = m_const


############################ Functions operating with the objects ##################

def extract_data_from_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Extract TT from the header
    tt = None
    for line in lines:
        if "PARAM: TT=" in line:
            tt = float(line.split("TT=")[-1].split(",")[0].strip())
        elif line.startswith("DATA_:"):
            # Stop parsing the header when reaching the data section
            break

    # Find the start of the data section (after DATA_:)
    data_start = next(i for i, line in enumerate(lines) if line.startswith("DATA_:")) + 1

    # Use numpy genfromtxt to load data while ignoring the header and first column (PNT)
    data = np.genfromtxt(filepath, skip_header=data_start, invalid_raise=False)

    # Ensure the data rows that contain header information are filtered out
    # Only take rows where valid numerical data exists
    data = data[~np.isnan(data).any(axis=1)]

    # Extract relevant columns from the data
    qh = data[:, 1]    # QH is in the 2nd column
    qk = data[:, 2]    # QK is in the 3rd column
    mn = data[:, 5]    # MN is in the 6th column (formerly M1)
    cnt = data[:, 8]   # CNT is in the 9th column
    cnt_err = np.sqrt(cnt)  # Assuming Poisson error for counts (CNT_err = sqrt(CNT))
    en = np.mean(data[:, 4])    # EN is in the 5th column

    # Create Dataset object
    dataset = Dataset(qh, qk, cnt, cnt_err, mn, en, tt)

    return dataset



def combine_datasets(dataset1, dataset2, qh_tolerance=1e-3, qk_tolerance=1e-3):
    # Create lists to hold the combined data
    combined_qh = []
    combined_qk = []
    combined_cnt = []
    combined_cnt_err = []
    combined_mn = []
    combined_tt = []
    combined_en = []

    # Convert the datasets into structured arrays for easier matching
    d1 = np.array(list(zip(dataset1.QH, dataset1.QK, dataset1.CNT, dataset1.CNT_err, dataset1.MN)))
    d2 = np.array(list(zip(dataset2.QH, dataset2.QK, dataset2.CNT, dataset2.CNT_err, dataset2.MN)))

    # Keep track of matched indices in dataset2
    matched_indices = []

    # Loop through dataset1 and find matches in dataset2
    for qh1, qk1, cnt1, cnt1_err, mn_1 in d1:
        # Find points in dataset2 that are close in QH and QK
        match = np.where((np.abs(d2[:, 0] - qh1) < qh_tolerance) & (np.abs(d2[:, 1] - qk1) < qk_tolerance))[0]

        if len(match) > 0:
            # There is a matching point in dataset2
            idx = match[0]  # Assuming first match for simplicity
            matched_indices.append(idx)

            qh2, qk2, cnt2, cnt2_err, mn_2 = d2[idx]

            # Combine CNT by summing
            combined_cnt_val = cnt1 + cnt2

            # Combine CNT errors in quadrature
            combined_cnt_err_val = np.sqrt(cnt1_err**2 + cnt2_err**2)

            # Combine M1 by summing
            combined_m1_val = mn_1 + mn_2

            # Use the average of TT and EN (assuming single values for TT and EN in each dataset)
            combined_tt_val = (dataset1.TT + dataset2.TT) / 2
            combined_en_val = (dataset1.EN + dataset2.EN) / 2

            # Append combined results
            combined_qh.append(qh1)
            combined_qk.append(qk1)
            combined_cnt.append(combined_cnt_val)
            combined_cnt_err.append(combined_cnt_err_val)
            combined_mn.append(combined_m1_val)
            combined_tt.append(combined_tt_val)
            combined_en.append(combined_en_val)

        else:
            # No match found in dataset2, use dataset1 values directly
            combined_qh.append(qh1)
            combined_qk.append(qk1)
            combined_cnt.append(cnt1)
            combined_cnt_err.append(cnt1_err)
            combined_mn.append(mn_1)
            combined_tt.append(dataset1.TT)
            combined_en.append(dataset1.EN)

    # Now add the points from dataset2 that were not matched
    for i, (qh2, qk2, cnt2, cnt2_err, mn_2) in enumerate(d2):
        if i not in matched_indices:
            combined_qh.append(qh2)
            combined_qk.append(qk2)
            combined_cnt.append(cnt2)
            combined_cnt_err.append(cnt2_err)
            combined_mn.append(mn_2)
            combined_tt.append(dataset2.TT)
            combined_en.append(dataset2.EN)

    # Create a new Dataset object with the combined data
    combined_dataset = Dataset(
        qh=np.array(combined_qh),
        qk=np.array(combined_qk),
        cnt=np.array(combined_cnt),
        cnt_err=np.array(combined_cnt_err),
        mn=np.array(combined_mn),
        en=np.mean(combined_en),  # Mean of all EN values
        tt=np.mean(combined_tt)   # Mean of all TT values
    )

    return combined_dataset



def plot_fits(data_objects):
    """
    Plots a grid of subplots showing the data points with error bars and the fitted curve.
    
    Parameters:
        data_objects (list): List of data objects, where each object contains x, y, yerr data and a best_fit_obj attribute.
    """
    num_plots = len(data_objects)
    
    # Determine the layout of subplots based on the number of plots
    cols = int(np.ceil(np.sqrt(num_plots)))
    rows = int(np.ceil(num_plots / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))
    axes = axes.flatten() if num_plots > 1 else [axes]  # Flatten in case of multiple axes
    
    for i, data_obj in enumerate(data_objects):
        ax = axes[i]

        # Extract data from each object
        x = data_obj.QK
        y = data_obj.Chi
        yerr = data_obj.Chi_err

        # Plot the data points with error bars
        ax.errorbar(x, y, yerr=yerr, fmt='o')

        # Generate x-values for plotting the fit
        x_fit = np.linspace(min(x), max(x), 500)

        # Plot the best fit based on the fit_type
        if data_obj.fit_type == "gauss":
            # Assuming best_fit_obj contains the parameters A, mu, sigma, and C for the Gaussian+constant fit
            A, mu, sigma, C = data_obj.best_fit_obj.values
            y_fit = (A / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((x_fit - mu)**2) / (2 * sigma**2)) + C
            ax.plot(x_fit, y_fit, label='Gaussian Fit', color='red')

        elif data_obj.fit_type == "const":
            # Assuming best_fit_obj contains just the constant C
            C = data_obj.best_fit_obj.values
            y_fit = np.full_like(x_fit, C)
            ax.plot(x_fit, y_fit, label='Constant Fit', color='b')

        # Add title and labels
        ax.set_title(f"{data_obj.EN:.1f} meV, {data_obj.TT:.1f} K: {data_obj.fit_type.capitalize()} Fit")
        ax.set_xlabel("qk")
        ax.set_ylabel("chi")
        ax.legend()

    # Remove any unused subplots if there are extra grid spaces
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()
