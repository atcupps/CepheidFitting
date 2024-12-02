import numpy as np
import matplotlib.pyplot as plt
import math

'''
This file is a program for fitting Cepheid variable stellar magnitude data
with a purposely-underfit Fourier series. Output is used qualitatively,
since the magnitudes are normalized and centered at zero to improve fitting.
Written by Andrew Cupps.

To run this fitting program:
- Replace `oid` and `period`
- Rename the data `{oid}.csv`, and place in the `data` folder
- Run: `python3 fitting.py` from the same folder as this script. Otherwise,
    you'll need to change the `path_to_data` string
'''
oid = '512115300010506'
period = 60.96616374
path_to_data = 'data/'

def fourier_series(x_arr, f_arr, n, L):
    '''
    Given data points `x_arr` and `f_arr`, such that the x_arr[i] and f_arr[i]
    together are one 2-d data point, as well as period length `L` and a series
    value `n`, calculate the Fourier series coefficients `a0, an, and bn`.
    
    Parameters
    ----------
    x_arr : array_like
        x-axis values of a function to fit a Fourier series to
    
    f_arr : array_like
        y-axis values of a function to fit a Fourier series to, such that
        `f_arr[i]` corresponds to `x_arr[i]`
    
    n : int
        Number of sin and cos terms in the series
    
    L : float
        Period of the function

    Returns
    -------
    (a0, a_coeffs, b_coeffs)
        These values can be used to construct a Fourier series of the data
        passed into this function; such a series would look like:

        ```
        a0 + 
            a_coeffs[0] * cos(2*1*pi*x/L) + 
            a_coeffs[1] * cos(2*2*pi*x/L) +
            ... +
            a_coeffs[n] * cos(2*n*pi*x/L) +
            
            b_coeffs[0] * sin(2*1*pi*x/L) +
            b_coeffs[1] * sin(2*2*pi*x/L) +
            ... +
            b_coeffs[n] * sin(2*n*pi*x/L)
        ```
    '''

    # Calculating a0
    prev_x = 0
    a0_integral = 0
    for (x, f) in zip(x_arr, f_arr):
        a0_integral += (x - prev_x) * f
        prev_x = x
    a0_integral += (L - prev_x) * (f_arr[0] + f_arr[-1]) / 2
    a0 = (1 / L) * a0_integral

    # Calculating a_n for each n
    a_coeffs = np.array([])
    for i in range(1, n + 1):
        an_func = np.array([])
        for (x, f) in zip(x_arr, f_arr):
            an_func = np.append(an_func, f * np.cos((2 * math.pi * i * x) / L))
        prev_x = 0
        an_integral = 0
        for (x, f) in zip(x_arr, an_func):
            an_integral += (x - prev_x) * f
            prev_x = x
        an_integral += (L - prev_x) * (an_func[0] + an_func[-1]) / 2
        a_coeffs = np.append(a_coeffs, (2 / L) * an_integral)
    
    # Calculating b_n for each n
    b_coeffs = np.array([])
    for i in range(1, n + 1):
        bn_func = np.array([])
        for (x, f) in zip(x_arr, f_arr):
            bn_func = np.append(bn_func, f * np.sin((2 * math.pi * i * x) / L))
        prev_x = 0
        bn_integral = 0
        for (x, f) in zip(x_arr, bn_func):
            bn_integral += (x - prev_x) * f
            prev_x = x
        bn_integral += (L - prev_x) * (bn_func[0] + bn_func[-1]) / 2
        b_coeffs = np.append(b_coeffs, (2 / L) * bn_integral)
    
    return (a0, a_coeffs, b_coeffs)

# Importing data
data = np.genfromtxt(f'{path_to_data}{oid}.csv', delimiter=',')[1:]
date = data[:,0]
mag = data[:,1]
err = data[:,2]

# Calculating phase and adjusted magnitudes + error values
phase = ((date - min(date)) / period) - np.floor((date - min(date)) / period)
err_adj = [f for _, f in sorted(zip(phase, err))] / np.median(mag)
mag_adj = [f for _, f in sorted(zip(phase, mag))] / np.median(mag) - 1
phase_sorted = sorted(phase)

# Find the `n` from 1 to 5 which best fits phase plot data
L = 1
best_fit_n = 1
best_fit_err = 100000000000000
for n in range(1, 5):
    (a0, an, bn) = fourier_series(phase_sorted, mag_adj, n, L)

    fitted_f = np.array([])
    for x_val in phase_sorted:
        f_val = a0
        for i in range(1, n+1):
            f_val += an[i-1]*math.cos((2*math.pi*i*x_val) / L) + bn[i-1]*math.sin((2*math.pi*i*x_val) / L)
        fitted_f = np.append(fitted_f, f_val)

    # Finding MSE
    num_points = len(phase_sorted)
    total_squared_error = 0
    for (y_pred, y) in zip(mag_adj, fitted_f):
        total_squared_error += (y_pred - y) ** 2
    mean_squared_error = total_squared_error / num_points
    if mean_squared_error < best_fit_err:
        best_fit_n = n
        best_fit_err = mean_squared_error

# Getting best-fit Fourier coefficients
(a0, an, bn) = fourier_series(phase_sorted, mag_adj, n, L)

# Using the best-fit Fourier coefficients, generate phase plot fit
fitted_f = np.array([])
for x_val in phase_sorted:
    f_val = a0
    for i in range(1, best_fit_n+1):
        f_val += an[i-1]*math.cos((2*math.pi*i*x_val) / L) + bn[i-1]*math.sin((2*math.pi*i*x_val) / L)
    fitted_f = np.append(fitted_f, f_val)

print(f"N: {best_fit_n}")
print(f"MSE: {best_fit_err}")

# Generate phase plot with Fourier series fit
plt.figure()
plt.errorbar(phase_sorted, mag_adj, yerr=err_adj, color='orange', alpha=0.7, fmt='o', zorder=1, label='Measured Data')
plt.plot(phase_sorted, fitted_f, zorder=2, label='Fit')
plt.gca().invert_yaxis()
plt.title(f'Adjusted Magnitude Data with Fit of {oid}\nPeriod of {period} days')
plt.xlabel('Phase')
plt.ylabel('Adjusted Magnitude')
plt.legend()
plt.show()