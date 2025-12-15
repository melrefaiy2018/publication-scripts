import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit

PH2KCAL = 1.364

# Sigmoid function
def sigmoid(x, x0, k):
    e = np.exp(k*(x-x0))
    y = e/(1.0+e)
    return y

# Function to plot results with different colors and simpler legend
def plot_res(resname, charge, t_type, t_points, color, lw):
    deltaY = charge[-1] - charge[0]
    a = math.ceil(abs(deltaY))
    lowest = min(charge)
    b = math.floor(lowest)

    x = t_points
    y = [(c - b) / a for c in charge]
    xdata = np.array(x)
    ydata = np.array(y)

    try:
        (popt, _) = curve_fit(sigmoid, xdata, ydata, bounds=([x[0], -4], [x[-1], 4]))
    except (RuntimeError, ValueError):
        print(f"Error fitting data for residue {resname}")
        return

    Npoints = 500
    xfit = np.linspace(xdata[0], xdata[-1], Npoints)
    yfit = [a * sigmoid(x, popt[0], popt[1]) + b for x in xfit]

    # Save the data to a file
    data_file = f"{resname}_fit_data.csv"
    with open(data_file, 'w') as f:
        f.write("xdata,ydata,xfit,yfit\n")
        for xd, yd, xf, yf in zip(xdata, ydata, xfit, yfit):
            f.write(f"{xd},{a * yd + b},{xf},{yf}\n")
    print(f"Saved data for {resname} to {data_file}")

    # Plotting
    plt.plot(xfit, yfit, '-', color=color, linewidth=lw,
             label=f'{resname}: pKa={popt[0]:.2f}')


# Function to fit and plot pKa for residues
def fitpka(resnames, input_file):
    lines = open(input_file).readlines()  # Input file provided as an argument
    headline = lines.pop(0)
    fields = headline.split()
    t_type = fields[0].strip()
    t_points = [float(x) for x in fields[1:]]

    # Create figure with enhanced settings
    plt.figure(figsize=(10, 6))
    plt.rc('font', family='serif')  # Set font to serif for publication quality
    plt.rcParams.update({'font.size': 18})  # Increase font size for readability

    # Define a set of unique publication-quality colors (Color Universal Design palette)
    pub_colors = [
        "#0072B2",  # Blue
        "#D55E00",  # Orange
        "#009E73",  # Green
        "#CC79A7",  # Pink/Purple
        "#F0E442",  # Yellow
        "#56B4E9",  # Light blue
        "#E69F00",  # Orange-brown
        "#000000",  # Black
        "#999999",  # Gray
        "#F4A582"   # Light Red
    ]

    # Ensure enough colors for the number of residues
    color_cycle = pub_colors * (len(resnames) // len(pub_colors) + 1)

    lw = 4  # Set desired line width here

    for i, resname in enumerate(resnames):
        found = False
        charge = []
        for line in lines:
            fields = line.split()
            if fields[0] == resname:
                found = True
                charge = [float(x) for x in fields[1:]]
                break

        if found:
            color = color_cycle[i]  # Get unique color for each residue
            plot_res(resname, charge, t_type, t_points, color, lw)
            plt.xlabel(f'{t_type}', fontsize=22)
            plt.ylabel('Charge', fontsize=22)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)

    # Adjust legend to be simpler
    plt.legend(loc='best', fontsize=14, frameon=False)

    # Save the enhanced figure
    plt.tight_layout()
    plt.savefig(f'{resnames[0]}_titration_publication.png', dpi=600)  # Save as high-res PNG
    plt.show()

    return

# if __name__ == "__main__":
#     # Input file and residues provided directly in the script
#     input_file = "sum_crg.out"  # Replace with your input file path
#     resnames = ["ASP-A0057_", "ASP-D0057_"]  # Replace with your residue names

#     fitpka(resnames, input_file)
