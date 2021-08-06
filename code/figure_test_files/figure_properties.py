# Importing matplotlib
import matplotlib.pyplot as plt
import matplotlib

# Defining scaling factor
b = 1.25

# Defining font-sizes
SMALL_SIZE = 8 * b
MEDIUM_SIZE = 10 * b
BIGGER_SIZE = 12 * b

# Defining figure properties
lwidth = 6.202
markersize = 4

# Defining 1/quantized conductance
R_K = 25812.8

# Initializing the standard figure parameters
fig_params = {
            'font.size'        : SMALL_SIZE,   # controls default text sizes
            'xtick.labelsize'  : SMALL_SIZE,   # fontsize of the tick labels
            'ytick.labelsize'  : SMALL_SIZE,   # fontsize of the tick labels
            'axes.labelsize'   : MEDIUM_SIZE,  # fontsize of the x and y labels
            'axes.titlesize'   : SMALL_SIZE,   # fontsize of the axes title
            'legend.fontsize'  : SMALL_SIZE,   # legend fontsize
            'figure.titlesize' : BIGGER_SIZE,  # fontsize of the figure title
            'text.usetex'      : True,         # Latex font
}

# Updating aixs properties
plt.rcParams.update(fig_params)

# Setting directory for figures
save_dir = '../data/figures/'
