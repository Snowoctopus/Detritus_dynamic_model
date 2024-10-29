
from pylab import *
from af001_energy_calculation import *
from Solver_2D import *
from Solver_3D import *
from Qc_essential_computation import *
from FigSetting2 import *
from matplotlib.pyplot import figure, stackplot, plot, legend, xlabel, ylabel, xticks, yticks, xlim, ylim, show, margins

rcParams['axes.xmargin'] = 0
rcParams.update({'font.size': 24})

# Original single-cell model parameters
m = 3.79146798299876E-19  # (mol C s-1 cell-1) maintenance carbohydrate consumption
Pmax = 0.00320513285659728
OT = 0.00863364097132997
Ynphoto_chl = 3.56099164557551  # Photosynthetic enzyme nitrogen to chlorophyll ratio
Cnbiosynth = 4.34728279914354E-10  # Biosynthesis protein nitrogen
Nconst_protein = 4.45336898828389E-15  # Constant protein pool in nitrogen
Nstore_max = 2.91679384515998E-15  # Nitrogen store max
Cnrna_variable = 6212.59249917364  # Variable part of RNA
Ypthylakoid_chl = 0.0281633095303638  # Phosphorus in thylakoid membrane to chlorophyll ratio
Pconst_other = 5.44534485638617E-17  # Constant part of phosphorus
Qp_max = 25.26/(3.097e16)  # Total phosphorus content in the cell
Cessential = 1.51786753491048E-15  # Essential carbon

E3 = evalue()
E = E3.E
Qc = 1.00*10**(-12)/12  # Biomass C per cell
YchlN_C = 4/55

CNprotein = 4.49  # C/N ratio in protein
YcyanoC_N = 2  # C/N molar ratio of cyanophycin
YpgC_P = 40  # C/P molar ratio of PG
Nunit = 1/Qc  # Unit conversion term
Punit = 1/Qc*30.97*10**6/(12*10**3)  # Phosphorus unit conversion

#==============================
# Original single-cell growth rate model (solver for D[])
#==============================
def calculate_single_cell_growth_rate():
    # Here D[] is the array calculated using the original solver for different NO3 concentrations
    NO3_concentrations = arange(0.01, 20+1e-10, 0.01)
    D = zeros(len(NO3_concentrations))
    
    # Example calculation for D[] based on NO3 concentration
    for i, NO3 in enumerate(NO3_concentrations):
        aNO3 = 0.05 * Qc / 86400
        Vn = aNO3 * NO3
        # Additional calculations from original single-cell model here...
        # Placeholder for actual solver:
        D[i] = 0.001 * NO3  # Replace with actual calculation based on solver_3D etc.
    
    return NO3_concentrations, D

#==============================
# Population dynamics model
#==============================
def get_single_cell_growth_rate(NO3, NO3_array, growth_rates):
    """根据当前NO3浓度获取单细胞生长速率D"""
    return np.interp(NO3, NO3_array, growth_rates)

def population_dynamics(time, NO3, NO3_array, growth_rates, m, initial_cells):
    """基于生长速率和死亡率计算种群数量变化"""
    N = np.zeros(len(time))
    N[0] = initial_cells
    for t in range(1, len(time)):
        mu = get_single_cell_growth_rate(NO3, NO3_array, growth_rates)  # 获取当前NO3对应的单细胞生长速率
        dN = (mu - m) * N[t-1]  # 生长和死亡的平衡
        N[t] = N[t-1] + dN * (time[t] - time[t-1])  # 使用欧拉法迭代
    return N

def detritus_dynamics(N, m, time):
    """基于死亡率计算碎屑生成"""
    D = np.zeros(len(time))
    for t in range(1, len(time)):
        D[t] = D[t-1] + m * N[t-1] * (time[t] - time[t-1])
    return D

def macromolecule_degradation(initial_content, decay_rate, time):
    """大分子降解模型"""
    return initial_content * np.exp(-decay_rate * time)

#==============================
# Main simulation
#==============================
# Get single cell growth rate data (D[])
NO3_concentrations, D = calculate_single_cell_growth_rate()

# Population dynamics simulation parameters
m = 0.1  # Death rate
initial_cells = 1  # Initial cell number
time = np.linspace(0, 50, 500)  # Time (days)
NO3 = 10  # Current NO3 concentration

# 1. Simulate population growth
N = population_dynamics(time, NO3, NO3_concentrations, D, m, initial_cells)

# 2. Simulate detritus formation
detritus = detritus_dynamics(N, m, time)

# Assume initial macromolecule composition in detritus (protein, RNA, etc.)
initial_protein = detritus[0] * 0.4  # Assume 40% protein in detritus
initial_rna = detritus[0] * 0.2  # Assume 20% RNA in detritus
protein_decay_rate = 0.05  # Protein decay rate
rna_decay_rate = 0.1  # RNA decay rate

# Simulate macromolecule degradation in detritus
protein_content = macromolecule_degradation(initial_protein, protein_decay_rate, time)
rna_content = macromolecule_degradation(initial_rna, rna_decay_rate, time)

# Plot results
plt.figure(figsize=(12, 6))

# Population size
plt.subplot(2, 2, 1)
plt.plot(time, N, label='Population Size')
plt.xlabel('Time (days)')
plt.ylabel('Number of Cells')
plt.title('Population Dynamics')

# Detritus formation
plt.subplot(2, 2, 2)
plt.plot(time, detritus, label='Detritus')
plt.xlabel('Time (days)')
plt.ylabel('Detritus (mol)')
plt.title('Detritus Formation')

# Macromolecule degradation
plt.subplot(2, 2, 3)
plt.plot(time, protein_content, label='Protein in Detritus')
plt.plot(time, rna_content, label='RNA in Detritus')
plt.xlabel('Time (days)')
plt.ylabel('Macromolecule Content (mol)')
plt.title('Macromolecule Degradation in Detritus')
plt.legend()

plt.tight_layout()
plt.show()
