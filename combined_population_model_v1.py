from FigSetting2 import *
from Qc_essential_computation import *
from Solver_3D import *
from af001_energy_calculation import *


def single_cell_model(NO3_input):
    '''
    From a800_05_12_07_SameQc.py
    '''

    rcParams['axes.xmargin'] = 0
    # rcParams['axes.ymargin'] = 0
    rcParams.update({'font.size': 24})
    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # Parameter sets
    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    m = 3.79146798299876E-19    # (mol C s-1 cell-1) maintenance carbohydrate consumption (idea from 172-7)
    Pmax = 0.00320513285659728  # (g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990)

    OT = 0.00863364097132997
    Ynphoto_chl = 3.56099164557551         # ((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
    Cnbiosynth = 4.34728279914354E-10      # (molN cell-1 s) Constant for variable part of biosynthesis protein nitrogen (193-37)
    Nconst_protein = 4.45336898828389E-15  # (molN cell-1) Constant protein pool in nitrogen (193-25)
    Nstore_max = 2.91679384515998E-15      # (molN cell-1) Constant protein pool in nitrogen (193-25)
    Cnrna_variable = 6212.59249917364      # (s) Constant for Variable part of RNA (193-26)
    Ypthylakoid_chl = 0.0281633095303638   # ((molP cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
    Pconst_other = 5.44534485638617E-17    # (molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
    Qp_max = 25.26 / (3.097e16)            # (molP cell-1) total phosphorus content in the cell (193-26)
    Cessential = 1.51786753491048E-15      # (molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%

    # ==============================
    E3 = evalue()
    E = E3.E
    Qc = 1.00 * 10 ** (-12) / 12  # (molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
    YchlN_C = 4 / 55

    # Conversion parameters================
    CNprotein = 3.82  # (molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
    YcyanoC_N = 2     # (molC molN) C/N molar ratio of cyanophycin
    YpgC_P = 40       # (molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
    Nunit = 1 / Qc    # ((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
    Punit = 1 / Qc * 30.97 * 10 ** 6 / (12 * 10 ** 3)  # ((ug P / mgC)/(molP cell-1) unit conversion term (164-20)

    # Photosynthesis================
    I = 64                            # (umolE m-2 s-1) light intensity
    Pchl = Pmax * (1 - exp(-OT * I))  # (C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)

    # Nutrient======================
    # NO3 = arange(0.01, 20 + 1e-10, 0.01)  # (umol L-1)
    NO3 = NO3_input
    Loop_array = arange(0, size(NO3), 1)  # array for loops
    Numbertoarray = ones(size(NO3))       # (dimensionless) Number to array converter

    # =========================================
    # DNA and RNA-> Kei 193-28
    # =========================================
    Molar_mass_DNA_AT_average = 307.47  # (g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_DNA_CG_average = 307.97  # (g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_AT_average = 316.47  # (g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_CG_average = 323.97  # (g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")

    # ------------------------------------
    # E coli
    # ------------------------------------
    CG_Ecoli = 0.506         # (dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli = 1 - CG_Ecoli  # (dimensionless)

    Molar_mass_DNA_Ecoli = Molar_mass_DNA_AT_average * CG_Ecoli + Molar_mass_DNA_CG_average * AT_Ecoli  # (g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli = Molar_mass_RNA_AT_average * CG_Ecoli + Molar_mass_RNA_CG_average * AT_Ecoli  # (g mol-1) Molar mass of RNA unit

    RNA_DNA_mass_ratio = 17.844 / 6.5239                                                    # (ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    RNA_DNA_molar_ratio = RNA_DNA_mass_ratio / Molar_mass_RNA_Ecoli * Molar_mass_DNA_Ecoli  # (mol mol-1)

    # ---------------------------------------------
    # Stoichiometric parameters for DNA and RNA
    # ---------------------------------------------

    CG = 0.563  # GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    YnucacidP_N = 1 / (3.5 * (1 - CG) + 4 * CG)  # (molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"

    YdnaC_N = 3.5 * (1 - CG) + 2.5 * CG           # (molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N = 3.25 * (1 - CG) + 2.5 * CG          # (molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)

    DNAmb = 2.1269                                 # (Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    Avogadro = 6.022 * 10 ** 23                    # (molecules mol-1) Avogadro constant
    Pdna_const = DNAmb * 2 * 10 ** 6 / Avogadro    # (molP cell-1) Constant part of DNA in phosphorus
    Prna_const = Pdna_const * RNA_DNA_molar_ratio  # (molP cell-1) Constant part of RNA in phosphorus
    # * Make sure to multiply by 2 as they are base PAIRs"
    Ndna_const = Pdna_const / YnucacidP_N          # (molN cell-1) Constant part of DNA in nitrogen
    Nrna_const = Ndna_const * RNA_DNA_molar_ratio  # (molN cell-1) Constatn part of RNA in phosphorus
    Ndna = Ndna_const                              # (molN cell-1) DNA in nitrogen (here assuming constant)

    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # Calculation of carbon usage (195-16)
    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # ===============================
    # Setting arrays
    # ===============================
    def o():
        return zeros(size(NO3)) * nan

    Vn = o()
    D = o()

    # ======================
    # Uptake parameter
    # ======================
    aNO3 = 0.05 * Qc / 86400  # (mol N cell-1 s-1 / umolN L-1) affinity of NO3 (initial value from d007_06_00)

    for i in Loop_array:
        # --------------------------------------------------------
        # Main calculation 199-21~
        # --------------------------------------------------------
        # Vn = aNO3 * NO3[i]  # (mol N cell-1 s-1) nitrogen uptake per cell
        Vn = aNO3 * NO3

        # ++++++++++++++++++++++++++++++++++++
        # For obtaining D
        # ++++++++++++++++++++++++++++++++++++
        A = ((1 + E) * Qc * Ynphoto_chl) / Pchl + Cnbiosynth
        B = Nconst_protein + (m * Ynphoto_chl) / Pchl
        G = ((1 + E) * Qc * YchlN_C) / Pchl
        H = (m * YchlN_C) / Pchl
        I = A * Cnrna_variable
        J = A + B * Cnrna_variable + G
        K = B + Nrna_const + Ndna + H
        L = ((1 + E) * Qc * Ypthylakoid_chl) / Pchl
        M = (m * Ypthylakoid_chl) / Pchl
        N = A * Cnrna_variable * YnucacidP_N
        O = L + B * Cnrna_variable * YnucacidP_N
        P = M + Nrna_const * YnucacidP_N + Ndna * YnucacidP_N + Pconst_other

        aN = I
        bN = J
        cN = K
        dN = -Vn

        aNf = float(aN)
        bNf = float(bN)
        cNf = float(cN)
        dNf = float(dN)

        DN = solver_3D(aNf, bNf, cNf, dNf)
        # DN.X is the solving result of the 3D function
        D[i] = DN.X
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # Obtaining output values (Similar to the previous steady state one)
    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Chl = ((1 + E) * D * Qc + m) / Pchl  # (molC chl cell-1) cN[i]chlorophyll concentration (193-25)
    # =========================
    # Nitrogen related
    # =========================
    Nchl = Chl * YchlN_C                            # (molN chl cell-1) Chlorophyll N concentration
    Nphoto = Chl * Ynphoto_chl                      # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth = D * Cnbiosynth                      # (molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein = Nphoto + Nconst_protein + Nbiosynth  # (molN cell-1) All the proteins in N (193-26)
    Nrna_variable = Nprotein * D * Cnrna_variable   # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Nrna = Nrna_const + Nrna_variable
    Qn = Nchl + Nconst_protein + Nphoto + Nbiosynth + Nrna + Ndna  # + Nstore

    # =========================
    # Phosphorus related
    # =========================
    Pdna = Ndna * YnucacidP_N           # (mol P cell-1) DNA in phosphorus
    Pthylakoid = Chl * Ypthylakoid_chl  # (molP cell-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna = Nrna * YnucacidP_N           # (molP cell-1) Phosphorus in RNA

    # ==============================
    # Carbon related: Solver setting up (for C) (Kei 200-33)
    # ==============================
    # Preparation---------------------------------------
    # Chlorophyll related
    Chl_const = m / Pchl         # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
    Chl_D = (1 + E) * Qc / Pchl  # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)

    # Nitrogen related
    Nchl_D = Chl_D * YchlN_C                        # (molN chl cell-1) Chlorophyll N concentration
    Nphoto_D = Chl_D * Ynphoto_chl                  # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nchl_const = Chl_const * YchlN_C                # (molN chl cell-1) Chlorophyll N concentration
    Nphoto_const = Chl_const * Ynphoto_chl          # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth_D = Cnbiosynth                        # (molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein_D = Nphoto_D + Nbiosynth_D             # (molN cell-1) All the proteins in N (193-26)
    Nprotein_const = Nphoto_const + Nconst_protein  # (molN cell-1) All the proteins in N (193-26)
    Nrna_D = Nprotein_const * Cnrna_variable        # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Nrna_D2 = Nprotein_D * Cnrna_variable           # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)

    # Constant carbon parameters------------------
    Cconst_protein = Nconst_protein * CNprotein  # (molC cell-1) carbon in other protein assumed constant (195-16)
    Cdna_const = Ndna_const * YdnaC_N            # (molC cell-1) carbon in constant part of DNA (195-16)

    # Calculating factors and Qc_essential-------------
    Qc_D2 = A * Cnrna_variable * YrnaC_N
    Qc_D = (1 + E) * Qc / Pchl + A * CNprotein + L * YpgC_P + B * Cnrna_variable * YrnaC_N
    Qc_const = m / Pchl + B * CNprotein + M * YpgC_P + Nrna_const * YrnaC_N + Ndna_const * YdnaC_N + Cessential

    Qc_essential = Qc_essential_computation(D, Chl_const, Chl_D, Nchl_const, Nchl_D, Nphoto_const, Nphoto_D,
                                            Nbiosynth_D, \
                                            Nprotein_const, Nprotein_D, Nrna_const, Nrna_D, Nrna_D2, Ypthylakoid_chl,
                                            YrnaC_N, \
                                            CNprotein, YpgC_P, Cdna_const, Cconst_protein, Cessential)

    i = Qc_essential >= Qc  # conditions where Mu max applies
    D[i] = 2 * (Qc - Qc_const) / (Qc_D + sqrt(Qc_D * Qc_D + 4 * Qc_D2 * (Qc - Qc_const)))

    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # for plotting
    # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Chl = Chl_const + Chl_D * D
    # =================
    # Nitrogen related
    # =================
    Nchl = Nchl_const + Nchl_D * D        # (molN chl cell-1) Chlorophyll N concentration
    Nphoto = Nphoto_const + Nphoto_D * D  # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth = Nbiosynth_D * D
    Nprotein = Nprotein_const + Nprotein_D * D  # (mol N cell-1) all the protein
    Nrna = Nrna_const + Nrna_D * D + Nrna_D2 * D * D
    Nessential = Nchl + Nphoto + Nbiosynth + Nconst_protein + Nrna + Ndna
    Qn_max = Nessential + Nstore_max
    # =========================
    # Phosphorus related
    # =========================
    Pthylakoid = Chl * Ypthylakoid_chl                    # (molP cell-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna = Nrna * YnucacidP_N                             # (molP cell-1) Phosphorus in RNA
    Pessential = Pthylakoid + Prna + Pdna + Pconst_other  # (molP cell-1)

    # ==================================
    # Unite conversion
    # ==================================
    # ----------------------------------
    # For N/C
    # ----------------------------------
    # Nstore_plot=Nstore*Nunit                                    # (ug N/ mgC) Nitrogen in storage
    Nconst_protein_plot = Nconst_protein * Nunit * Numbertoarray  # (ug N/ mgC) constant protein pool in nitrogen (193-25)(193-33)
    Nphoto_plot = Nphoto * Nunit                                  # (ug N/ mgC) Photosynthesis related protein nitrogen (193-25)(193-33)
    Nbiosynth_plot = Nbiosynth * Nunit                            # (ug N/ mgC) biosynthesis related protein in N (193-37)
    Ndna_const_plot = Ndna_const * Nunit * Numbertoarray          # (ug N/ mgC) Nitrogen in constant part of DNA
    Nrna_plot = Nrna * Nunit                                      # (ug N/ mgC) Nitrogen in constant part of RNA
    Nchl_plot = Nchl * Nunit                                      # (ug N/ mgC) Chlorophyll nitrogen (actually almost negligiable) (193-33)

    # ==============================
    # Carbon related
    # ==============================
    Cchl = Chl                          # (molC cell-1) carbon in chlorophyll (195-16)
    Cphoto = Nphoto * CNprotein         # (molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth = Nbiosynth * CNprotein   # (molC cell-1) carbon in biosynthesis protein (195-16)
    Crna = Nrna * YrnaC_N               # (molC cell-1) carbon RNA (195-16)
    CthylakoidPG = Pthylakoid * YpgC_P  # (molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    Cnstore = zeros(size(NO3))          # Nstore*YcyanoC_N

    Cother = Qc - Cphoto - Cbiosynth - Cconst_protein - Cchl \
             - Crna - Cdna_const \
             - Cessential - CthylakoidPG - Cnstore

    # =======================================
    # Unit conversions for C
    # =======================================

    percentorratio = 100  # 100: percent, 1:ratio
    Cphoto_plot = Cphoto / Qc * percentorratio
    Cbiosynth_plot = Cbiosynth / Qc * percentorratio
    Cconst_protein_plot = Cconst_protein / Qc * percentorratio * Numbertoarray
    Cchl_plot = Cchl / Qc * percentorratio
    Crna_plot = Crna / Qc * percentorratio * Numbertoarray
    Cdna_const_plot = Cdna_const / Qc * percentorratio * Numbertoarray
    Cother_plot = Cother / Qc * percentorratio
    Cessential_plot = Cessential / Qc * percentorratio * Numbertoarray
    CthylakoidPG_plot = CthylakoidPG / Qc * percentorratio
    Cnstore_plot = Cnstore / Qc * percentorratio

    D[i]=D[i]*86400 #s-1 to day-1

    # return the NO3 concentration and calculated growth rate 返回 NO3浓度 和 计算出的生长速率 D[]
    results_single_cell_model = {"NO3_array": NO3_input,
              "growth_rates": D,
              "Nchl": Nchl,
              "Nphoto": Nphoto,
              "Nbiosynth": Nbiosynth,
              "Nprotein": Nprotein,
              "Nrna": Nrna,
              "Ndna_const": Ndna_const,
              "Nessential": Nessential,
              "Nprotein_const": Nprotein_const,
              "Cchl": Cchl,
              "Cphoto": Cphoto,
              "Cbiosynth": Cbiosynth,
              "Crna": Crna,
              "CthylakoidPG": CthylakoidPG,
              "Cnstore": Cnstore,
              "Cessential": Cessential,
              }
    return results_single_cell_model

    # print("NO3 array:", NO3)
    # print("D array:", D)

#=============================================================
# Population Dynamic Module
def get_single_cell_growth_rate(NO3_input):
    """obtain the single cell growth rate based on the corresponding NO3 根据当前NO3浓度获取单细胞生长速率D"""
    current_NO3 = NO3_input
    results = single_cell_model(current_NO3)
    current_growth_rates = results["growth_rates"]
    print("current_growth_rates:", current_growth_rates)
    current_Nchl = results["Nchl"]
    current_Nphoto = results["Nphoto"]
    current_Nbiosynth = results["Nbiosynth"]
    current_Nprotein = results["Nprotein"]
    current_Nrna = results["Nrna"]
    current_Ndna_const = results["Ndna_const"]
    current_Nessential = results["Nessential"]
    current_Nprotein_const = results["Nprotein_const"]
    current_Cchl = results["Cchl"]
    current_Cphoto = results["Cphoto"]
    current_Cbiosynth = results["Cbiosynth"]
    current_Crna = results["Crna"]
    current_CthylakoidPG = results["CthylakoidPG"]
    current_Cnstore = results["Cnstore"]
    current_Cessential = results["Cessential"]

    current_results = {"current_growth_rates": current_growth_rates,
                       "current_Nchl": current_Nchl,
                       "current_Nphoto": current_Nphoto,
                       "current_Nbiosynth": current_Nbiosynth,
                       "current_Nprotein": current_Nprotein,
                       "current_Nrna": current_Nrna,
                       "current_Ndna_const": current_Ndna_const,
                       "current_Nessential": current_Nessential,
                       "current_Nprotein_const": current_Nprotein_const,
                       "current_Cchl": current_Cchl,
                       "current_Cphoto": current_Cphoto,
                       "current_Cbiosynth": current_Cbiosynth,
                       "current_Crna": current_Crna,
                       "current_CthylakoidPG": current_CthylakoidPG,
                       "current_Cnstore": current_Cnstore,
                       "current_Cessential": current_Cessential,
                       }

    return current_results

def population_dynamics(time, NO3_input, mr, K, initial_cell_abundance):
    """Calculate the population biomass based on the growth rate and the mortality rate基于生长速率和死亡率计算种群数量变化"""
    X = np.zeros(len(time)) # define the dimension of cell abundance
    C_cell = 23e-15         # (mol C cell-1) #Synechococcus C quota median value. ref: Baer et al. 2017 Stoichiometry of Prochlorococcus, Synechococcus, and small eukaryotic populations in the western North Atlantic Ocean
    N_cell = 2.4e-15        # (mol N cell-1)
    X[0] = initial_cell_abundance
    for t in range(1, len(time)):
        # 获取当前NO3对应的单细胞生长速率
        current_results = get_single_cell_growth_rate(NO3_input)
        mu = current_results["current_growth_rates"]
        print("Current growth rate: mu = ")
        growth_rate = mu * X[t-1]
        death_rate = mr * X[t - 1] * (1 + X[t - 1] / K)
        X[t] = X[t - 1] + (growth_rate - death_rate) * dt  # 使用欧拉法更新细胞丰度
        # 打印变量的维度
        print(f"X[t - 1] ndim: {np.ndim(X[t - 1])}")
        print(f"growth_rate ndim: {np.ndim(growth_rate)}")
        print(f"death_rate ndim: {np.ndim(death_rate)}")
        # 确保丰度不为负数
        if X[t] < 0:
            X[t] = 0
        #print("Cell abundance:", X)
        #macromolecular allocation in population
        #Carbon and nitrogen related
        population_results = get_single_cell_growth_rate(NO3_input)
        Nchl = population_results["current_Nchl"]
        Nphoto = population_results["current_Nphoto"]
        Nbiosynth = population_results["current_Nbiosynth"]
        Nprotein = population_results["current_Nprotein"]
        Nrna = population_results["current_Nrna"]
        Ndna_const = population_results["current_Ndna_const"]
        Nprotein_const = population_results["current_Nprotein_const"]
        Nessential = population_results["current_Nessential"]
        Cchl = population_results["current_Cchl"]
        Cphoto = population_results["current_Cphoto"]
        Cbiosynth = population_results["current_Cbiosynth"]
        Crna = population_results["current_Crna"]
        CthylakoidPG = population_results["current_CthylakoidPG"]
        Cnstore = population_results["current_Cnstore"]
        Cessential = population_results["current_Cessential"]

        #population macromolecular composition
        Nchl_population = Nchl * X[t] * C_cell                            # (molN chl cell-1) Chlorophyll N concentration
        Nphoto_population = Nphoto * X[t] * C_cell                        # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
        Nbiosynth_population = Nbiosynth * X[t] * C_cell                  # (molN cell-1) various part of biosynthesis related protein in N (193-37)
        Nprotein_population = Nphoto_population + Nbiosynth_population    # (molN cell-1) All the proteins in N (193-26)
        Nessential_population = Nessential * X[t] * C_cell
        C_phy_population = C_cell * X[t]                                  # Phytoplankon biomass in C content
        Ndna_population = Ndna_const * X[t] * C_cell
        Nrna_population = Nrna * X[t] * C_cell
        #Nprotein_const = Nphoto_const + Nconst_protein  # (molN cell-1) All the proteins in N (193-26)
        #Nrna = Nprotein_const * Cnrna_variable          # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
        #Nrna_D2 = Nprotein_D * Cnrna_variable           # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
        population_results = {"population_Nchl": Nchl_population,
                              "population_Nprotein": Nprotein_population,
                              "population_Nrna": Nrna_population,
                              "population_Ndna": Ndna_population,
                              "population_Nessential": Nessential_population,
                              "cell_abundance": X
                              #"population_Cchl": Cchl_population,
                              #"population_Crna": Crna,
                              #"population_CthylakoidPG": CthylakoidPG,
                              #"population_Cnstore": Cnstore,
                              #"population_Cphoto": Cphoto,
                              #"population_Cessential": Cessential
                              }
    return population_results

def detritus_dynamics(X, mr, NO3_input, time):
    C_cell = 23e-15  # (mol C cell-1) #Synechococcus C quota median value. ref: Baer et al. 2017 Stoichiometry of Prochlorococcus, Synechococcus, and small eukaryotic populations in the western North Atlantic Ocean
    degradation_rate = {
        'protein': 0.015/86400,
        'lipid': 0.015/86400,
        'chlorophyll': 0.015/86400,
        'carbohydrate': 0.015/86400,
        'rna': 0.015/86400,
        'dna': 0.015/86400,
    }

    detritus_macromolecular = population_dynamics(time, NO3_input, mr, K, initial_cell_abundance)
    population_Nprotein = detritus_macromolecular["population_Nprotein"]
    population_Nchl = detritus_macromolecular["population_Nchl"]
    population_Ndna = detritus_macromolecular["population_Ndna"]
    population_Nrna = detritus_macromolecular["population_Nrna"]
    population_Nessential = detritus_macromolecular["population_Nessential"]
    # capacity = K
    # 假设 Nchl 是你要检查的变量
    if isinstance(population_Nchl, np.ndarray):
        print("population_Nchl 是一个 NumPy 数组")
    else:
        print("population_Nchl 是一个单独的数值")

    """基于死亡率计算碎屑生成"""
    Detritus = np.zeros(len(time))
    Detritus_protein = np.zeros(len(time))
    Detritus_chl = np.zeros(len(time))
    Detritus_dna = np.zeros(len(time))
    Detritus_rna = np.zeros(len(time))
    Detritus_protein[0] = population_Nprotein[0]
    Detritus_chl[0] = population_Nchl[0]
    Detritus_dna[0] = population_Ndna
    Detritus_rna[0] = population_Nrna[0]
    for t in range(1, len(time)):
        Detritus[t] = Detritus[t-1] + mr * X[t-1] * (1 + X[t - 1] / K) * C_cell * (time[t] - time[t-1])
        print(f"Step {t}: Detritus_protein = {Detritus_protein[t - 1]}, X = {X[t - 1]}")

        # 蛋白质碎屑的生成和降解
        Detritus_protein[t] = Detritus_protein[t - 1] + mr * X[t - 1] * (1 + X[t - 1] / K) * Detritus_protein[t - 1] * dt \
                              - degradation_rate['protein'] * Detritus_protein[t - 1] * dt

        # 叶绿素碎屑的生成和降解
        Detritus_chl[t] = Detritus_chl[t - 1] + mr * X[t - 1] * (1 + X[t - 1] / K) * Detritus_chl[t - 1] * dt \
                          - degradation_rate['chlorophyll'] * Detritus_chl[t - 1] * dt

        # DNA 碎屑的生成和降解
        Detritus_dna[t] = Detritus_dna[t - 1] + mr * X[t - 1] *  (1 + X[t - 1] / K) * Detritus_dna[t - 1] * dt \
                          - degradation_rate['dna'] * Detritus_dna[t - 1] * dt

        # RNA 碎屑的生成和降解
        Detritus_rna[t] = Detritus_rna[t - 1] + mr * X[t - 1] *  (1 + X[t - 1] / K) * Detritus_rna[t - 1] * dt \
                          - degradation_rate['rna'] * Detritus_rna[t - 1] * dt

        detritus_results = {"detritus_chl": Detritus_chl,
                            "detritus_protein": Detritus_protein,
                            "detritus_dna": Detritus_dna,
                            "detritus_rna": Detritus_rna,
                            "detritus": Detritus,
                            }
    return detritus_results

# 主流程
# Population dynamics simulation parameters
import numpy as np
mr = 0.05 # Death rate (day-1)
initial_cell_abundance = 1000  # Initial cell number: cells/L
time = np.linspace(0, 1, 500)  # Time (days)
dt = time[1] - time[0]
NO3_input = 10  # Current NO3 concentration
K = 1000000

# 1. Simulate population growth
results_population = population_dynamics(time, NO3_input, mr, K, initial_cell_abundance)
X = results_population["cell_abundance"]
results_detritus = detritus_dynamics(X, mr, NO3_input, time)
detritus = results_detritus["detritus"]
detritus_chl = results_detritus["detritus_chl"]
detritus_protein = results_detritus["detritus_protein"]
detritus_dna = results_detritus["detritus_dna"]
detritus_rna = results_detritus["detritus_rna"]

# 可视化种群增长和碎屑生成
plt.figure(figsize=(12, 6))

# 种群数量变化
plt.subplot(2, 2, 1)
plt.plot(time, X, label='Population Size')
plt.xlabel('Time (days)')
plt.ylabel('Cell Abundance')
plt.title('Population Dynamics')

# 碎屑生成
plt.subplot(2, 2, 2)
plt.plot(time, detritus, label='Detritus')
plt.xlabel('Time (days)')
plt.ylabel('Detritus (mol)')
plt.title('Detritus Protein')

plt.tight_layout()
plt.show()
