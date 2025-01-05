import numpy as np
import math
import streamlit as st
import warnings
import matplotlib.pyplot as plt

# Define input parameters - use st.session_state to persist between runs
st.session_state.Tc = 20
st.session_state.Tk = st.session_state.Tc + 273.15  # Temperature in Kelvin
st.session_state.Pco2 = 1.0  # Partial Pressure CO2
st.session_state.R = 8.3145  # Gas constant, J/(mol*K)
st.session_state.initialpH = 5.0  # value to start newton iteration from
st.session_state.tol = 0.00001
st.session_state.iters = 0

# Parameters from Table 4 from Jones et al 2024 
params = {
    "H_CO2": [14.0, -1.3341e-2, -558.98, -422580.0],
    "K_w": [-4.098, -3.2452e3, 2.2362e5, -3.9840e7,13.957,-1262.3,856410],
    "K_hyd": [6.6330e-2, 9526.0],
    "K_ca_star": [233.52, 0, -11974, 0, -36.506],
    "K_b1": [-151.18, -0.088696, -1362.3, 0, 27.798],
}

###########################################################
# Calculate all equilibrium constants from Jones et al 2024
########################################################### 
def calculate_constants():
    
    Tk = st.session_state.Tk
    R = st.session_state.R
    # Calculate equilibrium constants based on temperature and pressure.
    
    # Water density calculation (from Table 7) 
    density_water = 753.596 + 1.877748 * Tk - 0.003562 * Tk**2

    # Step 1: Calculate ln(H_CO2) and H_CO2*
    a1, a2, a3, a4 = params["H_CO2"]
    ln_H_CO2 = a1 + a2 * Tk + a3 / Tk + a4 / Tk**2
    #log_H_CO2_star = (-ln_H_CO2 / np.log(10))  # Convert ln to log base 10
    H_CO2_star = np.exp(-ln_H_CO2)

    # Step 2: Calculate log(K_w) and K_w
    a1, a2, a3, a4 , a5, a6, a7= params["K_w"]
    log_K_w = a1 + a2/Tk + a3 / Tk**2 + a4/Tk**3 +(a5+ a6/Tk+ a7/Tk**2)*math.log10(density_water/1000)
    K_w = 10**log_K_w

    # Step 3: Calculate K_hyd
    a1, a2 = params["K_hyd"]
    K_hyd = a1 * np.exp(-a2 / (R * Tk))

    # Step 4: Calculate ln(K_ca_star) and K_ca*
    a1, a2, a3, a4, a5 = params["K_ca_star"]
    ln_K_ca_star = a1 + (a2*Tk)+(a3/Tk)+(a4/Tk**2)+a5*math.log(Tk)
    #log_K_ca_star = ln_K_ca_star / np.log(10)  # Convert ln to log base 10
    K_ca_star = np.exp(ln_K_ca_star)

    # Step 5: Calculate ln(K_b1) and K_b1
    a1, a2, a3, a4, a5 = params["K_b1"]
    ln_K_b1 = a1 + (a2*Tk)+(a3/Tk)+(a4/Tk**2)+a5*math.log(Tk)
    #log_K_b1 = ln_K_b1 / np.log(10)  # Convert ln to log base 10
    K_b1 = np.exp(ln_K_b1)

    KH=H_CO2_star*(1+K_hyd)
    Kca=K_ca_star*(1+(1/K_hyd))
    
    st.session_state.KH = KH
    st.session_state.K_w = K_w
    st.session_state.K_hyd = K_hyd
    st.session_state.Kca = Kca
    st.session_state.K_b1 = K_b1

# this function to produce same column in the excel sheet for all species.
def calculate_columns(pH, P):
    """Calculate concentrations for each species based on pH and constants."""
    
    # Convert pH to [H+]
    H_plus = 10 ** (-pH)

    # Calculate species concentrations
    OH_minus = st.session_state.K_w / H_plus
    CO2=st.session_state.KH*P
    H2CO3_star = st.session_state.K_hyd*CO2 # H2CO3*
    HCO3_minus = (st.session_state.Kca * H2CO3_star) / H_plus
    CO3_2_minus = (st.session_state.K_b1 * HCO3_minus) / H_plus

    return H_plus, OH_minus,CO2, H2CO3_star, HCO3_minus, CO3_2_minus

# function to calculate Ksp
def calc_Ksp():
    Tk = st.session_state.Tk
    logKsp = -59.2385 - 0.041377 * Tk - 2.1963 / Tk + (24.5724 * math.log10(Tk))
    st.session_state.Ksp = 10 ** logKsp

def newtoniteration(pH):
    
    KH = st.session_state.KH
    K_w = st.session_state.K_w
    K_hyd = st.session_state.K_hyd
    Kca = st.session_state.Kca
    K_b1 = st.session_state.K_b1
    Pco2 = st.session_state.Pco2
    Ksp = st.session_state.Ksp

    log10 = math.log(10)
    cH = 10**(-pH)
    dcHdpH = -log10*10**(-pH)
    
    cOH = K_w*10**(pH)
    dcOHdpH = K_w*log10*10**(pH)
    
    cHCO3 = Kca*KH*K_hyd*Pco2*10**(pH)
    dcHCO3dpH = log10*Kca*KH*K_hyd*Pco2*10**(pH)
    
    cCO3 = K_b1*cHCO3*10**(pH)
    dcCO3dpH = K_b1*(dcHCO3dpH*10**(pH)+cHCO3*log10*10**(pH))
    
    # Iron concentration
    cFe = Ksp/cCO3
    dcFedpH = -(Ksp/(cCO3**2))*dcCO3dpH  # Proper derivative using chain rule
    
    # Charge balance equation
    f = cH + 2*cFe - cOH - cHCO3 - 2*cCO3
    dfdpH = dcHdpH + 2*dcFedpH - dcOHdpH - dcHCO3dpH - 2*dcCO3dpH
    dpH = -f/dfdpH
    
    return dpH

# Determine pH using Newton-Raphson iteration
def calc_pH():
    
    pH = st.session_state.initialpH   # initial estimate of pH to start Newton iteration
    iters = 0
    iterating = 1
    while (iterating):
        iters += 1
        if (iters > 100):
            raise ValueError("Newton iterations have not converged")            
        dpH=newtoniteration(pH)
        pH = pH + dpH
        st.write('Iteration {0:4d} pH increment = {1:6.3e}'.format(iters,dpH))
        if (abs(dpH) < st.session_state.tol):
            iterating = 0    
            st.session_state.pH = pH
            st.session_state.iters = iters
         
    st.write('Solution pH = {0:6.3f}'.format(st.session_state.pH))

# Plot concentrations as a function of pH
def plot_concentrations():
    
    # Generate pH range
    pH = np.arange(0, 14, 0.1)
    P = st.session_state.Pco2
    Ksp = st.session_state.Ksp
     
    # Plot out RBF approximation
    plt.figure()
    
    # Calculate species concentrations
    H_plus, OH_minus,Co2, H2CO3_star, HCO3_minus, CO3_2_minus = calculate_columns(pH,P)
    #charge_balance=abs(H_plus-OH_minus-HCO3_minus-2*CO3_2_minus)
    Fe_plus2=Ksp/CO3_2_minus
    charge_balance=abs(H_plus+2*Fe_plus2-OH_minus-HCO3_minus-2*CO3_2_minus)

    #find the min charge balance and its corresponding pH value
    min_index=np.argmin(charge_balance)# to find the min value in charge balance and return the index in array
    min_charge_balance=charge_balance[min_index]# to find the value of change balance in the index 
    min_pH=pH[min_index]# to find the corresbonding value of the pH at the same index in charge balance

    # Plot the results
    plt.figure(figsize=(10, 7))
    plt.plot(pH, H_plus, linestyle='--', color='blue', label='[H+]')
    plt.plot(pH, OH_minus, linestyle='--', color='orange', label='[OH-]')
    plt.plot(pH, np.full_like(pH, Co2), linestyle='-', color='orange', label='[CO2]')
    plt.plot(pH, np.full_like(pH, P), linestyle='-', color='blue', label='PCo2')
    #plt.plot(pH, np.full_like(pH, H2CO3_star), linestyle='-', color='green', label='[H2CO3]')
    plt.plot(pH, HCO3_minus, linestyle='-', color='purple', label='[HCO3-]')
    plt.plot(pH, CO3_2_minus, linestyle='-', color='red', label='[CO3-2]')
    plt.plot(pH, Fe_plus2, linestyle='-', color='yellow', label='[Fe+2]')
    plt.plot(pH, charge_balance , linestyle='--', color='black', label='[Charge balance]')

    #Add the vertical line at the min charge balance
    plt.axvline(x=min_pH, color='red', linestyle=':')
    # Annotate the pH value on the graph
    plt.text(min_pH + 0.2, #x-location 
             min_charge_balance+2.5,#y-axis location
             f'pH = {min_pH:.3f}',# to print the value of pH
             color='black',
             #use bbox function to put the text inside  
              bbox=dict(facecolor='white', edgecolor='black'))

    plt.yscale('log')  # Set y-axis to logarithmic scale
    plt.ylim(10**-9, 10**2)  # Restrict y-axis range

    plt.xlabel('pH')
    plt.ylabel('Concentration (M)')
    tcval = '{0:.2f}'.format(st.session_state.Tc)
    tcstring = str(tcval)
    plt.title('Species Concentrations vs pH for temperature Tc='+tcstring+' degrees Celsius')
    plt.legend(loc='best')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    st.pyplot(plt)    
    
# Main program
st.title("Open Chemistry Calculator")
st.write('First tab shows newton iteration calculation of pH, second tab plots concentrations as a function of pH')

tab1, tab2 = st.tabs(["pH calculations", "Plot of concentrations vs pH"])
with tab1:
    
     # Create the number inputs
     row1 = st.columns([1,1])

     xval1 = row1[0].number_input("Temperature (degrees celsius)",0,50,20)
     st.session_state.Tc = xval1
     st.session_state.Tk = st.session_state.Tc + 273.15  # Temperature in Kelvin
     xval2 = row1[1].number_input("Initial pH",0,14,5)
     st.session_state.initialpH = xval2
     calculate_constants()
     calc_Ksp()
     calc_pH()

with tab2:
     calculate_constants()
     calc_Ksp()
     plot_concentrations()
