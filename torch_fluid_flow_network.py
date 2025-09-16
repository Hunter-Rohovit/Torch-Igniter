import CoolProp.CoolProp as CP
import cantera as ct
import math
from rocketcea.cea_obj import CEA_Obj
import scipy.optimize as opt
from pathlib import Path
import matplotlib.pyplot as plt

from custom_fuels import add_custom_fuels
#from custom_fuels import add_custom_fuels

#Inputs
fuel_type = 'C3H8:1'
oxidizer_type = 'O2:1'

pressure_fuel = 120; #psia
pressure_oxidizer = 120; #psia
pressure_ambient = 14.7; #psia
temp_fuel = 298.15; #K
temp_oxidizer = 298.15; #K


orifice_dia_fuel = 0.075; #in
orifice_dia_oxidizer = 0.106; #in

orifice_dia_outlet = 0.35; #in
cd = 0.6 #sharp-edged orifice assumption

#load custom fuels from custom_fuels.py library & initialize CEA
add_custom_fuels()
cea = CEA_Obj(oxName = 'GOX', fuelName = 'GC3H8')

#Define fluid properties from Cantera Library
fuel = ct.Solution('gri30.yaml')
oxidizer = ct.Solution('gri30.yaml')
fuel.TPX = temp_fuel, pressure_fuel*6894.76, fuel_type
oxidizer.TPX = temp_oxidizer, pressure_oxidizer*6894.76, oxidizer_type

rho_fuel = fuel.density #kg/m^3
rho_oxidizer = oxidizer.density #kg/m^3
gamma_fuel = fuel.cp/fuel.cv
gamma_oxidizer = oxidizer.cp/oxidizer.cv

print("Fuel Density [kg/m^3]:", rho_fuel)
print("Oxidizer Density [kg/m^3]:", rho_oxidizer)
print("Fuel Gamma:", gamma_fuel)
print("Oxidizer Gamma:", gamma_oxidizer)

#derived quantities
orifice_area_fuel = math.pi/4*(orifice_dia_fuel*0.0254)**2
orifice_area_oxidizer = math.pi/4*(orifice_dia_oxidizer*0.0254)**2
orifice_area_out = math.pi/4*(orifice_dia_outlet*0.0254)**2
pressure_fuel = pressure_fuel*6894.76 #Pa
pressure_oxidizer = pressure_oxidizer*6894.76 #Pa
pressure_ambient = pressure_ambient*6894.76 #Pa

#guess a chamber pressure
pressure_guess = 50 #psia
pressure_chamber = pressure_guess*6894.76 #Pa

def mass_balance(pressure_chamber):
    #use the critical pressure ratio to decide whether flow is choked or not
    pressure_critical_fuel = pressure_fuel*(2/(gamma_fuel+1))**(gamma_fuel/(gamma_fuel-1))
    pressure_critical_oxidizer = pressure_oxidizer*(2/(gamma_oxidizer+1))**(gamma_oxidizer/(gamma_oxidizer-1))

    if pressure_chamber <= pressure_critical_fuel:
        mdot_fuel = cd*orifice_area_fuel*math.sqrt(gamma_fuel*rho_fuel*pressure_fuel*(2/(gamma_fuel+1))**((gamma_fuel+1)/(gamma_fuel-1)))
    else:     
        mdot_fuel = cd*orifice_area_fuel*math.sqrt(2*rho_fuel*pressure_fuel*(gamma_fuel/(gamma_fuel-1))*((pressure_chamber/pressure_fuel)**(2/gamma_fuel)-(pressure_chamber/pressure_fuel)**((gamma_fuel+1)/gamma_fuel)))
    
    if pressure_chamber <= pressure_critical_oxidizer:
        mdot_oxidizer = cd*orifice_area_oxidizer*math.sqrt(gamma_oxidizer*rho_oxidizer*pressure_oxidizer*(2/(gamma_oxidizer+1))**((gamma_oxidizer+1)/(gamma_oxidizer-1)))
    else:
        mdot_oxidizer = cd*orifice_area_oxidizer*math.sqrt(2*rho_oxidizer*pressure_oxidizer*(gamma_oxidizer/(gamma_oxidizer-1))*((pressure_chamber/pressure_oxidizer)**(2/gamma_oxidizer)-(pressure_chamber/pressure_oxidizer)**((gamma_oxidizer+1)/gamma_oxidizer)))
   
    MR = mdot_oxidizer/mdot_fuel

    #get CEA values
    temp_chamber = cea.get_Tcomb(pressure_chamber/6894.76, MR)
    temp_chamber =temp_chamber*0.555556 #convert from Rankine to Kelvin
    rho_chamber = cea.get_Chamber_Density(pressure_chamber/6894.76, MR)
    rho_chamber = rho_chamber*16.018463 #convert from lbm/ft^3 to kg/m^3
    MW, gamma_chamber = cea.get_Chamber_MolWt_gamma(pressure_chamber/6894.76, MR)
    MW = MW*0.45359237 #convert from lbm/lbmol to kg/kmol

    pressure_critical_outlet = pressure_chamber*(2/(gamma_chamber+1))**(gamma_chamber/(gamma_chamber-1))

    if pressure_ambient < pressure_critical_outlet:
        mdot_out = cd*orifice_area_out*math.sqrt(gamma_chamber*rho_chamber*pressure_chamber*(2/(gamma_chamber+1))**((gamma_chamber+1)/(gamma_chamber-1)))
    else:
        mdot_out = cd*orifice_area_out*math.sqrt(2*pressure_chamber*rho_chamber*(gamma_chamber/(gamma_chamber-1))*((pressure_ambient/pressure_chamber)**(2/gamma_chamber)-(pressure_ambient/pressure_chamber)**((gamma_chamber+1)/gamma_chamber)))

    #Continuity Equation
    residue = mdot_fuel + mdot_oxidizer - mdot_out
    return residue

sol = opt.root_scalar(mass_balance, bracket=[pressure_ambient+1, min(pressure_fuel, pressure_oxidizer)-1], method='bisect')

#Get the mdots of the converged solution
def get_mdot(Pc):
        # same logic as inside mass_balance, but returns mdot_fuel, mdot_oxidizer
        pressure_critical_fuel =  pressure_fuel*(2/(gamma_fuel+1))**(gamma_fuel/(gamma_fuel-1))
        pressure_critical_oxidizer = pressure_oxidizer*(2/(gamma_oxidizer+1))**(gamma_oxidizer/(gamma_oxidizer-1))

        if Pc <= pressure_critical_fuel:
            mdot_fuel = cd*orifice_area_fuel*math.sqrt(gamma_fuel*rho_fuel*pressure_fuel*(2/(gamma_fuel+1))**((gamma_fuel+1)/(gamma_fuel-1)))
        else:
            mdot_fuel = cd*orifice_area_fuel*math.sqrt(2*rho_fuel*pressure_fuel*(gamma_fuel/(gamma_fuel-1))*((Pc/pressure_fuel)**(2/gamma_fuel)-(Pc/pressure_fuel)**((gamma_fuel+1)/gamma_fuel)))
        if Pc <= pressure_critical_oxidizer:
            mdot_oxidizer = cd*orifice_area_oxidizer*math.sqrt(gamma_oxidizer*rho_oxidizer*pressure_oxidizer*(2/(gamma_oxidizer+1))**((gamma_oxidizer+1)/(gamma_oxidizer-1)))
        else:
            mdot_oxidizer = cd*orifice_area_oxidizer*math.sqrt(2*rho_oxidizer*pressure_oxidizer*(gamma_oxidizer/(gamma_oxidizer-1))*((Pc/pressure_oxidizer)**(2/gamma_oxidizer)-(Pc/pressure_oxidizer)**((gamma_oxidizer+1)/gamma_oxidizer)))
        return mdot_fuel, mdot_oxidizer

if sol.converged:
    Pc_solution = sol.root
    print(f"Chamber Pressure: {Pc_solution/6894.76:.3f} psia")
    actual_mdot_fuel, actual_mdot_oxidizer = get_mdot(Pc_solution)
    print(f"Fuel Mass Flow Rate: {actual_mdot_fuel} kg/s")
    print(f"Oxidizer Mass Flow Rate: {actual_mdot_oxidizer:.6f} kg/s")
    print(f"Mixture Ratio: {actual_mdot_oxidizer/actual_mdot_fuel:.5f}")
else:
    print("Solver failed.")

# Counter diffusion flame
gas = ct.Solution('gri30.yaml')
gas.TP = gas.T, Pc_solution

width = 0.75 # distance between fuel and oxidizer outlets, in 
width = width*0.0254 #convert to meters
loglevel = 1 # 0 suppresses output, 5 produces very detailed output

flame = ct.CounterflowDiffusionFlame(gas, width=width)

flame.fuel_inlet.mdot = actual_mdot_fuel
flame.fuel_inlet.X = fuel_type
flame.fuel_inlet.T = temp_fuel

flame.oxidizer_inlet.mdot = actual_mdot_oxidizer
flame.oxidizer_inlet.X = oxidizer_type
flame.oxidizer_inlet.T = temp_oxidizer

flame.boundary_emissivities = 0.0, 0.0
flame.radiation_enabled = False

flame.solve(loglevel, auto=True)

flame.show() # shows current solution

fig, ax = plt.subplots()
plt.plot(flame.grid, flame.T)
ax.set_title('Temperature of the flame')
ax.set(ylim=(0,5000), xlim=(0.000, width))
# fig.savefig('./diffusion_flame.pdf')
plt.show()

# Extinction Strain Rate

# Define a limit for the maximum temperature below which the flame is considered as extinguished and the computation is aborted
temperature_limit_extinction = max(flame.oxidizer_inlet.T, flame.fuel_inlet.T)

# Initialize and solve
print('Creating the initial solution')
flame.solve(loglevel=0, auto=True)

# Define output locations
output_path = Path() / "diffusion_flame_extinction_data"
output_path.mkdir(parents=True, exist_ok=True)

hdf_output = "native" in ct.hdf_support()
if hdf_output:
    file_name = output_path / "flame_data.h5"
    file_name.unlink(missing_ok=True)

def names(test):
    if hdf_output:
        # use internal container structure for HDF
        file_name = output_path / "flame_data.h5"
        return file_name, test
    # use separate files for YAML
    file_name = output_path / f"{test}.yaml".replace("-", "_").replace("/", "_")
    return file_name, "solution"

file_name, entry = names("initial-solution")
flame.save(file_name, name=entry, description="Initial solution", overwrite=True)

exp_d_a = - 1. / 2.
exp_u_a = 1. / 2.
exp_V_a = 1.
exp_lam_a = 2.
exp_mdot_a = 1. / 2.

alpha = [35000.] # Set initial strain rate
delta_alpha = 1000. # Initial relative strain rate increase
delta_alpha_factor = 25. # Factor of refinement of the strain rate increase
delta_alpha_min = .001 # Limit of the refinement: Minimum normalized strain rate increase
delta_T_min = 1  # K, limit of the Temperature decrease

n = 0
n_last_burning = 0

T_max = [np.max(flame.T)]
a_max = [np.max(np.abs(np.gradient(flame.velocity) / np.gradient(flame.grid)))]

while True:
    n += 1
    # Update relative strain rates
    alpha.append(alpha[n_last_burning] + delta_alpha)
    strain_factor = alpha[-1] / alpha[n_last_burning]
    
    # Create an initial guess based on the previous solution
    # Update grid
    # Note that grid scaling changes the diffusion flame width
    flame.flame.grid *= strain_factor ** exp_d_a
    normalized_grid = flame.grid / (flame.grid[-1] - flame.grid[0])
    
    # Update mass fluxes
    flame.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
    flame.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
    
    # Update velocities
    flame.set_profile('velocity', normalized_grid, flame.velocity * strain_factor ** exp_u_a)
    flame.set_profile('spread_rate', normalized_grid, flame.spread_rate * strain_factor ** exp_V_a)
    
    # Update pressure curvature
    flame.set_profile('lambda', normalized_grid, flame.L * strain_factor ** exp_lam_a)
    
    try:
        flame.solve(loglevel=0)
    except ct.CanteraError as e:
        print('Error: Did not converge at n =', n, e)

    T_max.append(np.max(flame.T))
    a_max.append(np.max(np.abs(np.gradient(flame.velocity) / np.gradient(flame.grid))))
    
    if not np.isclose(np.max(flame.T), temperature_limit_extinction):
        n_last_burning = n # Proceed to next strain rate
        file_name, entry = names(f"extinction/{n:04d}")
        flame.save(file_name, name=entry, description=f"Solution at alpha = {alpha[-1]}", overwrite=True)

        print('Flame burning at alpha = {:8.4F}. Proceeding to the next iteration, with delta_alpha = {}'.format(alpha[-1], delta_alpha))
    elif ((T_max[-2] - T_max[-1] < delta_T_min) and (delta_alpha < delta_alpha_min)):
        # If the temperature difference is too small and the minimum relative strain rate increase is reached, save the last, non-burning, solution to the output file and break the loop
        file_name, entry = names(f"extinction/{n:04d}")
        flame.save(file_name, name=entry, overwrite=True, description=f"Flame extinguished at alpha={alpha[-1]}")

        print('Flame extinguished at alpha = {0:8.4F}.'.format(alpha[-1]),
              'Abortion criterion satisfied.')
        break
    else:
        # Procedure if flame extinguished but abortion criterion is not satisfied
        # Reduce relative strain rate increase
        delta_alpha = delta_alpha / delta_alpha_factor

        print('Flame extinguished at alpha = {0:8.4F}. Restoring alpha = {1:8.4F} and '
              'trying delta_alpha = {2}'.format(
                  alpha[-1], alpha[n_last_burning], delta_alpha))

        # Restore last burning solution
        file_name, entry = names(f"extinction/{n_last_burning:04d}")
        flame.restore(file_name, entry)

# Results
file_name, entry = names(f"extinction/{n_last_burning:04d}")
flame.restore(file_name, entry)

print('----------------------------------------------------------------------')
print('Parameters at the extinction point:')
print('Pressure p={0} bar'.format(flame.P / 1e5))
print('Peak temperature T={0:4.0f} K'.format(np.max(flame.T)))
print('Mean axial strain rate a_mean={0:.2e} 1/s'.format(flame.strain_rate('mean')))
print('Maximum axial strain rate a_max={0:.2e} 1/s'.format(flame.strain_rate('max')))
print('Fuel inlet potential flow axial strain rate a_fuel={0:.2e} 1/s'.format(flame.strain_rate('potential_flow_fuel')))
print('Oxidizer inlet potential flow axial strain rate a_ox={0:.2e} 1/s'.format(flame.strain_rate('potential_flow_oxidizer')))
print('Axial strain rate at stoichiometric surface a_stoich={0:.2e} 1/s'.format(flame.strain_rate('stoichiometric', fuel='C3H8')))

plt.figure()
plt.semilogx(a_max, T_max)
plt.xlabel(r'$a_{max}$ [1/s]')
plt.ylabel(r'$T_{max}$ [K]')
plt.show()
