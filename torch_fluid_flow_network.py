import CoolProp.CoolProp as CP
import cantera as ct
import math
from rocketcea.cea_obj import CEA_Obj
import scipy.optimize as opt
from pathlib import Path
import matplotlib.pyplot as plt


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
orifice_dia_outlet = 0.18; #in
cd = 0.6 #sharp-edged orifice assumption


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

print("Chamber Pressure: ", pressure_chamber/6894.76, "psia")
cea = CEA_Obj(oxName = 'GOX', fuelName = 'C3H8')
temp_chamber = cea.get_Tcomb(Pc = pressure_chamber/6894.76, MR = 3.64)
MW, gamma_chamber = cea.get_Chamber_MolWt_gamma(Pc = pressure_chamber/6894.76, MR = 3.64)

s = cea.get_full_cea_output(pressure_chamber/6894.76, MR = 3.64, eps =1)
print(s)

cea = CEA_Obj(oxName="LOX", fuelName="LH2")
#print(cea.get_Tcomb(Pc=100, MR=6.0))

'''def mass_balance(pressure_chamber):
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

    #bring in NASA CEA to get thermophysical properties of the combustion products
    cea = CEA_Obj(oxName = 'GOX', fuelName = 'Propane')
    Isp, cstar, temp_chamber, MW, gamma_chamber = cea.get_IvacCstrTc_ChmMwGam(pressure_chamber/6894.76, MR)
    rho_chamber = pressure_chamber*MW/(8314.5*temp_chamber) #kg/m^3
    print(mdot_fuel, mdot_oxidizer, MR, temp_chamber)
    pressure_critical_outlet = pressure_chamber*(2/(gamma_chamber+1))**(gamma_chamber/(gamma_chamber-1))

    if pressure_ambient < pressure_critical_outlet:
        mdot_out = cd*orifice_area_out*math.sqrt(gamma_chamber*rho_chamber*pressure_chamber*(2/(gamma_chamber+1))**((gamma_chamber+1)/(gamma_chamber-1)))
    else:
        mdot_out = cd*orifice_area_out*math.sqrt(2*pressure_chamber*rho_chamber*(gamma_chamber/(gamma_chamber-1))*((pressure_ambient/pressure_chamber)**(2/gamma_chamber)-(pressure_ambient/pressure_chamber)**((gamma_chamber+1)/gamma_chamber)))

    #Continuity Equation
    residue = mdot_fuel + mdot_oxidizer - mdot_out
    return residue

sol = opt.root_scalar(mass_balance, bracket=[pressure_ambient+1, min(pressure_fuel, pressure_oxidizer)-1], method='bisect')

if sol.converged:
    Pc_solution = sol.root
    print(f"Chamber Pressure: {Pc_solution/6894.76:.2f} psia")
else:
    print("Solver failed.")'''

# Counter diffusion flame

gas = ct.Solution('gri30.yaml')
gas.TP = gas.T, Pc_solution

width = 0.5 # distance between fuel and oxidizer outlets
loglevel = 1 # 0 suppresses output, 5 produces very detailed output

flame = ct.CounterFlowDiffusionFlame(gas, width=width)

flame.fuel_inlet.mdot = mdot_fuel
flame.fuel_inlet.X = fuel_type
flame.fuel_inlet.T = temp_fuel

flame.oxidizer_inlet.mdot = mdot_oxidizer
flame.oxidizer_inlet.X = oxidizer_type
flame.oxidizer_inlet.T = temp_oxidizer

flame.boundary_emissivities = 0.0, 0.0
flame.radiation_enabled = False

flame.solve(loglevel, auto=True)

flame.show() # shoes current solution

fig, ax = plt.subplots()
plt.plot(flame.grid, f.T)
ax.set_title('Temperature of the flame')
ax.set(ylim=(0,2500), xlim=(0.000, 0.020))
# fig.savefig('./diffusion_flame.pdf')
plt.show()
