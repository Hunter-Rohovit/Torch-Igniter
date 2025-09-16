from rocketcea.cea_obj import add_new_fuel

def add_custom_fuels():
    #the enthalpy of formation and density values are from NIST
    card_str_C3H8_G = """
    fuel C3H8  C 3  H 8  wt%=100.0
        h,cal=-24970.0  t(k)=298.15  rho=0.493
    """
    add_new_fuel("GC3H8", card_str_C3H8_G)

    card_str_H2_G = """
    fuel H2  H 2  wt%=100.0
        h,cal=0.0   t(k)=298.15   rho=8.99E-5
    """
    add_new_fuel("GH2", card_str_H2_G)

    card_str_CH4_G = """
    fuel CH4  C 1  H 4  wt%=100.0
        h,cal=-7480.0   t(k)=298.15   rho=0.668
    """
    add_new_fuel("GCH4", card_str_CH4_G)
