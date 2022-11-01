from CoolProp.CoolProp import PropsSI, PhaseSI
import numpy as np


# Temp conversions
def KtoC(temp): # K to °C
    return temp-273.15
def CtoK(temp): # °C to K
    return temp+273.15

# Dew temp https://en.citizendium.org/wiki/Water_dew_point
def dew(pressure): #partial pressure in kPa returns °C
    return KtoC((1668.21/(7.09171-np.log10(pressure)))+45.15)
pass


# Coolprop
def Cp(material,pressure,temp): # Returns specific heat at const pressure J/kg/K
    return PropsSI('C','P',pressure,'T',temp,material)
def Rho(material, pressure, temp): # Returns density at temp and pressure kg/m3
    return PropsSI('D','T',temp,'P',pressure,material) 
def Viscosity(material, pressure, temp): # Returns viscosity at temp and pressure Pa.s
    return PropsSI('V','T',temp,'P',pressure,material)
def Conductivity(material, pressure, temp): # Returns conductivity W/m/K
    return PropsSI('L','T',temp,'P',pressure,material)
def Phase(material, pressure, temp): # Returns phase of material at spec. temp an press.
    return PhaseSI("P", pressure, "T", temp, material)

# Area 
def circle(diameter): # Circle area in mm2
    return ((np.pi)*pow(diameter,2))/4 
def cylinder(diameter,length): # Cylinder wrap area in mm2
    return circle(diameter)*length
def cylinderLength(diameter,area): # Cylinder length from area in mm 
    return (area*4)/(np.pi*pow(diameter,2))


# Numbers
def reynolds(speed,diameter,k_viscosity): # Returns Reynolds number
    return (speed*diameter)/k_viscosity
def prandtl(capacity,conductivity,d_viscosity): # Returns Prandtl number
    return (capacity*d_viscosity)/conductivity
def graetz(diameter,length,reynolds,prandtl): # Returns Graetz number
    return (diameter/length)*reynolds*prandtl 

    
# Chimney effect pressure
def chimPress(density_in, density_out, length): # Returns chimney effect pressure
    return length*9.81*(-density_in+density_out) #  input kg/m3, m - output Pa


# Debug 
def DtoKvisc(viscosity, density): # Dynamic to Kinematic viscosity
    return viscosity/density


# Eqn choince
def choice_flowintube(reynolds, prandtl): # Equation choice for flow in pipe
    if reynolds > pow(10,4) and prandtl < 17000 and prandtl > 0.5:
        return 'colburn'
    elif reynolds > 2300 and reynolds < pow(10,5) and prandtl > 0.48 and prandtl < 592:
        return 'whitaker'
    elif reynolds < 2100 and prandtl > 0.5 and prandtl < 17000:
        return 'leveque'
    elif reynolds < 2300:
        return 'hausen'
def choice_aroundtube(reynolds, prandtl):
    if reynolds < 2300:
        return 'gnielinski'
    else:
        return 'flow_turbulent'


# Heat equations for flow in tube
def colburn(reynolds, prandtl): # Colurn equation for heat transfer in pipe
    return 0.023*pow(reynolds,0.8)*pow(prandtl,1/3)
def whitaker(reynolds, prandtl): # Whitaker -||-
    return 0.015*pow(reynolds,0.83)*pow(prandtl,0.42)
def leveque(gr): # Leveque -||-
    return 1.86*pow(gr,1/3)
def hausen(gr): # Hausen -||-
    return 3.66*((0.0668*graetz)/(1+(0.04*pow(gr,2/3))))


# Heat equation in between tubes
def gni(gr, diameter_in, diameter_out): # Gnielinski eqn. for heat transfer
    k = diameter_in/diameter_out
    return 3.66+1.2*pow(k,0.8)+(0.19*(1+0.14*pow(k,0.5))*
                                pow(gr,0.8))/(1+0.117*pow(gr,0.467))


## Heat transfer 
def alfa(nuselt,conductivity,diameter): # Heat transfer W/m2K
    return (nuselt*conductivity)/diameter
def heat_kwall(alfa1, alfa2, wall_thickness,conductivity): # Heat transfer for wall W/m2K
    return 1/((1/alfa1)+(wall_thickness/conductivity)+(1/alfa2))
def heat_ktube(alfa1, alfa2, d1, d2, conductivity): # Heat transfer for pipe W/m2K
    r1 = d1/2
    r2 = d2/2
    return 1/((1/alfa1)+((r2/conductivity)*np.log(r2/r1))+(r2/(alfa2*r1)))      