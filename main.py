import matplotlib.pyplot as plt
#import func.antoine as at
import func.comp as comp
import config as cf


## Base heat calculation

density = 0
capacity = 0
conductivity = 0
gas_viscosity = 0

dewtemp = comp.dew((cf.amb_pressure/1000)*cf.composition[1])#Calculate dew point temperature
#at.dewpoint(250, cf.composition[1]) 
calctemp = (cf.gas_ti+dewtemp)/2 #Average temp 

calctempK = comp.CtoK(calctemp)

print('INPUT')
print('Dew point temperature:', round(dewtemp,2), '[°C]')
print('Average temperature of gas:', round(calctemp,2), '[°C]', 
      round(comp.CtoK(calctemp),2), '[K]')
print('Pressure:', cf.amb_pressure, '[Pa]\n')

print('\nGAS PROPERTIES')
for i in range(len(cf.composition_materials)):
    idensity = comp.Rho(cf.composition_materials[i],cf.amb_pressure,
                                 calctempK) # kg/m3
    icapacity = comp.Cp(cf.composition_materials[i],cf.amb_pressure,
                                  calctempK) # J/kg/K
    iconductivity = comp.Conductivity(cf.composition_materials[i],
                                  cf.amb_pressure, calctempK) # W/m/K
    iviscosity = comp.Viscosity(cf.composition_materials[i], cf.amb_pressure,
                                  calctempK) # Visosity in Pa.s
    
    print(cf.composition_materials[i])
    print('Phase is',comp.getPhase(cf.composition_materials[i],
                                   cf.amb_pressure,comp.CtoK(calctemp)))
    print('Density:', round(idensity,2), '[kg/m3]' )
    print('Heat capacity:', round(icapacity,2), '[J/kg/K]')
    print('Volumetric capacity:', round(icapacity*idensity,2), '[J/m3/K]')
    print('Dynamic viscosity:', iviscosity, 'Pa.s')
    print('Thermal conductivity:', round(iconductivity,2), '[W/m/K]\n')
    
    density = density + idensity*cf.composition[i]
    capacity = capacity + icapacity*cf.composition[i]
    conductivity = conductivity + iconductivity*cf.composition[i]
    gas_viscosity = gas_viscosity + iviscosity*cf.composition[i]


massflow = cf.volume_flow*density # kg/h
heatflow = massflow*capacity*(cf.gas_ti-dewtemp) # J/h
gas_output = heatflow/3600 # W

water_density = comp.Rho('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_to-cf.water_ti)) # Get water density 
water_capacity = comp.Cp('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_to-cf.water_ti)) # Get water heat capacity 
water_conductivity = comp.Conductivity('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_to-cf.water_ti)) # Get water heat conductivity
water_viscosity = comp.Viscosity('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_to-cf.water_ti)) # Get water viscosity


water_massflow = heatflow/(water_capacity*(cf.water_to-cf.water_ti))
water_flow = water_massflow/water_density*1000 

print('\nCALCULATED PROPERTIES')
print('Gas density:', round(density,2),'[kg/m3]')
print('Gas capacity:', round(capacity,2), '[J/kg/K]')
print('Gas volumetric capacity:', round(density*capacity,2), '[J/m3/K]')
print('Gas massflow:', round(massflow,2), '[kg/h]')
print('Gas heatflow:', round(heatflow/1000,2), '[kJ/h]')
print('Total gas output:', round(gas_output,2), '[W]\n')

print('Water thermal conductivity:', round(water_conductivity,2), '[W/m/K]')
print('Water viscosity:', water_viscosity, '[Pa.s]')
print('Water density:', round(water_density,2), '[Kg/m3]')
print('Water capacity:',  round(water_capacity,2), '[J/kg/K]')
print('Water massflow:', round(water_massflow,2), '[kg/h]')
print('Water flow:', round(water_flow,2), '[l/h]\n')


## Flow and heat transfer calculation
#Asuming tube-in-tube design

gas_kviscosity = gas_viscosity/density
water_kviscosity = water_viscosity/water_density

gas_speed = cf.volume_flow/comp.circle(
            cf.diameter/1000)/3600 # Gas speed inside tube in m/s
gas_diffusivity = conductivity/(density*capacity) # m2/s
gas_reynolds = comp.reynolds(gas_speed, cf.diameter/1000, gas_kviscosity)
gas_prandtl = comp.prandtl(capacity, conductivity, gas_viscosity)

water_speed = (water_flow/(3600*1000))/(comp.circle((
    cf.shell_diameter)/1000)-comp.circle(cf.diameter/1000)) # Water speed around tube in m/s
water_diffusivity = water_conductivity/(water_density*water_capacity) # m2/s
water_reynolds = comp.reynolds(water_speed, cf.shell_diameter/1000, water_kviscosity) 
water_prandtl = comp.prandtl(water_capacity, water_conductivity, water_viscosity)

print('\nFLOW AND HEAT TRANSFER CALCULATION')
print('Gas speed:', round(gas_speed,2), '[m/s]')
print('Gas Reynolds:', round(gas_reynolds,0), '[-]')
print('Gas Prandtl', round(gas_prandtl,5), '[-]\n')

print('Water speed:', round(water_speed,5), '[m/s]')
print('Water Reynolds:', round(water_reynolds,0),'[-]')
print('Water Prandtl:', round(water_prandtl,5), '[-]\n')

intube_eqn = comp.choice_flowintube(gas_reynolds, gas_prandtl)
# For tube in tube design only one equation for water flow exists
aroundtube_eqn = comp.choice_aroundtube(water_reynolds, water_prandtl)
print('Equation chosen for gas:', intube_eqn)
print('Equation chosen for water:', aroundtube_eqn)

chimneypress = comp.chimPress(density, 1.204, cf.chimneylength/1000)

gas_nu = []
water_nu = []

gas_alfa = []
water_alfa = []

heat_transfer = []
calculated_length = []
pressure_exit = []

#Estimate tube length
if intube_eqn == 'leveque' or intube_eqn == 'hausen':
    print('\nCalculating heat transfer for water and gas as a function of length.')
    for l in range(cf.min_length,cf.max_length):
        gas_graetz = comp.graetz(cf.diameter, l, gas_reynolds, gas_prandtl)
        water_graetz = comp.graetz(cf.shell_diameter, l, water_reynolds, water_prandtl)
        
        if intube_eqn == 'leveque':
            inu = comp.leveque(gas_graetz)
        else:
            inu = comp.hausen(gas_graetz)
        
        lnu = comp.gni(water_graetz, cf.diameter+2*cf.tube_thickness, cf.shell_diameter)
        
        gas_nu.append(inu)
        water_nu.append(lnu) 
        
        ialfa = comp.alfa(inu, conductivity,
                          cf.diameter/1000)
        lalfa = comp.alfa(lnu, water_conductivity,
                          (cf.diameter+2*cf.tube_thickness)/1000)
        gas_alfa.append(ialfa)
        water_alfa.append(lalfa)
        
        k = 0.8*comp.heat_kwall(ialfa, lalfa, cf.tube_thickness/1000, cf.conductivity) #W/m2/K
        heat_transfer.append(k)
        
        needed_l = gas_output/(k*(calctemp-(cf.water_to+
                   cf.water_ti)/2)*comp.circle((cf.diameter+2*cf.tube_thickness))/1000)
        calculated_length.append(needed_l*1000) # back to mm
        if gas_reynolds < 2300:
            lamb_coef = 64/gas_reynolds   
            ipressloss = lamb_coef*((pow(gas_speed,2)*l)/(2*cf.diameter/1000))
        pressure_exit.append(chimneypress-ipressloss)
    print('Done.')      
else: #other eqns
    pass

print('\nEstimating calculated length.')
maxdif = 0.00001
limit = 5
mult = 2
done = False

while not done:
    print('Max difference is now:', maxdif, '[mm]')
    if maxdif > limit:
        break 
    k = 0
    for x in range(cf.min_length,cf.max_length):
        if abs(x-calculated_length[k]) <= maxdif:
            print('\nMatch found at L:', x, 'and calc L:', 
                  round(calculated_length[k],3), '[mm]')
            print('Max difference was:', maxdif, '[mm]')
            print('Actual difference:', round(abs(x-calculated_length[k]),3), '[mm]\n')
            done = True
            break
        k = k+1
    maxdif = maxdif*mult

if done:
    print('\nPROPERTIES FOR L:', round(calculated_length[k],3), '[mm]')
    print('Alfa internal', round(gas_alfa[k], 3),'[W/m2/K]')
    print('Alfa external:', round(water_alfa[k], 3), '[W/m2/K]')
    print('Heat transfer coeficient k:', round(heat_transfer[k], 3), '[W/m2/K]')
    print('Pressure at the exit:', round(pressure_exit[k], 2),'[Pa]\n')
    
    print('\nFINAL DIMENSIONS')
    print('Inner tube diameter_i:', cf.diameter, '[mm]')
    print('Inner tube diameter_e:', cf.diameter+2*cf.tube_thickness, '[mm]')
    print('Outer tube diameter_i:', cf.shell_diameter, '[mm]')
    print('Outer tube diameter_e:', cf.shell_diameter+2*cf.shell_thickness, '[mm]')
    print('Length:', round(calculated_length[k],3), '[mm]')
else:
    print('\nNo match found for selected min and max tube length.')
    print('Consider adjusting values of min and max length.')
    print('Min length:', cf.min_length)
    print('Max length:', cf.max_length)
    
if cf.plot:
    fig, ax = plt.subplots(3,1,figsize=(8, 8))
    fig.tight_layout(h_pad=5)
    x = list(range(cf.min_length,cf.max_length))
    
    ax[0].plot(x,x,label = 'L')
    ax[0].plot(x,calculated_length,label = 'Calculated L')
    ax[0].plot([calculated_length[k], calculated_length[k]], [cf.min_length,cf.max_length],
               label='Selected length')
    ax[0].legend(loc = 'upper left')
    ax[0].set_xlabel('l [mm]')
    ax[0].set_ylabel('l_expected [mm]')
    ax[0].set_title('Tube length')
    
    
    ax[1].plot([cf.min_length, cf.max_length],[cf.press_min,cf.press_min]
               , label= 'min p')
    ax[1].plot(x,pressure_exit, label = 'p at exit')
    ax[1].plot([calculated_length[k], calculated_length[k]], [5,chimneypress],
               label='Selected length')
    ax[1].legend(loc = 'upper right')
    ax[1].set_xlabel('l [mm]')
    ax[1].set_ylabel('p [Pa]')
    ax[1].set_title('Pressure loss')
    
    
    ax[2].plot(x,heat_transfer, label = 'k')
    ax[2].plot([calculated_length[k], calculated_length[k]], [heat_transfer[0],2],
               label='Selected length', color='green')
    ax[2].legend(loc = 'upper right')
    ax[2].set_xlabel('l [mm]')
    ax[2].set_ylabel('k [W/m2/K]')
    ax[2].set_title('Heat transfer coeficient')
    
    
    fig.suptitle('Spalinový výměník')
    plt.subplots_adjust(top=0.9)
    plt.show()