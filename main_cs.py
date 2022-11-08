#import matplotlib.pyplot as plt
import func.comp as comp
import config as cf
import math

## Base heat calculation

density = 0
capacity = 0
conductivity = 0
gas_viscosity = 0

dewtemp = comp.dew((cf.amb_pressure/1000)*cf.composition[1])#Calculate dew point temperature
calctemp = (cf.gas_ti+dewtemp+cf.dew_temp_add)/2 #Average temp 
calctempK = comp.CtoK(calctemp)

print('VSTUP')
print('Spaliny')
print('Vstupní teplota:', cf.gas_ti, '[°C]')
print('Objemový tok:', cf.volume_flow, '[Nm3/h]')
print('Kompozice spalin:', cf.composition_materials)
print('Kompozice spalin v obejmových zlomcích:', cf.composition)
print('Tlak:', cf.amb_pressure, '[Pa]')
print('Voda')
print('Vstupní teplota:', cf.water_ti, '[°C]')
print('Výstupní teplota:', cf.water_to, '[°C]')
print('Ostatní')
print('Minimální podtlak na výstupu:', cf.press_min, '[Pa]')
print('Vstupní průměr:', cf.input_diameter, '[mm]')
print('Minimální délka trubky:', cf.min_length, '[mm]')
print('Maximální délka trubky:', cf.max_length, '[mm]')

print('\nVÝPOČET')
print('Teplota rosného bodu spalin:', round(dewtemp,2), '[°C]')
print('Přídavek proti kondenzaci:', cf.dew_temp_add, '[°C]')
print('Střední teplota spalin:', round(calctemp,2), '[°C]', 
      round(comp.CtoK(calctemp),2), '[K]')

print('\nVLASTNOSTI SLOŽEK SPALIN')
print('MODUL COOLPROP')
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
    print('Fáze: ',comp.Phase(cf.composition_materials[i],
                                   cf.amb_pressure,comp.CtoK(calctemp)))
    print('Hustota:', round(idensity,2), '[kg/m3]' )
    print('Tepelná kapacita:', round(icapacity,2), '[J/kg.K]')
    print('Tepelná kapacita objemová:', round(icapacity*idensity,2), '[J/m3.K]')
    print('Dynamický viskozita:', iviscosity, '[Pa.s]\n')
    
    density = density + idensity*cf.composition[i]
    capacity = capacity + icapacity*cf.composition[i]
    conductivity = conductivity + iconductivity*cf.composition[i]
    gas_viscosity = gas_viscosity + iviscosity*cf.composition[i]


massflow = cf.volume_flow*density # kg/h
heatflow = massflow*capacity*(cf.gas_ti-(dewtemp+cf.dew_temp_add)) # J/h
gas_output = heatflow/3600 # W

water_density = comp.Rho('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water density 
water_capacity = comp.Cp('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water heat capacity 
water_conductivity = comp.Conductivity('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water heat conductivity
water_viscosity = comp.Viscosity('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water viscosity

air_density = comp.Rho('air', cf.amb_pressure, comp.CtoK(cf.air_temp))
chimneypress = comp.chimPress(density, air_density, cf.chimneylength/1000)

gas_kviscosity = gas_viscosity/density
water_kviscosity = water_viscosity/water_density

print('\nVYPOČTENÉ VLASTNOSTI')
print('Hustota spalin:', round(density,2),'[kg/m3]')
print('Tepelná kapacita spalin:', round(capacity,2), '[J/kg.K]')
print('Tepelná kapacita spalin objemová:', round(density*capacity,2), '[J/m3.K]')
print('Hmotnostní tok spalin:', round(massflow,2), '[kg/h]')
print('Výkon:', round(heatflow/1000,2), '[kJ/h]')
print('Výkon:', round(gas_output,2), '[W]\n')

print('Tepelná vodivost vody:', round(water_conductivity,2), '[W/m.K]')
print('Dynamická viskozita vody:', water_viscosity, '[Pa.s]')
print('Hustota vody:', round(water_density,2), '[Kg/m3]')
print('Tepelná kapacita vody:',  round(water_capacity,2), '[J/kg.K]')

print('Podtlak způsobený komínovým efektem:', round(chimneypress, 2), '[Pa]\n')

## Flow and heat transfer calculation

gas_prandtl = comp.prandtl(capacity, conductivity, gas_viscosity)
water_prandtl = comp.prandtl(water_capacity, water_conductivity, water_viscosity)
        
#Asuming multiple pipes design
if cf.multiple_pipes:
    results = []
    n_iterations = math.floor((cf.max_length-cf.min_length)/cf.step_size)
    for i,pipe_diameter in enumerate(cf.pipe_size): # For each pipe size
        for k,pipe_thickness in enumerate(cf.pipe_thickness): # For each pipe thickness    
            pipe_margin = cf.pipe_margin # Create local pipe_margin variable
            square_side = (cf.input_diameter)/math.sqrt(2) # Side of a square inside the circle
            # Maximum n pipes in col/row
            n_col = math.floor((square_side-cf.pipe_margin)/(pipe_diameter+cf.pipe_margin)); 
            n_row = n_col
            

            if n_col == 1: # If pipees dont fit on side, fit on diagonal
                n_col = 2              
                pipe_margin = (cf.input_diameter-2*pipe_diameter)/3 # Change margin to reflect the change
            
            for n in range(2,n_col+1): # For each number of pipes
                if n_row!= 1: # Calculate actual margin for each set of pipes
                    pipe_margin = (square_side-((n)*pipe_diameter))/(n+1)
                    n_row = n
                
                n_col = n
                n_pipes = n_col*n_row
                
                gas_speed = cf.volume_flow/(n_pipes*comp.circle((pipe_diameter-2*pipe_thickness)/1000))/3600 
                gas_reynolds = comp.reynolds(gas_speed,pipe_diameter/1000,gas_kviscosity)
                for l in range(cf.min_length, cf.max_length+1,cf.step_size): # For each length
                    eqn_choice = comp.choice_flowintube(gas_reynolds, gas_prandtl)
                    # Calculate Nusselt number based on flow characteristics
                    if eqn_choice == 'leveque':
                        gas_graetz = comp.graetz((pipe_diameter-2*pipe_thickness)/1000,
                                                 l/1000, gas_reynolds, gas_prandtl)
                        gas_nu = comp.leveque(gas_graetz)
                    elif eqn_choice == 'hausen':
                        gas_graetz = comp.graetz((pipe_diameter-2*pipe_thickness)/1000,
                                                 l/1000, gas_reynolds, gas_prandtl) 
                        gas_nu = comp.hausen(gas_graetz)
                    elif eqn_choice == 'whitaker':
                        gas_nu = comp.whitaker(gas_reynolds, gas_prandtl)
                    else:
                        gas_nu = comp.colburn(gas_reynolds, gas_prandtl)
                    
                    #Calculate alfa_i
                    alfa_i = comp.alfa(gas_nu, conductivity, (pipe_diameter-2*pipe_thickness)/1000) # W/m2.K
                    
                    #water_speed_min - Minimal water speed inbetween pipes for Re > 1500
                    #water_speed_max - Maximal water speed inbetween pipes for Re < 10^5
                    
                    #water_flow_min - Minimal volumetric waterflow - m3/s
                    #water_flow_max - Maximal volumetric waterflow - m3/s
                                                       
                    # Using a rectangle as an area of flow to estimate speed
                    # Crossectional area of all pipes in a row/col
                    pipe_crossection_area = ((l/(cf.n_partitions+1))*pipe_diameter)*n # mm2
                    # Final waterflow area
                    water_flow_area = (cf.input_diameter*(
                                       l/(cf.n_partitions+1)))-pipe_crossection_area # mm2
                  
                    rey_circumference = ((n_col+1)*((2*(cf.input_diameter-(
                                        n_col*pipe_diameter))/(n_col+1))+2*l/(cf.n_partitions+1)))/1000 #m
                    rey_diameter = ((4*comp.mm2tom2(water_flow_area))/(rey_circumference)) # m 

                    water_speed_min = comp.reySpeed(1500, rey_diameter, water_kviscosity) #m/s
                    water_speed_max = comp.reySpeed(pow(10,5), rey_diameter, water_kviscosity) #m/s
                    
                    water_flow_max = water_speed_max*comp.mm2tom2(water_flow_area) #m3/s
                    water_flow_min = water_speed_min*comp.mm2tom2(water_flow_area) #m3/s

                    # Water speed selection - affects the ration of massflow/dT
                    if water_speed_min<0.2:
                        water_speed_min = 0.01
                    if water_speed_max>0.8:
                        water_speed_max = 0.8
                    
                    calc_speed = water_speed_min
                    
                    #Calculate alfa_e - Dlouhý, T: Výpočty kotlů a spalinových výměníků - str. 98
                    c_z = 0.91+0.0125*(n_col-2)
                    
                    sigma_1 = (pipe_margin+pipe_diameter)/pipe_diameter
                    sigma_2 = sigma_1
                    c_s = 1/pow(1+(2*sigma_1-3)*pow(1-sigma_2/2,3),2)
                    
                    alfa_e = 0.2*c_z*c_s*(cf.conductivity/(pipe_diameter/1000))*pow(
                            (calc_speed*(pipe_diameter/1000))/water_kviscosity,0.65)*pow(
                                water_prandtl,0.33)
                    
                    heat_k = comp.heat_ktube(alfa_i,alfa_e,pipe_diameter-2*pipe_thickness, pipe_diameter,cf.conductivity) # W/m2K
                    pipes_area = comp.mm2tom2(comp.cylinder(pipe_diameter, l)*n_col*n_row) # m2 - Area of pipes 
                    
                    water_massflow = ((comp.mm2tom2(water_flow_area)
                                      *calc_speed)*water_density)*3600 # kg/h - Water massflow
                    water_deltaT = heatflow/(water_massflow*water_capacity) # K - Water heating dT
                    temp_change = (calctemp - (cf.water_to+water_deltaT)/2) # K - Calculated dT between pipe and water
                    
                    exchanger_power = pipes_area*heat_k*temp_change # W - Calculated exchanger power
                    
                    if gas_reynolds < 2300: # Pressure loss coeficient in pipe
                        lamb_coef = 64/gas_reynolds   
                    else:
                        lamb_coef = 0.06
                            
                    press_loss_pipe = lamb_coef*((density*pow(gas_speed,2)*(
                                 l/1000))/(2*(pipe_diameter-2*pipe_thickness)/1000)) # Pressure loss in pipe
                    press_loss_local = cf.local_press_mult*density*(pow(gas_speed,2)/2) # Local pressure loss in pipe entrance 
                    press_loss = (press_loss_pipe + press_loss_local)
                    
                    results.append([pipe_diameter,pipe_thickness,n_col,n_row,l,
                                    pipes_area,
                                    calc_speed,water_speed_min,water_speed_max,
                                    alfa_i,alfa_e,heat_k,comp.mm2tom2(water_flow_area),
                                    water_massflow,water_deltaT,temp_change,exchanger_power,press_loss])
                    results_notes=["Průměr trubky (mm)","Velikost stěny (mm)","N trubek ve sloupci",
                                   "N trubek v řadě","Délka trubek (mm)","Plocha trubek (mm2)",
                                   "Výpočtová rychlost proudění vody (m/s)",
                                   "Minimální možná rychlost vody (m/s)","Maximální možná rychlost (m/s)",
                                   "alfa_i (W/m2K)","alfa_e (W/m2K)","k (W/m2K)","Plocha prouodění vody (m2)"
                                   ,"Hmotnostní took vody (kg/h)","Voda dT (K)","Teplotní rozdíl (K)",
                                   "Výkon výměníku (W)","Tlaková ztráta (Pa)"]
                    
        selected_results = []
        for exchanger in results: #filter results
            if abs(exchanger[16]-gas_output) <= 25: #power filter
                if chimneypress-exchanger[17] >= 12: #pressure loss filter
                    selected_results.append(exchanger)   
