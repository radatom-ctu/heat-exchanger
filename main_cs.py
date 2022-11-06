import matplotlib.pyplot as plt
import func.comp as comp
import config as cf
import numpy as np
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
heatflow = massflow*capacity*(cf.gas_ti-dewtemp) # J/h
gas_output = heatflow/3600 # W

water_density = comp.Rho('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water density 
water_capacity = comp.Cp('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water heat capacity 
water_conductivity = comp.Conductivity('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water heat conductivity
water_viscosity = comp.Viscosity('Water', cf.amb_pressure+cf.water_pressure,
                         comp.CtoK(cf.water_ti)) # Get water viscosity

#Given output temperature
water_massflow = heatflow/(water_capacity*(cf.water_to-cf.water_ti))
water_flow = water_massflow/water_density*1000 

#Without output temperature
water_delta = heatflow/water_capacity #kg.K/h

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
print('Hmotnostní tok vody:', round(water_massflow,2), '[kg/h]')
print('Objemový tok vody:', round(water_flow,2), '[l/h]\n')

print('Podtlak způsobený komínovým efektem:', round(chimneypress, 2), '[Pa]\n')

## Flow and heat transfer calculation

gas_prandtl = comp.prandtl(capacity, conductivity, gas_viscosity)
water_prandtl = comp.prandtl(water_capacity, water_conductivity, water_viscosity)

#Asuming tube-in-tube design
if cf.one_pipe:   
    gas_speed = cf.volume_flow/comp.circle(
                cf.input_diameter/1000)/3600 # Gas speed inside tube in m/s
    gas_diffusivity = conductivity/(density*capacity) # m2/s
    gas_reynolds = comp.reynolds(gas_speed, cf.input_diameter/1000, gas_kviscosity)
    gas_prandtl = comp.prandtl(capacity, conductivity, gas_viscosity)
    
    water_speed = (water_flow/(3600*1000))/(comp.circle((
        cf.shell_diameter)/1000)-comp.circle(cf.input_diameter/1000)) # Water speed around tube in m/s
    water_diffusivity = water_conductivity/(water_density*water_capacity) # m2/s
    water_reynolds = comp.reynolds(water_speed, cf.shell_diameter/1000, water_kviscosity) 
    water_prandtl = comp.prandtl(water_capacity, water_conductivity, water_viscosity)
    
    print('\nVÝPOČET TOKU A PŘESTUPU TEPLA')
    print('Rychlost spalin:', round(gas_speed,2), '[m/s]')
    print('Reynoldsovo číslo spalin:', round(gas_reynolds,0), '[-]')
    print('Prandtlovo číslo spalin', round(gas_prandtl,5), '[-]\n')
    
    print('Rychlost vody:', round(water_speed,5), '[m/s]')
    print('Reynoldsovo číslo vody:', round(water_reynolds,0),'[-]')
    print('Prandtlovo číslo vody:', round(water_prandtl,5), '[-]\n')
    
    intube_eqn = comp.choice_flowintube(gas_reynolds, gas_prandtl)
    
    aroundtube_eqn = comp.choice_aroundtube(water_reynolds, water_prandtl)
    print('\nVÝPOČET NUSSELTOVA ČÍSLA')
    print('Předpoklad výměníku typ "trubka v trubce"')
    print('Rovnice pro spaliny:', intube_eqn)
    print('Rovnice pro vodu:', aroundtube_eqn)
    
    gas_nu = []
    water_nu = []
    
    gas_alfa = []
    water_alfa = []
    
    heat_transfer = []
    calculated_length = []
    pressure_exit = []
    pressure_loss = []
    
    #Estimate pipe length
    if intube_eqn == 'leveque' or intube_eqn == 'hausen':
        print('\nIterační výpočet Nusseltova čísla jako funkce délky výměníku.')
        for l in range(cf.min_length,cf.max_length):
            gas_graetz = comp.graetz(cf.input_diameter, l, gas_reynolds, gas_prandtl)
            water_graetz = comp.graetz(cf.shell_diameter, l, water_reynolds, water_prandtl)
            
            if intube_eqn == 'leveque':
                inu = comp.leveque(gas_graetz)
            else:
                inu = comp.hausen(gas_graetz)
            
            lnu = comp.gni(water_graetz, cf.input_diameter+2*cf.tube_thickness, cf.shell_diameter)
            
            gas_nu.append(inu)
            water_nu.append(lnu) 
            
            ialfa = comp.alfa(inu, conductivity,
                              cf.input_diameter/1000)
            lalfa = comp.alfa(lnu, water_conductivity,
                              (cf.input_diameter+2*cf.tube_thickness)/1000)
            gas_alfa.append(ialfa)
            water_alfa.append(lalfa)
            
            heat_k = comp.heat_ktube(ialfa, lalfa, cf.input_diameter,
                                     cf.outer_diameter, cf.conductivity) #W/m2K
            
            heat_transfer.append(heat_k)
            
            needed_l = gas_output/(heat_k*(calctemp-(cf.water_to+
                       cf.water_ti)/2)*comp.circleLen(cf.input_diameter/1000))
            
            calculated_length.append(needed_l*1000) # back to mm
            
            if gas_reynolds < 2300:
                lamb_coef = 64/gas_reynolds   
                ipressloss = lamb_coef*((pow(gas_speed,2)*l)/(2*cf.input_diameter/1000))
            pressure_loss.append(ipressloss)
            pressure_exit.append(chimneypress-ipressloss)      
    else: #other eqns
        pass
    
    # Finding calculated value match
    maxdif = 0.00001
    limit = 5
    mult = 1.01
    done = False
    
    while not done:
        if maxdif > limit:
            break 
        k = 0
        for x in range(cf.min_length,cf.max_length):
            if abs(x-calculated_length[k]) <= maxdif:
                print('\nShoda nalezena pro výpočetní L:', x, 'a vypočtené L:', 
                      round(calculated_length[k],3), '[mm]')
                print('Maximální rozdíl:', round(maxdif, 4), '[mm]')
                print('Sktuečný rozdíl:', round(abs(x-calculated_length[k]),3), '[mm]\n')
                done = True
                break
            k = k+1
        maxdif = maxdif*mult
    
    if done:
        print('\nVLASTNOSTI PRO VYPOČTENÉ L:', round(calculated_length[k],3), '[mm]')
        print('Spaliny Nu:', round(gas_nu[k],3), '[-]')
        print('Voda Nu:', round(water_nu[k],3), '[-]')
        print('Alfa vnitřní', round(gas_alfa[k], 3),'[W/m2.K]')
        print('Alfa vnější:', round(water_alfa[k], 3), '[W/m2.K]')
        print('Součinitel prostupu tepla:', round(heat_transfer[k], 3), '[W/m2.K]')
        print('Tlaková ztráta:', round(pressure_loss[k],2), '[Pa]')
        print('Podtlak na výstupu:', round(pressure_exit[k], 2),'[Pa]\n')
        
        # Re-run the calculation again using the iterated length 
        # the length calculated here should be the same as the iterated length
        if cf.test_length_value:
            gas_graetz = comp.graetz(cf.input_diameter, calculated_length[k],
                                     gas_reynolds, gas_prandtl)
            water_graetz = comp.graetz(cf.shell_diameter, calculated_length[k],
                                     water_reynolds, water_prandtl)
            if intube_eqn == 'leveque':
                inu = comp.leveque(gas_graetz)
            else:
                inu = comp.hausen(gas_graetz)
            
            lnu = comp.gni(water_graetz, cf.input_diameter+2*cf.tube_thickness, cf.shell_diameter)
            
            ialfa = comp.alfa(inu, conductivity,
                              cf.input_diameter/1000)
            lalfa = comp.alfa(lnu, water_conductivity,
                              (cf.input_diameter+2*cf.tube_thickness)/1000)
            
            heat_k = comp.heat_ktube(ialfa, lalfa, cf.diameter,
                                     cf.outer_diameter, cf.conductivity) #W/m2K
            
            needed_l = gas_output/(heat_k*(calctemp-(cf.water_to+
                       cf.water_ti)/2)*comp.circleLen(cf.input_diameter/1000))
            
            print('Přepočtená hodnota délky:', round(needed_l*1000,3),'[mm]')
            print('Rozdíl vůči iterované hodnotě:',
                  abs(round(needed_l*1000-calculated_length[k],3)),'[mm]')
            
        
        print('\nKONEČNÉ ROZMĚRY')
        print('Vnitřní trubka d_i:', cf.input_diameter, '[mm]')
        print('Vnitřní trubka d_e:', cf.input_diameter+2*cf.tube_thickness, '[mm]')
        print('Vnější trubka d_i:', cf.shell_diameter, '[mm]')
        print('Vnější trubka d_e:', cf.shell_diameter+2*cf.shell_thickness, '[mm]')
        print('Délka:', round(calculated_length[k],3), '[mm]')
    else:
        print('\nNenalezena shoda mezi výpočetní a vypočtenou délkou.')
        print('Zkuste změnit minimální a maximální výpočetní délku.')
        print('Min délka:', cf.min_length)
        print('Max délka:', cf.max_length)
           
       
    if cf.plot:
        fig, ax = plt.subplots(3,1,figsize=(8, 8))
        fig.tight_layout(h_pad=5)
        x = list(range(cf.min_length,cf.max_length))
        
        ax[0].plot(x,x,label = 'Výpočetní délka')
        ax[0].plot(x,calculated_length,label = 'Vypočtená délka')
        #ax[0].plot([calculated_length[k], calculated_length[k]], [cf.min_length,cf.max_length],
        #           label='Vybraná délka')
        ax[0].legend(loc = 'upper left')
        ax[0].set_xlabel('L [mm]')
        ax[0].set_ylabel('L [mm]')
        ax[0].set_title('Délka trubky')
        
        
        ax[1].plot([cf.min_length, cf.max_length],[cf.press_min,cf.press_min]
                   , label= 'Min tlak')
        ax[1].plot(x,pressure_exit, label = 'Tlak na výstupu')
        #ax[1].plot([calculated_length[k], calculated_length[k]], [5,chimneypress],
        #          label='Vybraná délka')
        ax[1].legend(loc = 'upper right')
        ax[1].set_xlabel('l [mm]')
        ax[1].set_ylabel('p [Pa]')
        ax[1].set_title('Podtlak na výstupu výměníku')
        
        
        ax[2].plot(x,heat_transfer, label = 'k')
        #ax[2].plot([calculated_length[k], calculated_length[k]], [heat_transfer[0],2],
        #           label='Vybraná délka', color='green')
        ax[2].legend(loc = 'upper right')
        ax[2].set_xlabel('l [mm]')
        ax[2].set_ylabel('k [W/m2.K]')
        ax[2].set_title('Součinitel prostupu tepůa')
        
        
        fig.suptitle('Spalinový výměník')
        plt.subplots_adjust(top=0.9)
        plt.show()
        
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
            

            if n_col == 1: # If not fit on side, fit on diagonal
                n_col = 2              
                pipe_margin = (cf.input_diameter-2*pipe_diameter)/3 # Change margin to reflect the change
            
            for n in range(2,n_col+1): # For each number of pipes
                if n_row!= 1: # Calculate actual margin for each set of pipes
                    pipe_margin = (square_side-((n)*pipe_diameter))/(n+1)
                    n_row = n
                
                n_col = n

                gas_speed = cf.volume_flow/(n*comp.circle((pipe_diameter-2*pipe_thickness)/1000))/3600 
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
                    # Crossectional area of all pipes
                    pipe_crossection_area = ((l/cf.n_partitions)*pipe_diameter)*n # mm2
                    # Final waterflow area
                    water_flow_area = (cf.input_diameter*(
                                       l/cf.n_partitions))-pipe_crossection_area # mm2
                  
                    rey_circumference = ((n_col+1)*((2*(cf.input_diameter-(
                                        n_col*pipe_diameter))/(n_col+1))+2*l/cf.n_partitions))/1000 #m
                    rey_diameter = ((4*comp.mm2tom2(water_flow_area))/(rey_circumference)) # m 

                    water_speed_min = comp.reySpeed(1500, rey_diameter, water_kviscosity) #m/s
                    water_speed_max = comp.reySpeed(pow(10,5), rey_diameter, water_kviscosity) #m/s
                    
                    water_flow_max = water_speed_max*comp.mm2tom2(water_flow_area) #m3/s
                    water_flow_min = water_speed_min*comp.mm2tom2(water_flow_area) #m3/s

                    if water_speed_min<0.2:
                        water_speed_min = 0.2
                    if water_speed_max>0.8:
                        water_speed_max = 0.8
                    
                    calc_speed = (water_speed_max+water_speed_min)/2 # average
                    
                    #Calculate alfa_e
                    c_z = 0.91+0.0125*(n_col-2)
                    
                    sigma_1 = (pipe_margin+pipe_diameter)/pipe_diameter
                    sigma_2 = sigma_1
                    c_s = 1/pow(1+(2*sigma_1-3)*pow(1-sigma_2/2,3),2)
                    
                    alfa_e = 0.2*c_z*c_s*(cf.conductivity/(pipe_diameter/1000))*pow(
                            (calc_speed*(pipe_diameter/1000))/water_kviscosity,0.65)*pow(
                                water_prandtl,0.33)
                    heat_k = comp.heat_ktube(alfa_i,alfa_e,pipe_diameter-2*pipe_thickness, pipe_diameter,cf.conductivity)
                    
                    pipe_area = comp.mm2tom2(comp.cylinder(pipe_diameter, l)*n_col*n_row) #m2
                    deltaT = 100 #k
                    
                    power = pipe_area*deltaT*heat_k
                    
                    results.append([pipe_diameter,pipe_thickness,n_col,n_row,l,
                                    calc_speed,water_speed_min,water_speed_max,alfa_i,alfa_e,heat_k,eqn_choice,power])
                    
    x = list(map(lambda x:x[4], results))[:n_iterations+1] 
    y = list(map(lambda y:y[12], results))[:n_iterations+1]
    
    fig, ax = plt.subplots(3,1,figsize=(8, 8))
    fig.tight_layout(h_pad=5)
    ax[0].plot(x,y,label = '')
    
