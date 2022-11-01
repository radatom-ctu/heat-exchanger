#INPUT
##Combustion gas information 
composition_materials = ['CarbonDioxide','Water','Oxygen','Nitrogen','Argon'] # Flue gas composition
composition = [0.087,0.097,0.096,0.712,0.0084] # Flue gass composition in volume fractions
gas_ti = 250 #°C - Flue gas input temperature
volume_flow = 20 #Nm3/h - Flue gas volume flow

##Cooling medium information
water_ti = 30 #°C - Water input temperature
water_to = 50 #°C - Water output temperature
water_pressure = 5*pow(10,5) #Pa - Water system gauge pressure

##Other 
press_min = 12 #Pa - Minimal underpressure at the exit of heat exchanger
chimneylength = 5000 #mm - Chimney length
air_temp = -15 #°C - Air temperature
dew_temp_add = 15 #°C or K - Margin against water condensation

#Constants 
amb_pressure = 101325 #Pa - Ambient temperature

#Materials
diameter = 130 #mm - Inner diameter of the heat exchanger pipe
outer_diameter = 140 #mm - Outer diameter of the heat exchanger pipe
tube_thickness = (outer_diameter-diameter)/2 #mm - Thickness of the heat exchanger pipe 
conductivity = 45 #W/mK - Thermal conductivity of the pipe material (carbon steel)

#Iteration variables
min_length = 50 #mm - Minimal pipe length for iteration 
max_length = 500 #mm - Maximal pipe length for iteration

## Tube in tube
shell_diameter = 142 #mm - Outer pipe inner diameter
shell_thickness = 5 #mm - Outer pipe wall thickness

#Options
plot = False # Plot yes/no
test_length_value = True # Rechecks the values after iteration by rerunning the calculation