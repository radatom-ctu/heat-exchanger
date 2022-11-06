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
input_diameter = 130 #mm - Inner diameter of the heat exchanger pipe
outer_diameter = 140 #mm - Outer diameter of the heat exchanger pipe
tube_thickness = (outer_diameter-input_diameter)/2 #mm - Thickness of the heat exchanger pipe 
conductivity = 35 #W/mK - Thermal conductivity of the pipe material (carbon steel)


#Iteration variables
min_length = 50 #mm - Minimal pipe length for iteration 
max_length = 500 #mm - Maximal pipe length for iteration
step_size = 10 #mm - Iteration step size (must be type int)

## Multiple pipes
pipe_size = [20,21.3,22.0,25.0,26.9,28,31.8,33.7,35.0,38.0,40,42.4,44.5,45.3,51,54] # Pipe sizes
pipe_thickness = [2.6]#,2.9,3.2] # Pipe thickness for size
water_tdelta = 20 # K - Maximal water temperature change wanted range(1,water_tdelta)


#Options
plot = False # Plot yes/no
output_temperature = False # Use water output temperature

## Multiple pipes
multiple_pipes = True # Multiple pipes calculation 
pipe_margin = 5 # mm - Minimal distance between each pipe, equal in x and y
n_partitions = 3 # Number of partitions 

## Tube in tube
one_pipe = False # Tube in tube calculation
test_length_value = True # Rechecks the values after iteration by rerunning the calculation
shell_diameter = 142 #mm - Outer pipe inner diameter
shell_thickness = 5 #mm - Outer pipe wall thickness


