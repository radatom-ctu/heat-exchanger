#INPUT
##Combustion gas information 
composition_materials = ['CarbonDioxide','Water','Oxygen','Nitrogen','Argon']
composition = [0.087,0.097,0.096,0.712,0.0084]
gas_ti = 250 #°C
volume_flow = 20 #Nm3/h

##Cooling medium information
water_ti = 30 #°C
water_to = 50 #°C
water_pressure = 5*pow(10,5) #Pa - overpressure/underpressure

##Other 
press_min = 12 #Pa
diameter = 130 #mm 
chimneylength = 5000 #mm

#Constants 
amb_pressure = 101325 #Pa

#Materials
outer_diameter = 140 #mm
tube_thickness = 5 #mm
conductivity = 45 #W/m/K

#Variables
min_length = 50 #mm
max_length = 500 #mm

## Tube in tube
shell_diameter = 142 #mm
shell_thickness = 5 #mm

#Options
plot = True