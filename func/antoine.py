#Antoinova rovnice pro vypocet teploty a tlaku 
import numpy as np

def pressure(T, A, B, C): #Returns pressure in kPa, imput in °C 
    p = 10 ** (A-(B/(C+T)))
    return p*0.1333223684 #Conversion from mmHg to kPa

def temperature(p, A, B, C): #Returns temp in °C, pressure in kPa 
    p = p/0.1333223684 #Conversion from kPa to mmHg
    T = (B/(A-np.log10(p)))-C
    return T

def dewpoint(temp, rel_humidity): #Temp input in °C, returns °C 
    T = ((4030*(temp+235))/((4030-(temp+235)*np.log(rel_humidity))))-235
    return T