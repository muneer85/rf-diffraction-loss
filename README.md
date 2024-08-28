# Rounded Obstacle Diffraction Loss Model

This MATLAB function provides the diffraction loss calculated using the rounded obstacle diffraction loss model by adding a correction factor to the Knife Edge diffraction model by considering the diffractive obstacle as a rounded obstacle having a radius of curvature instead of the ideal knife-edge assumption. It is built based on the Recommendation ITU-R P.526-15 [1]. This function uses a terrain data model to obtain the elevation profile between the transmitter and the receiver. 

For more details, the reader can refer to the following article [2]. If you use this code in your work, please cite the following article: 

M. M. Al-Zubi and M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for Internet Connectivity of Rural Areas," IEEE Open Journal of the Communications Society, 2023.
https://ieeexplore.ieee.org/document/10306284

## Output Parameters
An example illustrating how to use this function is attached here, i.e., “Example_Main.m”. This function provides the following output parameters: 

|Parameter| Description|
|:--|:--|
|L_fs|free space loss in dB|
|L_KE|Knife edge diffraction loss in dB|
|L_D|rounded obstacle diffraction loss in dB|
|L|total pathloss [L= L_fs + L_D] in dB|
|d1|the distance from the TX antenna to the top of the obstacle in m|
|d2|the distance from the RX antenna to the top of the obstacle in m|
|h|the distance from the top of the obstacle to the LOS line that connects the two antennas in m (h is negative if the obstacle is below this line)|
|R|the obstacle radius of curvature in m|
|F1|1st Fresnel zone radius at obstacle location in m|
|F1_max|maximum 1st Fresnel zone radius at the mid of path in m, i.e., d1 = d2| 
|d|distance between TX and RX in Km|
|Theta|diffraction angle in degree|

## Requirements
-	MATLAB
-	Mapping Toolbox                  
-	Symbolic Math Toolbox
-	Statistics and Machine Learning Toolbox
-	Antenna Toolbox
-	Communications Toolbox
-	Phased Array System Toolbox

## Referncies 
[1] “Propagation by diffraction,” Int. Telecommun. Union (ITU), Geneva, Switzerland, ITU-Recommendation P. 526–15, 2019.

[2] M. M. Al-Zubi and M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for Internet Connectivity of Rural Areas,"IEEE Open Journal of the Communications Society, 2023.

  

