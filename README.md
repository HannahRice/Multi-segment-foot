# Multi-segment-foot
Matlab code function that uses plantar pressure data to obtain the relative force acting under three segments (rearfoot, midfoot and forefoot) of the foot  during running stance and the location of the centre of pressure of each segment in a global motion capture coordinate system


% Hannah Rice August 2025

% Function that uses plantar pressure data to obtain the relative force acting under three segments of the foot 
% during running stance and the location of the centre of pressure of each segment in a QTM (Qualisys) coordinate system

% The function effectively treats the pressure plate as three force plates that represent the
% rearfoot, midfoot and forefoot for each trial, based on the relative location of
% each participant's foot markers during a static trial.

% Function is specific to a (Tekscan Mat 7101E) pressure plate with
% the following dimensions 
% Matrix Height 447.0 mm
% Matrix Width 487.7 mm
% Thickness 0.102 mm



%% File TYPE INFORMATION 
% Input pressure files are .csv files from Tekscan software and represent a matrix of force values for each pressure sensor where
% each row the location along the anterior-posterior plate axis in the direction of running; each column
% represents the location along the medial-lateral plate axis in the direction of running. Each data frame is represented in a new matrix
% separated by empty rows 

% Input static file provides the marker locations from the static trial
% from OpenSim, and presented in the QTM Coordinate system, saved as a .m file. 
% Marker names are called in 80 - 86 of the
% function . 


%% OUTPUT is a structure which contains for each of the three foot segments an n by 3 matrix, where n represents the number of frames 
% during stance and the columns are the force % under that foot region, the
% ML and the AP coordinate of the centre of pressure of that segment in the
% QTM coordinate system

The sample file S117_M.csv demonstrates the file type this code reads and provides users with an example.
