%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code for Chan, Anahtar, et al., Nature Nanotechnology, March 2020
% Last revised: May 4, 2020

function [t,final,finalc]=vABN_model()

%Based on the code structure from Kwong et al., PNAS (2015), 112(41):12627. 

clear all;
close all;
clc;

options = odeset('AbsTol',1e-10,'RelTol',1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call to the function that models vABN activity in mice 
% Models the first 2 hours after ABN administration. This can be changed by
% replacing 120 with any number of minutes. 
% Modify the dose by changing the first value in the backets. The default
% vABN dose is 10 uM (e.g. [10, 0, 0, 0, 0]). If you run the code in its
% original form, it will produce a graph of the exhaled breath
% concentration in an infected mouse.
[t,y] = ode15s(@(t,y) varying_parameters(t,y),0:0.1:120, [10, 0, 0, 0, 0], options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call to the function that models predicted vABN activity in humans.
% Uncomment in order to run. 
%[t,y] = ode15s(@(t,y) human(t,y),0:0.1:120, [27.03, 0, 0, 0, 0], options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uM of volatile reporter in the breath collection chamber
Rc=y(:,5);  

% Conversion from uM to ppb 
ExhPPM = Rc*1e-6*1000*24450;
ExhPPB = ExhPPM*1000;

% Graphs the predicted exhaled reporter concentration 
figure()
plot(t, ExhPPB)
xlabel('Time (min)')
ylabel('Parts per billion (ppb)')

% This outputs the exhaled breath predictions from the figure above
% to allow for graphing with an outside program.
[val, index] = max(ExhPPB);
save('ExhPPB') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = varying_parameters(t,y)
% Models vABN activity in MICE 
% This function allows for the user to vary all of the parameters in the
% model. It was used to generate all of graphs in Figure 3 of the main text. 

%%%% Body constants %%%%
Qm = 0.037; % Minute ventilation in mice (L/min)
Vl = 0.131e-3; % Tidal volume in mice (L)
Qmc = Qm / Vl; % Corrected minute volume aka. breathing rate

%%%% Transport rates %%%%
k_np_tissue = 0.05; % Diffusion rate of the vABNs into tissue (1/min), manually fit to in vivo data 
k_np_phago = 6e-4; % Clearance rate of vABNs via macrophages (1/min), from the literature 
k_reporter_tissue = 30.8; %  Diffusion rate of VOC reporters into tissue (1/min), computationally fit to in vivo data 
k_reporter_clear = 28.1; % Clearance rate of VOC reporters into blood (1/min), computationally fit to in vivo data  

%%%% NE Concentration %%%%
NE = 0.0035; % Measured concentration of NE in BALF from infected mice (uM)
% NE = 0; % Assumed concentration of NE in BALF from infected mice (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinetic parameters and partition coefficients for each modeled reporter 
% To toggle between the different reporters/properties, simply comment out 
% the code accordingly. By default, the predicted signal is calculated using 
% PFC1 as the reporter. 

%%%% PFC1 Cleavage kinetics and partition coefficients %%%%
k_cat = 186; % Turnover number for NE (1/min)  
Km = 10.91; % Michaelis constant for NE (uM) 
H_blood_air = 51.26; % VOC reporter tissue:blood partition coefficient
H_tissue_air = 34.49; % VOC reporter tissue:air partition coefficient

%%%% PFC3 Cleavage kinetics and partition coefficients %%%%
% k_cat = 53.4; % Turnover number for NE (1/min)
% Km = 4.42; % Michaelis constant for NE (uM) 
% H_blood_air = 31.50; % VOC reporter tissue:blood partition coefficient
% H_tissue_air = 36.64; % VOC reporter tissue:air partition coefficient

%%%% PFC5 Cleavage kinetics and partition coefficients %%%%
% k_cat = 142.2; % Turnover number for NE (1/min)  
% Km = 56.13; % Michaelis constant for NE (uM) 
% H_blood_air = 18.73; % VOC reporter tissue:blood partition coefficient
% H_tissue_air = 30.82; % VOC reporter tissue:air partition coefficient

%%%% PFC7 Cleavage kinetics and partition coefficients %%%%
% k_cat = 4.2; % Turnover number for NE (1/min) 
% Km = 22.34; % Michaelis constant for NE (uM) 
% H_blood_air = 18.73; % From PFC5 
% H_tissue_air = 30.82; % From PFC5 

%%%% Calculation of the tissue:blood partition coefficient %%%%
% Stays the same regardless of reporter. 
H_tissue_blood = H_tissue_air/H_blood_air; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Cleavage by nonspecific enzymes %%%%
NS_E = 2.62; % Concentration of nonspecific enzymes in the respiratory tissue (uM), computationally fit to in vivo data 
NS_k_cat = k_cat/60; % Turnover number for nonspecific enzymes (1/min), computationally fit to in vivo data 
NS_Km = Km * 35; % Michaelis constant for nonspecific enzymes (uM), computationally fit to in vivo data 

%%%% State Conditions %%%%
C_np_lumen = y(1); % Concentration of nanoparticles in the lumen (uM)
C_np_tissue = y(2); % Concentration of nanoparticles in the tissue (uM)
C_reporter_tissue = y(3); % Concentration of reporters in the tissue (uM)
C_reporter_lumen = y(4); % Concentration of reporters in the lumen (uM)
C_reporter_chamber = y(5); % Concentration of reporters in the chamber (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Differential equations %%%%

%%%% Equation S1 %%%%
dC_np_lumendt = -k_np_tissue*(C_np_lumen-C_np_tissue);

%%%% Equation S2 %%%%
dC_np_tissuedt = k_np_tissue*(C_np_lumen-C_np_tissue)-k_np_phago*C_np_tissue- (k_cat*NE*C_np_tissue)/(Km+C_np_tissue) ... 
- (NS_k_cat*NS_E*C_np_tissue)/(NS_Km+C_np_tissue);

%%%% Equation S6 %%%%
dC_reporter_tissuedt = -k_reporter_tissue*(C_reporter_tissue/H_tissue_air-C_reporter_lumen) -  k_reporter_clear*C_reporter_tissue/H_tissue_blood ...
+ (k_cat*NE*C_np_tissue)/(Km+C_np_tissue)+ (NS_k_cat*NS_E*C_np_tissue)/(NS_Km+C_np_tissue) ;

%%%% Equation S7 %%%%
dC_reporter_lumendt = k_reporter_tissue*(C_reporter_tissue/H_tissue_air - C_reporter_lumen)-Qmc*(C_reporter_lumen-C_reporter_chamber);

%%%% Equation S8 %%%%
dC_reporter_chamberdt = Qmc*(C_reporter_lumen-C_reporter_chamber);

dydt = [dC_np_lumendt; dC_np_tissuedt; dC_reporter_tissuedt;  dC_reporter_lumendt; dC_reporter_chamberdt];
       
return  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = human(t,y)
% Models vABN activity in HUMANS  

%%%% Body constants %%%%
Qm = 7.5; % Minute ventilation in humans (L/min)
Vl = 543.5e-3; % Tidal volume in humans (L)
Qmc = Qm / Vl; % Corrected minute volume aka. breathing rate

%%%% Transport rates %%%%
k_np_tissue = 0.05; % Diffusion rate of the vABNs into tissue (1/min), manually fit to in vivo data 
k_np_phago = 6e-4; % Clearance rate of vABNs via macrophages (1/min), from the literature 
k_reporter_tissue = 30.8; %  Diffusion rate of VOC reporters into tissue (1/min), computationally fit to in vivo data 
k_reporter_clear = 28.1; % Clearance rate of VOC reporters into blood (1/min), computationally fit to in vivo data  

%%%% NE Concentration %%%%
NE = 2.64; % Predicted concentration of human NE in BALF (uM)
% NE = 0.0035; % uM, mouse concentration
% NE = 0 ; 

%%%% PFC1 Cleavage kinetics and partition coefficients %%%%
k_cat = 186; % Turnover number for NE (1/min)  
Km = 10.91; % Michaelis constant for NE (uM) 
H_blood_air = 51.26; % VOC reporter tissue:blood partition coefficient
H_tissue_air = 34.49; % VOC reporter tissue:air partition coefficient

%%%% Calculation of the tissue:blood partition coefficient %%%%
H_tissue_blood = H_tissue_air/H_blood_air; 

%%%% Cleavage by nonspecific enzymes %%%%
NS_E = 2.62; % Concentration of nonspecific enzymes in the respiratory tissue (uM), computationally fit to in vivo data 
NS_k_cat = k_cat/60; % Turnover number for nonspecific enzymes (1/min), computationally fit to in vivo data 
NS_Km = Km * 35; % Michaelis constant for nonspecific enzymes (uM), computationally fit to in vivo data 

%%%% State Conditions %%%%
C_np_lumen = y(1); % Concentration of nanoparticles in the lumen (uM)
C_np_tissue = y(2); % Concentration of nanoparticles in the tissue (uM)
C_reporter_tissue = y(3); % Concentration of reporters in the tissue (uM)
C_reporter_lumen = y(4); % Concentration of reporters in the lumen (uM)
C_reporter_chamber = y(5); % Concentration of reporters in the chamber (uM)

%%%% Equation S1 %%%%
dC_np_lumendt = -k_np_tissue*(C_np_lumen-C_np_tissue);

%%%% Equation S2 %%%%
dC_np_tissuedt = k_np_tissue*(C_np_lumen-C_np_tissue)-k_np_phago*C_np_tissue- (k_cat*NE*C_np_tissue)/(Km+C_np_tissue) ... 
- (NS_k_cat*NS_E*C_np_tissue)/(NS_Km+C_np_tissue);

%%%% Equation S6 %%%%
dC_reporter_tissuedt = -k_reporter_tissue*(C_reporter_tissue/H_tissue_air-C_reporter_lumen) -  k_reporter_clear*C_reporter_tissue/H_tissue_blood ...
+ (k_cat*NE*C_np_tissue)/(Km+C_np_tissue)+ (NS_k_cat*NS_E*C_np_tissue)/(NS_Km+C_np_tissue) ;

%%%% Equation S7 %%%%
dC_reporter_lumendt = k_reporter_tissue*(C_reporter_tissue/H_tissue_air - C_reporter_lumen)-Qmc*(C_reporter_lumen-C_reporter_chamber);

%%%% Equation S8 %%%%
dC_reporter_chamberdt = Qmc*(C_reporter_lumen-C_reporter_chamber);

dydt = [dC_np_lumendt; dC_np_tissuedt; dC_reporter_tissuedt;  dC_reporter_lumendt; dC_reporter_chamberdt];
       
return  
