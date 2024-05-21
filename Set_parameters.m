%% Simulation and Model of Microgrid AC/DC (Model 14 Bus)

close all; clearvars; clc; 
%% Input:
    set_param('Model_14Bus_Microgrid_R2018','IgnoredZcDiagnostic','none');  
% Demand Function
    %P_nom_gen = (6.034e+05*3 +  2.303e+05*3 + 2.326e05 *3 + 1.446e04*3 + 1.0007e03 *3*2)/1e6;  % Generation Power(MVA)
    %           Main Grid       Diesel G        PV            BESS         Bus DC 
    
    %%%%generacion modificada
    P_nom_gen = (6.034e+05*3 +  2.303e+05*3 + 2.326e05 *3 + 1.446e04*3 + 1.0007e03 *3*2 + 1.446e04*3 + 2.326e05 *3)/1e6;  % Generation Power(MVA)
    %           Main Grid       Diesel G        PV            BESS           Bus DC       BESS2          Aero
    P_demand = [0.37 0.33 0.3 0.35 0.46 0.52 0.58 0.65 0.77 0.89 1 0.96 0.91 0.84 0.72 0.79 0.72 0.81 0.94 0.85 0.73 0.62 0.49 0.41]*P_nom_gen; t_hours = 1:1:24;
    Demand_case = input('Select the case to simulate: (1) Max demand; or (other) Min demand: '); 
    
%% Initialize: 
% ****** SIMULINK - POWERGUI and Microgrid ******

    Fnom=60;  % :System frequency (Hz)
    Phase = 3; % a, b, c
    Ts_Power=1/(33*Fnom)/100; Ts = Ts_Power; %%% Ts_Power=1/(33*Fnom)/100; Ts = Ts_Power; 
    Vnom_prim = 13.8e3;  Vnom_sec = 220; Vnom_ter_AC_DC = 150; Vnom_DC = 300; 
  
% ****** POWER AND ELECTRONIC PV SYSTEM ******
    Ts_Control=10*Ts_Power;
    % POWER PARAMETERS (PV)
        Pnom = 1e6;      % Inverter nominal 3-phase power (VA)
                Vnom_dc = 480;     % Nominal DC link voltage (V)
        Vnom_sec_PV= 0.85*Vnom_dc/2/sqrt(2)*sqrt(3); 
    % INVERTER CHOKE RL [Rpu Lpu]
        RLchoke=[ 0.15/100  0.15 ];  % in pu
        Pbase_sec=Vnom_sec_PV^2/Pnom;
        RL(1)=RLchoke(1)*Pbase_sec; RL(2)=RLchoke(2)*Pbase_sec/(2*pi*Fnom);
    % FILTER C PARAMETERS
        Qc=0.1*Pnom;       
        Pc=Qc/50;           
    % DC LINK ENERGY FOR 3/4 CYCLE OF PNOM
        Ceq= 3/4 * (Pnom/Fnom*2/Vnom_dc^2);
        Clink=Ceq*2;       
    % IGBT BRIDGE PARAMETERS
        Rs=1e6; Cs=inf; Ron=1e-3; Vf=0; Vfd=0;             

%****** DIESEL GENERATOR (DG) ********        
   
 Vnom_Diesel = 2.4e3; % (V) 
 Pgen_Diesel = 0; %3e6;%0; % (W)

%****** POWER MICROGRID PARAMETERS ******

% AC/DC LOADs: 

switch Demand_case
    case 1 % Max demand
        disp('Max demand'); max_min = 1;
%    AC LOADs (Primary GRID)
        S_Load_14 = 1600e3; % (kVA)  
        PF_Load_14 = 0.8; % Power factor 
        P_Load_14 = S_Load_14* PF_Load_14; % (W)
        Ql_Load_14 = S_Load_14 * sin(acos(PF_Load_14));  
        
        S_Load_12 = 800e3; % (kVA)
        PF_Load_12 = 0.8; % Power factor 
        P_Load_12 = S_Load_12* PF_Load_12; % (W)
        Ql_Load_12 = S_Load_12 * sin(acos(PF_Load_12)); 
        
        S_Load_11 = 400e3; % (kVA)
        PF_Load_11 = 0.8; % Power factor 
        P_Load_11 = S_Load_11* PF_Load_11; % (W)
        Ql_Load_11 = S_Load_11 * sin(acos(PF_Load_11)); 
        
        S_Load_10 = 800e3; % (kVA)
        PF_Load_10 = 0.8; % Power factor 
        P_Load_10 = S_Load_10* PF_Load_10; % (W)
        Ql_Load_10 = S_Load_10 * sin(acos(PF_Load_10));
         
        S_Load_9 = 320e3; % (kVA) NONLINEAL LOAD
        PF_Load_9 = 1; % Power factor
        P_Load_9 = S_Load_9* PF_Load_9; % (W)
        R_Load_9 = (Vnom_prim * 1.35)^2 / P_Load_9;
        Ql_Load_9 = S_Load_9 * sin(acos(PF_Load_9));        


%    AC LOADs (Secondary GRID) 
        S_Load_4 = 50e3; % (kVA)
        PF_Load_4 = 0.9; % Power factor 
        P_Load_4 = S_Load_4* PF_Load_4; % (W)
        Ql_Load_4 = S_Load_4 * sin(acos(PF_Load_4)); 
       
        S_Load_3 = 30e3; % (kVA)
        PF_Load_3 = 0.85; % Power factor 
        P_Load_3 = S_Load_3* PF_Load_3; % (W)
        Ql_Load_3 = S_Load_3 * sin(acos(PF_Load_3)); 
        
        S_Load_2 = 40e3; % (kVA)
        PF_Load_2 = 0.9; % Power factor 
        P_Load_2 = S_Load_2* PF_Load_2; % (W)
        Ql_Load_2 = S_Load_2 * sin(acos(PF_Load_2)); 
        
%    DC LOADs
        P_Load_DC = 2e3; % Nominal power (W)
        R_Load_DC = Vnom_DC^2/ P_Load_DC; % Ohms

    otherwise % Min demand
        disp('Min demand'); max_min = 0;
        Min_demand_percent = input('Enter the percentage of the Max demand: ');
%    AC LOADs (Primary GRID)
        S_Load_14 = 1600e3 * Min_demand_percent/100; % (kVA) 
        PF_Load_14 = 0.8; % Power factor 
        P_Load_14 = S_Load_14* PF_Load_14; % (W)
        Ql_Load_14 = S_Load_14 * sin(acos(PF_Load_14)); 
        
        S_Load_12 = 800e3 * Min_demand_percent/100; % (kVA)
        PF_Load_12 = 0.8; % Power factor 
        P_Load_12 = S_Load_12* PF_Load_12; % (W)
        Ql_Load_12 = S_Load_12 * sin(acos(PF_Load_12)); 
        
        S_Load_11 = 400e3 * Min_demand_percent/100; % (kVA)
        PF_Load_11 = 0.8; % Power factor 
        P_Load_11 = S_Load_11* PF_Load_11; % (W)
        Ql_Load_11 = S_Load_11 * sin(acos(PF_Load_11)); 
        
        S_Load_10 = 800e3 * Min_demand_percent/100; % (kVA)
        PF_Load_10 = 0.8; % Power factor 
        P_Load_10 = S_Load_10* PF_Load_10; % (W)
        Ql_Load_10 = S_Load_10 * sin(acos(PF_Load_10)); 
        
        S_Load_9 = 320e3 * Min_demand_percent/100; % (kVA) NONLINEAL LOAD
        PF_Load_9 = 1; % Power factor 
        P_Load_9 = S_Load_9* PF_Load_9; % (W)
        R_Load_9 = (Vnom_prim * 1.35)^2 / P_Load_9;
        Ql_Load_9 = S_Load_9 * sin(acos(PF_Load_9));   
        
%    AC LOADs (Secondary GRID) 
        S_Load_4 = 50e3 * Min_demand_percent/100; % (kVA)
        PF_Load_4 = 0.9; % Power factor 
        P_Load_4 = S_Load_4* PF_Load_4; % (W)
        Ql_Load_4 = S_Load_4 * sin(acos(PF_Load_4)); 
       
        S_Load_3 = 30e3 * Min_demand_percent/100; % (kVA)
        PF_Load_3 = 0.85; % Power factor 
        P_Load_3 = S_Load_3* PF_Load_3; % (W)
        Ql_Load_3 = S_Load_3 * sin(acos(PF_Load_3)); 
        
        S_Load_2 = 40e3 * Min_demand_percent/100; % (kVA)
        PF_Load_2 = 0.9; % Power factor 
        P_Load_2 = S_Load_2* PF_Load_2; % (W)
        Ql_Load_2 = S_Load_2 * sin(acos(PF_Load_2)); 
        
%    DC LOADs
        P_Load_DC = 2e3 * Min_demand_percent/100; % Nominal power (VA)
        R_Load_DC = Vnom_DC^2/ P_Load_DC; % Ohms
end

% TRANSFORMER PARAMETERS:

    Rcc_at = 0.015; Rcc_b = 0.03; 

%       TRANSFORMER (T3) - SLACK Bus (MAIN GRID) 
          Pnom_T3 = 4e6; % (VA) 
          V1_LL_T3 = 69e3; % (V)

%       TRANSFORMER (T1) - Bus6  - Bus5 (Micro-grid)
          Pnom_T1 = 1.5e6 ; % (VA) 

%       TRANSFORMER (T2)- Bus7 - Bus4 (Micro-grid)
          Pnom_T2 = 1.5e6 ; % (VA) 

%       TRANSFORMER - Bus6 (PV)                      
          Pnom_xfo=Pnom;       % (VA)
          TotalLeakage=0.06;   % (pu)
          W1_xfo=  [Vnom_prim TotalLeakage/25/2  TotalLeakage/2];
          W2_xfo=  [Vnom_sec_PV  TotalLeakage/25/2  TotalLeakage/2];
          Rm_xfo=200;          % (pu)
          Lm_xfo=200;          % (pu)

          
%% Lines
%    line length (Km) (Primary GRID)
    Rl = 0.394; Xl = 0.1168;
    dis8 = 2; dis9 = 6; dis10 = 6; dis11 = 3; dis12 = 6; dis13 = 3; dis14 = 2; 

%    line length (Km) (Secondary GRID)
    Rls = 0.198; Xls = 0.1089;
    dis1 = 0.15; dis2 = 0.2; dis3 = 0.15; dis4 = 0.4; dis5 = 0.4; dis6 = 0.4; dis7 = 0.1;  

%% ALGORITM OF PRIMARY CONTROL SYSTEMS 
%****** PANEL PHOTOVOLTAIC (PV) CONTROL SYSTEM ******

% CONTROL PARAMETERS
    % MPPT Control (Perturb & Observe Algorithm)
        Increment_MPPT = 0.01;      
        Limits_MPPT = [ 583 357 ]; 
    % VDC regulator (VDCreg)
        Kp_VDCreg =2; Ki_VDCreg = 400; 
        LimitU_VDCreg = 1.5; LimitL_VDCreg= -1.5;       
    % Current regulator (Ireg)
        RLff(1)= W1_xfo(2) + W2_xfo(2) + RLchoke(1); 
        RLff(2)= W1_xfo(3) + W2_xfo(3) + RLchoke(2); 
            Kp_Ireg= 0.3;Ki_Ireg= 20;           
            LimitU_Ireg= 1.5; LimitL_Ireg= -1.5;     
    % PWM Modulator Parameters 
        Fc= 33 * Fnom ;