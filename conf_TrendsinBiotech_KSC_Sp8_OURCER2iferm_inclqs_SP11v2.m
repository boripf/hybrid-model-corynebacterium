% Configuration file for definition of the SAN model

clc
clear EQ model p u P x y
%% Define model parameters
p.qsmax.Value = 0.50976% max substrate uptake rate [gS/gX/h]
p.KS.Value = 0.079651 % Ks value for glucose [g/L]
p.Yxs.Value = 0.37405 % Yield coefficent for Glucose [gX/gS]
p.Ypx.Value = 1369.8477 % Yield coefficient for Product [gP/gX]
% p.kd.Value = 0.001 % constant - cell death[1/h]?
p.lag.Value = 0.0040574

p.maxpO2.Value = 0.2095   % for %dO2 calculation 0.2095   0.605

%p.kLaO2.Value = 300 % h-1
p.alpha1.Value = 0.8303 %0.7117
p.alpha2.Value = 0.2956 %0.0657
p.beta.Value = 0.4575 %0.4298
p.K.Value = 0.0123 %0.0215
p.stirrerdiameter.Value = 0.046 % m
p.crossarea.Value = 0.01 %m2
%p.cO2max.Value = 0.000239 % mol/L
p.YO2X.Value = 0.038 %mol/g

%p.cCO2max.Value = 0.00116 %mol/L
p.YCO2X.Value = 0.037 %mol/g

p.y_CO2_in.Value = 0.032 % %0.032        0.020   0.016
p.y_O2_in.Value = 20.95 % %20.95        50.03   60.05
p.y_O2_O2in.Value = 99.9

p.epsilon.Value = 0.2
p.T0.Value = 298.15 %K
p.H_CO2_standard.Value = 3.4*10^(-2)  %mol/(L atm) - sp�ter dann 1/x um auf L atm/ mol zu kommen
p.deltaH_R_CO2.Value = 2400 % von Sven �bernommen (calc_offgas_corrSDA)

p.H_O2_standard.Value = 1.3*10^(-3) 
p.deltaH_R_O2.Value = 1700 % von Sven �bernommen (calc_offgas_corrSDA

p.R.Value = 8.314462 % kg*m2/s2/mol/K

p.kLafactor.Value = 0.80

%p.kplus.Value = 211.32      %[1/h];Kappelmann
%p.kminus.Value = 19.152       %[1/h];Kappelmann
p.k1.Value = 129.6
p.k2.Value = 41760000
p.Kw.Value = 10^(-13.8)

% Define model constant
p.c_glu.Value = 350 % conc. of glucose in Feed [g/L]

%p.a.Value = 0.0123 % exponential correlation exponent for viscosity (data driven!) 

%% Model equations
% Feed
EQ.F_in = 'EQ.F_in'
% Volume
EQ.Volume = 'EQ.Volume'
% Substrate uptake rate qs
EQ.qs = 'p.qsmax.Value * ( x(3) / (x(3) + p.KS.Value)) * (1/(1 + EQ.IPTG *(exp(x(10) * p.lag.Value))))'
% Biomass growth rate
EQ.qx = 'EQ.qs * p.Yxs.Value'
%IPTG on/off
EQ.IPTG = 'EQ.IPTG'
% Glucose in
% EQ.Glc_in = 'EQ.Glc_in'
% Glucose metabolized
%EQ.Smet = 'EQ.Smet'

%Temperature --> Kelvin
EQ.T = 'EQ.T + 273.15'
% Temperature dependence of Henry constant - O2
EQ.H_O2 = '1./[p.H_O2_standard.Value*exp(p.deltaH_R_O2.Value*(1./EQ.T -1/p.T0.Value))]'
% Temperature dependence of Henry constant - CO2
EQ.H_CO2 = '1./[p.H_CO2_standard.Value*exp(p.deltaH_R_CO2.Value*(1./EQ.T -1/p.T0.Value))]'

% Airflow import + convert from NL/min --> m3/h
EQ.Fair = 'EQ.Fair * 60/1000'
% Oxygenflow import + convert from NL/min --> m3/h
EQ.FO2 = 'EQ.FO2 * 60/1000'
%sum of Gasflow
EQ.Fgas = 'EQ.FO2 + EQ.Fair'

%pressure --> UMRECHNEN Pascal
EQ.p = 'EQ.p *100000'

%Gas hold up - Vair
EQ.Vair = 'p.epsilon.Value * x(1)'

EQ.yO2_out = '(x(6) * 101325 / EQ.p)*100'
EQ.yCO2_out = '(x(8) * 101326 / EQ.p)*100'

EQ.pH = 'EQ.pH'
EQ.K1 = 'exp(-11.582 - (918.9/EQ.T))'
EQ.K2 = 'EQ.K1/p.Kw.Value'

% Stirrer speed
EQ.stirrerspeed = 'EQ.stirrerspeed'
%kLa O2
EQ.kLaO2 = 'p.K.Value * ((EQ.stirrerspeed*60)^p.alpha1.Value) * (p.stirrerdiameter.Value^p.alpha2.Value) * ((EQ.Fgas/p.crossarea.Value)^p.beta.Value)'
% OTR mol/(Lh)
EQ.OTR = 'EQ.kLaO2 * (x(6)/EQ.H_O2 - x(5))'
% OUR mol/(Lh)
EQ.OUR = 'EQ.qx * p.YO2X.Value * x(2)'
% OUR mol/h
EQ.OUR2 = 'EQ.OUR * EQ.Volume *(-1)'
% dO2 %
EQ.dO2 = 'x(5)* 100 / (p.maxpO2.Value/EQ.H_O2)'   

%kLa CO2
EQ.kLaCO2 = 'EQ.kLaO2 *p.kLafactor.Value'
% CTR mol/(Lh)
EQ.CTR = '- EQ.kLaCO2 * (x(8)/EQ.H_CO2 - x(7))'
%CER mol/(Lh)
EQ.CER = 'EQ.qx * p.YCO2X.Value * x(2)'
% CER mol/h
EQ.CER2 = 'EQ.CER * EQ.Volume'



%% Differential equations
% dV/dt --> Volume
x.V_L.dxdt = 'EQ.F_in' 
% dX/dt --> Biomass
x.X_g.dxdt = 'x(2) *  EQ.qx - x(2) * EQ.F_in / x(1) ' 
% dS/dt --> Substrate
x.S_g.dxdt = 'EQ.F_in / x(1) * ( p.c_glu.Value - x(3) ) - EQ.qs * x(2)' 
% dP/dt --> Product (growth related)
x.P_g.dxdt = 'EQ.qx * p.Ypx.Value * EQ.IPTG - x(4) * EQ.F_in / x(1)'


%ddO2/dt --> dissolved oxygen
x.dO2_molperL.dxdt = 'EQ.OTR - EQ.OUR' 
%dpO2/dt
x.pO2.dxdt = '(((-EQ.OTR*x(1))-(EQ.Fgas * x(6)*101325/p.R.Value/EQ.T) + (EQ.p * p.y_O2_in.Value/100 * EQ.Fair/p.R.Value/EQ.T) + (EQ.p * p.y_O2_O2in.Value/100 * EQ.FO2/p.R.Value/EQ.T)) * p.R.Value * EQ.T/(EQ.Vair/1000))/101325'

%ddCO2/dt --> dissolved carbondioxide
x.dCO2_molperL.dxdt = 'EQ.CER - EQ.CTR - ((p.k1.Value + p.k2.Value*10^(EQ.pH-14)) *x(7)) + (((p.k2.Value/EQ.K2) + (p.k1.Value/EQ.K1)*10^(-EQ.pH)) *x(9))'
%dpCO2/dt - partial pressure of CO2 in gas phase
x.pCO2.dxdt = '(((EQ.CTR*x(1))-(EQ.Fgas * x(8)*101325/p.R.Value/EQ.T) + (EQ.p * p.y_CO2_in.Value/100 * EQ.Fair/p.R.Value/EQ.T)) * p.R.Value * EQ.T/(EQ.Vair/1000))/101325'
%ddHCO3/dt - dissolved HCO3
x.dHCO3.dxdt = '((p.k1.Value + + p.k2.Value*10^(EQ.pH-14) ) * x(7)) - (((p.k2.Value/EQ.K2) + (p.k1.Value/EQ.K1)*10^(-EQ.pH)) * x(9))'

% dSmet/dt --> metabolized Substrate (qs related)
x.Smet_g.dxdt = 'EQ.qs * x(2) * x(1)'
%% Model outputs
y.V.eq = 'EQ.Volume'
y.V.Unit = 'L'

y.X.eq = 'x(2)'
y.X.Unit = 'g/L'

y.S.eq = 'x(3)'
y.S.Unit = 'g/L'

y.P.eq = 'x(4)'
y.P.Unit = 'g/L'

y.qs.eq = 'EQ.qs'
y.qs.Unit = 'g/[g/L]'

y.qx.eq = 'EQ.qx'
y.qx.Unit = '1/h'

y.IPTG.eq = 'EQ.IPTG'
y.IPTG.Unit = '-'

y.Smet.eq = 'x(10)'
y.Smet.Unit = 'g'

y.OUR2.eq = 'EQ.OUR2'
y.OUR2.Unit = 'mol/h'

y.OTR.eq = 'EQ.OTR'
y.OTR.Unit = 'mol/(Lh)'

y.dO2.eq = 'EQ.dO2' % equation --> in %; x bzw. dxdt in mol/L
y.dO2.Unit = '%'

y.CER2.eq ='EQ.CER2'
y.CER2.Unit = 'mol/h'

y.dCO2.eq = 'x(7)' 
y.dCO2.Unit = 'mol/L'

y.CTR.eq = 'EQ.CTR'
y.CTR.Unit = 'mol/(Lh)'

y.kLaO2.eq = 'EQ.kLaO2'
y.kLaO2.Unit = 'h-1'

y.CER.eq ='EQ.CER'
y.CER.Unit = 'mol/(Lh)'

y.dO2mol.eq = 'x(5)' 
y.dO2mol.Unit = 'mol/L'

y.pCO2mol.eq = 'x(8)' 
y.pCO2mol.Unit = 'atm'

y.pO2.eq = 'x(6)' 
y.pO2.Unit = 'atm'

y.H_O2.eq = 'EQ.H_O2' 
y.H_O2.Unit = 'Latm/mol'

y.H_CO2.eq = 'EQ.H_CO2' 
y.H_CO2.Unit = 'Latm/mol'

y.OUR.eq = 'EQ.OUR'
y.OUR.Unit = 'mol/Lh'

y.yO2_out.eq = 'EQ.yO2_out'
y.yO2_out.Unit = '%'

y.yCO2_out.eq = 'EQ.yCO2_out'
y.yCO2_out.Unit = '%'

y.dHCO3.eq = 'x(9)' 
y.dHCO3.Unit = 'mol/L'

%y.app.eq = 'x(1)^p.a.Value'
%y.app.Unit = 'cm/s'
%%

model.Name = 'TrendinBiotech_KSC_IPTG';          % name of the model
model.Type = 'ode';                     % only ode is supported
model.p = p;                            % model parameters
model.EQ = EQ;                          % equations of the model
model.x = x;                            % states & differential equations & information
model.y = y;                            % output equations & information

% define the inputs required by the model
model.u_model = {'F_in', 'Volume','IPTG'}; 

% names of the input data timevariables
model.u_data = { 'FeedGlc_Lh', 'ReactorVolume_hampel','IPTG'};           


model.recompile = 0;                   % 0 or 1 - recompile model every time?
model.solvertype = '15s';              % solver type for matlab solving or failover of C, supports: 23, 15s
model.C_filename = '';                 % currently not used
model.simmode = 'memory';              % can be Matlab or memory
model.interp ='interpfast';             % selct fast or slow matlab buildin interpolation


%% Testsimulation %% get the inputs out of B data file
u_names = model.u_data;
u = B.getvars(u_names);

% convert units of inputs from days to hours
for i=1:length(u_names)
   u.(u_names{i}).Time = u.(u_names{i}).Time .* 24;
end

%% simulation

t = 0:1:27;
xinit = [1.6 0.3 22.8 0];

tic
    [t,y1,x1]=BPTT.tools.model_sim_conf(u,xinit,t,model,1);
toc
%%

subplot(1,2,1)
plot(t,x1(:,1))

subplot(1,2,2)
plot(t,x1(:,2))

% legend(fieldnames(model.x))
%%

subplot(2,3,1)
plot(t,y1(5,:))
title('qs trajectory')
ylabel(y.qs.Unit);
xlabel('Time [h]');

subplot(2,3,2)
plot(t,y1(6,:))
title('� trajectory')
ylabel(y.qx.Unit);
xlabel('Time [h]');

subplot(2,3,3)
plot(t,y1(1,:))
title('Volume')
ylabel(y.V.Unit);
xlabel('Time [h]');

subplot(2,3,4)
hold on
plot(t,y1(2,:))
B.plot('DCW','*')
hold off

legend('CDW Model', 'CDW Experiment')
title('Biomass CDW')
ylabel(y.X.Unit);
xlabel('Time [h]');

subplot(2,3,5)
hold on
plot(t,y1(3,:))
B.plot('cGlc','*')
hold off

legend('S Model', 'S Experiment')
title('S conc.')
ylabel(y.S.Unit);
xlabel('Time [h]');

subplot(2,3,6)
title('Feed rate')
ylabel('L/h');
xlabel('Time [h]');

hold on
plot(t,y1(1,:),'*')
B.plot('WeightFeedGlc_feedrate_hampel')
hold off

legend('Feed model', 'Feed experiment')