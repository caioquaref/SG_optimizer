function [ RMSSE ] = GOFAST_run( x,Hz,A,B,C,Me,gc,DVPAW,clb,cub,wt_fact,driver )
% Launches GOFAST and returns the optimization objective, RMSSE.
%==========================================================================
%   x       = the design variables (the parameters changed in the driver file)
%   Hz      = frequency definition of trace
%   A       = road load coefficient A
%   B       = road load coefficient B
%   C       = road load coefficient C
%   Me      = effective mass, lbm
%   gc      = gravitational constant (lbm*ft) / (lbf*s^2)
%   DVPAW   = voltage at maximum pedal
%   clb     = lower bound on the constraints
%   cub     = upper bounds on the constraints
%   wt_fact = weight factors
%   driver  = driver algorithm number
%==========================================================================
global VTi % from main program
global output_data % to main program
global aDi WI_Di DD DDw CED CEDw DELPVS % sent out
persistent niter

if isempty(niter); niter = 0; end
niter = niter + 1;

% Spaces are acceptable in the names of the model_path and lnc_path
GF_path = 'C:\PWT-Tools\GOFAST\GOFAST_2021_07';
model_path = 'C:\PWT-Tools\GOFAST\Data\10766 DT TRX T6HO Glamis Cycle';
lnc_file = '10766 DT TRX T6HO GlamisOptim.lnc';
drv_file = '- DT T6HO Glamis Optimum 0.4 Driver (TEST).drv';

% Write the design variables to the driver file
ModifyDriverFile(x, model_path, drv_file, driver);

% Delete files from the previous run in the _output folder  
out_dir = strrep(lnc_file,'.lnc','_output');
log_file = strcat(model_path,'\',out_dir,'\','*.log');
zip_file = strcat(model_path,'\',out_dir,'\','*.zip');
delete(log_file)
delete(zip_file)

% To avoid a GF license error, the call from MATLAB must be made from the 
% same folder as GOFAST.EXE, so change the directory temporarily and revert
% upon completion.

current_folder = pwd;
cd ..
cd (GF_path)
full_path = char(strcat(GF_path,'\GOFAST.exe'," ",'"',model_path,'\',lnc_file,'"'));
system (full_path) % launch GOFAST
cd (current_folder)

% Extract solutions from the _output folder
out_file = strcat(model_path,'\',out_dir,'\','TH_Veh1_Task5_Iter1.mat');
load (out_file, 'tTH')

time_sim = tTH.time.v;  % sec
%accel_sim = tTH.dvdt_veh.v; % m/sec/sec
vel_sim = tTH.v_veh.v;  % kph
%dist_sim = tTH.s_veh.v; % meters

% DELPVS = DELta Position Voltage Sensor
% DVPAW converts GF pedal position to equivalent DELPVS. This is used
% to calculate the pedal metric.
DELPVS_sim = tTH.pos_acc.v/100 * DVPAW; 

% Save MPG for output at the end of optimization
unadj_MPG = tTH.FEcum_comb.v * 0.6213712 * 3.785412; % from km/l to MPG
output_data(niter,1) = unadj_MPG(end);

% To be consistent with the SAE spec, resample velocity to 10Hz, double 
% smooth it, then re-derive the acceleration and displacement increments 
% accordingly.

sec = 0:1/Hz:time_sim(end);
mph = interp1(time_sim,vel_sim,sec) * 0.6213712;
DELPVS = interp1(time_sim,DELPVS_sim,sec);

% Driven Vehicle Speed (mph)
VDi = movmean(movmean(mph,5),5);	% apply 5 sample moving average twice
VDi(VDi < 0.03/1609.344*3600) = 0;	% zero out velocity below 0.03 m/s

% Driven Acceleration and Displacement Increment
aDi = gradient(VDi*5280/3600,1/Hz);	% feet/sec^2
dDi = VDi*5280/3600 ./ Hz;          % feet

% Driven Accumulated Distance, miles
DD = sum(dDi)/5280;
% NOTE: If you get an error on the next line it is probably because the 
%       end time of the trace definition in the mss file is not the same as 
%       what this program is using. The end times need to be the same.
DDw = sum(dDi .* wt_fact)/5280;

% Driven Road Load
FI_Di = Me*aDi/gc;                  % Driven Inertial Force, lbf
FRL_Di = (A + B*VDi + C*VDi.^2);	% Driven Road Load Force, lbf
FENG_Di = FI_Di + FRL_Di;           % Driven Engine Force, lbf
FENG_Di(FENG_Di<0) = 0;     % zero out negative force, per SAE J2951

% Driven Engine Work Increment (used to calculate the IWR)
WI_Di = (FI_Di .* dDi)*1.355818; % convert ft-lbf to Joules
WI_Di(WI_Di<0)=0; % want only positive values

% Driven Cycle Energy (VDE)
CED = sum(FENG_Di .* dDi)*1.355818; % convert ft-lbf to Joules
CEDw = sum(FENG_Di .* dDi .* wt_fact)*1.355818; % convert ft-lbf to Joules
MJperMile_Dw = CEDw/DDw/1e6;

output_data(niter,2) = MJperMile_Dw; % MJ/mile
nx = size(x,2);
output_data(niter,3:(nx+2)) = x';

% Return the objective
RMSSE = sqrt(sum((VDi .* wt_fact - VTi .* wt_fact).^2)/numel(VDi)); % mph

disp('----- Simulation Update -----')
display(niter)
display(x)
display(RMSSE)
[constraints,ceq] = DriverMetrics(x,clb,cub,Hz,wt_fact);
% NOTE: These represent the most recent simulated values and may not reflect
%       the settings the optimizer will choose for the next simulation.
Sceq = ['Equality constraint = ', num2str(ceq)];
disp(Sceq)
Slc = ['Upper constraint (c - cub <0) = ', num2str(constraints(1:5))];
disp(Slc)
Suc = ['Lower constraint (clb - c <0) = ', num2str(constraints(6:10))];
disp(Suc)
disp('-----------------------------')

end