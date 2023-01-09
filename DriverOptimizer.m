% GOFAST Driver Optimizer
% ========================================================================= 
% =========================================================================
%                  Finds the Best Driver Settings for the
%                           GOFAST Driver Model
%                     Using Constrained Minimization
% =========================================================================
%           Calculations are made in accordance with SAE J2951
% =========================================================================
%
%   Milton Hubbard
%   November 2019
%   Auburn Hills, MI USA
% =========================================================================
clear all
format shortEng
format compact

% Before running this program it's a good idea to clear the cache folder:
% C:\Users\<user>\AppData\Local\Temp\<user>

global aTi VTi WI_Ti DT DTw CET CETw % sent out
global DR ER EER ASCR IWR CSI PEDAL % from frunction DriverMetrics
global output_data DD DDw CED CEDw % from function GOFAST_run

% Specify the drive cycle to be optimized:
%    UDDS, EPA75, UDDS_4bag, HWFET, US06, WLTCb, RDE38, CRCITY, generic

% NOTE: Only one fuel economy maneuver in the launch file is allowed.

cycle = "Glamis";
DVPAW = 3.84; % = "Delta Voltage Pedal At WOT"

% Set the fuel economy ABCs and mass

%A = 40;             % lbf, SAE J2951 test case at 4000 TW
%B = 0.4;            % lbf/mph, SAE J2951 test case at 4000 TW
%C = 0.02;           % lbf/mph^2, SAE J2951 test case at 4000 TW

% WS T6SO SWB 4x4
A=51.98;
B=0.3753;
C=0.02925;
ETW = 6000; % vehicle test weight, lbm

A=244.86;
B=0.3876;
C=0.21006;
ETW = 6300; % vehicle test weight, lbm

gc = 32.174;        % (lbm*ft) / (lbf*s^2)
Me = ETW * 1.015;	% effective mass

% Set the lower and upper bound constraints on the driver metrics

clb(1) = - 0.4;	cub(1) = - clb(1);	% EER (symetrical)
clb(2) = - 1.5;	cub(2) = 1.0;		% ASCR (asymetrical)
clb(3) = - 1.5;	cub(3) = - clb(3);	% IWR (symetrical)

if cycle == "EPA75" || cycle == "UDDS" || cycle == "UDDS_4bag"
    clb(4) = -15;   cub(4) = 10;        % CSI
    clb(5) = 100;   cub(5) = 170;       % PEDAL
elseif cycle == "CRCITY"
    clb(4) = -30;   cub(4) = 30;        % CSI
    clb(5) = 50;   cub(5) = 280;       % PEDAL
else
    clb(4) = -15;   cub(4) = 15;        % CSI
    clb(5) = 100;   cub(5) = 250;       % PEDAL        
end

% Specify the driver algorithm

driver = 0.4; % BE SURE TO USE THE SAME SETTING IN THE LAUNCH FILE!

% Specify the initial values (x0), lower bounds (lb), and upper bounds (ub)
% for the design variables. 

if driver == 0.4
  
    % HO settings
    x0(1) = 0.1;	lb(1)=0.01;	ub(1) = 0.15;	% Driver resolution
    x0(2) = 0.12;	lb(2)=0.01;	ub(2) = 1;	% Low pass filter
    x0(3) = 0.20;	lb(3)= 0.;	ub(3) = 2;	% Min accel pedal position [%]
    x0(4) = 0.20;	lb(4)= 0.;	ub(4) = 2.;	% Integral Fading Rate [dec%]    
    x0(5) = 3.7;    lb(5)= 1.;	ub(5)= 13.; % Integral gain for delta speed control (gas pedal)
    x0(6) = 3.6;	lb(6)= 1.;	ub(6)= 13.; % Proportional gain for delta speed control (gas pedal)
    x0(7) = 4.8;    lb(7)= 1.;	ub(7)= 13.; % Integral gain for delta speed control (brake pedal)
    x0(8) = 7.4;    lb(8)= 1.;	ub(8)= 18.; % Proportional gain for delta speed control (brake pedal)
     
    % CR settings
    x0(1) = 0.3;	lb(1)=0.01;	ub(1) = 1;	% Driver resolution
    x0(2) = 0.1;	lb(2)=0.01;	ub(2) = 1;	% Low pass filter
    x0(3) = 1;      lb(3)= 0;	ub(3) = 5;	% Min accel pedal position [%]
    x0(4) = 0.1;    lb(4)= 0;	ub(4) = 2.;	% Integral Fading Rate [dec%]    
    x0(5) = 7;      lb(5)= 1.;	ub(5)= 14.; % Integral gain for delta speed control (gas pedal)
    x0(6) = 7;      lb(6)= 1.;	ub(6)= 14.; % Proportional gain for delta speed control (gas pedal)
    x0(7) = 7;      lb(7)= 1.;	ub(7)= 14.; % Integral gain for delta speed control (brake pedal)
    x0(8) = 7;      lb(8)= 1.;	ub(8)= 18.; % Proportional gain for delta speed control (brake pedal)

    % Gen-2 SO settings US06
    x0(1) = 0.05;	lb(1)=0.01;	ub(1) = 0.1;% Driver resolution
    x0(2) = 0.3;	lb(2)=0.01;	ub(2) = 1;	% Low pass filter
    x0(3) = 0.3;    lb(3)= 0;	ub(3) = 2;	% Min accel pedal position [%]
    x0(4) = 0.3;    lb(4)= 0;	ub(4) = 2.;	% Integral Fading Rate [dec%]    
    x0(5) = 4;      lb(5)= 1.;	ub(5)= 10.; % Integral gain for delta speed control (gas pedal)
    x0(6) = 4;      lb(6)= 1.;	ub(6)= 11.; % Proportional gain for delta speed control (gas pedal)
    x0(7) = 5;      lb(7)= 1.;	ub(7)= 11.; % Integral gain for delta speed control (brake pedal)
    x0(8) = 6;      lb(8)= 1.;	ub(8)= 14.; % Proportional gain for delta speed control (brake pedal)

    % Gen-4 SO settings US06
    x0(1) = 0.05;	lb(1)=0.01;	ub(1) = 0.1;% Driver resolution
    x0(2) = 0.35;	lb(2)=0.01;	ub(2) = 1;	% Low pass filter
    x0(3) = 1.1;	lb(3)= 0.;	ub(3) = 2;	% Min accel pedal position [%]
    x0(4) = 0.30;	lb(4)= 0.;	ub(4) = 2.;	% Integral Fading Rate [dec%]    
    x0(5) = 3.6;    lb(5)= 1.;	ub(5)= 10.; % Integral gain for delta speed control (gas pedal)
    x0(6) = 1.5;    lb(6)= 1.;	ub(6)= 11.; % Proportional gain for delta speed control (gas pedal)
    x0(7) = 5;      lb(7)= 1.;	ub(7)= 11.; % Integral gain for delta speed control (brake pedal)
    x0(8) = 7;      lb(8)= 1.;	ub(8)= 14.; % Proportional gain for delta speed control (brake pedal)
    
elseif driver == 0.5
 
    x0(1) = 0.06;	lb(1)= 0.01;ub(1) = 0.9;% Proportional Gain [-]
    x0(2) = 0.06;	lb(2) = 0.;	ub(2) = 0.9;% Integral Gain [-]
    x0(3) = 0.3;	lb(3) = 0.;	ub(3) = 1;	% Feed forward contribution [kph]
    x0(4) = 0.7;    lb(4) = 0.;	ub(4) = 1;	% Integral preservation [%/100]
    x0(5) = 0.65;   lb(5) = 0.;	ub(5) = 1;	% Predictive contribution [%/100]
    x0(6) = 1.5;    lb(6) = 1.;	ub(6) = 5;	% Prediction time [sec]    

elseif driver == 0.6
 
    x0(1) = 0.06;	lb(1) = 0.01;	ub(1) = 0.9;	% Proportional Gain [-]
    x0(2) = 0.06;	lb(2) = 0.0;	ub(2) = 0.9;	% Integral Gain [-]
    x0(3) = 0.3;	lb(3) = 0.;     ub(3) = 1;	% Feed forward contribution [kph]
    x0(4) = 0.7;	lb(4) = 0.;     ub(4) = 1;	% Integral preservation [%/100]
    x0(5) = 0.65;	lb(5) = 0.;     ub(5) = 1;	% Predictive contribution [%/100]
    x0(6) = 1.5;	lb(6) = 1.;     ub(6) = 5;	% Prediction time [sec]
    
    x0(7) = 0.1;	lb(7) = 0.01;	ub(7) = 1;      % Driver resolution
    x0(8) = 0.1;	lb(8) = 0.01;	ub(8) = 1;      % Low pass filter
    x0(9) = 0.0;	lb(9) = 0.0;	ub(9) = 0.7;	% Min accel pedal position [%]    
    
    %x0(10) = 0.2;	lb(10) = 0.01;  ub(10) = 0.9;	% Proportional gain (pedal)
    %x0(11) = 0.7;	lb(11) = 0.01;	ub(11) = 0.9;	% Integral gain (pedal)
    %x0(12) = 0.2;	lb(12) = 0.01;	ub(12) = 0.9;	% Proportional gain (brake)
    %x0(13) = 0.2;	lb(13) = 0.01;	ub(13) = 0.9;	% Integral gain (brake) 

elseif driver == 0.7 % "Alpha Driver"
        
    x0(1) = 0.27;	lb(1) = 0.01;	ub(1) = 0.9;	% Proportional Gain [-]
    x0(2) = 0.03;	lb(2) = 0.01;	ub(2) = 0.9;	% Integral Gain [-]
    x0(3) = 0.5;    lb(3) = 0.;     ub(3) = 1.0;    % Pedal gain map
    x0(4) = 0.1;    lb(4) = 0.;     ub(4) = 0.5;	% Integral reset speed [kph]    
    x0(5) = 0.4;	lb(5) = 0.;     ub(5) = 1.;     % Brake gain map    

elseif driver == 0.8 % "PC Toolbox"
    
    x0(1) = 0.35;	lb(1) = 0.01;	ub(1) = 0.9;	% Proportional Gain [-]
    x0(2) = 0.03;	lb(2) = 0.01;	ub(2) = 0.9;	% Integral Gain [-]
    x0(3) = 0.35;   lb(4) = 0.01;	ub(4) = 1.0;    % Pedal gain map
    x0(4) = 0.18;   lb(5) = 0.;     ub(5) = 0.5;	% Integral reset speed [kph]    
    x0(5) = 0.007;	lb(6) = 0.;     ub(6) = 0.1;	% Feed forward contribution
    x0(6) = 0.36;	lb(7) = 0.;     ub(7) = 0.5;	% Anti-windup gain

end    

% Load the official time-speed definitions of the target cycle traces

load Drive_CyclesFull_v4.mat
Hz = 10;
switch cycle
    case "UDDS"
        sec = UDDS_10Hz(:,1)';
        mph = UDDS_10Hz(:,2)';
    case "EPA75"
        sec = EPA75_10Hz(:,1)';
        mph = EPA75_10Hz(:,2)';
    case "UDDS_4bag"
        sec = UDDS_4bag_10Hz(:,1)';
        mph = UDDS_4bag_10Hz(:,2)';
    case "HWFET"
        sec = HWFET_10Hz(:,1)';
        mph = HWFET_10Hz(:,2)';        
    case "US06"
        sec_1Hz = US06_1Hz(:,1)';
        mph_in_1Hz = US06_1Hz(:,2)';
        sec = 0:1/Hz:sec_1Hz(end); % convert to 10Hz
        mph = interp1(sec_1Hz,mph_in_1Hz,sec);
    case "WLTCb"
        sec_1Hz = WLTCb_1Hz(:,1)';
        kph_1Hz = WLTCb_1Hz(:,2)';
        checksum=sum(kph_1Hz);
        assert(checksum ~= 83758.6);
        mph_in_1Hz = kph_1Hz*0.6213712;
        sec = 0:1/Hz:sec_1Hz(end); % convert to 10Hz
        mph = interp1(sec_1Hz,mph_in_1Hz,sec);
    case "CRCITY"
        sec_1Hz = CRCITY_1Hz(:,1)';
        mph_in_1Hz = CRCITY_1Hz(:,2)';
        sec = 0:1/Hz:sec_1Hz(end); % convert to 10Hz
        mph = interp1(sec_1Hz,mph_in_1Hz,sec);
    case "RDE38"
        sec_1Hz = RDE_38sec_1Hz(:,1)';
        kph_1Hz = RDE_38sec_1Hz(:,2)';
        mph_in_1Hz = kph_1Hz*0.6213712;
        sec = 0:1/Hz:sec_1Hz(end); % convert to 10Hz
        mph = interp1(sec_1Hz,mph_in_1Hz,sec);
    case "Glamis"
        sec_1Hz = Glamis_1Hz(:,1)';
        mph_in_1Hz = Glamis_1Hz(:,2)';
        sec = 0:1/Hz:sec_1Hz(end); % convert to 10Hz
        mph = interp1(sec_1Hz,mph_in_1Hz,sec);        
    otherwise
        disp("*** Unrecognized drive cycle ***");
        stop
end

% To be consistent with Emissions Summary Report, apply the 3-bag fuel 
% economy weighting factors to the driver metrics for the EPA75 3-bag cycle
% only. Application of the weights will cause the reported ratio of cycle 
% energy to distance to equal that of the unweighted 2-bag cycle. The GOFAST 
% simulation should be run using the 3-bag trace definition in this case.

wt_fact = ones(size(sec));
if cycle == "EPA75"
    wt_fact(sec<=505)=0.43;
    wt_fact(sec>1369)=0.57;
end

% When the US06 is being optimized, if desired the optimization can be 
% limited to the highway segment only by setting the weighting for the 
% city bags to zero.
if cycle == "US06"
    wt_fact(sec<130)=1.;  % affects bag 1 (the first city portion)
    wt_fact(sec>495)=1.;  % affects bag 3 (the second city portion)
end
    
% Target Vehicle Speed (mph)
VTi = movmean(movmean(mph,5),5);    % apply 5 sample moving average twice
VTi(VTi < 0.03/1609.344*3600) = 0;  % zero out velocity below 0.03 m/s

% Target Acceleration and Displacement Increment
aTi = gradient(VTi*5280/3600,1/Hz);	% feet/sec^2
dTi = VTi*5280/3600 ./ Hz;          % feet

% Target Accumulated Distance, miles
DT = sum(dTi)/5280;   
DTw = sum(dTi .* wt_fact)/5280;

% Target Road Load
FI_Ti = Me*aTi/gc;                  % Target Inertial Force, lbf
FRL_Ti = (A + B*VTi + C*VTi.^2);	% Target Road Load Force, lbf
FENG_Ti = FI_Ti + FRL_Ti;           % Target Engine Force, lbf
FENG_Ti(FENG_Ti<0) = 0;     % zero out negative force, per SAE J2951

% Target Engine Work Increment (used to calculate the IWR)
WI_Ti = (FI_Ti .* dTi)*1.355818; % convert ft-lbf to Joules
WI_Ti(WI_Ti<0)=0; % want only positive values

% Target Cycle Energy (VDE)
CET = sum(FENG_Ti .* dTi)*1.355818; % convert ft-lbf to Joules
MJperMile_T = CET/DT/1e6;
CETw = sum(FENG_Ti .* dTi .* wt_fact)*1.355818; % convert ft-lbf to Joules
MJperMile_Tw = CETw/DTw/1e6; % this should agree with the STATS VDE report

% Set optimizer options

AA = [];	% linear inequality constraints <none used>
bb = [];    % linear inequality constraints <none used>
Aeq = [];   % linear equality constraints <none used>
beq = [];   % linear equality constraints <none used>

options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',20,...
    'FiniteDifferenceStepSize',0.015,'OptimalityTolerance',4.e-4,...
    'ConstraintTolerance',0.001,'ScaleProblem','obj-and-constr',...
    'StepTolerance', 4.e-4,'Display','iter-detailed',...
     'TypicalX',[0.05 0.35 1.1 0.3 3.6 1.5 5 7],...  
    'PlotFcn',{@optimplotfval,@optimplotconstrviolation,@optimplotx});

% Objective function
objfunc = @(x)GOFAST_run(x,Hz,A,B,C,Me,gc,DVPAW,clb,cub,wt_fact,driver);

% Nonlinear constraint function
nlcfunc = @(x)DriverMetrics(x,clb,cub,Hz,wt_fact); 

% Launch the optimizer
[x,fval,exitflag,output] = ...
    fmincon(objfunc,x0,AA,bb,Aeq,beq,lb,ub,nlcfunc,options);

RMSSE = fval;
display(output)

% Because the values of the driver metrics in global memory represent the 
% results of the LAST solution, they may not be the same as the BEST  
% solution determined by the optimizer. Thus it is necessary to make one  
% last run using the final values of the design variables. These will  
% produce the correct driver metrics, which are passed back to the main
% program through the global variables. 

RMSSE_final = GOFAST_run(x,Hz,A,B,C,Me,gc,DVPAW,clb,cub,wt_fact,driver);
MJperMile_Dw = CEDw/DDw/1e6;

format
fid = fopen('Optimization Results Summary.txt','a');
fprintf(fid,'\n========= %s =========\n',datetime);
fprintf(fid,'Drive cycle = %s\n',cycle);
fprintf(fid,'Driver model = %2.1f\n',driver);
fprintf(fid,'Iterations = %u\n',output.iterations);
fprintf(fid,'GOFAST runs (function evals) = %u\n',output.funcCount);
fprintf(fid,'Max constraint violation = %9.3f\n',output.constrviolation);
fprintf(fid,'Message: %s \n\n',output.message);

fprintf(fid,'--- SAE J2951 Drive Metrics (10Hz, unweighted) ---\n');
fprintf(fid,'D = Driven     T = Target\n');
fprintf(fid,'CED (MJ) = %8.5f\n',CED/1e6);
fprintf(fid,'CET (MJ) = %8.5f\n',CET/1e6);
fprintf(fid,'ER (%%) = (CED-CET)/CET*100 = %8.4f\n\n',ER);

fprintf(fid,'DD (meters) = %9.1f\n',DD*1609.344);
fprintf(fid,'DT (meters) = %9.1f\n',DT*1609.344);
fprintf(fid,'DR (%%) = (DD-DT)/DT*100 = %8.4f\n\n',DR);

fprintf(fid,'--- The following are weighted ---\n');
fprintf(fid,'Objective:    RMSSE (mph) = %7.3f  (Target: < 0.4)\n',RMSSE_final);
fprintf(fid,'Constraint 1) EER (%%)     = %7.3f  (Target: |EER| < 0.5 %%)\n',EER);
fprintf(fid,'Constraint 2) ASCR (%%)    = %7.3f  (Target: |ASCR| < 1.5 %%)\n',ASCR);
fprintf(fid,'Constraint 3) IWR (%%)     = %7.3f  (Target: |IWR| < 1.5 %%)\n\n',IWR);

fprintf(fid,'--- FCA Drive Metrics (10Hz, weighted) ---\n');
fprintf(fid,'Constraint 4) CSI (mph/s)   = %8.3f  (Target for EPA75: -15 to +10)\n',CSI);
fprintf(fid,'Constraint 5) PEDAL (%%/100) = %8.3f  (Target for EPA75: 100 to 175)\n\n',PEDAL);

fprintf(fid,'--- VDE/Distance Ratio (weighted) ---\n');
if cycle == "WLTCb"
    fprintf(fid,'CED/DD (MJ/km)          = %8.5f\n',MJperMile_Dw/1.609344);
    fprintf(fid,'CET/DT (MJ/km)          = %8.5f\n',MJperMile_Tw/1.609344);
else
    fprintf(fid,'CED/DD (MJ/mile)        = %8.5f\n',MJperMile_Dw);
    fprintf(fid,'CET/DT (MJ/mile)        = %8.5f\n',MJperMile_Tw);
end
fprintf(fid,'Constraint 6) (D-T)/T %% = %8.5f  (Target = 0)\n\n',...
    (MJperMile_Dw - MJperMile_Tw)/MJperMile_Tw*100);

fid = fopen('Optimization Results Summary.txt','a');
n = 1:numel(x);
fprintf(fid,'--- Design Variable Summary ---\n');
fprintf(fid,'%s\n','  x     Variable(x)');
fprintf(fid,'%3u %12.8f\n',[n;x]);
fprintf(fid,'--- End of Optimization ---\n');
fclose(fid);

% Write miscellaneous outputs to a text file
dlmwrite('history_output.txt',output_data,' ')

msgbox({'Optimization has finished. For details refer' ...
        'to the Optimization Results Summary document.'})