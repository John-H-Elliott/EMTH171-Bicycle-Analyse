%{
Case Study #2 EMTH171  2020
E-Bike

Task 3 Description:


Orginal code by: P.J. Bones   UCECE
Code edited by: John Elliott and Samuel Vallance
Last modified:	09/10/2020
%}


clear
clc
close all


%{
===========================================================================
                            | Variable Setup | 
===========================================================================
%}


% Constants
g = 9.81;		% Gravitational acceleration in m.s^-2

% E-bike constants
M = 95;		    % Mass of e-bike + rider in kg
Jw = 0.5;		% Motor/wheels moment of intertia in kg.m^2
fd0 = 2;		% Static (& low speed) drag force in N
beta = 0.25;	% Wind resistance factor in N.m^-2.s^2
km = 3.2;		% Motor factor in N.m per amp
rw = 0.35;		% Wheel radius in m
rp = 0.09;      % Pedal radius in m
rs1 = 0.045;    % (Low Gear) radius in m
rs2 = 0.030;    % (High Gear) radius in m
pmax = 500;     % Max-Power output in W

% Simulation parameters
tdelta = 0.01;	% Time step in seconds
Im = 0;         % CHANGE
Tp = 6.48;        % CHANGE
alpha = 0;      % Angle of incline in rad for Task 2

Tpf = 11.6;       % CHANGE
Imf = 6.3;        % CHANGE
alphaf = 0.05;  % Final angle of incline in rad for Task 2

% Calculate derived quantities
Tm = Im * km;   % Fixed motor torque in N.m
Tmf = Imf * km; % Final fixed motor torque in N.m
Jef = Jw + rw^2 * M;
                % Effective moment of intertia in kg.m^2

% Initial conditions (standing start)
tarray = 0;     % Time taken
rsc = rs2;      % Current gear radius in m
Tpc = Tp;
alphac = alpha;
Tmc = Tm;       % Current motor torque in Nm
Tmh = Tm;       % Holder motor torque in Nm (For when power > 500 W)
om = 0;         % Initial omega_m ('om')
u = 0;          % Speed along road in m/sec
s = 0;          % Distance along road from start in m

% Distances
flatr = 5000;   % Flat road distance in m
slopem = 5000;  % Slope mapped distance in m
sloper = slopem/cos(alphaf);
                % Slope road distance in m
                
% Battery parameters
Qtotal = 18;
Qmax = Qtotal*0.05;
Qused = 0;
tf = 1800;
p = 0;
earray = 0;

%{
===========================================================================
                 | Euler's Method and Speed Calculations | 
===========================================================================
%}


% Continually compute u and s as the e-bike proceeds, advancing time  
% by tdelta sec each step until required distance is traveled
index = 1;
while (s <= (flatr + sloper)) && (Qused <= Qtotal)
   % Compute the motor acceleration at the start of the step
   windf = beta * u^2;
   gravityf =  M * g * sin(alphac);
   domdt = (Tmc + (rsc / rp) * Tpc - rw * (fd0 + windf + gravityf)) / Jef;

   % Estimate the state at the end of the time step (n = 'next')
   omn = om + tdelta * domdt;   	% Estimate om by Euler's method
   un = omn * rw;                   % Estimate u at end of step
   sn = s + tdelta * (un + u) / 2;	% Estimate s at end of step (Trap.)
   p = (Tpc * ((omn * rsc) / rp));    % Estimate power used
   
   domdtarray(index) = domdt;       % Store current values
   uarray(index) = u;
   sarray(index) = s;
   earray(index + 1) = earray(index) + (p * tdelta);
                                    % Stores new total energy used
   tarray(index + 1) = tarray(index) + tdelta;
                                    % Stores new time step taken
   index = index + 1;
   Qused = Qused + (Tmc / km) * tdelta / 3600;
                                    % Stores total battery used 
   
   % Advance to next time step
   om = omn;
   u = un; 
   s = sn;
   
   % Checks to ensure power is below maximum output.
   if (Tmc * omn) > pmax
       Tmc = (pmax / omn);       % Power is at max, so torque is limited
   else
       Tmc = Tmh;                % Returns torque to fixed value
   end
   
   % Checks if next step is at the slope so that variables can be updated. 
   % This process could have comparisons reduced if the while loop is made 
   % as a function.
   if (s >= flatr)
       Tmc = Tmf;
       Tmh = Tmf;
       Tpc = Tpf;
       alphac = alphaf;
       rsc = rs1; 
   end
end


%{
===========================================================================
                      | Display Results and Graphs | 
===========================================================================
%}

if s <= (flatr + sloper)
    fprintf("BATTERY DEAD\n")
else
    % Store final values (As the caculated values only get asigned in the 
    % arrays the next loop).
    domdtarray(index) = domdt;   
    uarray(index) = un;
    sarray(index) = s;
    earray(index + 1) = earray(index) + (p * tdelta);
    
    if (Qused > Qmax)
        fprintf("ERROR: BATTERY OVER 5%%\n")
        fprintf("Battery Used: %.4f Ah\n", Qused)
    end
    if (tarray(end) > tf)
        fprintf("ERROR: OVER 30 MINUTES\n")
    end
    fprintf("\nFinal Time: %.4f s\n", tarray(end))
    fprintf("Final Angular Speed (Motor): %.4f rads^-1\n", omn)
    fprintf("Final Linear Speed (Motor): %.4f ms^-1\n", un)
    fprintf("Final Distance: %.4f m\n", sn)


    figure(1)
    % Task 0 using Euler's method:  d\omega_m/dt
    subplot (3,1,1), plot (tarray, domdtarray)
    ylabel ('d\omega_m/dt for e-bike (rad.s^{-2})')
    xlabel ('Time (sec)')
    % Task 0 using Euler's method:  e-bike speed u
    subplot (3,1,2), plot (tarray, uarray)
    ylabel ('speed (m.s^{-1})')
    xlabel ('Time (sec)')
    % Task 0 using Trapezium method:  e-bike distance s
    subplot (3,1,3), plot (tarray, sarray)
    ylabel ('distance (m)')
    xlabel ('Time (sec)')
    
    figure(2)
    plot (tarray, earray(:,2:end))
    ylabel ('Energy J')
    xlabel ('Time (sec)')
end