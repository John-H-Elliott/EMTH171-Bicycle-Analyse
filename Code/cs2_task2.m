%{
Case Study #2 EMTH171  2020
E-Bike

Task 2 Description:
The e-bike starts from a standing start on a level road with fully charged 
batteries. A motor current of Im = 2 A and a pedal torque of Tp = 8 Nm 
(using high gear) is applied until the e-bike reaches the start of a hill 
after 2 km. The road now slopes upward with Alpha = 0.07 rad, so the rider 
engages low gear and applies p T = 15 N.m and m I = 6.5 A. The map distance 
from the start of the rise to the top is 4 km. Note that the maximum power 
that can be applied by the e-bike motor is 500 W.

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

t=0;
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
Im = 2;		    % Inital motor current in A for Task 2
Tp = 8;         % Inital pedal torque in Nm for Task 2
alpha = 0;      % Angle of incline in rad for Task 2

Tpf = 15;       % Final pedal torque in Nm for Task 2
Imf = 6.5;      % Final motor current in A for Task 2
alphaf = 0.07;  % Final angle of incline in rad for Task 2

% Calculate derived quantities
Tm = Im * km;   % Fixed motor torque in N.m
Tmf = Imf * km; % Final fixed motor torque in N.m
Jef = Jw + rw^2 * M;
                % Effective moment of intertia in kg.m^2

% Initial conditions (standing start)
tarray = 0;     % Time taken
rsc = rs2;      % Current gear radius in m
Tmc = Tm;       % Current motor torque in Nm
Tmh = Tm;       % Holder motor torque in Nm (For when power > 500 W)
om = 0;         % Initial omega_m ('om')
u = 0;          % Speed along road in m/sec
s = 0;          % Distance along road from start in m

% Distances
flatr = 2000;   % Flat road distance in m
slopem = 4000;  % Slope mapped distance in m
sloper = slopem/cos(alphaf);
                % Slope road distance in m
                
% Battery parameters
Qtotal = 18;
Qused = 0;


%{
===========================================================================
                 | Euler's Method and Speed Calculations | 
===========================================================================
%}
% Continually compute u and s as the e-bike proceeds, advancing time  
% by tdelta sec each step until required distance is traveled
index = 1;
ha = 0;
while (s <= (flatr + sloper))
   % Compute the motor acceleration at the start of the step
   windf = beta * u^2;
   gravityf =  M * g * sin(alpha);
   domdt = (Tmc + (rsc / rp) * Tp - rw * (fd0 + windf + gravityf)) / Jef;

   % Estimate the state at the end of the time step (n = 'next')
   omn = om + tdelta * domdt;   	% Estimate om by Euler's method
   un = omn * rw;                   % Estimate u at end of step
   sn = s + tdelta * (un + u) / 2;	% Estimate s at end of step (Trap.)
   
   domdtarray(index) = domdt;       % Store current values
   uarray(index) = u;
   sarray(index) = s;
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
       Tmc = (pmax / omn);          % Power is at max, so torque is limited
   else
       Tmc = Tmh;                   % Returns torque to fixed value
   end
   
   % Checks if next step is at the slope so that variables can be updated. 
   % This process could have comparisons reduced if the while loop is made 
   % as a function.
   if (s >= flatr)          
       Tmc = Tmf;           % Changes motor torque to sloped value        
       Tmh = Tmf;           % Holder motor torque, for when power > 500 W
       Tp = Tpf;            % Changes pedal torque to sloped value
       alpha = alphaf;      % Changes angle to sloped angle
       rsc = rs1;           % Changes gear from low to high
       ha(end+1) = omn;
   end
end
% Store final values (As the caculated values only get asigned in the 
% arrays the next loop).
domdtarray(index) = domdt;   
uarray(index) = un;
sarray(index) = s;

Qpercent = Qused/Qtotal * 100; 


%{
===========================================================================
                      | Display Results and Graphs | 
===========================================================================
%}

fprintf("Battery Used: %.4f %%\n", Qpercent)
fprintf("Final Time: %.4f s\n", tarray(end))
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