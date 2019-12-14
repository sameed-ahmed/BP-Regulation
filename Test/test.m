function test

% close all

% Initial condition
x0 = 5;

% Time at which to keep steady state, change a parameter, etc.
tchange = 1;
% Initial time; Final time; % Time vector;
t0 = 0; tend = tchange + 9; tspan = [t0, tend];

% ode options
options = odeset();

[t,x] = ode45(@(t,x) odefun(t,x,tchange), tspan, x0, options);

figure
plot(t,x)

end

function dx = odefun(t,x,tchange)

a = 2;
b = 2;

if     t <  tchange
    
elseif tchange <= t && t < 4*tchange
    a = 3;
    b = 1/3 * t + 5/3;
elseif 4*tchange <= t
    a = 3;
    b = 3;
end

dx = a - b;

end