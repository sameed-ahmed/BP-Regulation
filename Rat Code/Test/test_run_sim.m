% Run simulation

function test_run_sim

close all

% Parameters
a = 3; b = 5;
pars = [a; b];

% Initial value
x0  = [1  ; 0; 3      ];
xp0 = [a*b; b; a*b + b];

% Initial time; Final time; Time vector;
t0 = 0; tf = 1; tspan = [t0, tf];

% Solve dae
[t,x] = ode15i(@(t,x,xp) test_func(t,x,xp,pars), tspan, x0, xp0);

% True solution
xtrue = [exp(a*b*t)          , ...
         b*t                 , ...
         exp(a*b*t) + b*t + 2];

figure
plot(t,x(:,1),t,xtrue(:,1))
legend('Numerical Solution','True Solution')
title('\fontsize{16}x_1')

figure
plot(t,x(:,2),t,xtrue(:,2))
legend('Numerical Solution','True Solution')
title('\fontsize{16}x_2')

figure
plot(t,x(:,3),t,xtrue(:,3))
legend('Numerical Solution','True Solution')
title('\fontsize{16}x_3')

end


































