% Differential algebraic equation system f(t,x(t),x'(t);theta) = 0.

function f = test_func(t,x,xp,pars)

% Retrieve parameters by name.
a = pars(1);
b = pars(2);

% Retrieve variables by name.
x1 = x(1); x1p = xp(1); 
x2 = x(2); x2p = xp(2); 
x3 = x(3); x3p = xp(3); 

% Differential algebraic equation system f(t,x,x') = 0.
f = zeros(length(x),1);

% x1  = exp(a*x2)
f(1)  = x1 - exp(a*x2);
% x2  = x3 - x1 - 2
f(2)  = x2 - x3 + x1 + 2;
% x3' = a*b*x1 + x2'
f(3)  = x3p - a*b*x1 - x2p;

end






























