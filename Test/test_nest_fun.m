function test_nest_fun

x0 = 2;

% Last place ode solver was called
xLast = [];
% xLast = 2;
% ODE solution structure
sol = []; 
% % Test var defiined in nested fun.
% test_inner = []; 

% xLast
% x1 = constr(x0);
% xLast
% x  = objfun(x1);
% xLast

for i = 1:10
    tic
    x = objfun(x0);
    x0 = x;
    toc
end

x
% test_inner

    function y = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            x0 = 2;
            sol = 5*x0;
            xLast = x;
            pause(1);
            disp('obj not equal')
        else
            disp('obj equal')
        end
        y = sol;
%         test_inner = 1;
    end

    function c = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            x0 = 2;
            sol = 5*x0;
            xLast = x;
            disp('con not equal')
        else
            disp('con equal')
        end
        c = sol;
    end

end





















