function test_fun_handle

% a = 1;
% b = 2;

% test = @(pars) myfun(a,b,pars);
test = @myfun;

test(4)

end

% function c = myfun(a,b,pars)
function c = myfun(pars)

a = 1;
b = 2;
c = (a + b) * pars;

end