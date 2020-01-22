function test_sim_ann

options = [];
% options = optimoptions('simulannealbnd', 'Display','iter');

fun = @dejong5fcn;

x0 = [0,0];
lb = [-64,-64];
ub = [64,64];

x1 = simulannealbnd(fun,x0,            lb,ub,options)

x2 = fmincon       (fun,x0,[],[],[],[],lb,ub,options)

end