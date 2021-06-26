function dRq = lyap_obsv(dRq,dQ,Rq,rhoperms,listParameter)
%Find dRp/drho
%   Derivative of Cholesky factor Rq, Analytical optimization, perturbation
%   method

sdpoptions=sdpsettings('solver','sdpt3','debug',1);
constraints = [];
dRq = triu(dRq);
eps = sdpvar(1);
n = size(dRq,1);
for i = 1:size(rhoperms,1)
    rhos = num2cell(rhoperms(i, :));
    sum = 0;
    for j = 1:length(listParameter)
        sum = sum + dQ{j}(rhos{j});
    end
    
    constraint1 = transpose(dRq)*value(Rq(rhos{:})) + transpose(value(Rq(rhos{:})))*dRq - value(sum);
    constraints = [constraints, constraint1<=eps*eye(n), constraint1>=-eps*eye(n)];

%     constraint1 = Rp_dot*Rp' + Rp*Rp_dot'-dP;
    
end
    lyap_obsv = optimize(constraints,eps,sdpoptions)
    solution = value(dRq);
end

