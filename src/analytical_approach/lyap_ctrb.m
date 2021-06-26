function dRp = lyap_ctrb(dRp,dP,Rp,rhoperms,listParameter)
%Find dRp/drho
%   Derivative of Cholesky factor Rp, Analytical optimization, perturbation
%   method

sdpoptions=sdpsettings('solver','sdpt3','debug',1);
constraints = [];
dRp = tril(dRp);
eps = sdpvar(1);
n = size(dRp,1);
for i = 1:size(rhoperms,1)
    rhos = num2cell(rhoperms(i, :));
    sum = 0;
    for j = 1:length(listParameter)
        sum = sum + dP{j}(rhos{j});
    end
    
    constraint1 = dRp*transpose(value(Rp(rhos{:}))) + value(Rp(rhos{:}))*transpose(dRp) - value(sum);
    constraints = [constraints, constraint1<=eps*eye(n), constraint1>=-eps*eye(n)];

%     constraint1 = Rp_dot*Rp' + Rp*Rp_dot'-dP;
    
end
    lyap_ctrb = optimize(constraints,eps,sdpoptions)
    solution = value(dRp);
end

