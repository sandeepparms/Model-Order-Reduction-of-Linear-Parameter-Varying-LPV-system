function P = ctrb_gram(sys,rho_grid,rho_dot,P,dP,Q,listParameter,rhoperms)
[A,B,~,~]=ssdata(sys);
LMIConstraints_ctrb=[];
sdpoptions=sdpsettings('solver','sdpt3','debug',1);               % Solver settings
namelist=sys.A.Domain.IVName;                                     % List containing names of parameters
%valuelist=ones(size(sys.A.Domain.IVData,1),1);

for ridx=1:size(rhoperms,1)
    rhos = num2cell(rhoperms(ridx, :));
    
        stringpart1 = 'A_rho=lpvsubs(A,namelist,[rhos{1';
        for cost= 2:length(listParameter)
            stringpart1 =[stringpart1 '};rhos{' num2str(cost)]; 
        end
        eval([stringpart1 '}]);'])
%     A_rho=lpvsubs(A,namelist,[rhos{1};rhos{2}]);

        stringpart1 = 'B_rho=lpvsubs(B,namelist,[rhos{1';
        for cost= 2:length(listParameter)
            stringpart1 =[stringpart1 '};rhos{' num2str(cost)]; 
        end
        eval([stringpart1 '}]);'])
%     B_rho=lpvsubs(B,namelist,[rhos{1};rhos{2}]);
    
        % Positive definiteness of the Controllability gramian at every grid point
        LMI01=P(rhos{:});    
        LMIConstraints_ctrb=[LMIConstraints_ctrb,LMI01>=0];
        
        for j=1:length(rho_dot)
            rate=rho_dot(j);
                sum=0;
                for l=1:size(dP,2)
                    sum=sum+dP{l}(rhos{1});
                end
             
            % LMI should be satisfied at rho_dot(min) and rho_dot(max)    
            LMI04=-(rate*sum)+(A_rho*P(rhos{:}))+(P(rhos{:})*A_rho')+(B_rho*B_rho');
            LMIConstraints_ctrb=[LMIConstraints_ctrb,LMI04<=0];        
        end
end

diagnostic_ctrb = optimize(LMIConstraints_ctrb,cost_ctrb(P,Q,rhoperms),sdpoptions)

end

