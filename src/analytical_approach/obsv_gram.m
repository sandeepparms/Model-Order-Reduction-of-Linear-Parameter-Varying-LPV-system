function Q = obsv_gram(sys,rho_grid,rho_dot,Q,dQ,k,P,listParameter,rhoperms)

[A,~,C,~]=ssdata(sys);
no=size(A,1);
LMIConstraints_obsv=[];
sdpoptions=sdpsettings('solver','sdpt3','debug',1);            % Solver settings
namelist=sys.A.Domain.IVName;                                  % List containing names of parameters
% valuelist=ones(size(sys.A.Domain.IVData,1),1);
    for ridx=1:size(rhoperms,1)
        rhos = num2cell(rhoperms(ridx, :));

        stringpart1 = 'A_rho=lpvsubs(A,namelist,[rhos{1';
        for cost= 2:length(listParameter)
            stringpart1 =[stringpart1 '};rhos{' num2str(cost)]; 
        end
        eval([stringpart1 '}]);'])
%         A_rho=lpvsubs(A,namelist,[rhos{1};rhos{2}]);

        stringpart1 = 'C_rho=lpvsubs(C,namelist,[rhos{1';
        for cost= 2:length(listParameter)
            stringpart1 =[stringpart1 '};rhos{' num2str(cost)]; 
        end
        eval([stringpart1 '}]);'])

%         C_rho=lpvsubs(C,namelist,[rhos{1};rhos{2}]);
            
        % Positive definiteness of the Observability gramian at every grid point
        LMI02=Q(rhos{:});
        LMIConstraints_obsv=[LMIConstraints_obsv,LMI02>=0];
        for j = 1:length(rho_dot)
        rate=rho_dot(j);
        sum=0;
            for l=1:length(dQ)
            sum=sum+(dQ{l}(rhos{l}));
            end    
            
            % LMI should be satisfied at rho_dot(min) and rho_dot(max)
            LMI03=((rate*sum)+(A_rho'*Q(rhos{:}))+(Q(rhos{:})*A_rho)+C_rho'*C_rho);            
            LMIConstraints_obsv=[LMIConstraints_obsv, LMI03<=0];          
                   
        end 
    end
   diagnostic_obsv = optimize(LMIConstraints_obsv,cost_obsv(Q,k,no,P,rhoperms),sdpoptions)
end

