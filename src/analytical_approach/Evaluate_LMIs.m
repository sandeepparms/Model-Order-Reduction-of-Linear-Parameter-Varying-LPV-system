% Evaluate Lyapuniv inequalities and positive definitness of Lyapunov
% matrices on a denser grid.
P0=double(P0);
P_mat{1}=double(P_mat{1});
P_mat{2}=double(P_mat{2});
% P_mat{3}=double(P_mat{3});

Q0=double(Q0);
Q_mat{1}=double(Q_mat{1});
Q_mat{2}=double(Q_mat{2});
% Q_mat{3}=double(Q_mat{3});

namelist=sys.A.Domain.IVName; 
rho_grid = linspace(-1,1,100);
rhoperms = cell(1,no_params);
[rhoperms{:}] = ndgrid(rho_grid);
rhoperms = reshape(cat(no_params + 1, rhoperms{:}), [], no_params);

for ridx=1:size(rhoperms,1)
        rhos = num2cell(rhoperms(ridx, :));

        stringpart1 = 'A_rho=lpvsubs(A,namelist,[rhos{1';
        for cost= 2:length(listParameter)
            stringpart1 =[stringpart1 '};rhos{' num2str(cost)]; 
        end
        eval([stringpart1 '}]);'])

        stringpart1 = 'C_rho=lpvsubs(C,namelist,[rhos{1';
        for cost= 2:length(listParameter)
            stringpart1 =[stringpart1 '};rhos{' num2str(cost)]; 
        end
        eval([stringpart1 '}]);'])
        
        stringpart1 = 'B_rho=lpvsubs(B,namelist,[rhos{1';
        for cost= 2:length(listParameter)
            stringpart1 =[stringpart1 '};rhos{' num2str(cost)]; 
        end
        eval([stringpart1 '}]);'])
            
        % Positive definiteness of the Observability gramian at every grid point
        LMI01=value(Q(rhos{:}));
        LMI02=value(P(rhos{:}));
        
        for j = 1:length(rho_dot)
        rate=rho_dot(j);
        sum1 = 0;
        sum2 = 0;
            for l=1:size(dP,2)
                sum2=sum2+dP{l}(rhos{1});
            end
            for l=1:length(dQ)
                sum1=sum1+(dQ{l}(rhos{l}));
            end                
            % LMI should be satisfied at rho_dot(min) and rho_dot(max)
            LMI03=value(((rate*sum1)+(A_rho'*Q(rhos{:}))+(Q(rhos{:})*A_rho)+C_rho'*C_rho));
            LMI04=value(-(rate*sum2)+(A_rho*P(rhos{:}))+(P(rhos{:})*A_rho')+(B_rho*B_rho'));
            
            eigTest_current1 = eig(LMI01);
            eigTest_current1 = sum(real(eigTest_current1) < 0);
            eigTest1(ridx,j) = eigTest_current1;
            
            eigTest_current2 = eig(LMI02);
            eigTest_current2 = sum(real(eigTest_current2) < 0);
            eigTest2(ridx,j) = eigTest_current2;
            
            eigTest_current3 = eig(LMI03);
            eigTest_current3 = sum(real(eigTest_current3) > 0);
            eigTest3(ridx,j) = eigTest_current3;
            
            eigTest_current4 = eig(LMI04);
            eigTest_current4 = sum(real(eigTest_current4) > 0);
            eigTest4(ridx,j) = eigTest_current4;
        end 
end

eigTest1 = sum(sum(eigTest1));
if (eigTest1 ==0)
    fprintf('exp stability guaranteed on denser grid!\n')
else
    fprintf('exp stability not guaranteed on denser grid!\n')
end

eigTest2 = sum(sum(eigTest2));
if (eigTest2 ==0)
    fprintf('exp stability guaranteed on denser grid!\n')
else
    fprintf('exp stability not guaranteed on denser grid!\n')
end

eigTest3 = sum(sum(eigTest3));
if (eigTest3 ==0)
    fprintf('exp stability guaranteed on denser grid!\n')
else
    fprintf('exp stability not guaranteed on denser grid!\n')
end

eigTest4 = sum(sum(eigTest4));
if (eigTest4 ==0)
    fprintf('exp stability guaranteed on denser grid!\n')
else
    fprintf('exp stability not guaranteed on denser grid!\n')
end