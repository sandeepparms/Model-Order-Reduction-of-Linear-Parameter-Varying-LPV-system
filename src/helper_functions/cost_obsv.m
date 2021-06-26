function sum = cost_obsv(Q,k,no,P,rhoperms)

sum=0;
P_1=eye(no);             % Uses Identity matrix in the first iteration
    if k==1
        for ridx=1:size(rhoperms,1)
            rhos = num2cell(rhoperms(ridx, :));
            sum=sum+trace((Q(rhos{:}))*(P_1));            
        end               
    else
        for ridx=1:size(rhoperms,1)
            rhos = num2cell(rhoperms(ridx, :));
            sum=sum+trace(Q(rhos{:})*value(P(rhos{:})));
        end

    end

end

