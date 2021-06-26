function sum = cost_ctrb(P,Q,rhoperms)

sum=0;
    for ridx=1:size(rhoperms,1)
        rhos = num2cell(rhoperms(ridx, :));        
        sum=sum+trace(P(rhos{:})*value(Q(rhos{:})));        
    end

end

