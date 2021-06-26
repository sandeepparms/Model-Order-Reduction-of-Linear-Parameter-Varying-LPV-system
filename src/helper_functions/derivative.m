function T_dot = derivative(no_params,T)
%DERIVATIVE of transformation matrix
%   Central difference method upto 4 dimension. Each diemnsion is
%   discretized with 9 points [-1.02, -1, -0.98, -0.02, 0, 0.02, 0.98, 1, 1.02]
init=3;
n = no_params;
ord=size(T,1);
del=0.04; % 2*length of interval for FDM
T_dot=zeros(ord,ord,3^n);
if n == 1               % 1 parameter
    for i = 1:3
        ind=i + 2*i;
        T_dot(:,:,i) = (T(:,:,ind) - T(:,:,ind-2))/del;
    end
elseif n==2             % 2 parameter. Since every dimension is discretized in a same manner, there are fixed relations between the indices.        
    for i = 1:3^n
        if i==1
            ind = init + sum(9.^(1:n-1)); % Geometric series
        elseif mod(i-1,3)==0
            ind = ind + 21;
        else
            ind = ind + 3;                
        end
        T_dot(:,:,i) = (T(:,:,ind) - T(:,:,ind-2))/del + (T(:,:,ind+8) - T(:,:,ind-10))/del;
    end
elseif n==3
    for i = 1:3^n
        if i==1
            ind = init + sum(9.^(1:n-1)); % Geometric series
        elseif mod(i-1,3)==0
            ind = ind + 21;
        elseif mod(i-1,9)==0
            ind = ind + 183;
        else
            ind = ind + 3;                
        end
        T_dot(:,:,i) = (T(:,:,ind) - T(:,:,ind-2))/del...
            + (T(:,:,ind+8) - T(:,:,ind-10))/del...
            +(T(:,:,ind+80) - T(:,:,ind-82))/del;
    end 
elseif n==4
    for i = 1:3^n
        if i==1
            ind = init + sum(9.^(1:n-1)); % Geometric series
        elseif mod(i-1,3)==0
            ind = ind + 21;
        elseif mod(i-1,9)==0
            ind = ind + 183;
        elseif mod(i-1,27)==0
            ind = ind + 1641;
        else
            ind = ind + 3;                
        end
        T_dot(:,:,i) = (T(:,:,ind) - T(:,:,ind-2))/del...
            + (T(:,:,ind+8) - T(:,:,ind-10))/del...
            +(T(:,:,ind+80) - T(:,:,ind-82))/del...
            +(T(:,:,ind+728) - T(:,:,ind-730))/del;
    end
elseif n==5
    for i = 1:3^n
        if i==1
            ind = init + sum(9.^(1:n-1)); % Geometric series
        elseif mod(i-1,3)==0
            ind = ind + 21;
        else
            ind = ind + 3;                
        end
        T_dot(:,:,i) = (T(:,:,ind) - T(:,:,ind-2))/del...
            + (T(:,:,ind+8) - T(:,:,ind-10))/del...
            +(T(:,:,ind+80) - T(:,:,ind-82))/del...
            +(T(:,:,ind+728) - T(:,:,ind-730))/del...
            +(T(:,:,ind+6560) - T(:,:,ind-6562))/del;
    end
elseif n==6
    for i = 1:3^n
        if i==1
            ind = init + sum(9.^(1:n-1)); % Geometric series
        elseif mod(i-1,3)==0
            ind = ind + 21;
        else
            ind = ind + 3;                
        end
        T_dot(:,:,i) = (T(:,:,ind) - T(:,:,ind-2))/del...
            + (T(:,:,ind+8) - T(:,:,ind-10))/del...
            +(T(:,:,ind+80) - T(:,:,ind-82))/del...
            +(T(:,:,ind+728) - T(:,:,ind-730))/del...
            +(T(:,:,ind+6560) - T(:,:,ind-6562))/del...
            +(T(:,:,ind+59048) - T(:,:,ind-59050))/del;
    end    
end


end

