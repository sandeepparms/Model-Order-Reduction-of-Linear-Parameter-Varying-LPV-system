%% Project arbeit "Model Order Reduction of LPV systems"
% Author:Sandeep Parameshwara, MSc. Mechatronics
% Mentor: Lennart Heeren, Institute of Control Systems, TUHH

clear;
thisdate = datestr(now, 'yyyy-mm-dd_HH-MM');
%% define system parameters and form the PSS object.
k0=0.5; krho=0.3;

% Assumption: all parameters have identical rate of parameter variation.
rho_grid = [-1,0,1];
rho_dot=[-0.1 0.1];
no_params=1;

%no. of masses
n=3;

% for i=1:no_params
%     rho_dot{i} = [-1 1];
% end

rho=cell(1,no_params);
for i=1:no_params
    rho{i}=pgrid(strcat('rho',num2str(i)),rho_grid,rho_dot);
end

r_grid=rgrid(rho{:});         % Create Rectangular grid

% Define dependence of spring stiffness on parameters.

for i=1:no_params
    k(i)=k0+rho{i}*krho;
end
% for i=2:n
%     k(i)=k0+rho{1}*krho;
% end

k(2)=k(1);
k(3)=k(1);
% k(4)=k(1);
% k(5)=k(1);
% k(6)=k(1);

d=0.75;                                      %Damping
m=1;                                         %Mass
F=1;                                         %external force on n'th mass in newtons
sys=sys_eqns_alternate(n,k,d,F);
tic;
[SYSB,G,Tg,Tig] = lpvbalreal(sys);  % Using LPVBALREAL command
toc;
%% extract PSS data
[A,B,C,D]=ssdata(sys);
nX=size(A,1);
nU=size(B,2);                               % Read number of states, inputs and outputs
nY=size(C,1);
rho_grid = [-1,0,1];
prod_size = [2*n 2*n];                      %For preallocation purposes
rhoperms = cell(1,no_params);
[rhoperms{:}] = ndgrid(rho_grid);           % Creates a table with all possible combination values of rho
rhoperms = reshape(cat(no_params + 1, rhoperms{:}), [], no_params);
%% Define Basis functions
yalmip('clear');
sdpoptions=sdpsettings('solver','sdpt3','debug',0,'verbose',0);

listParameter = fieldnames(sys.Parameter);  %List containing names of parameters

Q_mat = cell(1,no_params);
for i = 1:no_params
   Q_mat{i} =  sdpvar(2*n, 2*n,'symmetric');
end

% stringPart1 = 'Q = @(rho1';
% stringPart2 = ') Q0 + Q1 * rho1';

% for i = 2:length(listParameter)
%    
%     stringPart1 = [stringPart1 ',' listParameter{i}];
%     stringPart2 = [stringPart2 '+ Q' num2str(i) '*' listParameter{i}];
% end
% eval([stringPart1 stringPart2 ';'])
Q0 = sdpvar(2*n, 2*n,'symmetric');
Q= @ (rho1) Q0+Q_mat{1}*rho1;
dQ1=@(rho1) Q_mat{1};                            %This part can be automated if the basis function is decided
% dQ2=@(rho2) Q_mat{2};
% dQ3=@(rho3) Q_mat{3};

stringpart = 'dQ = {dQ1';
for i = 2:length(listParameter)
    stringpart = [stringpart ',' strcat('dQ',num2str(i))];
end
eval([stringpart '};'])                        % Creates cell array containing derivative function handles of the basis function
% dQ={dQ1,dQ2};

P_mat = cell(1,no_params);
for i = 1:no_params
   P_mat{i} =  sdpvar(2*n, 2*n,'symmetric');
end

% stringPart1 = 'P = @(rho1';
% stringPart2 = ') P0 + P1 * rho1';
% 
% for i = 2:length(listParameter)
%    
%     stringPart1 = [stringPart1 ',' listParameter{i}];
%     stringPart2 = [stringPart2 '+ P' num2str(i) '*' listParameter{i}];
% end
% 
% eval([stringPart1 stringPart2 ';'])
P0 = sdpvar(2*n, 2*n,'symmetric');
P= @ (rho1) P0 + P_mat{1}*rho1;
dP1=@(rho1) P_mat{1};                             %This part can be automated if the basis function is decided
% dP2=@(rho2) P_mat{2};
% dP3=@(rho3) 0;

stringpart = 'dP = {dP1';
for i = 2:length(listParameter)
    stringpart = [stringpart ',' strcat('dP',num2str(i))];
end
eval([stringpart '};'])                          % Creates cell array containing derivative function handles of the basis function

% dP={dP1,dP2};

%% Algorithm initial and stop conditions
maxiter=10;           % Max no. of iterations
reltol=1e-2;          % Convergence tolerance
go=1;                
k=1;
J=zeros(maxiter,1);   % Preallocation of cost Array
%% Algorithm
tic
while go==1
    
    Q  = obsv_gram(sys,rho_grid,rho_dot,Q,dQ,k,P,listParameter,rhoperms);
    J(k) = cost_obsv(Q,k,nX,P,rhoperms);             %Optimization of Q w.r.t fixed P
    k=k+1;
    
    P= ctrb_gram(sys,rho_grid,rho_dot,P,dP,Q,listParameter,rhoperms);
    J(k)= cost_ctrb(P,Q,rhoperms);                   %Optimization of P w.r.t fixed Q
    
    % Check if the convergence is achieved. If not, continue until it reaches max no. of iterations
    if k==1
        reldiff=2*reltol;
    else
        reldiff=abs((J(k)-J(k-1))/J(k));
    end
    
    if (k>=maxiter || (reldiff<reltol))
        go=0;
    end
        k=k+1;
end
%% Cholesky Factorization
% rho_grid = [-1,0,1];
Nrho=length(rho_grid);

% rhoperms = cell(1,no_params);
% [rhoperms{:}] = ndgrid(rho_grid);           % Creates a table with all possible combination values of rho
% rhoperms = reshape(cat(no_params + 1, rhoperms{:}), [], no_params);
    
    stringpart1 = 'Rq = @(rho1';
    stringpart2 = ') chol(value(Q(rho1';
    stringpart3 = ',''upper'')';
        for z = 2:length(listParameter)
            stringpart1 = [stringpart1 ',' listParameter{z}];
            stringpart2 = [stringpart2 ',' listParameter{z} ];
        end
    eval([stringpart1 stringpart2 '))' stringpart3 ';']);
    
    
%     Rq= @(rho1,rho2) chol(value(Q(rho1,rho2)),'upper');
     
    stringpart1 = 'Rp = @(rho1';
    stringpart2 = ') chol(value(P(rho1';
    stringpart3 = ',''lower'')';
        for z = 2:length(listParameter)
            stringpart1 = [stringpart1 ',' listParameter{z}];
            stringpart2 = [stringpart2 ',' listParameter{z} ];
        end
    eval([stringpart1 stringpart2 '))' stringpart3 ';']);
%     Rp= @(rho1,rho2) chol(value(P(rho1,rho2)),'lower');
    
    
    stringpart1 = 'pro = @(rho1';
    stringpart2 = ') Rq(rho1';
    stringpart3 = ')*Rp(rho1';
        for z = 2:length(listParameter)
            stringpart1 = [stringpart1 ',' listParameter{z}];
            stringpart2 = [stringpart2 ',' listParameter{z}];
            stringpart3 = [stringpart3 ',' listParameter{z}];
        end
    eval([stringpart1 stringpart2 stringpart3 ');'])
    
%     pro=@(rho1,rho2) Rq(rho1,rho2)*Rp(rho1,rho2);

%% Singular Value Decomposition

%Lennart
% rhoperms = linspace(-1,1,9)';


%preallocation of results as 3D matrices. Can be reshape to ND matrices at the end
U_vals = zeros([prod_size,size(rhoperms,1)]);
S_vals = U_vals;
V_vals = U_vals;
systrans.A = zeros(2*n,2*n,size(rhoperms, 1));
systrans.B = zeros(2*n,1,size(rhoperms, 1));
systrans.C = zeros(1,2*n,size(rhoperms, 1));
systrans.D = zeros(1,1,size(rhoperms, 1));
systrans=ss(systrans.A,systrans.B,systrans.C,systrans.D);


%reshape sys.Data into 3D matrix for easier indexing later on
sysdata = reshape(sys.Data, size(sys.Data, 1), size(sys.Data, 2), []);

%loop over each unique combination
for ridx = 1:size(rhoperms, 1)
    rhos = num2cell(rhoperms(ridx, :));
    [U_vals(:, :, ridx), S_vals(:, :, ridx), V_vals(:, :, ridx)] = svd(pro(rhos{:}));
    Ti(:,:,ridx) = (sqrt(S_vals(:, :, ridx)))\  U_vals(:, :, ridx)' * Rq(rhos{:});
    T(:,:,ridx) = Rp(rhos{:}) * V_vals(:, :, ridx) / (sqrt(S_vals(:, :, ridx)));
%     systrans(:, :, ridx) = ss2ss(sysdata(:, :, ridx), Ti(:,:,ridx));
end
%% Optimization to find out derivative of cholesky factors and computation of T_dot

dRp = sdpvar(2*n,2*n,'full');
dRp = lyap_ctrb(dRp,dP,Rp,rhoperms,listParameter);
dRp=double(dRp);
                                                          % Solving equation 7.69 from thesis using convex optimization
dRq = sdpvar(2*n,2*n,'full');
dRq = lyap_obsv(dRq,dQ,Rq,rhoperms,listParameter);
dRq=double(dRq);


Z_p1 = zeros(2*n,2*n,size(rhoperms,1));                      
Z_p2 = zeros(2*n,2*n,size(rhoperms,1));
dS = zeros(2*n, 2*n, size(rhoperms,1));
dV = zeros(2*n, 2*n, size(rhoperms,1));
for i = 1:size(rhoperms,1)
    rhos = num2cell(rhoperms(i,:));
    Z_p1(:,:,i) = Rp(rhos{:})' * Rq(rhos{:})' * Rq(rhos{:}) * Rp(rhos{:}); % Equation 7.70
    
    Z_p2(:,:,i) = (dRp' * Rq(rhos{:})' * Rp(rhos{:}) * Rq(rhos{:}))... 
                + (Rp(rhos{:})'*dRq'*Rq(rhos{:})*Rp(rhos{:}))...
                + (dRp'*Rq(rhos{:})'*Rp(rhos{:})*Rq(rhos{:}))'...     % calculating Z_p2 term from thesis
                + (Rp(rhos{:})'*dRq'*Rq(rhos{:})*Rp(rhos{:}))';
              
    for pert = 1:(2*n)
        % Derivative of sigma terms
         dS(pert,pert,i) = ((V_vals(:,pert,i))' * Z_p2(:,:,i) * V_vals(:,pert,i))/(2*S_vals(pert,pert,i));
         
         % Derivative of columns of V
         vi =1:2*n;
         vi(vi==pert) = []; %Removes one index to keep the denominator non-zero(see the below equation)
         sum=0;
         for pert2 = 1:length(vi)
             j=vi(pert2);            
             sum = sum+(((V_vals(:,j,i))' * Z_p2(:,:,i) * V_vals(:,pert,i))...
                        /((S_vals(pert,pert,i))^2-(S_vals(j,j,i))^2))*V_vals(:,j,i);
         end
         dV(:,pert,i) = sum;      
    end
end

% Calculate equation 7.68
T_dot = zeros(2*n,2*n,size(rhoperms,1));
for i = 1:size(rhoperms,1)
    rhos = num2cell(rhoperms(i,:));
    T_dot(:,:,i) = (dRp * V_vals(:,:,i) * S_vals(:,:,i)^(-0.5))...
                 + Rp(rhos{:}) * ((dV(:,:,i)* S_vals(:,:,i)^(-0.5)) - (0.5 * V_vals(:,:,i) * (S_vals(:,:,i)^(-1.5)) * dS(:,:,i)));             
end

%% Compute Balancing realization
rho_dot = linspace(-0.1,0.1,3);

%Lennart
% rhoperms = linspace(-1,1,3)';

for i=1:size(rhoperms,1)
    for j=1:3
    A_dot(:,:,j,i) = Ti(:,:,i)*sysdata.A(:,:,1,1,i) * T(:,:,i) - Ti(:,:,i) * T_dot(:,:,i)*rho_dot(j);       
    B_dot(:,:,j,i) = Ti(:,:,i)*sysdata.B(:,:,1,1,i);
    C_dot(:,:,j,i) = sysdata.C(:,:,1,1,i)*T(:,:,i);         % Computing matrix in definition 7.3.12 in thesis
    D_dot(:,:,j,i) = sysdata.D(:,:,1,1,i);
    end    
end
systrans=ss(A_dot,B_dot,C_dot,D_dot);


%reshape into N-D matrices, where N = num_rhos + 2
% newshape = [prod_size,3, repelem(numel(rho_grid), no_params)];
% systrans=reshape(systrans,repelem(numel(rho_grid), no_params));
% U_vals = U_vals(:,:,2:(size(rhoperms,1)-1));
% S_vals = S_vals(:,:,2:(size(rhoperms,1)-1));
% V_vals = V_vals(:,:,2:(size(rhoperms,1)-1));
% U_vals = reshape(U_vals(:,:,2:4), newshape);
% S_vals = reshape(S_vals(:,:,2:4), newshape);
% V_vals = reshape(V_vals(:,:,1:21), newshape);
rho_dot = linspace(-0.1,0.1,2);
rho = pgrid('rho',[-1,0,1],rho_dot);

rhoDot = pgrid('rhoDot',[-0.1 0 0.1]);
% Balance the system and truncate the model
sysbal=pss(systrans,rgrid(rhoDot,rho));
toc;
order=3;
elim=order+1:2*n;
systrunc1=modred(sysbal,elim,'Truncate');
systrunc2=modred(SYSB,elim,'Truncate');
% lpvnorm(sys)
% lpvnorm(sysbal-systrunc1)
% lpvnorm(SYSB-systrunc2)
%% Computing LPV system norm
% Construct basis
F = rho+rhoDot;
bf = basis(F,'rhoDot',1,'rho',1);
% xb = [1;bf;bf^2;bf^3];
lpvnorm(sysbal,bf)
%% Time domain response
%Step response
t=linspace(0,30,100)';
% u=zeros(1,length(t))';
u = ones(size(t));
ptraj.time=t;
straj.time=t;
ptraj.rhoDot=0.1*cos(t);
ptraj.rho=sin(0.1*t);

straj.rho1=sin(t);

lpvlsim(systrunc1,ptraj,u,t);
% lpvlsim(sysbal,ptraj,u,t)
hold on;
lpvstep(sys,straj);
% lpvstep(systrunc2,ptraj)
% hold off
% legend('original','varying','constant')
% lpvstep(systrunc,ptraj);