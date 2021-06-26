%% Project arbeit "Model Order Reduction of LPV systems"
% Author:Sandeep Parameshwara, MSc. Mechatronics
% Mentor: Lennart Heeren, Institute of Control Systems, TUHH

clear;
thisdate = datestr(now, 'yyyy-mm-dd_HH-MM');
%% define system parameters and form the PSS object.
k0=0.5; krho=0.3;

% Assumption: all parameters have identical rate of parameter variation.
rho_dot=[-1 1];
no_params=2;

%no. of masses
n=4;

% for i=1:no_params
%     rho_dot{i} = [-1 1];
% end

rho=cell(1,no_params);
for i=1:no_params
    rho{i}=pgrid(strcat('rho',num2str(i)),-1:1:1,rho_dot);
end

r_grid=rgrid(rho{:});         % Create Rectangular grid

% Define dependence of spring stiffness on parameters.

for i=1:no_params
    k(i)=k0+rho{i}*krho;
end
% for i=2:n
%     k(i)=k0+rho{1}*krho;
% end

% k(2)=k(1);
k(3)=k(1);
k(4)=k(2);
% k(6)=k(2);
% k(7)=k(1);
% k(8)=k(2);
% k(9)=k(1);
% k(10)=k(2);

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
rho_grid=-1:1:1;
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
Q= @ (rho1) Q0+Q_mat{1}*(rho1);
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
P= @ (rho1) P0 + P_mat{1}*(rho1);
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
Nrho=length(rho_grid);
    
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

%preallocation of results as 3D matrices. Can be reshape to ND matrices at the end
U_vals = zeros([prod_size, size(rhoperms, 1)]);
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

newshape = [prod_size, repelem(numel(rho_grid), no_params)];
%%
Ti_dot = gradient(Ti);
T_dot = gradient(T);

for i=2:size(rhoperms,1)-1
    T_dot(:,:,i)=(T(:,:,i+1)-T(:,:,i-1))/(2*0.1);
    Ti_dot(:,:,i)=(Ti(:,:,i+1)-Ti(:,:,i-1))/(2*0.1);
end

[U_vals(:, :, ridx+1), S_vals(:, :, ridx+1), V_vals(:, :, ridx+1)] = svd(pro(-1.1));
Ti(:,:,ridx+1) = (sqrt(S_vals(:, :, ridx+1)))\  U_vals(:, :, ridx+1)' * Rq(-1.1);
T(:,:,ridx+1) = Rp(-1.1) * V_vals(:, :, ridx+1) / (sqrt(S_vals(:, :, ridx+1)));

[U_vals(:, :, ridx+2), S_vals(:, :, ridx+2), V_vals(:, :, ridx+2)] = svd(pro(1.1));
Ti(:,:,ridx+2) = (sqrt(S_vals(:, :, ridx+2)))\  U_vals(:, :, ridx+2)' * Rq(1.1);
T(:,:,ridx+2) = Rp(1.1) * V_vals(:, :, ridx+2) / (sqrt(S_vals(:, :, ridx+2)));

Ti_dot(:,:,size(rhoperms,1)) = (Ti(:,:,ridx+2)-Ti(:,:,ridx))/(2*0.1);
T_dot(:,:,size(rhoperms,1)) =  (T(:,:,ridx+2)-T(:,:,ridx))/(2*0.1);

Ti_dot(:,:,1) = (Ti(:,:,2)-Ti(:,:,ridx+1))/(2*0.1);
T_dot(:,:,1) =  (T(:,:,2)-T(:,:,ridx+1))/(2*0.1);

% Ti_dot = reshape(Ti_dot,newshape);
% T_dot = reshape(T_dot,newshape);
%% State transformation
for i=1:size(rhoperms,1)
    A_dot(:,:,i) = Ti(:,:,i)*sysdata.A(:,:,1,1,i)*T(:,:,i)-Ti(:,:,i)*T_dot(:,:,i);
    B_dot(:,:,i) = Ti(:,:,i)*sysdata.B(:,:,1,1,i);
    C_dot = sysdata.C(:,:,1,1,i)*T(:,:,i);
    D_dot = sysdata.D(:,:,1,1,i);
end
systrans=ss(A_dot,B_dot,C_dot,D_dot);


%reshape into N-D matrices, where N = num_rhos + 2
% newshape = [prod_size, repelem(numel(rho_grid), no_params)];
% systrans=reshape(systrans,repelem(numel(rho_grid), no_params));
U_vals = reshape(U_vals(:,:,1:21), newshape);
S_vals = reshape(S_vals(:,:,1:21), newshape);
V_vals = reshape(V_vals(:,:,1:21), newshape);

% Balance the system and truncate the model
sysbal=pss(systrans,r_grid);
toc;
order=2;
elim=order+1:2*n;
systrunc1=modred(sysbal,elim,'Truncate');
systrunc2=modred(SYSB,elim,'Truncate');
% lpvnorm(sys)
% lpvnorm(sysbal-systrunc1)
% lpvnorm(SYSB-systrunc2)

%% Plot
%rho=pgrid('rho',-1:1:1,rho_dot);
% rho_rgrid=rgrid(rho{1},rho{2});
% for i=1:length(rho_grid)
%     for j=1:length(rho_grid)
%         systrans(:,:,j,i)= ss2ss(sys.Data(:,:,j,i),Ti(rho_grid(j),rho_grid(i)));
%     end
% end
% sysbal=pss(systrans,rho_rgrid);
% Plot Singular Values
% for i=1:2*n
%     plot(rho_grid,reshape(Svals(i,i,:),[1,length(rho_grid)]));
%     title('Singular Values');
%     hold on;
% end
%% Time domain response
%Step response
t=linspace(0,30,100);
ptraj.time=t;
ptraj.rho1=sin(t);
% ptraj.rho2=sin(t);
% ptraj.rho2=sin(t);
% ptraj.rho3=sin(t);
% ptraj.rho4=sin(t);
lpvstep(sys,ptraj);
hold on;
lpvstep(systrunc1,ptraj)
lpvstep(systrunc2,ptraj)
hold off
legend('original','varying','constant')
% lpvstep(systrunc,ptraj);