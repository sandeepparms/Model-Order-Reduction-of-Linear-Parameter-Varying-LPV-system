%% Project Arbeit "Comparison of Variants of the Balanced Truncation Method for Linear Parameter-Varying Systems"
% Author:Sandeep Parameshwara, MSc. Mechatronics
% Mentor: Lennart Heeren, Institute of Control Systems, TUHH

clear;
thisdate = datestr(now, 'yyyy-mm-dd_HH-MM');
%% define system parameters and form the PSS object.
k0=0.5; krho=0.3;

% Assumption: all parameters have identical rate of parameter variation. 
% We can't form grid with lesser grid points than this becuase we need system
% matrices at [-1 0 1] and grid has to be extended [-1-delta, -1, 0, 1,
% 1+delta]. Consider using a denser grid if required.
rho_grid = [-1.02, -1, -0.9:0.1:1, 1.02];

% Rate of parameter variation
rho_dot=[-1 1];

no_params=2;

%no. of masses
n=4;

rho=cell(1,no_params);
for i=1:no_params
    rho{i}=pgrid(strcat('rho',num2str(i)),rho_grid,rho_dot);
end

% Define dependence of spring stiffness on parameters.
for i=1:no_params
    k(i)=k0+rho{i}*krho;
end
% for i=3:2:n
%     k(i)=k(1);
% end
% for i=4:2:n
%     k(i)=k(2);
% end
k(3)=k(1);
k(4)=k(2);


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

prod_size = [2*n 2*n];                      %For preallocation purposes
rhoperms = cell(1,no_params);
[rhoperms{:}] = ndgrid(rho_grid);           % Creates a table with all possible combination values of rho
rhoperms = reshape(cat(no_params + 1, rhoperms{:}), [], no_params);

rho_grid_actual=[-1,0,1];
rhoperms_actual = cell(1,no_params);
[rhoperms_actual{:}] = ndgrid(rho_grid_actual);           % Creates a table with all possible combination values of rho
rhoperms_actual = reshape(cat(no_params + 1, rhoperms_actual{:}), [], no_params);
%% Define Basis functions
yalmip('clear');
sdpoptions=sdpsettings('solver','sdpt3','debug',0,'verbose',0);

listParameter = fieldnames(sys.Parameter);  %List containing names of parameters

Q_mat = cell(1,no_params);
for i = 1:no_params
   Q_mat{i} =  sdpvar(2*n, 2*n,'symmetric');
end
Q0 = sdpvar(2*n, 2*n,'symmetric');

stringPart1 = 'Q = @(rho1';
stringPart2 = ') Q0 + Q_mat{1}*(rho1)';

for i = 2:length(listParameter)
   
    stringPart1 = [stringPart1 ',' listParameter{i}];
    stringPart2 = [stringPart2 ' +Q_mat{' num2str(i) '}*(' listParameter{i} ')'];
end
eval([stringPart1 stringPart2 ';'])

% Q= @ (rho1,rho2,rho3) Q0 + Q_mat{1}*(rho1) + Q_mat{2}*(rho2) + Q_mat{3}*(rho3);
dQ1=@(rho1) Q_mat{1};
dQ2=@(rho2) Q_mat{2};
% dQ3=@(rho3) Q_mat{3};
% dQ4=@(rho4) Q_mat{4};

stringpart = 'dQ = {dQ1';
for i = 2:length(listParameter)
    stringpart = [stringpart ',' strcat('dQ',num2str(i))];
end
eval([stringpart '};'])                        % Creates cell array containing derivative function handles of the basis function

P_mat = cell(1,no_params);
for i = 1:no_params
   P_mat{i} =  sdpvar(2*n, 2*n,'symmetric');
end
P0 = sdpvar(2*n, 2*n,'symmetric');

stringPart1 = 'P = @(rho1';
stringPart2 = ') P0 + P_mat{1}*(rho1)';

for i = 2:length(listParameter)
   
    stringPart1 = [stringPart1 ',' listParameter{i}];
    stringPart2 = [stringPart2 ' +P_mat{' num2str(i) '}*(' listParameter{i} ')'];
end
eval([stringPart1 stringPart2 ';'])


% P= @ (rho1,rho2,rho3) P0 + P_mat{1}*(rho1) + P_mat{2}*(rho2);
dP1=@(rho1) P_mat{1};
dP2=@(rho2) P_mat{2};
% dP3=@(rho3) P_mat{3};
% dP4=@(rho4) P_mat{4};

stringpart = 'dP = {dP1';
for i = 2:length(listParameter)
    stringpart = [stringpart ',' strcat('dP',num2str(i))];
end
eval([stringpart '};'])                          % Creates cell array containing derivative function handles of the basis function

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
% Redefine the grid for calculation.
rho_grid = [-1.02, -1, -0.98, -0.02, 0, 0.02, 0.98, 1, 1.02];
Nrho=length(rho_grid);
rhoperms2 = cell(1,no_params);
[rhoperms2{:}] = ndgrid(rho_grid);           % Creates a table with all possible combination values of rho
rhoperms2 = reshape(cat(no_params + 1, rhoperms2{:}), [], no_params);
    
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
U_vals = zeros([prod_size,size(rhoperms2,1)]);
S_vals = U_vals;
V_vals = U_vals;
systrans.A = zeros(2*n,2*n,size(rhoperms2, 1));
systrans.B = zeros(2*n,1,size(rhoperms2, 1));
systrans.C = zeros(1,2*n,size(rhoperms2, 1));
systrans.D = zeros(1,1,size(rhoperms2, 1));
systrans=ss(systrans.A,systrans.B,systrans.C,systrans.D);


%reshape sys.Data into 3D matrix for easier indexing later on
sysdata = reshape(sys.Data, size(sys.Data, 1), size(sys.Data, 2), []);

%loop over each unique combination
for ridx = 1:size(rhoperms2, 1)
    rhos = num2cell(rhoperms2(ridx, :));
    [U_vals(:, :, ridx), S_vals(:, :, ridx), V_vals(:, :, ridx)] = svd(pro(rhos{:}));
    Ti(:,:,ridx) = (sqrt(S_vals(:, :, ridx)))\  U_vals(:, :, ridx)' * Rq(rhos{:});
    T(:,:,ridx) = Rp(rhos{:}) * V_vals(:, :, ridx) / (sqrt(S_vals(:, :, ridx)));
end

% newshape = [prod_size, repelem(numel(rho_grid), no_params)];
%% Get derivatives of Transformation matrices for N-dimensional grid
% Central Difference method for N-dimensional grid
T_dot = derivative(no_params,T);

%% State transformation
% grid along rho_dot aswell
rho_dot = linspace(1,1,3);
vec1(:)=find(ismember(rhoperms,rhoperms_actual,'rows'));  %Find row indices of the actual grid w.r.t initial grid
vec2(:)=find(ismember(rhoperms2,rhoperms_actual,'rows')); %Find row indices of the actual grid w.r.t redefined grid

% Loop through possible combination of rho's in actual grid
for it=1:3^(no_params)
    x = vec1(it);
    y = vec2(it);
    for j=1:length(rho_dot)
        A_dot(:,:,it,j) = Ti(:,:,y)*sysdata.A(:,:,x)*T(:,:,y)-Ti(:,:,y)*T_dot(:,:,it)*rho_dot(j);       
        B_dot(:,:,it,j) = Ti(:,:,y)*sysdata.B(:,:,x);
        C_dot(:,:,it,j) = sysdata.C(:,:,x)*T(:,:,y);
        D_dot(:,:,it,j) = sysdata.D(:,:,x);
    end
end

systrans=ss(A_dot,B_dot,C_dot,D_dot);

%reshape into N-D matrices, where N = num_rhos + 2
newshape = [prod_size, repelem(numel(rho_grid_actual), no_params)];
systrans=reshape(systrans,repelem(numel(rho_grid_actual), no_params+1));

rho_dot = linspace(-1,1,2);

rho1 = pgrid('rho1',[-1,0,1],rho_dot);
rho2 = pgrid('rho2',[-1,0,1],rho_dot);
% rho3 = pgrid('rho3',[-1,0,1],rho_dot);
% rho4 = pgrid('rho4',[-1,0,1],rho_dot);
rhoDot = pgrid('rhoDot',[-1 0 1],rho_dot);
r_grid=rgrid(rho1,rho2,rhoDot);

% Balance the system and truncate the model
sysbal=pss(systrans,r_grid);
toc;
order=4;
elim=order+1:2*n;
systrunc1=modred(sysbal,elim,'Truncate');
systrunc2=modred(SYSB,elim,'Truncate');
%% LPV Norm
% Calculate lpv norm and error bounds
% Define basis for truncated and balanced system norm calculation

F=0.1*rhoDot+rho1+rho2+rho1*rho2 + rho2*rhoDot + rho1*rhoDot;
pf1=0.1+rho2+rho1;
pf2=1+rho2+rhoDot;
pf3=1+rho1+rhoDot;
b=basis(F,'rhoDot',pf1,'rho1',pf2,'rho2',pf3);
xb=[1,b,b^2];
lpvnorm(sysbal,xb)
lpvnorm(systrunc1,xb)

H=rho1+rho2+rho1*rho2;
h1=1+rho2;
h2=1+rho1;
c=basis(H,'rho1',h1,'rho2',h2);
xc=[1,c,c^2];
lpvnorm(sys,xc)
lpvnorm(SYSB,xc)
%% Plot

% Plot singular values
figure(1);
sigmas = zeros(2*n,length(rho_grid));
for i = 1:2*n
    sigmas(i,:)=S_vals(i,i,:,:);
    plot(rho_grid,sigmas(i,:));
    hold on;
    constG = G(i)*ones(1,9);
    plot(rho_grid,constG,'-.');
end
axis([-1.02 1.02 0 2.5]);
title('Singular values');

%% Bode Plot
% Use lpvsplit to extract data at frozen parameter values
figure(2);
rho_bode = [-1 0 1];
for i = 1:length(rho_bode)
    plpv = lpvsplit(sys,'rho1',rho_bode(i),'rho2',rho_bode(i));
    bode(plpv,'-b');
    hold on;
    trunclpv = lpvsplit(systrunc1,'rho1',rho_bode(i),'rho2',rho_bode(i),'rhoDot',0);
    bode(trunclpv,'--r');
    truncstab = lpvsplit(systrunc2,'rho1',rho_bode(i),'rho2',rho_bode(i));
    bode(trunclpv,'-.c');
end 
legend('full order system','truncated par Var','truncated constant')

%% Time domain response

% Define trajectory and plot the response
%Step response
t=linspace(0,20,100);
% u=ones(length(t),1);
ptraj.time=t;
ptraj.rhoDot=.01*ones(1,length(t));
ptraj.rho1=ones(1,length(t));
% ptraj.rho1=0.01*sin(t);
ptraj.rho2=ones(1,length(t));
% ptraj.rho3=ones(1,length(t));
% ptraj.rho4=ones(1,length(t));
straj.time=t;
straj.rho1=ones(1,length(t));
straj.rho2=ones(1,length(t));
% straj.rho3=ones(1,length(t));
% straj.rho4=ones(1,length(t));

figure(3);
lpvstep(sys,straj);
hold on;
lpvstep(systrunc1,ptraj);
lpvstep(systrunc2,straj);
legend('original','trunc-param var','trunc-const')

figure(4);
lpvimpulse(sys,straj);
hold on;
lpvimpulse(systrunc1,ptraj);
lpvimpulse(systrunc2,straj);
legend('original','trunc-param var','trunc-const')


%% Compute RMSE
% Compute and plot room mean square error between balanced and truncated
% systems for constant and parameter dependent models.
[x1,y1]=lpvimpulse(sys,straj);
[x2,y2]=lpvimpulse(systrunc1,ptraj);
[x3,y3]=lpvimpulse(SYSB,straj);
[x4,y4]=lpvimpulse(systrunc2,straj);
for i=1:length(t)
    error1(i)=sqrt((x1(i)-x2(i))^2);
    error2(i)=sqrt((x3(i)-x4(i))^2);
end
figure(5);
plot(t,error1,t,error2)
title('Root mean square error, Impulse response');
legend('paramter varying','constant');

% 2D singular values
% S_vals=reshape(S_vals,[20 20 9 9]);
% sigma2d_1=zeros(length(rho_grid),length(rho_grid));
% sigma2d_2=zeros(length(rho_grid),length(rho_grid));
% for i=1:length(rho_grid)
%     for j=1:length(rho_grid)
%         sigma2d_1(j,i)=S_vals(1,1,j,i);
%         sigma2d_2(j,i)=S_vals(2,2,j,i);
%     end
% end
% figure(5);
% surf(rho_grid,rho_grid,sigma2d_1);
% hold on;
% surf(rho_grid,rho_grid,sigma2d_2);
% axis([-1.02 1.02 -1.02 1.02 0 2.5]);
% title('Singular values in 2D parameter space');
%% Time domain simulation for custom signals
t = 0:0.01:10;
straj.time = t;
straj.rho1 = 0.1*sin(t);
straj.rho2 = 0.1*sin(t);
ptraj.rhoDot = 0.01*cos(t);
u = [zeros(size(0:0.01:3)) 0.1*ones(size(3.01:0.01:5)),...
-0.1*ones(size(5.01:0.01:7)) zeros(size(7.01:0.01:10))]';
lpvlsim(systrunc1,ptraj,u,t)