clear;
%% define system parameters and form the PSS object.
k0=0.5; krho=0.3;
rho_dot=[-.1 .1];
no_params=2;

rho={};
for i=1:no_params
    rho{i}=pgrid(strcat('rho',num2str(i)),-1:1:1,rho_dot);
end

%no. of masses
n=4;

% Define dependence of spring stiffness on parameters.
k(1)=k0+rho{1}*krho;
k(2)=k0+rho{2}*krho;
k(3)=k(1);
k(4)=k(2);
% for i=5:6
%     k(i)= k0+rho{1}*krho;
% end

d=0.75;
m=1;
F=1;                                         %external force on n'th mass in newtons
sys=sys_eqns_alternate(n,k,d,F);
[SYSB,G,T,Ti] = lpvbalreal(sys);            % Using LPVBALREAL command

%% extract PSS data
[A,B,C,D]=ssdata(sys);
nX=size(A,1);
nU=size(B,2);
nY=size(C,1);

%% LMI's
yalmip('clear');
sdpoptions=sdpsettings('solver','sdpt3','debug',1);

listParameter = fieldnames(sys.Parameter);
Q0=sdpvar(2*n, 2*n,'symmetric');
Q1=sdpvar(2*n, 2*n,'symmetric');
Q2=sdpvar(2*n, 2*n,'symmetric');
Q= @ (rho1,rho2) Q0 + Q1*(rho1);% + Q2*(rho2);
dQ1=@(rho1) 0 + Q1;                           %Automate this Part!
dQ2=@(rho2) 0;% + Q2;

stringpart = 'dQ = {dQ1';
for i = 2:length(listParameter)
    stringpart = [stringpart ',' strcat('dQ',num2str(i))];
end
eval([stringpart '};'])
% dQ={dQ1,dQ2};


% stringPart1 = 'Q = @(rho1';
% stringPart2 = ') Q0 + Q1 * rho1';
% 
% for i = 2:length(listParameter)
%    
%     stringPart1 = [stringPart1 ',' listParameter{i}];
%     stringPart2 = [stringPart2 '+ Q' num2str(i) '*' listParameter{i}];
% end
% eval([stringPart1 stringPart2 ';'])

P0=sdpvar(2*n, 2*n,'symmetric');
P1=sdpvar(2*n, 2*n,'symmetric');
P2=sdpvar(2*n, 2*n,'symmetric');
P=@(rho1,rho2) P0 + P1*(rho1);% + P2*cos(rho2);
dP1=@(rho1) 0 + P1;                             %Automate this Part!
dP2=@(rho2) 0;% + P2;

stringpart = 'dP = {dP1';
for i = 2:length(listParameter)
    stringpart = [stringpart ',' strcat('dP',num2str(i))];
end
eval([stringpart '};'])

% dP={dP1,dP2};

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
%%
rho_grid=-1:1:1;
maxiter=10;
reltol=1e-4;
go=1;
k=1;
J=zeros(maxiter,1);
%% 
tic
while go==1
    
    Q  = obsv_gram(sys,rho_grid,rho_dot,Q,dQ,k,P,listParameter);
    J(k) = cost_obsv(Q,rho_grid,k,nX,P,listParameter);
    k=k+1;
    
    P= ctrb_gram(sys,rho_grid,rho_dot,P,dP,Q,listParameter);
    J(k)= cost_ctrb(P,rho_grid,Q,listParameter);
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
%%
Nrho=length(rho_grid);
    
    stringpart1 = 'Rq = @(rho1';
    stringpart2 = ') chol(value(Q(rho1';
    stringpart3 = ',''upper'')';
        for z = 2:length(listParameter)
            stringpart1 = [stringpart1 ',' listParameter{z}];
            stringpart2 = [stringpart2 ',' listParameter{z} ];
        end
    eval([stringpart1 stringpart2 '))' stringpart3 ';']);
    
    
%     Rq= @(rho1,rho2) chol(value(Q(rho1,rho2)),'upper');                    %How to automate argument list in Function handle?
     
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
            stringpart3 = [stringpart3 ',' listParameter{z} ');'];
        end
    eval([stringpart1 stringpart2 stringpart3])
    
%     pro=@(rho1,rho2) Rq(rho1,rho2)*Rp(rho1,rho2);
%%       
for i=Nrho:-1:1
    rho1=rho_grid(i);
    for j=Nrho:-1:1
        rho2=rho_grid(j);
        [U,S,V]= svd(pro(rho1,rho2));
        [Uvals(:,:,j,i),Svals(:,:,j,i),Vvals(:,:,j,i)]= svd(pro(rho1,rho2));
        Ti= inv(sqrt(S))*U'*Rq(rho1,rho2);
        T = Rp(rho1,rho2)*V*inv(sqrt(S));
        systrans(:,:,j,i)= ss2ss(sys.Data(:,:,j,i),Ti);
    end
end
rho_rgrid=rgrid(rho{1},rho{2});
sysbal=pss(systrans,rho_rgrid);
%% Plot Singular Values
% for i=1:2*n
%     plot(rho_grid,reshape(Svals(i,i,:),[1,length(rho_grid)]));
%     title('Singular Values');
%     hold on;
% end
toc;
%% Truncate the model
%rho=pgrid('rho',-1:1:1,rho_dot);
% rho_rgrid=rgrid(rho{1},rho{2});
% for i=1:length(rho_grid)
%     for j=1:length(rho_grid)
%         systrans(:,:,j,i)= ss2ss(sys.Data(:,:,j,i),Ti(rho_grid(j),rho_grid(i)));
%     end
% end
% sysbal=pss(systrans,rho_rgrid);

order=2;
elim=order+1:2*n;
systrunc=modred(sysbal,elim,'Truncate');
% lpvnorm(sys)
% lpvnorm(sys-sysbal)