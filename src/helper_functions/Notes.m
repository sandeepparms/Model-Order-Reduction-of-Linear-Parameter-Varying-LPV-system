    A1= Ti(:,:,i)*sysdata.A(:,:,1,1,i)*T(:,:,i)-Ti(:,:,i)*T_dot(:,:,i-1)*rho_dot(1);
    A2= Ti(:,:,i)*sysdata.A(:,:,1,1,i)*T(:,:,i)-Ti(:,:,i)*T_dot(:,:,i-1)*rho_dot(2);
    V = [reshape(A1,1,[]);reshape(A2,1,[])];
    Vq = interp1([-1;1],V,0);
    A_dot(:,:,i-1)=reshape(Vq,size(A1));
    lpvinterp(pmat(A_dot),{'rhoDot';'rho1'},{0 ;-1})