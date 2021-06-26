function sys=sys_eqns_alternate(n,k,d,F)

A=pmat();

A(1,1)=-2*d;
A(1,2)= d;
A(1,n+1)=-(k(1)+k(2));
A(1,n+2)=k(2);
for i=2:n-1
    A(i,i-1)=d;
    A(i,i)=-3*d;
    A(i,i+1)=d;
    A(i,n+i-1)=k(i);
    A(i,n+i)=-(2*k(i)+k(i+1));
    A(i,n+i+1)=k(i+1);
end
A(n,:)=[zeros(1,n-2) d -2*d zeros(1,n-2) k(n) -2*k(n)];
for i=n+1:2*n
    A(i,i-n)=1;
end
B=[zeros(n-1,1); F; zeros(n,1)];
 
 %write output matrix C.Only displacement of n'th block is measured.
 C=[zeros(1,n) zeros(1,n-1) 1];
 
 %No feedthrough
 D=0;
 
 %Create PSS object
 sys=ss(A,B,C,D);
end
