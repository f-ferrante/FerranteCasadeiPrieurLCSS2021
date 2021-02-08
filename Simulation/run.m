%Run the solver and generates the solution vector [x1, x2, x3, x4, x5, x6, x7, x8] and the grid of
%points (t,z) over which the solution has been computed
[z1,t1, x1, x2, x3, x4, x5, x6, x7, x8]=hyp_static_bc();
z=z1(1,:);
t=t1(:,1);
%%
%%%Errors means
Ns=length(t);

e1(1,:)=x1(1,:)-x3(1,:);
e2(1,:)=x1(1,:)-x5(1,:);
e3(1,:)=x1(1,:)-x7(1,:);

e12(1,:)=x2(1,:)-x4(1,:);
e22(1,:)=x2(1,:)-x6(1,:);
e32(1,:)=x2(1,:)-x8(1,:);

E1=e1(1,:)+e2(1,:)+e3(1,:);
E12=e12(1,:)+e22(1,:)+e32(1,:);

t2=floor(Ns*0.1);
e1(1,:)=x1(t2,:)-x3(t2,:);
e2(1,:)=x1(t2,:)-x5(t2,:);
e3(1,:)=x1(t2,:)-x7(t2,:);

e12(1,:)=x2(t2,:)-x4(t2,:);
e22(1,:)=x2(t2,:)-x6(t2,:);
e32(1,:)=x2(t2,:)-x8(t2,:);

E2=e1(1,:)+e2(1,:)+e3(1,:);
E22=e12(1,:)+e22(1,:)+e32(1,:);


t3=floor(Ns*0.25);
e1(1,:)=x1(t3,:)-x3(t3,:);
e2(1,:)=x1(t3,:)-x5(t3,:);
e3(1,:)=x1(t3,:)-x7(t3,:);

e12(1,:)=x2(t3,:)-x4(t3,:);
e22(1,:)=x2(t3,:)-x6(t3,:);
e32(1,:)=x2(t3,:)-x8(t3,:);

E3=e1(1,:)+e2(1,:)+e3(1,:);
E32=e12(1,:)+e22(1,:)+e32(1,:);


t4=floor(Ns*1);
e1(1,:)=x1(t4,:)-x3(t4,:);
e2(1,:)=x1(t4,:)-x5(t4,:);
e3(1,:)=x1(t4,:)-x7(t4,:);

e12(1,:)=x2(t4,:)-x4(t4,:);
e22(1,:)=x2(t4,:)-x6(t4,:);
e32(1,:)=x2(t4,:)-x8(t4,:);


E4=e1(1,:)+e2(1,:)+e3(1,:);
E42=e12(1,:)+e22(1,:)+e32(1,:);
figure
hold on 
grid on
plot(z, E1(1,:)/3,'linewidth',2)
plot(z, E2(1,:)/3,'linewidth',2)
plot(z, E3(1,:)/3,'linewidth',2)
plot(z, E4(1,:)/3,'linewidth',2)
figure
hold on 
grid on
plot(z, E12(1,:)/3,'linewidth',2)
plot(z, E22(1,:)/3,'linewidth',2)
plot(z, E32(1,:)/3,'linewidth',2)
plot(z, E42(1,:)/3,'linewidth',2)
