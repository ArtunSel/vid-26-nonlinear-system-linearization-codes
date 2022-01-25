clear all,close all,clc;
%%
A1=[5,-4,3;-5,5,3;-1,-5,4];
Q1=[-5,3,-3;-1,-1,-3;-3,5,-4];
Q2=[-4,1,1;4,-4,-2;1,4,0];
Q3=[-1,-4,-1;-5,-3,-5;-3,-3,4];
% eq-pt
% x_eq=[1.443529522096883;1.770144022538879;-0.968037133334780];
% u_eq=[0.585423279694439;0.601236144699445];

% randi([imin,imax],___)
% A1=randi([-5,5],3)
% Q1=randi([-5,5],3)
% Q2=randi([-5,5],3)
% Q3=randi([-5,5],3)
param.A1=A1;
param.Q1=Q1;
param.Q2=Q2;
param.Q3=Q3;
% [t,y] = ode45(@(t,y) myode(t,y,ft,f,gt,g), tspan, ic, opts);
[z,fval] =  fmincon(@(z)cost_fun1(z,param),ones(5,1))
x_eq=z(1:3)
u_eq=z(4:5)
%%
%%
clearvars -except A1 Q1 Q2 Q3 x_eq u_eq
syms x1 x2 x3 u1 u2 real
syms x real
syms u real
x=[x1;x2;x3]
u=[u1;u2]
f_vec=A1*x+[x'*Q1*x;x'*Q2*x;x'*Q3*x]+[10;15;25]+[0;1;0]*u1+[0;0;1]*u2;
f_vec=simplify(f_vec)

Asym=jacobian(f_vec,x)
A=subs(Asym,{x1,x2,x3,u1,x2},{x_eq(1),x_eq(2),x_eq(3),u_eq(1),u_eq(2)})
A=double(A)

Bsym=jacobian(f_vec,u)
B=subs(Bsym,{x1,x2,x3,u1,x2},{x_eq(1),x_eq(2),x_eq(3),u_eq(1),u_eq(2)})
B=double(B)
%%
A
B
C=[1,0,0;0,1,0]
D=0*C*B
rank(ctrb(A,B))
rank(obsv(A,C))

plant=ss(A,B,C,D)
zpk(plant)
tzero(plant)
eig(plant.a)

% design a weight filters [LPF]
W1=makeweight(1,[3,0.1],0.01)
bodemag(W1)
evalfr(W1,0)    % unity dc-gain

% [K,CL,gamma,info] = mixsyn(plant,tracking-error,ctrl,output)
% [K,CL,gamma,info] = mixsyn(plant,W1,eye(2),eye(2))
[K,CL,gamma,info] = mixsyn(plant,W1,W1,eye(2))
% [K,CL,gamma,info] = mixsyn(plant,W1,[],[])
sigma(CL)
%%
function cost_val=cost_fun1(z,param)
A1=param.A1;
Q1=param.Q1;
Q2=param.Q2;
Q3=param.Q3;
x1=z(1);
x2=z(2);
x3=z(3);
x=[x1;x2;x3];
u1=z(4);
u2=z(5);
f_vec=A1*x+[x'*Q1*x;x'*Q2*x;x'*Q3*x]+...
    [10;15;25]+[0;1;0]*u1+[0;0;1]*u2;
% cost_val=norm([f_vec;abs(u1-5);abs([y2]-15)],2)
% cost_val=norm([f_vec;abs(u1-5);abs([x2^2]-15)],2)
cost_val=norm(f_vec,2);
end