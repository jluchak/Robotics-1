%% Constant Variables
T=2; % Total time duration [s]
%syms(a,b,c,d,ql,qr);
d=10; % Distance between the shoulders
a=20; % Length of upper arm [m]
b=20; % Length of forearm [m]
c=5; % Length of hand [m]
Sl=[0,0]; % origin/left shoulder
Sr=[d,0]; % right shoulder
ql=deg2rad(150); % 
qr=deg2rad(30);
theta=deg2rad(45);
%% Equations
El=[a*cos(ql),a*sin(ql)];
Er=[a*cos(qr)+d,a*sin(qr)];
C=[(Er(1)+El(1))/2,(Er(2)+El(2))/2];
vec_ErC=Er-C;
norm_ErC=norm(vec_ErC);
m=vec_ErC/norm_ErC; % identity vector ErC
n=[-1*m(2),m(1)]; % identity vector WC
norm_WC=sqrt(b^2-norm_ErC^2);
vec_WC=m*norm_WC;
W=C+vec_WC;
H=W+[c*cos(theta),c*sin(theta)];
%H=[sqrt((2*b^2-a^2*(1-cos(ql-qr)))/(a^2*(1-cos(ql-qr))))*(a*cos((ql+qr)/2)*sin((ql-qr)/2))+c*cos(theta)+d/2+(a*cos((ql+qr)/2)*cos((ql-qr)/2)),sqrt((2*b^2-a^2*(1-cos(ql-qr)))/(a^2*(1-cos(ql-qr))))*(a*sin((ql+qr)/2)*sin((ql-qr)/2)-d/2)+c*sin(theta)+(a*sin((ql+qr)/2)*cos((ql-qr)/2))];
%H=[-1*sqrt((2*b^2-a^2*(1-cos(ql-qr)))/(a^2*(1-cos(ql-qr))))*(a/2*sin(qr)-a/2*sin(ql))+c*cos(theta)+d/2+(a/2*cos(qr))+(a/2*cos(ql)),-1*sqrt((2*b^2-a^2*(1-cos(ql-qr)))/(a^2*(1-cos(ql-qr))))*(a/2*cos((ql))-a/2*cos(qr)-d/2)+c*sin(theta)+(a/2*sin(qr))+(a/2*sin(ql))];

%Hx=H(1);
%Hy=H(2);
plot([Sl(1),El(1)],[Sl(2),El(2)]); hold on;
plot([El(1),W(1)],[El(2),W(2)]);
plot([W(1),H(1)],[W(2),H(2)]);
plot([Sr(1),Er(1)],[Sr(2),Er(2)]);
plot([Er(1),W(1)],[Er(2),W(2)]);
plot([El(1),Er(1)],[El(2),Er(2)]);
