clear
close all 
clc

%Img=imread("C:\Users\Dell\OneDrive - Politecnico di Milano\4° anno\1° semestre\Image analysis and computer vision\HOMEWORK 2022-2023\stilizzata.png");
Img=imread("C:\Users\Dell\OneDrive - Politecnico di Milano\4° anno\1° semestre\Image analysis and computer vision\HOMEWORK 2022-2023\image 1.jpg");
Img=rgb2gray(Img);
Img=im2double(Img);
edges = edge(Img, 'canny', [0.025 0.06]);
figure(1), imshow(edges);
hold on;

%Computation of C1
[x, y] = getpts();
scatter(x,y,40,'filled');
A=[x.^2 x.*y y.^2 x y ones(size(x))];
N = null(A);
cc = N(:, 1);
[a1, b1, c1, d1, e1, f1] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));
C1=[a1 b1/2 d1/2; b1/2 c1 e1/2; d1/2 e1/2 f1];
p1=[x(1), y(1), 1];
p2=[x(5), y(5), 1];

%Computation of C2
[x, y] = getpts();
scatter(x,y,40,'filled');
A=[x.^2 x.*y y.^2 x y ones(size(x))];
N = null(A);
cc = N(:, 1);
[a2, b2, c2, d2, e2, f2] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));
C2=[a2 b2/2 d2/2; b2/2 c2 e2/2; d2/2 e2/2 f2];
p3=[x(1), y(1), 1];
p4=[x(5), y(5), 1];

%Computation of L1 and L2
L1=cross(p1,p3);
L2=cross(p2,p4);

%print contour of img
syms 'g';
syms 's';
eq1 = a1*g^2 + b1*g*s + c1*s^2 + d1*g + e1*s + f1;
eq2 = a2*g^2 + b2*g*s + c2*s^2 + d2*g + e2*s + f2;
eqns= [eq1==0, eq2==0];
fimplicit(eqns);
plot([p1(1), p3(1)], [p1(2), p3(2)], 'g');
plot([p2(1), p4(1)], [p2(2), p4(2)], 'g');

%Computation of image of circular points
S = solve(eqns, [g,s]);
s1 = [double(S.g(1));double(S.s(1));1];
s2 = [double(S.g(2));double(S.s(2));1];
s3 = [double(S.g(3));double(S.s(3));1];
s4 = [double(S.g(4));double(S.s(4));1];
II = s1;
JJ = s2;
imDCCP = II*JJ' + JJ*II';
imDCCP = imDCCP./norm(imDCCP);

%Computation and plot of h
h=cross(II, JJ);
h=h/h(3);
x = linspace(1,100000,1000000);
y=((JJ(2)-II(2))/(JJ(1)-II(1)))*(x-II(1))+II(2);
plot(x,y,'linewidth',2,'color','r')

%Computation of diameters parallel to x and y for both C1 and C2
d1=C1*II;
d2=C1*JJ;
d3=C2*II;
d4=C2*JJ;
d1= d1/d1(3);
d2= d2/d2(3);
d3= d3/d3(3);
d4= d4/d4(3);
v1=cross(d1,d3);
v2=cross(d2,d4);

%Computation and plot of centres of C1 and C2
O1=cross(d1,d2);
O1=O1./O1(3);
O2=cross(d3,d4);
O2=O2./O2(3);
scatter(O1(1), O1(2), 40,'filled');
scatter(O2(1), O2(2), 40,'filled');

%Computation of the image of the axis
a=cross(O1,O2);
x = linspace(1,100000,1000000);
y=((O2(2)-O1(2))/(O2(1)-O1(1)))*(x-O1(1))+O1(2);
plot(x,y,'linewidth',2,'color','r')

%computation of matrix K
syms aa U0 V0 f;
%omega = [aa^2 0 -U0*aa^2; 0 1 -V0; -U0*aa^2 -V0 (f^2+aa^2*U0^2+V0^2)];
omega = [1 0 -U0; 0 1 -V0; -U0 -V0 (f^2+U0^2+V0^2)];
eq1 = II.' * omega * II == 0;
eq2 = v1.' * omega * v2 == 0;
eq3 = JJ.' * omega * JJ ==0;
S = solve ([eq1 eq2 eq3],[aa U0 V0 f]);
f = double(S.f);
U0 = double(S.U0);
V0 = double(S.V0);
f = abs(f(1));
U0 = abs(U0(1));
V0 = abs(V0(1));
K = [f 0 U0; 0 f V0; 0 0 1];


%Computation of cone axis orientation
P=[K [0 0 0]'];
pih=P.'*h;
theta1 = 90 + rad2deg(atan(pih(1)/pih(2)));
theta2 = 90 + rad2deg(atan(pih(1)/pih(3)));

%Computation of alpha
piL=P.'*L1';
thetaL =  180 + rad2deg(atan(piL(1)/piL(3)));
alpha = thetaL - theta2;
