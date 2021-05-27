%% PARABOLIC PDEs SOLVER FOR THE COURSE "NUMERICAL METHODS" MAIN CODE
%  BY Alessandro Arduino

%% CLEAN UP

clear all
close all
clc

%% DATA INPUT

minx = -1;
stepx = 1;
maxx = 1;

miny = -1;
stepy = 1;
maxy = 1;

mint = 0;
stept = 0.05;
maxt = 50;

theta = 0.5;

S = 10;

u_fun = @(x,y,t) cos(2*x^2-y)+y^2; % choosen function for the 
f_longterm_fun = @(x,y) -(-4*sin(2*x^2-y)-(16*x^2+1)*cos(2*x^2-y)+2); % calculated f function
q_fun = @(x,y) -sin(2*x^2-y)*4*x; % dirichlet boundary condition
%% CALCULATIONS

x_val = minx:stepx:maxx; % calculate x coordinates
y_val = miny:stepy:maxy; % calculate y coordinates 
t_val = mint:stept:maxt; % calculate time vector

[nodes,topology,S_vect] = createSimpleMesh(x_val,y_val,S);

u = zeros(length(x_val),1);
f = zeros(length(x_val),1);

for i = 1:size(nodes,1)
    xi = nodes(i,1);
    yi = nodes(i,2);
    u(i) = u_fun(xi,yi,0);
    f(i) = f_longterm_fun(xi,yi);
end

u0 = ones(length(u),1);

[P,H] = PHcalc(nodes,topology,S_vect);
f = f_scalprod(f,nodes,topology);

%% BOUNDARY CONDITIONS
%
% The way the boundary conditions are created is automatic as it is due to
% the type decided for the experiment: the upper and lower portion of the
% boundary will have dirichlet while the left and right sides will have
% neumann
%
% Dirichlet 
full(H)
[fstationary,Hstationary] = dirichletBoundary(H,f,x_val,y_val,u); 

% Neumann

f = neumanntBoundary(f,x_val,y_val,q_fun,nodes)';
fstationary = neumanntBoundary(fstationary,x_val,y_val,q_fun,nodes)';

% SOLUTION OF THE STEADY STATE PROBLEM
epsilon = 10^(-50);
kmax = 10^90;
[ustationary,rf,kf,rk] = PCG_Jacobi(Hstationary,fstationary,epsilon,kmax,u0);

u_time = zeros(length(ustationary),length(t_val));
u_time(:,1) = ustationary;

for i = 1:(length(t_val)-1)
    
    u_time(:,i+1) = thetaMethod(H,P,u_time(:,i),f,theta,stept,epsilon,kmax,x_val,y_val,ustationary,q_fun,nodes);
    
end

