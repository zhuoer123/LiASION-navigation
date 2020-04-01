%% Compute a pole orbit around earth
clc;clear
close all

%% character creation

earth = celestial_body;
earth.radius = 6378.13643; % [km]
earth.mu = 3.986e5; % [km^3/s^2]
 
moon = celestial_body;
moon.radius = 1738; % [km]
moon.mu = 4902.799; % [km^3/s^2]

sun = celestial_body;
sun.radius = 696000; % [km]
sun.mu = 1.327e11; % [km^3/s^2]

mu_ME = moon.mu/(moon.mu + earth.mu);
mu_SE = earth.mu/(earth.mu + sun.mu);

L_EM = 384400; % [km], Earth-Moon distance
L_SE = 149598023; % [km], Sun-Earth distance

L_points = lagrangePoints(mu_ME);
x_L1 = L_points(1,1);
x_L2 = L_points(1,2);

G = 6.67408e-20;

%%
% Normalization stuff (divide by these to get normalized versions)
LU = L_EM;
T_SU = 2*pi*sqrt(LU^3/(earth.mu+moon.mu));
DU = L_EM; % [km]
TU = 1/(2*pi/T_SU); % [s]
VU = DU/TU; % [km/s]
AU = DU/TU^2; % [km/s^2]

normalizers = struct('time_norm',TU,'vel_norm',VU,'accel_norm',AU);


% Find ICs for L1 and L2 Lyapunov orbits with same Jacobi Constant
% Use bisection to find orbit with matching Jacobi constant

% Numerical parameters
eps = 5e-2; % epsilon bound for node placement event

%% First (reference) halo orbit
L_km = L_EM;
% Az_L1 = 200000; %[km]
Az_L1 = 30000; %[km]

halo_IC_L2_southern = halo_computeplot(mu_ME, Az_L1, L_km, "L2", "north", 0);
fprintf("First orbit found.\n")
%%%% Integrate and plot both orbits

ode_opts = odeset('RelTol',5e-14,'AbsTol',1e-15);

T_L2 = halo_IC_L2_southern{2};
X0_L2 = halo_IC_L2_southern{1};
[~, X_hist_L2] = ode113(@(t,X) CR3BP(t,X,mu_ME), [0 T_L2], X0_L2, ode_opts);

halo1_JC = jacobi_constant(halo_IC_L2_southern{1},mu_ME);
%% 6 elements

  
Rm = 1735.97; %km
coe = [1000+Rm;0;90;0;0;0];
%%%    [a;e;i;Omiga;w;theta]

[r,v] = coe2state(coe,moon.mu);

% Integrate
x0 = [r;v]';
T = 2*pi*sqrt((coe(1))^3/moon.mu);
epoch = 10;  %%²ÉÑù¼ä¸ô
[t,xp] = ode113('f',0:epoch:T*2,x0);

%%%% Compute satellitte's Position and Velocities in Rotation Frame 
P = zeros(round(T/epoch),3);
V = zeros(round(T/epoch),3);
for i = 1:round(T/epoch)
    time(i,:) = juliandate([2012,3,5,0,0,0+i]);
    [P(i,:),V(i,:)] = planetEphemeris(time(i,:),'Earth','Moon');  
end
Xp = state2coe(P,V,T,epoch,moon.mu,xp,L_km,mu_ME);



%% plot

scatter3(x_L2*L_km, 0, 0, 'd', 'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','L2'); hold on
plot3(Xp(:,1),Xp(:,2),Xp(:,3),'b','DisplayName','Moon pole orbit'); hold on
plot3((1-mu_ME)*L_km, 0, 0, 'ok', 'markerfacecolor', 'y', 'markersize', 10, 'DisplayName', 'Moon'); hold on % Smaller primary
plot3(X_hist_L2(:,1)*L_km, X_hist_L2(:,2)*L_km, X_hist_L2(:,3)*L_km, 'r-','DisplayName', '$$L_2$$ Halo Orbit'); hold on
legend()
xlabel('$$x [km]$$')
ylabel('$$y [km]$$')
zlabel('$$z [km]$$')
grid on
% axis equal

%title("Initial and target halo orbits")

%%%% Some codes refer to Saichikine's project. 
