
display("CASE 1")
global D_p C_d x_db p_relief K_relief ; 
global A_a A_b l_cyl mass1 c ; 
global Q_p Q_pa Q_pb Q_at Q_bt Q_pt Q_r ; 
global p_p p_a p_b p_t ;
global y ydot ;
global w_pump ; 
global x_s F_load ;

% Parameters
D_p = 0.0001 ; % Pump displacement volume per revolution 
C_d = 0.065 ; % Valve parameters
x_db = 10 ; % Valve deadband in percentage of total travel

% Relief valve parameters
p_relief = 20 * 10^6 ; % [N/m^2] = Pa 
K_relief = 0.01 * 10^-6 ;

% Cylinder parameters
A_a = 0.01 ;
A_b = 0.005 ;
l_cyl = 1.0 ; % [m] 
mass1 = 10000 ; % [kg] 
c = 0.0 ; % [N/(m/sec)]

% Input conditions
w_pump = 25 ; % Pump input shaft speed; [rev/sec] 
F_load = mass1 * 9.81 ; % kg m/sec^2 = N

% Initial conditions on head and rod-end volume of cylinder.
y = 0.1 ; % Initial cylinder position
ydot = 0.0 ; % Initial cylinder velocity 
p_p = p_relief ; % Pump pressure 
p_a = F_load/A_a ;
p_b = 0.0 ;
Q_r = D_p * w_pump ; 
p_t = 0.0 ; % Tank pressure

t_0 =0.0; 
t_f =5.0 ; 
t_sample = 0.001 ; 
z=zeros(1);
x =zeros(4) ;
z(1) = y ; 
x(1) = ydot ; 
x(2) = p_a ;
x(3) = p_b ;
x(4) = p_p;
z_out=[Q_p Q_pa Q_pb Q_at Q_bt Q_r x_s F_load/1000 y ydot p_a p_b p_p] ;

% Position of valve at different time as trapezoidal function
for (t=t_0: t_sample:t_f)
if  t<1.0
    x_s = 0.0 ; 
elseif (t>= 1.0 && t<=1.25) 
    x_s = (100/0.25) *(t-1.0) ; 
elseif (t> 1.10 && t<=3.0) 
    x_s = 100 ;
elseif (t> 3.0 && t<=3.25) 
    x_s = 100 - (100/0.25) * (t - 3.0) ;
else
    x_s = 0.0 ;
end

% Solve Algebraic Equations....
% From Mathworks Documents
options = optimset('Algorithm','levenberg-marquardt','Display','off','TolFun',1e-11,'TolX',1e-10); % Minimizing value of x using lm algorithem to get
                                                                                                   % accurate solution (Global Minimum Concept)
x = fsolve('cyl_ae1',x,options) ; % Using above options for finding optimum value for x
ydot = x(1) ;   % velocity of cylinder
p_a = x(2) ;   % pressure on A side
p_b = x(3) ;   % pressure on B side
p_p = x(4) ;   % pump pressure 

% Solve ODEs...
t_span=[t,t+t_sample] ; 
[T,z1] = ode45('cyl_dyn1',t_span, z); % Solving for output
[m,n]=size(z1);
z(:)=[z1(m,:)] ;
y = z(1) ;
z_out=[ z_out
        Q_p Q_pa Q_pb Q_at Q_bt Q_r x_s F_load/1000 y ydot p_a p_b p_p ] ;
    
end

[m,n]=size(z_out);
t_inc = (t_f-t_0)/m ;
tout=t_0:t_inc:t_f-t_inc;
tout = tout' ;

figure(1) ;
subplot(2,2,1) ; 
plot(tout, z_out(:,1), 'k',tout, z_out(:,2), 'b',tout, z_out(:,5), 'm',tout, z_out(:,6), 'c'); 
xlabel('Time (sec)') ; ylabel('Flow rate (m^3/sec)') ; 
legend('Q_p','Q_{pa}','Q_{bt}','Q_r');
subplot(2,2,2) ;
plot(tout, z_out(:,7), 'k',tout, z_out(:,8), 'b'); 
xlabel('Time (sec)') ; ylabel('Spool position and External Load') ; 
legend('x_s','F_{load}');
subplot(2,2,3) ;
plot(tout, z_out(:,9) , 'k',tout, z_out(:,10), 'b');
xlabel('Time (sec)') ; ylabel('Cylinder position and velocity') ; 
legend('y','ydot');
subplot(2,2,4) ; 
plot(tout, z_out(:,11), 'k',tout, z_out(:,12), 'b',tout, z_out(:,13), 'g'); 
xlabel('Time (sec)') ; ylabel('Pressure [Pa]') ; 
legend('p_a’,’p_b’,’p_p');

