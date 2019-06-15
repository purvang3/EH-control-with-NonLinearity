
display("Case-2") 
 % For CASE 2 %
global D_p C_d x_db p_relief tau_relief K_relief ; 
global p_max p_min p_pre V_disch K_acc C_acc; 
global A_a A_b l_cyl mass1 c ; 
global beta V_hose_pv V_hose_va V_hose_vb ; 
global Q_p Q_pa Q_pb Q_at Q_bt Q_pt Q_r Q_acc ; 
global p_p p_a p_b p_t p_acc V_acc ; global y ydot ; 
global w_pump ; 
global x_s F_load ;
% Parameters
D_p = 0.0001 ; % Pump displacement volume per revolution

C_d = 0.065 ; % Valve parameters
x_db = 10 ; % Valve deadband in percentage of total travel
% Relief valve parameters
p_relief = 20 * 10^6 ; % [N/m^2] = Pa 
tau_relief = 0.025 ; % time constant of relief valve dynamics
K_relief = 0.01*10^-6 ;
% Accumulator parameters
p_max = 20*10^6 ; % [N/m^2] = Pa 
p_min = 15*10^6 ; % [N/m^2] = Pa 
p_pre = 15*10^6 ; % [N/m^2] = Pa 
V_disch = 0.005 ; % [m^3] 
V_acc = 0.0 ; % initial fluid volume in the accumulator 
K_acc = 1.0*10^-6 ; 
C_acc = V_disch/(p_max - p_min) ;

% Cylinder parameters
A_a = 0.01 ; 
A_b = 0.005 ; 
l_cyl = 1.0 ; % [m] 
mass1 = 10000 ; % [kg] 
c = 10.0 ; % [N/(m/sec)]
% Fluid parameters and hose volumes
beta = 15.0*(10^8) ; % Bulk modulus 
V_hose_pv = 0.0001 ; % [m^3] 
V_hose_va = 0.0001 ; % [m^3] 
V_hose_vb = 0.0001 ; % [m^3]
% Input conditions
w_pump = 25 ; % Pump input shaft speed; [rev/sec]
F_load = mass1 * 9.81 ; % kg m/sec^2 = N
% Initial conditions on head and rod-end volume of cylinder.
y = 0.1 ; % Initial cylinder position 
ydot = 0.0 ; % Initial cylinder velocity
p_p = 0.0 ; % Pump pressure 
p_a = F_load/A_a ; 
p_b = 0.0 ;
p_acc = p_min ; 
V_acc = 0.0 ; 
Q_r = 0.0 ;
p_t = 0.0 ; % Tank pressure

t_0 =0.0; 
t_f =5.0 ; 
t_sample = 0.001 ; 
z=zeros(7,1);
z(1) = y ; % Position of cylinder
z(2) = ydot ; % Velocity of cylinder 
z(3) = p_a ; % Pressure at side A
z(4) = p_b ; % Pressure at side B
z(5) = p_p;  % Pump pressure
z(6) = p_acc ; % Accumulator Pressure
z(7) = V_acc ; % Accumulator Volume

z_out=[Q_p Q_pa Q_pb Q_at Q_bt Q_r Q_acc x_s F_load/1000 z(1) z(2) z(3) z(4) z(5) z(6)/10^7 z(7)*1000 ] ; % Output Variables
for (t=t_0: t_sample:t_f)
    % Valve Position as function of trapezoidal signal
if t<1.0 
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

% Solve ODEs...
t_span=[t,t+t_sample] ; 
[T,z1] = ode45('cyl_dyn2',t_span, z); 
[m,n]=size(z1);
z(:)=[z1(m,:)] ;
y = z(1) ;
ydot = z(2) ;
p_a = z(3) ; 
p_b = z(4) ;
p_p = z(5) ; 
p_acc = z(6) ;
V_acc = z(7) ;
z_out=[z_out;
       Q_p Q_pa Q_pb Q_at Q_bt Q_r Q_acc x_s F_load/1000 z(1) z(2) z(3) z(4) z(5) z(6)/10^7 z(7)*1000 ] ;
end
[m,n]=size(z_out);
t_inc = (t_f-t_0)/m ; 
tout=t_0:t_inc:t_f-t_inc; 
tout = tout' ;

figure(2) ; 
subplot(3,2,1) ;
plot(tout, z_out(:,1), 'k',tout, z_out(:,2), 'b',tout, z_out(:,5), 'm', tout, z_out(:,6), 'c',tout, z_out(:,7), 'g'); 
xlabel('Time (sec)') ; ylabel('Flow rate (m^3/sec)') ;
legend('Q_p','Q_{pa}','Q_{bt}','Q_r', 'Q_{acc}');
subplot(3,2,2) ; 
plot(tout, z_out(:,8), 'k',tout, z_out(:,9), 'b');
xlabel('Time (sec)') ; ylabel('Spool position and External Load') ; 
legend('x_s','F_{load}');
subplot(3,2,3) ; 
plot(tout, z_out(:,10) , 'k',tout, z_out(:,11), 'b'); 
xlabel('Time (sec)') ; ylabel('Cylinder position and velocity') ; 
legend('y','ydot');
subplot(3,2,4) ; 
plot(tout, z_out(:,12), 'k',tout, z_out(:,13), 'b',tout, z_out(:,14), 'g'); 
xlabel('Time (sec)') ; 
ylabel('Pressure [Pa]') ; 
legend('p_a','p_b','p_s');
subplot(3,2,5) ;
plot(tout, z_out(:,15), 'k',tout, z_out(:,16), 'b'); 
xlabel('Time (sec)') ;ylabel('Accumulator state') ;
legend('Pressure/10^7','Fluid Volume (liter)');

