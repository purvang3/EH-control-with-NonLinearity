
function zdot=cyl_dyn2(t,z)
global D_p C_d x_db p_relief tau_relief K_relief ; 
global p_max p_min p_pre V_disch K_acc C_acc; 
global A_a A_b l_cyl mass1 c ; 
global beta V_hose_pv V_hose_va V_hose_vb ; 
global Q_p Q_pa Q_pb Q_at Q_bt Q_pt Q_r Q_acc ;
global p_p p_a p_b p_t p_acc V_acc ; 
global y ydot ; 
global w_pump ; 
global x_s F_load ;
y = z(1) ;  % Cylinder Position
ydot = z(2) ; % Velocity
p_a = z(3) ;  % A side pressure
p_b = z(4) ;  % B side pressure
p_p = z(5) ;  % pump pressure  
p_acc = z(6) ; % Accumulator Pressure
V_acc = z(7) ; %Accumulator Volume


if (abs(x_s) >= x_db ) 
    A_pa = ((20*10^-6)/(100-x_db))* (abs(x_s) - x_db) ;
    A_pb = (( 10*10^-6)/(100-x_db))* (abs(x_s) - x_db) ;
    A_at = ((40*10^-6)/(100-x_db))* (abs(x_s) - x_db) ; 
    A_bt = ((10*10^-6)/(100-x_db))* (abs(x_s) - x_db) ;
else
    A_pa = 0.0 ; 
    A_pb = 0.0 ;
    A_at = 0.0 ; 
    A_bt = 0.0 ; 
end

% Flow rates
Q_p = D_p * w_pump ; 
Q_pt = 0.0 ; % closed center valve.
if(x_s >= 0.0 ) 
    Q_pa = C_d * A_pa * (abs(p_p - p_a))^0.5 ; 
    Q_bt = C_d * A_bt * (abs(p_b - p_t))^0.5 ;
    Q_pb = 0.0 ; 
    Q_at = 0.0 ; 
else
    Q_pb = C_d * A_pb * (abs(p_p - p_b))^0.5 ;
    Q_at = C_d * A_at * (abs(p_b - p_t))^0.5 ; 
    Q_pa = 0.0 ; 
    Q_bt = 0.0 ; 
end

% Case 2: with Accumulator
if  V_acc <= 0.0 && p_p < p_acc
    Q_acc = 0.0; 
else
    Q_acc = sign(p_p - p_acc) * K_acc * (abs(p_p - p_acc))^0.5 ; 
end

% Non-Ideal relief valve, no transient dynamics
if  p_p < p_relief
    Q_r = 0.0 ; 
else
    Q_r = K_relief * (p_p - p_relief) ;
end
zdot=zeros(7,1) ;

% ODEs...
zdot(1) = ydot ;

zdot(2) = (1/mass1)*(-c * ydot + p_a * A_a - p_b * A_b - F_load) ;
if x_s >= 0.0
    zdot(3) = (beta/(V_hose_va + y * A_a))*(Q_pa - ydot * A_a) ; 
    zdot(4) = (beta/(V_hose_vb + (l_cyl - y) * A_b))*(-Q_bt + ydot * A_b);
    zdot(5) = (beta/(V_hose_pv + V_acc))*(Q_p-(Q_pa + Q_pt + Q_r + Q_acc)); 
else
    zdot(3) = (beta/(V_hose_va + y * A_a))*(-Q_at - ydot * A_a) ;
    zdot(4) = (beta/(V_hose_vb + (l_cyl - y) * A_b))*(Q_pb + ydot * A_b); 
    zdot(5) = (beta/(V_hose_pv + V_acc))*(Q_p-(Q_pb + Q_pt + Q_r + Q_acc)); 
end
zdot(6) = (1/C_acc)* Q_acc ;
zdot(7) = Q_acc;
return; 