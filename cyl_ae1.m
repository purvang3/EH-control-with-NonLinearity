%%% cyl_ae1.m

function F=cyl_ae1(x)
global D_p C_d x_db p_relief K_relief ; 
global A_a A_b l_cyl mass1 c ; 
global Q_p Q_pa Q_pb Q_at Q_bt Q_pt Q_r ;
global p_p p_a p_b p_t ;
global y ydot ; 
global w_pump ; 
global x_s F_load ;
ydot = x(1);  % velocity of cylinder
p_a = x(2);   % pressure on A side
p_b = x(3);   % pressure on B side
p_p = x(4);   % pump pressure

% Valve dynamics while considerinng Deadband and respecteed Area
if (abs(x_s) >= x_db ) 
    A_pa = ((20*10^-6)/(100-x_db))* (abs(x_s) - x_db) ;
    A_pb = ((10*10^-6)/(100-x_db))* (abs(x_s) - x_db) ;
    A_at = ((40*10^-6)/(100-x_db))* (abs(x_s) - x_db) ; 
    A_bt = ((10*10^-6)/(100-x_db))* (abs(x_s) - x_db) ; 
else
    A_pa = 0.0 ; 
    A_pb = 0.0 ;
    A_at = 0.0 ;
    A_bt = 0.0 ; 
end

% Flow rates for individual Section
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
Q_acc = 0.0 ;

% Ideal relief valve, without transient dynamics 

if (p_p < p_relief)
    Q_r = 0.0 ; 
    if x_s >= 0.0
        F = [ (- c * ydot + p_a * A_a - p_b *A_b - F_load) ; 
            (Q_pa - ydot * A_a) *10^6 ; (-Q_bt + ydot * A_b)*10^6; (Q_p-(Q_pa + Q_pt + Q_r + Q_acc))* 10^6] ; 
    else
        F = [ (- c * ydot + p_a * A_a - p_b *A_b - F_load) ; 
            (-Q_at - ydot * A_a)* 10^6 ; (Q_pb + ydot * A_b)* 10^6 ;
            (Q_p-(Q_pb + Q_pt + Q_r + Q_acc))* 10^6 ] ; 
    end
else
    Q_r = K_relief * (p_p - p_relief) ; 
    if x_s >= 0.0 
        F = [(- c * ydot + p_a * A_a - p_b *A_b - F_load); 
            (Q_pa - ydot * A_a)* 10^6 ; (-Q_bt + ydot * A_b)* 10^6 ; 
            (Q_p-(Q_pa + Q_pt + Q_r + Q_acc))* 10^6 ] ; 
    else
        F = [ (- c * ydot + p_a * A_a - p_b *A_b - F_load) ;
            (-Q_at - ydot * A_a)* 10^6 ; (Q_pb + ydot * A_b)* 10^6 ; 
            (Q_p-(Q_pb + Q_pt + Q_r + Q_acc))* 10^6 ] ;
    end
end
return;