clear all

%% Original Transfer Function Matrix 
% Gain matrix
K = [-0.098 -0.036 -0.014;
     -0.043 -0.092 -0.011;
     -0.012 -0.016 -0.102];
R = (inv(K))';
T = [122 149 158;
     147 130 156;
     153 151 118];
L = [17 27 32;
     25 16 33;
     31 34 16];
K_N = K./(T + L)
A = K.* R
A_N = K_N.* inv(K_N')
I = A_N./ A

X = [ A_N(1,1)/A_N(1,1) A_N(1,2)/A_N(1,1) A_N(1,3)/A_N(1,1);
      A_N(2,1)/A_N(2,2) A_N(2,2)/A_N(2,2) A_N(2,3)/A_N(2,2);
      A_N(3,1)/A_N(3,3) A_N(3,2)/A_N(3,3) A_N(3,3)/A_N(3,3)]

%% Equivalent Transfer Function
Kp = K./ A; % Equivalent Gain
Tp = I.* T; % Equivalent Time Constant
Lp = I.* L; % Equivalent Time Delay


for k = 1:3
    for m = 1:3
        G(k,m) = tf([K(k,m)],[T(k,m), 1]);
        G(k,m).iodelay = L(k,m);
    end
end

allmargin_ = [allmargin(-G(1,1)),allmargin(-G(2,2)),allmargin(G(3,3))];
GM = [-allmargin_(1).GainMargin(1),-allmargin_(2).GainMargin(1),allmargin_(3).GainMargin(1)];
GMF = [allmargin_(1).GMFrequency(1),allmargin_(2).GMFrequency(1),allmargin_(3).GMFrequency(1)];

K_ZN = GM./2.2;
T_I_ZN = 2*pi./(1.2*GMF);
T_D_ZN = 2*pi./(8*GMF);
    
min_error=999;
for F=2.0:0.05:3
    K_C = K_ZN/F;
    T_I = F*T_I_ZN;
    K_I = K_C./T_I;
    max_Lc = 0;
    Gc=[tf([K_C(1) K_I(1)],[1 0]),0,0;
        0,tf([K_C(2) K_I(2)],[1 0]),0;
        0,0,tf([K_C(3) K_I(3)],[1 0])];
    for w=0.1:0.02:2.5
        W = -1+det(eye(3)+freqresp(G*Gc,w));
        Lc = 20*log10(abs(W/(1+W)));
        if(Lc>max_Lc)
             max_Lc=Lc;
        end
    end
    error = abs(max_Lc-6);
    if(error<min_error)
    min_error=error;
    F_match=F;
    end 
end
 
K_C = K_ZN./F_match;
T_I = F_match.*T_I_ZN;
K_I = K_C./T_I;

min_error_2=999;
for FD=2.0:0.5:20
    T_D = T_D_ZN/FD;
    K_D=K_C.*T_D;
    max_Lc_2 = 0;
    Gc=[tf([K_C(1) K_I(1)],[1 0]),0,0;
        0,tf([K_D(2) K_C(2) K_I(2)],[0 1 0]),0;
        0,0,tf([K_D(3) K_C(3) K_I(3)],[0 1 0])];
    for w=0.1:0.05:2.5
        W = -1+det(eye(3)+freqresp(G*Gc,w));
        Lc = 20*log10(abs(W/(1+W)));
        if(Lc>max_Lc_2)
             max_Lc_2=Lc;
        end
    end
    error = abs(max_Lc_2-6);
    if(error<min_error_2)
    min_error_2=error;
    F_match_D=FD;
    end 
end
 
K_C = K_ZN./F_match;
T_I = F_match.*T_I_ZN;
K_I = K_C./T_I;
T_D = T_D_ZN/F_match_D;
K_D = K_C.*T_D;
