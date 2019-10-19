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

%% Detuning factor F

for k = 1:3
    for m = 1:3
        Go(k,m) = tf([K(k,m)],[T(k,m), 1]);
        Go(k,m).iodelay = L(k,m);
    end
end

wmax = 1.5; dw = 0.05; q=1;  Lc_max = 1000;
margin_ = [allmargin(-Go(1,1)),allmargin(-Go(2,2)),allmargin(-Go(3,3))];
gaimM = [-margin_(1).GainMargin(1),-margin_(2).GainMargin(1),-margin_(3).GainMargin(1)];
phaseM = [margin_(1).GMFrequency(1),margin_(2).GMFrequency(1),margin_(3).GMFrequency(1)];

K_ZN = gaimM./2.2;
Tzn = 2*pi./(1.2*phaseM);

for F = 2:0.02:3
    p = 1;
    for w = 0:dw:wmax; 
        Kc = K_ZN/F;
        Ti = F * Tzn;
        ti = Kc./Ti;
        %Ku1 = 1/(10^(50.8/20)); Ku2 = 1/(10^(52.4/20)); Ku3 = 1/(10^(50.7/20));
        %Kc1 = Ku1/(2.2*F); Kc2 = Ku2/(2.2*F); Kc3 = Ku3/(2.2*F);
        %wu1 = 0.278; wu2 = 0.296; wu3 = 0.296;
        %ti1 = 2*pi*F/(1.2*wu1); ti2 = 2*pi*F/(1.2*wu2); ti3 = 2*pi*F/(1.2*wu3);
        Gc1 = tf([Kc(1,1), ti(1,1)],[1 0]);
        Gc2 = tf([Kc(1,2), ti(1,2)],[1 0]);
        Gc3 = tf([Kc(1,3), ti(1,3)],[1 0]);
        Gc = [Gc1 0 0;
              0 Gc2 0;
              0 0 Gc3];
        I3 = eye(3, 3); 
        W = -1 + det(I3 + freqresp(Go*Gc,w));
        Lc(1,p) = 20*log10(abs(W/(1+W)));
        Lc_total(1,q) = 20*log10(abs(W/(1+W)));
        p = p+1; q = q+1;
    end
    Lc_m = abs(6-max(Lc))
    if (Lc_m < Lc_max)
        F_m = F
        Lc_max = Lc_m;
    end
    % F = F % Display the progress during iteration
end
plot(Lc_total) % Display the changing process of the detuning

%% Obtain the final parameters for the controller with BLT

Kc_blt = K_ZN/F_m;
Ti_blt = F_m * Tzn;
ti_blt = Kc./Ti_blt;