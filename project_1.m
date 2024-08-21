%% clearvars; clc
time = 0.01;
%% Plot style options
bopts = bodeoptions('cstprefs');
bopts.PhaseUnits = 'deg';
bopts.grid = 'on';
bopts.Title.Interpreter = 'latex';
bopts.Title.FontSize=20;
bopts.XLabel.Interpreter = 'latex';
bopts.XLabel.FontSize = 14;
bopts.YLabel.Interpreter = 'latex';
bopts.YLabel.FontSize = 14;
bopts.PhaseMatching = 'on';        %Match phase
bopts.PhaseMatchingFreq = 0;       %at 0 [rad/sec]
bopts.PhaseMatchingValue = -90;    % to -90 degrees

nopts = nicholsoptions('cstprefs');
nopts.PhaseUnits = 'deg';
nopts.grid = 'on';
nopts.Title.Interpreter = 'latex';
nopts.Title.FontSize=20;
nopts.XLabel.Interpreter = 'latex';
nopts.XLabel.FontSize = 14;
nopts.YLabel.Interpreter = 'latex';
nopts.YLabel.FontSize = 14;
nopts.TickLabel.FontWeight = 'bold';
nopts.PhaseMatching = 'on';        %Match phase
nopts.PhaseMatchingFreq = 0;       %at 0 [rad/sec]
nopts.PhaseMatchingValue = -90;    % to -90 degrees

topts = timeoptions('cstprefs');
topts.grid = 'on';
topts.Title.Interpreter = 'latex';
topts.Title.FontSize=20;
topts.XLabel.Interpreter = 'latex';
topts.XLabel.FontSize = 14;
topts.YLabel.Interpreter = 'latex';
topts.YLabel.FontSize = 14;

newcolors = [62,94,168; 112,143,57; 120,5,5; 217,142,50; 114,34,130; 217,89,149]/256;
set(0,'DefaultAxesColorOrder',newcolors)
%% 2
% FINF THE "TF" FUNCTION-FROM U TO PHI AND U TO TETHA

% Define the symbolic variables
syms T TETHA PHI U S K_a  L_a R_a K_b m r l J_b J_p zeta g V_b f n_g

% Write the equations
eq1 = K_b*n_g/(L_a*S + R_a)*(K_a*U-V_b) == T; 
eq2 = K_b*PHI*S*n_g == V_b ; 
eq3 = (-m*r*l*S^2 /((J_p+m*l^2)*S^2 +zeta*S+m*l*g))*PHI == TETHA;
eq4 = TETHA *(m*r*l*S^2) + PHI*((J_b+m*r^2)*S+f)*S ==T;

% Solve the system of equations for T, TETHA, PHI, V_b
solutions = solve([eq1, eq2, eq3, eq4], [T, TETHA, PHI, V_b]);

% Display the solutions
sol_PHI_TO_TETHA=solutions.PHI/TETHA;
sol_T = solutions.T;
sol_TETHA = solutions.TETHA/U;
sol_PHI = solutions.PHI/U;
%sol_vb = solutions.vb

sol_PHI_1 = subs(sol_PHI, {K_a L_a R_a K_b m r l J_b J_p zeta g f n_g},{12 0 3 0.0243 0.04 0.33 0.22 0.01 0.0003 0.000375 9.81 0.01 5.2});
simplify(sol_PHI_1)

% Define the coefficients of sol_PHI_1

A = 706352400000;
B = 118462500000;
C = 27271015200000;
D = 14791750000;
E = 24777578272;
F = 778169140125;
G = 826713790056;

% Define the numerator and denominator coefficients
numerator_coeffs = [A, B, C];
denominator_coeffs = [D, E, F, G,0];

% Create the transfer function
boda_tf_PHI = tf(numerator_coeffs, denominator_coeffs);

% Plot the frequency response (Bode plot) of the transfer function
fig01 = figure(1); clf

P1 = boda_tf_PHI;


%find the tf from u to phi when zeta limit as inf
limit_inf = limit(sol_PHI, zeta, inf); % Evaluates the limit of f as x approaches infinity
limit_inf_1 = subs(limit_inf, {K_a L_a R_a K_b m r l J_b J_p  g f n_g},{12 0 3 0.0243 0.04 0.33 0.22 0.01 0.0003 9.81 0.01 5.2});
simplify(limit_inf_1);
disp(limit_inf_1) % Output should be 0
P2 = tf(9477/6250,[10767/250000,28729281/625000000,0]);

%%  3

%plot of both tf from u to phi and limit_tf from u to phi
figure;
bode(P1,P2,bopts);
title('Bode Plot of L ${\phi}$','interpreter','latex');
ylabel('${\phi}(t)$','interpreter','latex')
legend('$ P\_{1}(s)$','$ P\_{2}(s)$','fontsize',14,'interpreter','latex')
saveas(gcf, 'Bode Plot of L of ϕ.png');
grid on;


%transfer function from u to theta

sol_TETHA_1= subs(sol_TETHA, {K_a L_a R_a K_b m r l J_b J_p zeta g f n_g},{12 0 3 0.0243 0.04 0.33 0.22 0.01 0.0003 0.000375 9.81 0.01 5.2});
figure;
bode(tf_form(sol_TETHA_1));
grid on;
title('Bode Plot of transfer function from u to ${\theta}$','interpreter','latex');
print(gcf, 'Bode Plot of transfer function from u to theta', '-dpng');


%transfer function from phi to theta

sol_PHI_to_tetha =-m*r*l*S^2 /((J_p+m*l^2)*S^2 +zeta*S+m*l*g);
sol_PHI_to_tetha_subs= subs(sol_PHI_to_tetha, {K_a L_a R_a K_b m r l J_b J_p zeta g f n_g},{12 0 3 0.0243 0.04 0.33 0.22 0.01 0.0003 0.000375 9.81 0.01 5.2});
%simplify(sol_PHI_to_tetha_subs)

TF_sol_PHI_to_tetha_subs = tf_form(sol_PHI_to_tetha_subs);
figure;
bode(TF_sol_PHI_to_tetha_subs);
set(gca, 'XLim', [10^-2 10^2]);

title('Bode Plot of transfer function from ${\phi}$ to ${\theta}$','interpreter','latex');
print(gcf, 'Bode Plot of transfer function from phi to theta', '-dpng');
grid on;

%% 4-a,b,c
% THE CONTROL OF PHI 
% Define the function from laplas to forie using the coefficients

figure
margin(P1);grid on;
figure
nichols(P1);grid on;
legend('P');
print(gcf, 'nichols_P', '-dpng');

Wgk = 55;
tf_phi_u_forie = @(W) ((-A*W.^2 + B*W*1i + C)) ./  (D*W.^4 - 1i*E*W.^3 -F*W.^2 + G*W*1i);
k= 1/abs(tf_phi_u_forie(Wgk)); 
figure
nichols(k*P1);grid on;
legend('KP');
print(gcf, 'nichols_KP', '-dpng');

Wg2 = 12; %A low transition frequency was chosen so ...
% that there would be no amplification of frequencies in the graph of the 3 decibels
c2 = (10*S+Wg2)/(10*S);
tf_c2 =tf([10,Wg2],[10, 0]);
tf_c2_55=tf([10,55],[10, 0]);

figure
margin(k*P1*tf_c2);grid on;
figure
nichols(k*P1*tf_c2,k*P1*tf_c2_55);grid on;
legend('kP1c2','kP1c2*');
print(gcf, 'nichols_KPC2', '-dpng');


% Compute the phase margin 
[~, PM] = margin(k*P1*tf_c2);
% Display the phase margin
disp(['Phase Margin: ', num2str(PM), ' degrees'])

Wg3 = 52;
n=1;
c3_angle_deg = (45+PM+1.7); %wg=55  t under 3db 
c3_angle_n = deg2rad(c3_angle_deg/n);
c3_alpha = (1+sin(c3_angle_n))/(1-sin(c3_angle_n));
c3 = (sqrt(c3_alpha)*S+Wg3)/(S+sqrt(c3_alpha)*Wg3);

tf_c3 =tf([sqrt(c3_alpha),Wg3],[1, sqrt(c3_alpha)*Wg3]);

figure
margin(k*P1*tf_c2*tf_c3);grid on;
figure
nichols(k*P1*tf_c2*tf_c3,k*P1*tf_c2_55*tf_c3);grid on;
legend('kP1c2c3','kP1c2*c3');
print(gcf, 'nichols_KPC2C3', '-dpng');

%to chek the mg of L is 45 as we want

% Given symbolic expression for 'L_phi'


L_phi = k*c2*c3^n*sol_PHI_1;

% Create transfer function 'TF'
TF_L_phi = tf_form(L_phi);
figure
bode(TF_L_phi)
grid on
title ('bode diagram of L')
% Compute the phase margin 
[~, ~,freq_0dB_TF_L_phi1,freq_0dB_TF_L_phi2] = margin(TF_L_phi); 
print(gcf, 'bode of L', '-dpng');


% Display the phase margin
disp(['Phase Margin: ', num2str(PM), ' degrees'])
% Get the Bode plot data
[mag, phase, w] = bode(TF_L_phi);


% Find the frequency where the magnitude is 0 dB
mag_dB = squeeze(20*log10(mag));  % Squeeze to get rid of singleton dimensions
[~, idx] = min(abs(mag_dB));
freq_0dB_TF_L_phi3 = w(idx);  % Find frequency where magnitude is 0 dB

figure
nichols(TF_L_phi);
grid on;
[~, PM_TF] = margin(TF_L_phi);
% Display the phase margin
disp(['Phase Margin_L_PHI_FINAL: ', num2str(PM_TF), ' degrees'])
disp(['wg_TF_L_phi: ', num2str(freq_0dB_TF_L_phi1), ' [rad/sec]'])
disp(['wg_TF_L_phi: ', num2str(freq_0dB_TF_L_phi2), ' [rad/sec]'])
disp(['wg_TF_L_phi: ', num2str(freq_0dB_TF_L_phi3), ' [rad/sec]'])

%% 4-d%% 
% now tf T
T_phi = L_phi/(1+L_phi);

% Create transfer function 'TF'
TF_T_phi = tf_form(T_phi);
figure
bode(TF_T_phi)

title ('T(S)')
hold on 

h = findobj(gcf, 'Type', 'Axes');
axMag = h(2);
plot(axMag, axMag.XLim, [3 3], '--', 'LineWidth', 0.5,'Color',[0.5 0.5 0.5])
title('Bode Plot of T ${\phi}$','interpreter','latex');
%ylabel('${\phi}(t)$','interpreter','latex');
grid on
print(gcf, 'bode of T', '-dpng');
pause(time)
hold off
figure
nichols(TF_T_phi)
grid on

wb_TF_T_phi = bandwidth(TF_T_phi);
disp(['wb_TF_T_phi: ', num2str(wb_TF_T_phi), ' [rad/sec]'])

%%   4-e
% now tf Sd

S_phi = sol_PHI_1/(1+L_phi);

% Create transfer function 'TF'
TF_S_phi = tf_form(S_phi);
figure
bode(TF_S_phi)
hold on
h = findobj(gcf, 'Type', 'Axes');
axMag = h(2);
plot(axMag, axMag.XLim, [-20 -20], '--', 'LineWidth', 5,'Color',[0.5 0.5 0.5])
title('Bode Plot of Sd ${\phi}$','interpreter','latex');
ylabel('${\phi}(t)$','interpreter','latex');
grid on
saveas(gcf, 'Bode Plot of S_d  of ϕ.png');
hold off
figure
nichols(TF_S_phi)
grid on


%responde of the close loop to distrubance

disturbance = tf(1 ,[0.1 ,  1]);

% Define time vector for simulation
t = linspace(0,5,1001);  % Adjust the time vector as needed

Y_dist = TF_S_phi ;
tetha__dist = Y_dist*TF_sol_PHI_to_tetha_subs;
U_BEFORE_TF =T_phi/sol_PHI_1;

% Create transfer function 'TF'
U_dist = tf_form(U_BEFORE_TF);
 

% Define the disturbance input (step input)
dist_input = ones(size(t));

% Simulate the responses
[y_response, t_out] = lsim(Y_dist*disturbance, dist_input, t);
[u_response, ~] = lsim(U_dist*disturbance, dist_input, t);
[tehtha_response, ~] = lsim(tetha__dist*disturbance, dist_input, t);

% Plot the responses

figure
plot(t_out,step(disturbance,t_out))
hold on
plot(t_out, u_response);
title('Output Respons  $d {\rightarrow} u$  ','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('u(t) [V] ','interpreter','latex');
legend('d(t)','u(t)','interpreter','latex')
hold off
grid on
print(gcf, 'Time respone from d to u ', '-dpng');

figure
plot(t_out,step(disturbance,t_out))
hold on
plot(t_out, y_response);
title('Output Response $d {\rightarrow} {\phi}$  ','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('${\phi(t)}$ [deg]','interpreter','latex');
legend('d(t)','\phi(t)')
hold off
grid on
print(gcf, 'Time respone from d to phi', '-dpng');

figure
plot(t_out,step(disturbance,t_out))
hold on
plot(t_out, tehtha_response);
title(' output Response $d {\rightarrow} {\theta}$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('${\theta(t)} [deg] $','interpreter','latex');
legend('d(t)','\theta(t)')
hold off
grid on
print(gcf, 'Time respone from d to theta', '-dpng');

%% 5

t2 = linspace(0,2,1001);

% FIND THE TIME OPTIMAL CONTROL

Ra1=0.33;
Kb1=0.0243;
ng1=5.2;
f1=0.01;
m1=0.04;
r1=0.33;
%Jb1=0.003;
Jb1=0.01;
Ka1=12;

%the variable

k_subs = Kb1*ng1/f1;
teo_subs = (Jb1+m1*r1^2)/f1;

% Define the given parameters
phi_f = pi;
phi_0 = 0;
i_max = 1;


% Compute t_sw and t_f
a_t_sw = abs(phi_f - phi_0) / (k_subs * i_max);
t_sw = a_t_sw + teo_subs * log(1 + sqrt(1 - exp(-a_t_sw/teo_subs)));
t_f = a_t_sw + 2 * teo_subs * log(1 + sqrt(1 - exp(-a_t_sw/teo_subs)));


s_tf1 = tf('s');

T_i = i_max * (1 - 2 * exp(-t_sw * s_tf1) + exp(-t_f * s_tf1));
G_i_phi = T_i*k_subs/(teo_subs*s_tf1^2+s_tf1); % tetha=0 
G_v_i = T_i*((Kb1^2*ng1^2+Ra1*f1)+(Jb1+m1*r1^2)*Ra1*s_tf1)/(Ka1*(f1 + s_tf1*(m1*r1^2 + Jb1))); 



figure
plot(t2,step(T_i,t2))
title('Output Response i(t) opt','interpreter','latex');
xlabel('Time [s]');
ylabel('i opt [A]');
grid on
print(gcf, 'Output Response i(t) opt', '-dpng');

figure
plot(t2,step(G_v_i,t2))
title('Output Response u(t) opt','interpreter','latex');
xlabel('Time [s]');
ylabel('u opt [V]');
ylim([-1 1])
grid on
print(gcf, 'Output Response u(t) opt', '-dpng');

figure
plot(t2,rad2deg(step(G_i_phi,t2)))
title('Output Response ${\phi (t) }$  as i(t) current','interpreter','latex');
xlabel('Time [s]');
ylabel('\phi [deg]');
grid on
print(gcf, 'Output Response phi(t) opt', '-dpng');
%% 6
% Now return to the system with zeta as in Table 1


REAL_PHI = TF_T_phi ;
REAL_TETHA = REAL_PHI*TF_sol_PHI_to_tetha_subs;
REAL_U_BEFORE_TF =T_phi/sol_PHI_1;

% Create transfer function 'TF'
REAL_U = tf_form(REAL_U_BEFORE_TF);
 

% Define the disturbance input (step input)
input = rad2deg(step(G_i_phi,t));

% Simulate the responses
[REAL_PHI_response, t_out] = lsim(REAL_PHI, input, t);
[REAL_U_response, ~] = lsim(REAL_U, input, t);
[REAL_TETHA_response, ~] = lsim(REAL_TETHA, input, t);

% Plot the responses
figure
plot(t_out, REAL_U_response);
title('Control Signal $r {\rightarrow} u$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('U(t) [V]','interpreter','latex');
grid on
print(gcf, 'Control Signal r to u', '-dpng');


figure
plot(t,input,'--',LineWidth=1)
hold on
plot(t_out, REAL_PHI_response);
title('Output Response  $r {\rightarrow}{\phi}$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('${\phi (t) }$ [deg]','interpreter','latex');
legend('${r_\phi}$','${\phi(t)}$','interpreter','latex')
hold off
grid on
print(gcf, 'Output Response r to phi', '-dpng');

figure
plot(t_out, REAL_TETHA_response);
title('Output Response  $r {\rightarrow} {\theta}$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('\theta(t) [deg]');
grid on
print(gcf, 'Output Response r to theta', '-dpng');

%% 7
%  External loop for dampening tetha


G7 = sol_PHI_to_tetha_subs;
TF_G7 = TF_sol_PHI_to_tetha_subs;
figure
nichols(-TF_G7)
grid on
w_lpf1 = getGainCrossover(TF_G7,1);
lpf1 = 1/(S/w_lpf1+1);
print(gcf, 'nichols chart of p', '-dpng');



tf_g_lpf1 = tf_form(G7*lpf1);
figure
nichols(-tf_g_lpf1)
title ('nichols Chart PC1')
grid on
print(gcf, 'nichols chart of pc1', '-dpng');

w_lpf2_mean2 = getGainCrossover(tf_g_lpf1,1);
phase_left = find_phase(tf_g_lpf1, w_lpf2_mean2(1));
phase_right = find_phase(tf_g_lpf1, w_lpf2_mean2(2));


w_lpf2 = tan(deg2rad(40))*w_lpf2_mean2(2);
lpf2 = 1/(S/w_lpf2+1);
C_7 = lpf1*lpf2;
tf_g_lpf1_lpf_2 = tf_form(G7*C_7);

%w_lpf2_mean2 = getGainCrossover(tf_g_lpf1_lpf_2,1);
phase2_left = find_phase(tf_g_lpf1_lpf_2, w_lpf2_mean2(1));
phase2_right =find_phase(tf_g_lpf1_lpf_2, w_lpf2_mean2(2));

nichols(-tf_g_lpf1_lpf_2)
title ('nichols Chart PC1C2')
grid on
print(gcf, 'nichols chart of pc1c2', '-dpng');

% C(wj)<1 for all w,
figure
h= bodeplot(tf_form(C_7));
title ('bode diagram of controller C7')
grid on
hold on
print(gcf, 'bode diagram of controller C7', '-dpng');

ax = getaxes(h);
ylim(ax(1),[-100 10])
wb_c7 = bandwidth(tf_form(C_7));
disp(['wb_c7: ', num2str(wb_c7), ' [rad/sec]']);


L_PHI_TETHA = G7*C_7;
T_PHI_TETHA = G7/(1-L_PHI_TETHA );

figure
plot(t2,rad2deg(step(tf_form(T_PHI_TETHA),t2)))
hold on
plot (t2,rad2deg(step(TF_G7,t2)))
yline(1,"r--")
yline(-1,"r--")
legend ('close loop','{P_{\theta \phi}}')
grid on

hold off
title('Output Response  r {\rightarrow} {\theta}');
xlabel('Time [s]','interpreter','latex');
ylabel('\theta(t) [deg]');
print(gcf, 'time response of  r to theta', '-dpng');

%% 8
% simulate the response of the system in Fig. 2(c)

P_all = T_phi*sol_PHI_to_tetha_subs;
L_all = P_all*C_7;

G_all_U = REAL_U_BEFORE_TF/(1-L_all);
G_all_PHI = T_phi/(1-L_all);
G_all_TETHA = P_all/(1-L_all);


TF_G_all_U =tf_form(G_all_U);
TF_G_all_PHI =tf_form(G_all_PHI);
TF_G_all_TETHA =tf_form(G_all_TETHA);

figure
bode(-tf_g_lpf1_lpf_2)
grid on

figure
nichols(-tf_g_lpf1_lpf_2);
hold on
nichols(-tf_form(L_all));
grid on
hold off
legend ('P{\theta}{\phi}C{\theta}','P{\theta}rC{\theta}')
% Define the disturbance input (step input)
input = rad2deg(step(G_i_phi,t));
print(gcf, 'nicols chart compering to close loops', '-dpng');

% Simulate the responses

[TF_G_all_U_response, ~] = lsim(TF_G_all_U, input, t);
[TF_G_all_PHI_response, t_out2] = lsim(TF_G_all_PHI, input, t);
[TF_G_all_TETHA_response, ~] = lsim(TF_G_all_TETHA, input, t);

% Plot the responses
figure
plot(t_out2, REAL_U_response);
hold on
plot(t_out2, TF_G_all_U_response);
title('Output Response  $r {\rightarrow} u $','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('U(t) [V]','interpreter','latex');
%legend('u(t)','u(t) (C_{\theta}=0)')
legend('u(t) (C_{\theta}=0)','u(t)')
hold off
grid on
print(gcf, 'Output Response r to u_Q8', '-dpng');

figure
plot(t_out2,input,'--',LineWidth=1)
hold on
plot(t_out2, TF_G_all_PHI_response);
title('Output Response  $r {\rightarrow}{\phi}$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('${\phi (t) }$ [deg]','interpreter','latex');
legend('r_{\phi}(t)','\phi(t)')
hold off
grid on
print(gcf, 'Output Response r to phi_Q8', '-dpng');

figure
plot(t_out2, REAL_TETHA_response);
hold on
plot(t_out2, TF_G_all_TETHA_response);
title('Output Response  $r {\rightarrow}{\theta}$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
ylabel('\theta(t) [deg]');
legend('{\theta(t)} (C_{\theta}=0)','{\theta(t)}')
hold off
grid on
print(gcf, 'Output Response r to theta_Q8', '-dpng');

%% general function

function TF_syms = tf_form(syms)
% Extract coefficients from the symbolic expression 'sol_PHI_to_tetha_subs'
[num, den] = numden(syms);  % Separate numerator and denominator

num_coeffs = sym2poly(num);  % Convert symbolic numerator to polynomial coefficients
den_coeffs = sym2poly(den);  % Convert symbolic denominator to polynomial coefficients

% Create transfer function 'TF'
TF_syms = tf(double(num_coeffs), double(den_coeffs));
end

function phase = find_phase(TF, omega)
[~, phase] = bode(TF, omega);
%phase_deg = phase * 180 / pi; % Convert phase to degrees
%disp(['Phase at ', num2str(omega), ' rad/s: ', num2str(phase), ' degree']);
end


