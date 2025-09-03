clc
clear all
close all

fid = fopen('exp1_armazen.txt','r');
if fid < 0
  error('Não foi possível abrir data2.txt');
end

% lê as 12 colunas como números (%f), tratando 'X' como NaN
C = textscan(fid, ...
    repmat('%f',1,12), ...          % '%f%f…%f' 12 vezes
    'Delimiter', ',', ...           % separador vírgula
    'CollectOutput', true, ...      % junta tudo em C{1}
    'TreatAsEmpty', {'X'}, ...      % todo "X" vira vazio
    'EmptyValue', NaN);             % valor vazio = NaN

fclose(fid);

% separa a matriz [n×12] em vetores “plano”
dados   = C{1};
x_real       = dados(:,1);
y_real       = dados(:,2);
lv      = dados(:,3);
rv      = dados(:,4);
ID      = dados(:,5);
hamming = dados(:,6);
tx      = dados(:,7);
ty      = dados(:,8);
tz      = dados(:,9);
yaw     = dados(:,10);
pitch   = dados(:,11);
roll    = dados(:,12);


% Processamento de ruído e plot como antes
n = length(x_real);
sigma_left  = 0;
sigma_right = 0;

leftEncoder  = lv;
rightEncoder = rv;

ruido_left  = sigma_left  * randn(n,1);
ruido_right = sigma_right * randn(n,1);

leftEncoder  = ((leftEncoder  + ruido_left)  * 0.0975);
rightEncoder = ((rightEncoder + ruido_right) * 0.0975);

velocidadeEncoder = (leftEncoder+rightEncoder)/2;
angularEncoder = (-leftEncoder+rightEncoder)/(0.381);  % 0.381

%%figure
%%hold on
%%plot(time, leftEncoder,  'r', 'LineWidth', 2)
%%plot(time, rightEncoder, 'g', 'LineWidth', 2)
%%grid on
%%xlabel('Tempo (s)')
%%ylabel('Velocidade (m/s)')
%%legend('Left','Right')


dt = 0.05;
%%N = length(velocidadeEncoder);
%%
%%xMedido = zeros(N,1);
%%yMedido = zeros(N,1);
%%theta   = zeros(N,1);


% posição inicial
xMedido(1) = x_real(1);
yMedido(1) = y_real(1);
theta(1)   = pi;

%% Parâmetros do método RK4
a1  = 1/6;
a2  = 2/6;
a3  = 2/6;
a4  = 1/6;
q21 = 1/2;
q31 = 0;
q32 = 1/2;
q41 = 0;
q42 = 0;
q43 = 1;
dt  = 0.05;

N = length(velocidadeEncoder);

for i = 2:N
    % estados iniciais neste passo
%%    x0 = xMedido(i-1);
%%    y0 = yMedido(i-1);
%%    t0 = theta(i-1);
%%    v  = velocidadeEncoder(i-1);
%%    w  = angularEncoder(i-1);
%%
%%    % cálculo de F1
%%    F1x = v * cos(t0);
%%    F1y = v * sin(t0);
%%    F1t = w;
%%
%%    % cálculo de F2 usando q21
%%    t2   = t0 + q21*F1t*dt;
%%    x2   = x0 + q21*F1x*dt;
%%    y2   = y0 + q21*F1y*dt;
%%    % entrada v,w permanece a mesma (assumindo u constante no passo)
%%    F2x = v * cos(t2);
%%    F2y = v * sin(t2);
%%    F2t = w;
%%
%%    % cálculo de F3 usando q31, q32
%%    t3   = t0 + q31*F1t*dt + q32*F2t*dt;
%%    x3   = x0 + q31*F1x*dt + q32*F2x*dt;
%%    y3   = y0 + q31*F1y*dt + q32*F2y*dt;
%%    F3x = v * cos(t3);
%%    F3y = v * sin(t3);
%%    F3t = w;
%%
%%    % cálculo de F4 usando q41, q42, q43
%%    t4   = t0 + q41*F1t*dt + q42*F2t*dt + q43*F3t*dt;
%%    x4   = x0 + q41*F1x*dt + q42*F2x*dt + q43*F3x*dt;
%%    y4   = y0 + q41*F1y*dt + q42*F2y*dt + q43*F3y*dt;
%%    F4x = v * cos(t4);
%%    F4y = v * sin(t4);
%%    F4t = w;
%%
%%    % combinações finais
%%    xMedido(i) = x0 + dt*(a1*F1x + a2*F2x + a3*F3x + a4*F4x);
%%    yMedido(i) = y0 + dt*(a1*F1y + a2*F2y + a3*F3y + a4*F4y);
%%    theta(i)   = t0 + dt*(a1*F1t + a2*F2t + a3*F3t + a4*F4t);

    xMedido(i) = xMedido(i-1) + velocidadeEncoder(i)*cos(theta(i-1))*dt;
    yMedido(i) = yMedido(i-1) + velocidadeEncoder(i)*sin(theta(i-1))*dt;
    theta(i)   = theta(i-1) + angularEncoder(i)*dt;
end

%% var_erro_x = 2.3693
%% var_erro_y = 5.9983e-03


time = [0:0.05:0.05*(n-1)];

figure
plot(time,xMedido,'r','LineWidth', 2)
hold on
plot(time,x_real,'b','LineWidth', 2)
grid on
xlabel('Tempo (s)')
ylabel('Posição X (m)')
legend('Medido','Real')

figure
plot(time,yMedido,'r','LineWidth', 2)
hold on
plot(time,y_real,'b','LineWidth', 2)
grid on
xlabel('Tempo (s)')
ylabel('Posição Y (m)')
legend('Medido','Real')

figure
plot3(xMedido, yMedido, zeros(size(xMedido)), 'r','LineWidth', 2)
hold on
plot3(x_real, y_real, zeros(size(x_real)), 'b','LineWidth', 2)
grid on
xlabel('Y')
ylabel('X')
zlabel('Z')
title('Caminho realizado pelo robô', 'FontSize', 14, 'FontWeight', 'bold')
view(3)           % Visualização 3D
zlim([0 0.01])    % "Chão" bem fino, quase plano
legend('Odometria', 'Original', 'FontSize', 14, 'Location', 'southeast')





erro_x = xMedido' - x_real;
media_erro_x = mean(erro_x);
desvio_erro_x = std(erro_x);
var_erro_x = var(erro_x)


erro_y = yMedido' - y_real;
media_erro_y = mean(erro_y);
desvio_erro_y = std(erro_y);
var_erro_y = var(erro_y)


