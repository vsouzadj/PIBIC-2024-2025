clc
clear all
close all

%%nome = 'data2.txt';
%%%
%%
%%[fid, msg] = fopen(nome, 'r');
%%if fid < 0
%%    error('Não foi possível abrir "%s": %s', nome, msg);
%%end
%%
%%% Lê cabeçalho
%%headerLine = fgetl(fid);
%%header = strsplit(headerLine, ',');  % cell array 1×11
%%
%%% Pré-aloca M (célula) e extras (célula variável)
%%M = cell(1, numel(header));
%%M(1,:) = header;
%%extras = {};  % só os valores numéricos das colunas 5–11
%%
%%% Processa dados
%%row = 2;
%%while true
%%    linha = fgetl(fid);
%%    if ~ischar(linha), break; end
%%    campos = strsplit(linha, ',');
%%    if numel(campos) ~= numel(header)
%%        warning('Linha %d ignorada (esperava %d campos):\n%s', row, numel(header), linha);
%%        continue;
%%    end
%%
%%    linhaCell = cell(1, numel(header));
%%    extraRow = {};  % vai armazenar só os numéricos de 5–11
%%    for k = 1:numel(header)
%%        txt = campos{k};
%%        if k <= 4
%%            % sempre decimal
%%            val = str2double(txt);
%%            if isnan(val)
%%                error('Erro na linha %d, coluna %d: esperava número, encontrou "%s".', row, k, txt);
%%            end
%%            linhaCell{k} = val;
%%        else
%%            % coluna 5–11: pode ser 'X' ou número
%%            if ~strcmp(txt, 'X')
%%                val = str2double(txt);
%%                if isnan(val)
%%                    warning('Linha %d, coluna %d: valor não reconhecido e será ignorado: "%s".', row, k, txt);
%%                    linhaCell{k} = [];
%%                else
%%                    linhaCell{k} = val;
%%                    extraRow{end+1} = val;  %#ok<SAGROW>
%%                end
%%            else
%%                linhaCell{k} = [];
%%            end
%%        end
%%    end
%%
%%    M(row, :) = linhaCell;
%%    if ~isempty(extraRow)
%%        extras(end+1,1:length(extraRow)) = extraRow;  %#ok<SAGROW>
%%    end
%%    row = row + 1;
%%end
%%fclose(fid);
%%
%%% Converte as 4 colunas iniciais para matriz numérica
%%dados4 = cell2mat( M(2:end, 1:4) );
%%time = (0:0.01:0.01*(size(dados4,1)-1))';
%%
%%% Exemplos de extração das duas primeiras
%%data1 = dados4(:,3);
%%data2 = dados4(:,4);

function Rz = rotZ(deg)
  rad = deg2rad(deg);
  Rz = [cos(rad), -sin(rad), 0;
        sin(rad),  cos(rad), 0;
        0,         0,        1];
end

function Ry = rotY(deg)
  rad = deg2rad(deg);
  Ry = [cos(rad),  0, sin(rad);
        0,         1, 0;
       -sin(rad),  0, cos(rad)];
end

function Rx = rotX(deg)
  rad = deg2rad(deg);
  Rx = [1, 0,         0;
        0, cos(rad), -sin(rad);
        0, sin(rad),  cos(rad)];
end


% abre o arquivo
fid = fopen('odo_cam_armazen_6.txt','r');
% fid = fopen('odo_cam_armazen_5.txt','r'); esse cheguei a colocar no pibic
if fid < 0
  error('Não foi possível abrir cam_odo_armazen.txt');
end

% lê as 12 colunas como números (%f), tratando 'X' como NaN
C = textscan(fid, ...
    repmat('%f',1,13), ...          % '%f%f…%f' 12 vezes
    'Delimiter', ',', ...           % separador vírgula
    'CollectOutput', true, ...      % junta tudo em C{1}
    'TreatAsEmpty', {'X'}, ...      % todo "X" vira vazio
    'EmptyValue', NaN);             % valor vazio = NaN

fclose(fid);

% separa a matriz [n×12] em vetores “plano”
dados   = C{1};
x_real       = dados(:,1);
y_real       = dados(:,2);
yawRobo = dados(:,3);
lv      = dados(:,4);
rv      = dados(:,5);
ID      = dados(:,6);
hamming = dados(:,7);
tx      = dados(:,8);
ty      = dados(:,9);
tz      = dados(:,10);
yaw     = dados(:,11);
pitch   = dados(:,12);
roll    = dados(:,13);


% Processamento de ruído e plot como antes
n = length(x_real);
sigma_left  = 0;  % 1.3
sigma_right = 0; % 0.3;  % 1.3

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
%%a1  = 1/6;
%%a2  = 2/6;
%%a3  = 2/6;
%%a4  = 1/6;
%%q21 = 1/2;
%%q31 = 0;
%%q32 = 1/2;
%%q41 = 0;
%%q42 = 0;
%%q43 = 1;
dt  = 0.05;

N = length(velocidadeEncoder);

for i = 2:N
%%    % estados iniciais neste passo
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

%%
%%theta = theta * (180/pi);
%theta = atan2(sin(theta), cos(theta));

theta = mod(theta + pi, 2*pi) - pi;
%%theta(1:193) = -1 * theta(1:193);

% 4.6 rad  263

% -1.55

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
plot(xMedido,yMedido,'r')
hold on
plot(x_real,y_real,'b')


erro_x = xMedido' - x_real;
media_erro_x = mean(erro_x);
desvio_erro_x = std(erro_x);
var_erro_x = var(erro_x)


erro_y = yMedido' - y_real;
media_erro_y = mean(erro_y);
desvio_erro_y = std(erro_y);
var_erro_y = var(erro_y)


% tag 0 x 2m  y 4.95m
% tag 9 x 4m  y 4.95m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKF


% 0.030919   0.1062   % RK4
% 0.030948   0.1063   % euler

xVarEnc =  var(erro_x);
yVarEnc =  var(erro_y);
if(xVarEnc > yVarEnc)
  thethaVarEnc = xVarEnc + 0.1*xVarEnc;
else
  thethaVarEnc = yVarEnc + 0.1*yVarEnc;
end


% X: th = [0.093554; 0.0021035; -0.093582];  modelo = 0.093554*exp(0.0021035*x) - 0.093582;

% Z: th = [8.6174e-05; 1.9101];  modelo = (8.6174e-05*exp(1.9101*x));

%%xVarCam = 0.001;
%%yVarCam = 0.002;
%%thethaVarCam = 0.003;

Q = [xVarEnc 0 0; 0 yVarEnc 0; 0 0 thethaVarEnc];   % Variância do processo: encoder
%%R = [xVarCam 0 0; 0 yVarCam 0; 0 0 thethaVarCam];   % Variância da medição: câmera

x_est     = zeros(3, n);      % estimativa do estado
x_kalman  = zeros(3, n);      % saída filtrada

% Inicializar estimativas
x_est(:,1) = [x_real(1); y_real(1); pi];
% x_est(:,1) = [0.751210; 0.761793; 0];  % no mapa grande
P_est = Q;

n = length(xMedido);
x_kalman(:,1) = x_est(:,1);

tamos = 0.05;

velEncoder = velocidadeEncoder;
velAng = angularEncoder;

for k = 2:n

   if (velAng(k) == 0)
     velAng(k) = 4.3220e-06 * (randi([0,1])*2 - 1);  % gera -1 ou +1 aleatório
   end

   x_model = x_est(1,k-1) - ((velEncoder(k)/velAng(k))*sin(x_est(3,k-1))) + ((velEncoder(k)/velAng(k))*sin(x_est(3,k-1)+(velAng(k)*tamos)));
   y_model = x_est(2,k-1) + ((velEncoder(k)/velAng(k))*cos(x_est(3,k-1))) - ((velEncoder(k)/velAng(k))*cos(x_est(3,k-1)+(velAng(k)*tamos)));
   theta_model = x_est(3,k-1) + velAng(k)*tamos;

   x_pred = [x_model; y_model; theta_model];

   g1 = [1 0 ( ((velEncoder(k)/velAng(k))*-cos(x_est(3,k-1))) + ((velEncoder(k)/velAng(k))*cos(x_est(3,k-1)+(velAng(k)*tamos))))];
   g2 = [0 1 ( ((velEncoder(k)/velAng(k))*-sin(x_est(3,k-1))) + ((velEncoder(k)/velAng(k))*sin(x_est(3,k-1)+(velAng(k)*tamos))))];
   g3 = [0 0 1];

   G = [g1; g2; g3];

   P_pred = G * P_est * G' + Q;

   % x,y,yawRobo,lv,rv,ID,hamming,tx,ty,tz,yaw,pitch,roll
   if ~any(isnan(dados(k,:))) && abs(velAng(k)) < 0.1*max(velAng(:)) % medições disponíveis

       h_ = x_pred;

       H = eye(3);

       if(ID(k) == 0 || ID(k) == 5)

         % --- Ângulos de Euler para camera -> mundo (exemplo) ---
         yaw_camera_mundo   = 0;   % Rotação em Z
         pitch_camera_mundo = 0;     % Rotação em Y
         roll_camera_mundo  = -90;  % Rotação em X

         yaw_tag_cam   = yaw(k);
         pitch_tag_cam = pitch(k);
         roll_tag_cam  = roll(k);

         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);

         R_mapa_camera = R_tag_mundo * R_tag_cam;

         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa

         tx(k) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tx(k))) - 0.093582))*randn + tx(k);
         tz(k) = 0.01*sqrt((8.6174e-05*exp(1.9101*tz(k))))*randn + tz(k);

         aux = 7.8073e-07*exp(1.5452*yaw_camera_mapa_rad);
         yaw_camera_mapa_rad = aux*randn + yaw_camera_mapa_rad;

         yaw_camera_mapa_deg = (7.8073e-07*exp(1.5452*yaw_camera_mapa_deg))*randn + yaw_camera_mapa_deg;

         R = [10^-4*(0.093554*exp(0.0021035*abs(tx(k))) - 0.093582) 0 0; 0 10^-4*(8.6174e-05*exp(1.9101*tz(k))) 0; 0 0 aux];

         K = P_pred * H'*pinv(H * P_pred * H' + R);

         if(ID(k) == 0)
           x_est(:, k) = x_pred + K * ([2.95+tx(k); -3.625+tz(k); yaw_camera_mapa_rad] - h_);  % Atualizar x_est com a nova estimativa
         else
          x_est(:, k) = x_pred;
          P_est = P_pred;

%%         x_est(:, k) = x_pred + K * ([2.6+tx(k); -10.175+tz(k); yaw_camera_mapa_rad] - h_);  % Atualizar x_est com a nova estimativa
         end

         P_est = (eye(size(K * H)) - K * H) * P_pred;
       end

        % x,y,lv,rv,ID,hamming,tx,ty,tz,yaw,pitch,roll
       if(ID(k) == 2 || ID(k) == 7 || ID(k) == 8)

         % --- Ângulos de Euler para camera -> mundo (exemplo) ---
         yaw_camera_mundo   = 0;   % Rotação em Z
         pitch_camera_mundo = 0;     % Rotação em Y
         roll_camera_mundo  = 90;  % Rotação em X

         yaw_tag_cam   = yaw(k);
         pitch_tag_cam = pitch(k);
         roll_tag_cam  = roll(k);

         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);

         R_mapa_camera = R_tag_mundo * R_tag_cam;

         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa

         tx(k) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tx(k))) - 0.093582))*randn + tx(k);
         tz(k) = 0.01*sqrt((8.6174e-05*exp(1.9101*tz(k))))*randn + tz(k);

         aux = 7.8073e-07*exp(1.5452*yaw_camera_mapa_rad);
         yaw_camera_mapa_rad = aux*randn + yaw_camera_mapa_rad;

         R = [10^-4*(0.093554*exp(0.0021035*abs(tx(k))) - 0.093582) 0 0; 0 10^-4*(8.6174e-05*exp(1.9101*tz(k))) 0; 0 0 aux];

         K = P_pred * H'*pinv(H * P_pred * H' + R);

         if(ID(k) == 2)
           x_est(:, k) = x_pred + K * ([-1.725-tx(k); -1.875-tz(k); yaw_camera_mapa_rad] - h_);  % Atualizar x_est com a nova estimativa

%%           x_est(:, k) = x_pred;
%%           P_est = P_pred;
         elseif(ID(k) == 8)
      %     x_est(:, k) = x_pred + K * ([-3.725-tx(k); 10.3-tz(k); yaw_camera_mapa_rad] - h_);

           x_est(:, k) = x_pred;
           P_est = P_pred;
         else
 %          x_est(:, k) = x_pred + K * ([-2.36131-tx(k); 4.27557-tz(k); yaw_camera_mapa_rad] - h_);

            x_est(:, k) = x_pred;
           P_est = P_pred;
         end

         P_est = (eye(size(K * H)) - K * H) * P_pred;
       end




       if(ID(k) == 3)

%%        % --- Ângulos de Euler para camera -> mundo (exemplo) ---
         yaw_camera_mundo   = 0;   % Rotação em Z
         pitch_camera_mundo = 0;     % Rotação em Y
         roll_camera_mundo  = 90;  % Rotação em X

         yaw_tag_cam   = yaw(k);
         pitch_tag_cam = pitch(k);
         roll_tag_cam  = roll(k);

         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);

         R_mapa_camera = R_tag_mundo * R_tag_cam;

         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa

         tz(k) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tz(k))) - 0.093582))*randn + tz(k);
         tx(k) = 0.01*sqrt((8.6174e-05*exp(1.9101*tx(k))))*randn + tx(k);

         aux = 7.8073e-07*exp(1.5452*yaw_camera_mapa_rad);
         yaw_camera_mapa_rad = aux*randn + yaw_camera_mapa_rad;

         R = [10^-4*(0.093554*exp(0.0021035*abs(tz(k))) - 0.093582) 0 0; 0 10^-4*(8.6174e-05*exp(1.9101*tx(k))) 0; 0 0 aux];

         K = P_pred * H'*pinv(H * P_pred * H' + R);

         x_est(:, k) = x_pred + K * ([-1.95+tz(k); -8.45-tx(k); yaw_camera_mapa_rad] - h_);  % Atualizar x_est com a nova estimativa

         P_est = (eye(size(K * H)) - K * H) * P_pred;

%%     x_est(:, k) = x_pred;
%%    P_est = P_pred;
       end

       if(ID(k) == 4)

         % --- Ângulos de Euler para camera -> mundo (exemplo) ---
%%         yaw_camera_mundo   = 90;   % Rotação em Z
%%         pitch_camera_mundo = 0;     % Rotação em Y
%%         roll_camera_mundo  = 90;  % Rotação em X
%%
%%         yaw_tag_cam   = yaw(k);
%%         pitch_tag_cam = pitch(k);
%%         roll_tag_cam  = roll(k);
%%
%%         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
%%         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);
%%
%%         R_mapa_camera = R_tag_mundo * R_tag_cam;
%%
%%         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
%%         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa
%%
%%         tz(k) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tz(k))) - 0.093582))*randn + tz(k);
%%         tx(k) = 0.01*sqrt((8.6174e-05*exp(1.9101*tx(k))))*randn + tx(k);
%%
%%         aux = 7.8073e-07*exp(1.5452*yaw_camera_mapa_rad);
%%         yaw_camera_mapa_rad = aux*randn + yaw_camera_mapa_rad;
%%
%%         R = [10^-4*(0.093554*exp(0.0021035*abs(tz(k))) - 0.093582) 0 0; 0 10^-4*(8.6174e-05*exp(1.9101*tx(k))) 0; 0 0 aux];
%%
%%         K = P_pred * H'*pinv(H * P_pred * H' + R);
%%
%%         x_est(:, k) = x_pred + K * ([3.4+tz(k); 3.025-tx(k); yaw_camera_mapa_rad] - h_);  % Atualizar x_est com a nova estimativa
%%
%%         P_est = (eye(size(K * H)) - K * H) * P_pred;

          x_est(:, k) = x_pred;
          P_est = P_pred;
       end

   else
     x_est(:, k) = x_pred;
     P_est = P_pred;
   end

   x_kalman(:, k) = x_est(:, k);
end


tempo = time;

figure
subplot(3,1,1)
plot(tempo,x_kalman(1,:), 'ob')
hold on
plot(tempo,xMedido, 'g')
plot(tempo,x_real, 'k')
grid on
legend('EKF','Encoder','Original')
title("EKF para a posição X")
xlabel('Tempo (s)')
ylabel('Posição (m)')

%%figure
subplot(3,1,2)
plot(tempo,x_kalman(2,:), 'ob')
hold on
plot(tempo,yMedido, 'g')
plot(tempo,y_real, 'k')
grid on
legend('EKF', 'Encoder', 'Original')
title("EKF para a posição Y")
xlabel('Tempo (s)')
ylabel('Posição (m)')


%%figure
subplot(3,1,3)
plot(tempo,x_est(3,:), 'g')
hold on
plot(tempo,theta,'k')
plot(tempo,(yawRobo*pi)/180,'m')
grid on
legend('EKF','Encoder','Original')
title("EKF para a orientação")
xlabel('Tempo (s)')
ylabel('Posição (m)')



figure
hold on
plot(x_kalman(1,:),x_kalman(2,:),'og','LineWidth', 2)
plot(xMedido,yMedido,'b','LineWidth', 2)
plot(x_real,y_real,'r','LineWidth', 2)
grid on
%%legend('EKF','Encoder', 'Original', 'location','southeast')
legend('EKF','Encoder', 'Original')
title("Caminho do robô")
xlabel('Posição X (m)')
ylabel('Posição Y (m)')

figure
plot3(x_real, y_real, zeros(size(x_real)), 'r','LineWidth', 2)
hold on
plot3(x_kalman(1,:), x_kalman(2,:), zeros(size(x_real)), '-g','LineWidth', 2)
plot3(xMedido, yMedido, zeros(size(xMedido)), 'b','LineWidth', 2)
grid on
xlabel('Y')
ylabel('X')
zlabel('Z')
title('Caminho realizado pelo robô', 'FontSize', 14, 'FontWeight', 'bold')
view(3)           % Visualização 3D
zlim([0 0.01])    % "Chão" bem fino, quase plano
legend('Original', 'EKF', 'Odometria','FontSize', 14, 'Location', 'southeast')

figure
plot(x_real,y_real,'r','LineWidth', 2)
grid on
hold on
plot(x_kalman(1,:),x_kalman(2,:),'ob','LineWidth', 2)
legend('Original','EKF')
title("Caminho do robô")
xlabel('Posição X (m)')
ylabel('Posição Y (m)')

erro_x = x_real' - x_kalman(1,:);
rmse_x = sqrt(mean(erro_x.^2))   % 0.059281 0.059389 0.059349 0.059209 0.059260 0.059277 0.059444 0.059330 0.059329 0.059331

erro_y = y_real' - x_kalman(2,:);
rmse_y = sqrt(mean(erro_y.^2))   % 0.088378 0.088614 0.088353 0.088708 0.088327 0.088604 0.088369 0.088760  0.088688 0.088509

erro_theta = ((pi*yawRobo')/180) - x_kalman(3,:);
erro_theta = mod(erro_theta+ pi, 2*pi) - pi;
rmse_rad = sqrt(mean(erro_theta.^2));
rmse_deg = rmse_rad * 180/pi     % 6.3580 6.3580 6.3583  6.35725 6.3577 6.3579 6.3586 6.3586 6.3579 6.3584

vetor = [6.3580 6.3580 6.3583  6.35725 6.3577 6.3579 6.3586 6.3586 6.3579 6.3584];


