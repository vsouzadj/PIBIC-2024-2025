clc
clear all
close all


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


% Abre o arquivo para leitura linha a linha
fid = fopen('odo_cam_armazen_6.txt','r');
if fid < 0
  error('Não foi possível abrir odo_cam_armazen_2.txt');
end

% Inicializa os vetores de dados
x_real   = [];
y_real   = [];
yawRobo  = [];
lv       = [];
rv       = [];
ID       = [];
hamming  = [];
tx       = [];
ty       = [];
tz       = [];
yaw      = [];
pitch    = [];
roll     = [];

% Loop de leitura linha a linha
while ~feof(fid)
  linha = strtrim(fgetl(fid));

  % Se for linha "tamos", insere zeros em todos os campos
  if strcmpi(linha, 'tamos')
    x_real(end+1,1)   = 0;
    y_real(end+1,1)   = 0;
    yawRobo(end+1,1)  = 0;
    lv(end+1,1)       = 0;
    rv(end+1,1)       = 0;
    ID(end+1,1)       = 0;
    hamming(end+1,1)  = 0;
    tx(end+1,1)       = 0;
    ty(end+1,1)       = 0;
    tz(end+1,1)       = 0;
    yaw(end+1,1)      = 0;
    pitch(end+1,1)    = 0;
    roll(end+1,1)     = 0;
    continue;
  end

  % Divide a linha em 13 campos separados por vírgula
  partes = strsplit(linha, ',');

  % Valida número de colunas
  if numel(partes) ~= 13
    warning('Linha ignorada, colunas esperadas: 13, encontradas: %d', numel(partes));
    continue;
  end

  % Converte valores, tratando 'X' como NaN
  valores = nan(1,13);
  for i = 1:13
    if ~strcmpi(partes{i}, 'X')
      valores(i) = str2double(partes{i});
    end
  end

  % Armazena em cada vetor
  x_real(end+1,1)   = valores(1);
  y_real(end+1,1)   = valores(2);
  yawRobo(end+1,1)  = valores(3);
  lv(end+1,1)       = valores(4);
  rv(end+1,1)       = valores(5);
  ID(end+1,1)       = valores(6);
  hamming(end+1,1)  = valores(7);
  tx(end+1,1)       = valores(8);
  ty(end+1,1)       = valores(9);
  tz(end+1,1)       = valores(10);
  yaw(end+1,1)      = valores(11);
  pitch(end+1,1)    = valores(12);
  roll(end+1,1)     = valores(13);
end

fclose(fid);



% Processamento de ruído e plot como antes
%%n = length(x_real);
%%sigma_left  = 0;
%%sigma_right = 0;
%%
%%leftEncoder  = lv;
%%rightEncoder = rv;
%%
%%ruido_left  = sigma_left  * randn(n,1);
%%ruido_right = sigma_right * randn(n,1);
%%
%%leftEncoder  = ((leftEncoder  + ruido_left)  * 0.0975);
%%rightEncoder = ((rightEncoder + ruido_right) * 0.0975);
%%
%%velocidadeEncoder = (leftEncoder+rightEncoder)/2;
%%angularEncoder = (-leftEncoder+rightEncoder)/(0.381);  % 0.381

%%figure
%%hold on
%%plot(time, leftEncoder,  'r', 'LineWidth', 2)
%%plot(time, rightEncoder, 'g', 'LineWidth', 2)
%%grid on
%%xlabel('Tempo (s)')
%%ylabel('Velocidade (m/s)')
%%legend('Left','Right')


%%dt = 0.05;
%%N = length(velocidadeEncoder);
%%
%%xMedido = zeros(N,1);
%%yMedido = zeros(N,1);
%%theta   = zeros(N,1);


% posição inicial
%%xMedido(1) = x_real(1);
%%yMedido(1) = y_real(1);
%%theta(1)   = pi;
%%
%%%% Parâmetros do método RK4
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
%%dt  = 0.05;
%%
%%N = length(velocidadeEncoder);
%%
%%for i = 2:N
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

%%    xMedido(i) = xMedido(i-1) + velocidadeEncoder(i)*cos(theta(i-1))*dt;
%%    yMedido(i) = yMedido(i-1) + velocidadeEncoder(i)*sin(theta(i-1))*dt;
%%    theta(i)   = theta(i-1) + angularEncoder(i)*dt;
%%end

%% var_erro_x = 2.3693
%% var_erro_y = 5.9983e-03

%%
%%time = [0:0.05:0.05*(n-1)];
%%
%%figure
%%plot(time,xMedido,'r','LineWidth', 2)
%%hold on
%%plot(time,x_real,'b','LineWidth', 2)
%%grid on
%%xlabel('Tempo (s)')
%%ylabel('Posição X (m)')
%%legend('Medido','Real')
%%
%%figure
%%plot(time,yMedido,'r','LineWidth', 2)
%%hold on
%%plot(time,y_real,'b','LineWidth', 2)
%%grid on
%%xlabel('Tempo (s)')
%%ylabel('Posição Y (m)')
%%legend('Medido','Real')
%%
%%figure
%%plot3(xMedido, yMedido, zeros(size(xMedido)), 'r','LineWidth', 2)
%%hold on
%%plot3(x_real, y_real, zeros(size(x_real)), 'b','LineWidth', 2)
%%grid on
%%xlabel('Y')
%%ylabel('X')
%%zlabel('Z')
%%title('Caminho realizado pelo robô', 'FontSize', 14, 'FontWeight', 'bold')
%%view(3)           % Visualização 3D
%%zlim([0 0.01])    % "Chão" bem fino, quase plano
%%legend('Odometria', 'Original', 'FontSize', 14, 'Location', 'southeast')
%%
%%
%%
%%
%%
%%erro_x = xMedido' - x_real;
%%media_erro_x = mean(erro_x);
%%desvio_erro_x = std(erro_x);
%%var_erro_x = var(erro_x);
%%
%%
%%erro_y = yMedido' - y_real;
%%media_erro_y = mean(erro_y);
%%desvio_erro_y = std(erro_y);
%%var_erro_y = var(erro_y);


x_cam = [];
y_cam = [];
theta_cam = [];
%%figure
%%plot(x_real,y_real,'r','LineWidth', 2)
%%hold on

for m = 1:length(tx)
  x_cam(m) = NaN;
  y_cam(m) =  NaN;
  theta_cam(m) = NaN;

  if ~any(isnan(tx(m,:)))  % medições disponíveis
    if(ID(m) == 0 || ID(m) == 5) % ja foi
         % --- Ângulos de Euler para camera -> mundo (exemplo) ---
         yaw_camera_mundo   = 0;   % Rotação em Z
         pitch_camera_mundo = 0;     % Rotação em Y
         roll_camera_mundo  = -90;  % Rotação em X

         yaw_tag_cam   = yaw(m);
         pitch_tag_cam = pitch(m);
         roll_tag_cam  = roll(m);

         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);

         R_mapa_camera = R_tag_mundo * R_tag_cam;

         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa

         tx(m) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tx(m))) - 0.093582))*randn + tx(m);
         tz(m) = 0.01*sqrt((8.6174e-05*exp(1.9101*tz(m))))*randn + tz(m);
         %yaw_camera_mapa_rad = sqrt(varAng)*randn + yaw_camera_mapa_rad;

         yaw_camera_mapa_deg = (7.8073e-07*exp(1.5452*yaw_camera_mapa_deg))*randn + yaw_camera_mapa_deg;

         if yaw_camera_mapa_deg < 0
          yaw_camera_mapa_deg = yaw_camera_mapa_deg + 360;
         end


         if(ID(m) == 0)
           x_cam(m) = 2.95+tx(m);
           y_cam(m) = -3.625+tz(m);
           theta_cam(m) = yaw_camera_mapa_deg;
         else
           x_cam(m) = 2.6+tx(m);
           y_cam(m) = -10.175+tz(m);
           theta_cam(m) = yaw_camera_mapa_deg;
         end
       end

        % x,y,yawRobo,lv,rv,ID,hamming,tx,ty,tz,yaw,pitch,roll
       if(ID(m) == 8 || ID(m) == 2 || ID(m) == 7) % ja foi
         % --- Ângulos de Euler para camera -> mundo (exemplo) ---
         yaw_camera_mundo   = 0;   % Rotação em Z
         pitch_camera_mundo = 0;     % Rotação em Y
         roll_camera_mundo  = 90;  % Rotação em X

         yaw_tag_cam   = yaw(m);
         pitch_tag_cam = pitch(m);
         roll_tag_cam  = roll(m);

         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);

         R_mapa_camera = R_tag_mundo * R_tag_cam;

         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa

         tx(m) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tx(m))) - 0.093582))*randn + tx(m);
         tz(m) = 0.01*sqrt((8.6174e-05*exp(1.9101*tz(m))))*randn + tz(m);

         if yaw_camera_mapa_deg < 0
          yaw_camera_mapa_deg = yaw_camera_mapa_deg + 360;
         end

         yaw_camera_mapa_deg = (7.8073e-07*exp(1.5452*yaw_camera_mapa_deg))*randn + yaw_camera_mapa_deg;

         if(ID(m) == 8)
           x_cam(m) = -3.725-tx(m);
           y_cam(m) = 10.3-tz(m);
           theta_cam(m) = yaw_camera_mapa_deg;
         elseif(ID(m) == 2)
           x_cam(m) = -1.725-tx(m);
           y_cam(m) = -1.875-tz(m);
           theta_cam(m) = yaw_camera_mapa_deg;
         else
           x_cam(m) = -2.36131-tx(m);
           y_cam(m) = 4.27557-tz(m);
           theta_cam(m) = yaw_camera_mapa_deg;
         end

       end

       if(ID(m) == 3) % ja foi
         % --- Ângulos de Euler para camera -> mundo (exemplo) ---
         yaw_camera_mundo   = 0;   % Rotação em Z
         pitch_camera_mundo = 0;     % Rotação em Y
         roll_camera_mundo  = 90;  % Rotação em X

         yaw_tag_cam   = yaw(m);
         pitch_tag_cam = pitch(m);
         roll_tag_cam  = roll(m);

         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);

         R_mapa_camera = R_tag_mundo * R_tag_cam;

         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa

         tz(m) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tz(m))) - 0.093582))*randn + tz(m);
         tx(m) = 0.01*sqrt((8.6174e-05*exp(1.9101*tx(m))))*randn + tx(m);
         yaw_camera_mapa_deg = (7.8073e-07*exp(1.5452*yaw_camera_mapa_deg))*randn + yaw_camera_mapa_deg;

         if yaw_camera_mapa_deg < 0
          yaw_camera_mapa_deg = yaw_camera_mapa_deg + 360;
         end

         x_cam(m) = -1.95+tz(m);
         y_cam(m) = -8.45-tx(m);
         theta_cam(m) = yaw_camera_mapa_deg;
       end

       if(ID(m) == 4)
         % --- Ângulos de Euler para camera -> mundo (exemplo) ---
         yaw_camera_mundo   = -90;   % Rotação em Z
         pitch_camera_mundo = 0;     % Rotação em Y
         roll_camera_mundo  = -90;  % Rotação em X

         yaw_tag_cam   = yaw(m);
         pitch_tag_cam = pitch(m);
         roll_tag_cam  = roll(m);

         R_tag_mundo = rotZ(yaw_camera_mundo) * rotY(pitch_camera_mundo) * rotX(roll_camera_mundo);
         R_tag_cam   = rotZ(yaw_tag_cam)   * rotY(pitch_tag_cam)   * rotX(roll_tag_cam);

         R_mapa_camera = R_tag_mundo * R_tag_cam;

         yaw_camera_mapa_rad = atan2( R_mapa_camera(2,3), R_mapa_camera(1,3) );
         yaw_camera_mapa_deg = rad2deg(yaw_camera_mapa_rad);  % orientação da camera/robô em relação ao mapa

         tz(m) = 0.01*sqrt((0.093554*exp(0.0021035*abs(tz(m))) - 0.093582))*randn + tz(m);
         tx(m) = 0.01*sqrt((8.6174e-05*exp(1.9101*tx(m))))*randn + tx(m);
         yaw_camera_mapa_deg = (7.8073e-07*exp(1.5452*yaw_camera_mapa_deg))*randn + yaw_camera_mapa_deg;

         if yaw_camera_mapa_deg < 0
          yaw_camera_mapa_deg = yaw_camera_mapa_deg + 360;
         end

         x_cam(m) = 3.4+tz(m);
         y_cam(m) = 3.0255-tx(m);
         theta_cam(m) = yaw_camera_mapa_deg;
       end

  end
end



%%figure
%%plot(x_real,y_real,'r','LineWidth', 2)
%%grid on
%%hold on
%%plot(x_cam,y_cam,'-b','LineWidth', 2)
%%legend('Original','Câmera')
%%title("Caminho do robô")
%%xlabel('Posição X (m)')
%%ylabel('Posição Y (m)')



%%figure
%%plot3(x_real, y_real, zeros(size(x_real)), 'r','LineWidth', 2)
%%hold on
%%plot3(x_cam, y_cam, zeros(size(y_cam)), 'm','LineWidth', 2)
%%grid on
%%xlabel('Y')
%%ylabel('X')
%%zlabel('Z')
%%title('Caminho realizado pelo robô', 'FontSize', 14, 'FontWeight', 'bold')
%%view(3)           % Visualização 3D
%%zlim([0 0.01])    % "Chão" bem fino, quase plano
%%legend('Original', 'AprilTag', 'FontSize', 14, 'Location', 'southeast')



%%figure
%%plot(yawRobo,'r')
%%hold on
%%plot(theta_cam,'b')
%%legend('Original yaw','Câmera yaw')

