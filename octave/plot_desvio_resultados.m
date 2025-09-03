clc
clear all
close all

distancia = [60; 120; 180; 240];

% DESVIO PADRÃO DO ERRO Z

z_april = [6.9326e-03; 0.021595; 0.055000; 0.092544];

z_aruco = [0.049550; 0.3355; 0.6170; 3.2564];

figure
hold on
plot(distancia, z_april,'*-r' ,'LineWidth', 2)
plot(distancia, z_aruco, '*-b' ,'LineWidth', 2)
xlabel('Distância da câmera para a tag (cm)','FontSize', 14)
ylabel('Desvio padrão (cm)','FontSize', 14)
title('Desvio padrão do erro em Z', 'FontSize', 14, 'FontWeight', 'bold')
legend('AprilTag', 'ArUco', 'FontSize', 14, 'Location', 'east');


% DESVIO PADRÃO DO ERRO X

x_april = [5.4254e-03; 0.017303; 0.018249; 0.019899];

x_aruco = [0.016276; 0.048440; 0.045932; 0.078909];

figure
hold on
plot(distancia, x_april,'*-r' ,'LineWidth', 2)
plot(distancia, x_aruco, '*-b' ,'LineWidth', 2)
xlabel('Distância da câmera para a tag (cm)','FontSize', 14)
ylabel('Desvio padrão (cm)','FontSize', 14)
title('Desvio padrão do erro em X', 'FontSize', 14, 'FontWeight', 'bold')
legend('AprilTag', 'ArUco', 'FontSize', 14, 'Location', 'east');



% DESVIO PADRÃO DO ERRO YAW

yaw_april = [7.6165e-05; 1.1799e-03; 1.2402e-03; 1.2146e-03];

yaw_aruco = [0.057157; 0.20708; 98.937; 0.9791];

figure
hold on
plot(distancia, yaw_april,'*-r' ,'LineWidth', 2)
plot(distancia, yaw_aruco, '*-b' ,'LineWidth', 2)
xlabel('Distância da câmera para a tag (cm)','FontSize', 14)
ylabel('Desvio padrão (graus)','FontSize', 14)
title('Desvio padrão do erro em YAW', 'FontSize', 14, 'FontWeight', 'bold')
legend('AprilTag', 'ArUco', 'FontSize', 14, 'Location', 'east');




% DESVIO PADRÃO DO ERRO PITCH

pitch_april = [2.1663e-03; 0.068553; 0.064675; 0.095021];

pitch_aruco = [42.667; 114.84; 172.48; 149.59];

figure
hold on
plot(distancia, pitch_april,'*-r' ,'LineWidth', 2)
plot(distancia, pitch_aruco, '*-b' ,'LineWidth', 2)
xlabel('Distância da câmera para a tag (cm)','FontSize', 14)
ylabel('Desvio padrão (graus)','FontSize', 14)
title('Desvio padrão do erro em PITCH', 'FontSize', 14, 'FontWeight', 'bold')
legend('AprilTag', 'ArUco', 'FontSize', 14, 'Location', 'east');



% DESVIO PADRÃO DO ERRO ROLL

roll_april = [6.5171e-04; 0.010165; 0.018407; 8.8978e-03];

roll_aruco = [2.0126; 3.0994; 8.1400; 16.383];

figure
hold on
plot(distancia, roll_april,'*-r' ,'LineWidth', 2)
plot(distancia, roll_aruco, '*-b' ,'LineWidth', 2)
xlabel('Distância da câmera para a tag (cm)','FontSize', 14)
ylabel('Desvio padrão (graus)','FontSize', 14)
title('Desvio padrão do erro em ROLL', 'FontSize', 14, 'FontWeight', 'bold')
legend('AprilTag', 'ArUco', 'FontSize', 14, 'Location', 'east');



