clear; clc;
%Programa Pulso Simples - (Maio/2020) revisão

%------------------------------------------------------------------------%
% Este código calcula o comportamento da Magnetização quando aplicamos um
% pulso CP Spin Echo;        90y_t_180x_t_eco
%
% Moraes, T.B.; Colnago, L.A.; 
% Simulação de sinais de RMN através das equações de Bloch, 
% v. 37(8), p. 1410-1416, Química Nova, 2014.
% https://doi.org/10.5935/0100-4042.20140210 
%------------------------------------------------------------------------%

% === Parâmetros entrada ===
% Amostra  
T1 = 600;   % ms
T2 = 400;   % ms

% Espectrômetro
theta1 = 90;       % Ângulo 1° Pulso 90°  ( em graus) 
theta2 = 180;      % Ângulo 2° Pulso 180°
T      = 1000;     % (ms) Tempo total aquisição FID
Tn1    = 150;      % (ms) Tempo até primeiro pulso

% convertendo ângulo para radianos;
theta1 = theta1*pi/180;   
theta2 = theta2*pi/180;   

% Frequencia em [1/ms] % tudo está em ms, logo 25*10^-3 (1/ms) = 25 Hz !
x0 = 25 * 10^-3 ;      

% T2* = 1/((1/T2) + (1/T2inom));  
% Inomogeneidade de campo: T2inom = gamaH * delta(B0)
% T2* => T22;  decaimento natural do FID = exp(-tempo/T2*);
T22 = 20; %ms 

% Largura FWHM da distribuição de Isocromatas (Hz)
FWHM = 1/(pi*T22);

% Região de integração do pico (Isocromatas)
% 39*FWHM integra 99,5% da região isocromatas
% 76*FWHM integra 99,9% da região isocromatas
int = 39;
x01 = x0 - int*FWHM/2; 
x02 = x0 + int*FWHM/2;

% Offset da frequencia Omega = (2pi*gamma)*2pi*(fref - f0)
% df = Omega / 2pi         
% Niso: número de isocromatas de Frequências
% Caso surja Echos "falsos" na simulação -> aumentar Niso
Niso = 1000;
df = linspace(x01,x02,Niso);       

% Lorentziana Isocromatas - Campo Inomogêneo
fL = (T22)./( 1 + ((df - x0).^2).*(2*pi*T22).^2 );

% (ms)  tempo entre pontos
dT = 1;             

% ===== Matrizes Rotação =====
Rtheta1 = [cos(theta1) 0 sin(theta1);0 1 0 ; -sin(theta1) 0 cos(theta1)];  % Fase y
Rtheta2 = [1 0 0; 0 cos(theta2) -sin(theta2); 0 sin(theta2) cos(theta2)];  % Fase x

E1 = exp(-dT/T1);  E2 = exp(-dT/T2); 
E = [E2 0 0; 0 E2 0; 0 0 E1];       B = [0 0 1-E1].';

% ===== Contadores Sinal =====
N0 = round(T/dT);        % nº pontos total;    tempo total = N0*dT
N1 = round(Tn1/dT);      % tempo até primeiro pulso
N2 = N0 - N1;            % nº pontos depois do ultimo pulso

M  = zeros(3,N0);         % Cria vetor magnetização tamanho (3 x N0)
Ms = zeros(3,N0);         % Cria vetor magnetização Ms (armazenar sinal)
M0 = [0;0;1];             % Magnetização posição inicial no eixo z

for f=1:length(df)   
  % Phi = 2pi*(Fref - f0)*dT; onde df = (Fref - f0)
  phi = 2*pi*(df(f))*dT;     
  Rphi = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
       
  M(:,1)= Rtheta1*M0 + B;             % 1° Pulso 90°y
  for k=2:(N1+1)                     % evolução sinal
	  M(:,k) = E*Rphi*M(:,k-1)+B;
  end;
  
  M(:,N1+2)= Rtheta2*M(:,N1+1)+B;           % 2° Pulso 180ºx 
  for k=2:(N2-1)                            % evolução sinal
    M(:,k+N1+1) = E*Rphi*M(:,k+N1)+B;
  end;
  
  % ------- Calculo peso Inomogeneidade - T2* ----------
  % Pico Lorentziano
  g = (T22)./((1+((df(f)-x0).^2).*(2*pi*T22).^2));
   
  % Somando as componentes x y z das isocromata
  Ms = [g*M(1,:)+Ms(1,:); g*M(2,:)+Ms(2,:); g*M(3,:)+Ms(3,:)];
end;

% Normaliza sinal FID
Ms = Ms/max(Ms(:,1));

% Curvas de Referência
tempo = [0:N0-1]*dT;  
A = sqrt(Ms(1,1)^2+Ms(2,1)^2);
CurvaT2 = A*exp(-tempo/T2);  
CurvaT22 = A*exp(-tempo/T22);

% Módulo do Sinal real e imaginário 
 Mod = (Ms(1,:).^2+Ms(2,:).^2).^0.5;

% ===== Graficando Resultados ======
figure(1)
p = plot(df,fL,'k'); 
title('Distribuicao Lorentziana'); 
legend(strcat('T2*=',num2str(T22),'ms'))
xlabel('Frequência (x10^3 Hz)'); ylabel('Intensidade'); grid on;
set(p,'LineWidth',1.5);

figure(2)
p=plot(tempo,Ms(1,:),'b-',tempo,Ms(2,:),'k-',tempo,Ms(3,:),'r--',tempo, CurvaT2,'g:',tempo, Mod,'y-');
legend('M_x','M_y','M_z','T_2','Mod'); title('SpinEcho: 180°y - 90°x');
xlabel('Tempo (ms)'); ylabel('Intensidade'); 
grid on;  set(p,'LineWidth',1.5); 

% Comando exportar imagem 
%print -dpng PulsoSimples.png -r300   

%  Comando exportar dados .dat
%dlmwrite('Teste.dat', [tempo.'  Ms(1,:).'  Ms(2,:).'  Ms(3,:).'], '\t')  



