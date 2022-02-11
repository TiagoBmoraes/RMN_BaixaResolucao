clear; clc;
%Programa Invers�o Recupera��o - (Maio/2020) revis�o

%------------------------------------------------------------------------%
% Este c�digo calcula o comportamento da Magnetiza��o quando aplicamos um
% de Invers�o Recupera��o      180x_t_90x_Fid
%
% Moraes, T.B.; Colnago, L.A.; 
% Simula��o de sinais de RMN atrav�s das equa��es de Bloch, 
% v. 37(8), p. 1410-1416, Qu�mica Nova, 2014.
% https://doi.org/10.5935/0100-4042.20140210 
%------------------------------------------------------------------------%

% === Par�metros entrada ===
% Amostra  
T1 = 600;   % ms
T2 = 400;   % ms

% Espectr�metro
theta1 =  180;      % Angulo 1� Pulso 180�  (em graus)
theta2 =   90;      % Angulo 2� Pulso 90�   
T      = 1500;      % (ms)  tempo total aquisi��o FID
Tn1    =  100;      % tempo at� Pulso 90�

% convertendo �ngulo para radianos;
theta1 = theta1*pi/180;
theta2 = theta2*pi/180;

% Frequencia em [1/ms] % tudo est� em ms, logo 25*10^-3 (1/ms) = 25 Hz !
x0 = 25 * 10^-3 ;      

% T2* = 1/((1/T2) + (1/T2inom));  
% Inomogeneidade de campo: T2inom = gamaH * delta(B0)
% T2* => T22;  decaimento natural do FID = exp(-tempo/T2*);
T22 = 20; %ms 

% Largura FWHM da distribui��o de Isocromatas (Hz)
FWHM = 1/(pi*T22);

% Regi�o de integra��o do pico (Isocromatas)
% 39*FWHM integra 99,5% da regi�o isocromatas
% 76*FWHM integra 99,9% da regi�o isocromatas
int = 39;
x01 = x0 - int*FWHM/2; 
x02 = x0 + int*FWHM/2;

% Offset da frequencia Omega = (2pi*gamma)*2pi*(fref - f0)
% df = Omega / 2pi         
% Niso: n�mero de isocromatas de Frequ�ncias
% Caso surja Echos "falsos" na simula��o -> aumentar Niso
Niso = 1000;
df = linspace(x01,x02,Niso);       

% Lorentziana Campo Inomog�nio
fL = (T22)./( 1 + ((df - x0).^2).*(2*pi*T22).^2 );

% (ms)  tempo entre pontos
dT = 1;             

% ===== Matrizes Rota��o =====
%Rtheta1 = [cos(theta1) 0 sin(theta1);0 1 0 ; -sin(theta1) 0 cos(theta1)];  % Fase y
Rtheta1 = [1 0 0; 0 cos(theta1) sin(theta1); 0 -sin(theta1) cos(theta1)];  % Fase x
Rtheta2 = [1 0 0; 0 cos(theta2) sin(theta2); 0 -sin(theta2) cos(theta2)];  % Fase x

E1 = exp(-dT/T1);  E2 = exp(-dT/T2); 
E = [E2 0 0; 0 E2 0; 0 0 E1];       B = [0 0 1-E1]'; 

% ===== Contadores Sinal =====
N0 = round(T/dT);        % n� pontos total;    tempo total = N0*dT
N1 = round(Tn1/dT);      % tempo at� primeiro pulso
N2 = N0 - N1;            % n� pontos depois do ultimo pulso

M = zeros(3,N0);                % Criar vetor magnetiza��o Msignal
Ms = zeros(3,N0);                % Criar vetor magnetiza��o Msignal
M0 = [0;0;1];             % Magnetiza��o posi��o inicial no eixo z

for f=1:length(df)
  phi = 2*pi*df(f)*dT;     
  Rphi = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    
  M(:,1)= Rtheta1*M0 + B;             % Pulso 180�x
  for k=2:(N1+1)                       % evolu��o sinal
	  M(:,k) = E*Rphi*M(:,k-1)+B;
  end;
      
  M(:,N1+2)= Rtheta2*M(:,N1+1)+B;           % pulso de 90�x 
  for k=2:N2-1                              % evolu��o sinal
    M(:,k+N1+1) = E*Rphi*M(:,k+N1)+B;
  end;
 
  % ------- Calculo peso Inomogeneidade - T2* ----------
  % Pico Lorentziano
  g = T22./((1+((df(f)-x0).^2)*(2*pi*T22)^2));
   
  % Somando as componentes x y z das isocromata
  Ms = [g*M(1,:)+Ms(1,:); g*M(2,:)+Ms(2,:); g*M(3,:)+Ms(3,:)];
end;

% Normaliza sinal FID
Ms = Ms/max(Ms(:,1));

% Curvas
tempo = [0:N0-1]*dT;  

% M�dulo Sinal
 Mod = (Ms(1,:).^2+Ms(2,:).^2).^0.5;

% ===== Graficando Resultados ======
figure(1)
p = plot(df,fL,'k'); 
title('Distribuicao Lorentziana'); 
legend(strcat('T2*=',num2str(T22),'ms'))
xlabel('Frequ�ncia (x10^3 Hz)'); ylabel('Intensidade'); grid on;
set(p,'LineWidth',1.5);

figure(2)
p=plot(tempo,Ms(1,:),'b-',tempo,Ms(2,:),'k-',tempo,Ms(3,:),'r--');
legend('M_x','M_y','M_z'); title('Invers�o Recupera��o: 180x - 90x - Fid');
xlabel('Tempo (ms)'); ylabel('Intensidade'); 
grid on;  set(p,'LineWidth',1.5);

   
% Comando exportar imagem 
%print -dpng PulsoSimples.png -r300   

%  Comando exportar dados .dat
%dlmwrite('Teste.dat', [tempo.'  Ms(1,:).'  Ms(2,:).'  Ms(3,:).'], '\t')  


