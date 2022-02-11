clear; clc;
%Programa Pulso Simples - (Maio/2020) revis�o

%------------------------------------------------------------------------%
% Este c�digo calcula o comportamento da Magnetiza��o quando aplicamos um
% pulso simples de �ngulo theta;
% Tiago Bueno de Moraes  -  (tiagobuemoraes@gmail.com) - Qu�mica Nova 2013
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
theta = 90;             % (graus) Angulo do Pulso de Radio-Frequ�ncia
theta = theta*pi/180;   % convertendo �ngulo para radianos;
T = 1000;               % (ms)  Tempo total aquisi��o FID

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
Niso = 1500;
df = linspace(x01,x02,Niso);       

% Lorentziana Isocromatas - Campo Inomog�neo
fL = (T22)./( 1 + ((df - x0).^2).*(2*pi*T22).^2 );

% (ms)  tempo entre pontos
dT = 1;            

% ===== Matrizes Rota��o =====
Rtheta = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];  % Fase x
E1 = exp(-dT/T1);  E2 = exp(-dT/T2); 
E = [E2 0 0; 0 E2 0; 0 0 E1];       B = [0 0 1-E1].';

% ===== Contadores Sinal =====
N0 = round(T/dT);         % n� pontos total;    tempo total = N0*dT
M  = zeros(3,N0);         % Cria vetor magnetiza��o tamanho (3 x N0)
Ms = zeros(3,N0);         % Cria vetor magnetiza��o Ms (armazenar sinal)
M0 = [0;0;1];             % Magnetiza��o posi��o inicial no eixo z


for f=1:length(df)   
  phi = 2*pi*(df(f))*dT;     
  Rphi = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
       
  M(:,1)= Rtheta*M0 + B;             % Aplica��o pulso angulo theta
  for k=2:(N0)                       % evolu��o sinal
	  M(:,k) = E*Rphi*M(:,k-1)+B;
  end;
  
  % ------- Calculo peso Inomogeneidade - T2* ----------
  % Pico Lorentziano
  g = (T22)./((1+((df(f)-x0).^2).*(2*pi*T22).^2));
   
  % Somando as componentes x y z das isocromata
  Ms = [g*M(1,:)+Ms(1,:); g*M(2,:)+Ms(2,:); g*M(3,:)+Ms(3,:)];
end;

% Normaliza sinal FID
%Ms = Ms/max(Ms(:,1));

% Curvas de Refer�ncia
tempo = [0:N0-1]*dT;  
A = sqrt(Ms(1,1)^2+Ms(2,1)^2);
CurvaT2 = A*exp(-tempo/T2);  
CurvaT22 = A*exp(-tempo/T22);

% ===== Grafico Resultados ======
figure(1)
subplot(2,1,1);
p = plot(df,fL,'k'); 
title('Distribuicao Lorentziana'); 
legend(strcat('T2*=',num2str(T22),'ms'))
xlabel('Frequ�ncia (x10^3 Hz)'); ylabel('Intensidade'); grid on;
set(p,'LineWidth',1.5);


subplot(2,1,2);
p=plot(tempo,Ms(1,:),'b-',tempo,Ms(2,:),'k-',tempo,Ms(3,:),'r--',tempo, CurvaT2,'g:',tempo, CurvaT22,'y-');
legend('M_x','M_y','M_z','T_2'); title('Pulso Simples');
xlabel('Tempo (ms)'); ylabel('Intensidade'); 
grid on;  set(p,'LineWidth',1.5); 

% Comando exportar imagem 
%print -dpng PulsoSimples.png -r300   

%  Comando exportar dados .dat
%dlmwrite('Teste.dat', [tempo.'  Ms(1,:).'  Ms(2,:).'  Ms(3,:).'], '\t')  



