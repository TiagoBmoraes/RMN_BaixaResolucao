clear; clc;
%Programa Pulso Simples - (Maio/2020) revis�o

%------------------------------------------------------------------------%
% Este c�digo calcula o comportamento da Magnetiza��o quando aplicamos um
% CPMG (matricial) de NMR.      90y_t_[180x_t_eco]  
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
theta1 = 90;    % Angulo 1� Pulso 90� y  (em graus)
theta2 = 180;   % Angulo 2� Pulso 180� x
T      = 1000;   % (ms) Tempo total aquisi��o FID
Tn1    = 50;    % (ms) Tempo at� Bloco de 180�
Tn2    = 100;   % (ms) Tempo entre Pulsos de 180� (CPMG: Tn2 = 2*Tn1)
Np     = 8;     % Numero de pulsos de 180�

% Se colocar Np maior do que cabe dentro do tempo de aquisi��o vai dar
% erro, ent�o, diminuir Np ou aumentar T

% convertendo �ngulo para radianos;
theta1 = theta1*pi/180;
theta2 = theta2*pi/180;

% Frequencia em [1/ms] % tudo est� em ms, logo 25*10^-3 (1/ms) = 25 Hz !
x0 = 50 * 10^-3 ;      

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
RthetaCP    = [cos(theta1) 0 sin(theta1);0 1 0 ; -sin(theta1) 0 cos(theta1)];  % Fase y
RthetaBloco = [1 0 0; 0 cos(theta2) -sin(theta2); 0 sin(theta2) cos(theta2)];  % Fase x

E1 = exp(-dT/T1);  E2 = exp(-dT/T2); 
E = [E2 0 0; 0 E2 0; 0 0 E1];       B = [0 0 1-E1]'; 

% ===== Contadores Sinal =====
N0 = round(T/dT);              % n� pontos total;    tempo total = N0*dT
N1 = round(Tn1/dT);            % n� pontos CP at� Bloco  (tau/2)
N2 = round(Tn2/dT);            % n� pontos entre pulsos Bloco

% Caso queria N�mero de Echo que cabem no Tempo de aquisi��o:
%Np = floor((N0-(N1+1))/(N2+1));    % n� de pulsos no Bloco
N3 = N0 - (N1+1) - (Np*(N2+1));    % n� pontos depois do ultimo pulso

M = zeros(3,N0);          % Criar vetor magnetiza��o Msignal
Ms = zeros(3,N0);         % Criar vetor magnetiza��o Msignal
M0 = [0;0;1];             % Magnetiza��o posi��o inicial no eixo z


for f=1:length(df)
  phi  = 2*pi*(df(f))*dT;     
  Rphi = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    
  M(:,1)= RthetaCP*M0 + B;             % pulso de theta�x
  for k=2:N1+1                         % evolu��o 
	  M(:,k) = E*Rphi*M(:,k-1)+B;
  end;
  N1c=N1+1;  
   
  for n=1:Np
      M(:,N1c+1) = RthetaBloco*M(:,N1c)+B;     % pulso de 180�x
    for k=1:N2                                 % evolu��o
	  M(:,k+N1c+1) = E*Rphi*M(:,k+N1c)+B;
    end;
  N1c=N1c+N2+1;
  end;
    
  % evolu��o at� o final
  for k=1:N3-1                              
    M(:,k+N1c+1) = E*Rphi*M(:,k+N1c)+B;
  end;
  
% ------- Calculo peso Inomogeneidade - T2* ----------
  % Pico Lorentziano
  g = T22./((1+((df(f)-x0).^2)*(2*pi*T22)^2));
   
  % Somando as componentes x y z das isocromata
  Ms = [g*M(1,:)+Ms(1,:) ; g*M(2,:)+Ms(2,:) ; g*M(3,:)+Ms(3,:)];
  
end;

% Normaliza sinal FID
Ms = Ms/max(abs(Ms(:,1)));

% Curvas
tempo = [0:N0-1]*dT;  
CurvaT2 = exp(-tempo/T2);  
CurvaT22 = exp(-tempo/T22);

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
p=plot(tempo,Ms(1,:),'b-',tempo,Ms(2,:),'k-',tempo,Ms(3,:),'r--',tempo, CurvaT2,'g-',tempo, Mod,'y-'); 
legend('M_x','M_y','M_z','T_{2}','Mod'); title('CPMG: 90y [180x - eco]_{n} Echos [Positivos]_{n}');
xlabel('Tempo (ms)'); ylabel('Intensidade'); 
grid on;  set(p,'LineWidth',1.5);   

% Comando exportar imagem 
%print -dpng PulsoSimples.png -r300   

%  Comando exportar dados .dat
%dlmwrite('Teste.dat', [tempo.'  Ms(1,:).'  Ms(2,:).'  Ms(3,:).'], '\t')  

