% Igor Lima & Sabrina Sá 
% Orientador: Andrey Chaves
clc; clear all;

a0 = 0.5292; %raio de Bohr
Ry = 13.605*10^3; %energia de Rydberg

N = 500; %numero de interfaces
B = 60/a0; %largura das regioes de barreira
W = 50/a0; %largura das regioes de poÃ§o
V0 = 80/Ry; %altura do potencial das barreiras

WX = 200/a0; %tamanho do eixo x

mw = 0.5; %massas efetivas
mb = 0.5;

deltax = WX/N;

for iInt = 1:N;
 d(iInt) = (iInt-1)*deltax; %posiÃ§oes das interfaces 
end

for iV = 1:N
 if(d(iV) >= 15/a0 && d(iV) <= 15/a0 + B || d(iV) >= 15/a0 + B + W && d(iV) <= 15/a0 + B + W + B)    
  V(iV) = V0; 
  m(iV) = mb;
 else
  V(iV) = 0;
  m(iV) = mw;
 end
end

V(N+1) = 0;
m(N+1) = mw;

%plot(d*a0,V*Ry)
%ylim([-1 100])

file1 = fopen('Transmissao.dat','w');

dE = 0.1/Ry; %passo da energia
for E = 0.1/Ry:dE:100/Ry;
 
 for j = 1:N %indice da interface
  k(j) = sqrt((E-V(j))*m(j)); 
  kp(j) = sqrt((E-V(j+1))*m(j+1));   
  M = [exp(1i*k(j)*d(j)) exp(-1i*k(j)*d(j));exp(1i*k(j)*d(j))*(1i*k(j)/m(j)) exp(-1i*k(j)*d(j))*(-1i*k(j)/m(j))];
  Mp = [exp(1i*kp(j)*d(j)) exp(-1i*kp(j)*d(j));exp(1i*kp(j)*d(j))*(1i*kp(j)/m(j+1)) exp(-1i*kp(j)*d(j))*(-1i*kp(j)/m(j+1))];
  MB{j} = Mp\M; %matriz relativa Ã  interface j
 end
 
 MT = eye(2); %iniciar com uma matriz identidade 2x2
 for j = 1:N
  MT = MB{j}*MT;%multiplicacao das matrizes de cada interface, (T 0)^T = MT (1 R)^T
 end
 
 Ref = abs(-MT(2,1)/MT(2,2))^2; %coeficiente de reflexao    
 Tr = 1-Ref; %coeficiente de transmissao
 fprintf(file1,'%E \t %E \n',E*Ry,Tr);    
end

fclose(file1);
   
