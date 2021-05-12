%Sabrina da Silva Sá
%Orientador: Andrey Chaves

clc; clear all;

N = 10000; %numero de pontos no eixo da posição
nstates = 10;

r = zeros(10001,1);
m = zeros(10002,1);

m(10002) = 1;

Ry = 13.605*10^3; %energia de Rydberg em meV
a0 = 0.5292; %raio de Bohr em Angstron
wf = zeros(N,N);
E = zeros(N,N);
hbar = 0.6582; %meV*ps

dr(1) = 0.0001/a0; %tamanho do passo em x
dr(2) = 0.04/a0; %tamanho do passo em x
dr(3) = 0.05/a0; %tamanho do passo em x

rmin = 10^-8; %melhor nao usar o ponto r = 0, pode gerar casos do tipo "1/0"

rdim = rmin;
for ir = 1:N+1
 m(ir) = 0.2;
 if(ir <= 1000)
  rdim = rdim + dr(1); 
 elseif(ir >= 1000 && ir <= 5000)
  rdim = rdim + dr(2); 
 elseif(ir >= 5000)
  rdim = rdim + dr(3); 
 end
 r(ir) = rdim;
end

epsilon1 = 1.0; %constante dieletrica em cima do plano do semicondutor 2D 
epsilon2 = 1.0; %constante dieletrica em baixo do plano do semicondutor 2D 
epsilon = 13; %constante dieletrica do proprio semicondutor 2D - numero dado em comparacao com a do vacuo; nesse caso, por exemplo, epsilon1 = 13*epsilon0.
d = 5.0/a0; %espessura da camada do semicondutor 2D (em Angstron/a0)
rho0 = d*epsilon/(epsilon1+epsilon2); %comprimento de blindagem eletrostatica, ou "screening"

for ir = 1:N+1
   %V(ir) = -2/r(ir); %interacao de Coulomb	
   V(ir) = -(2*pi/((epsilon1+epsilon2)*rho0))*(struve(0,r(ir)/rho0)-bessely(0,r(ir)/rho0)); %interacao de Rytova-Keldysh, N.S. Rytova, Proc. MSU, Phys., Astron. 3, 30 (1967), L. V. Keldysh, JETP Lett. 29, 658 (1979)
end

file1 = fopen('potInt.dat','w');

for ir = 1:N+1
 fprintf(file1,'%8.6E \t %8.6E \t %8.6E\n',r(ir)*a0,V(ir)*Ry);
end

l = 0; %momentum angular
l2 = l^2;

k = 0;
for ir = 1:N+1
    
 if(ir ~= 1)
  rojm1=r(ir-1);
 else
  rojm1=rmin;
 end
 
 roj=r(ir);
 if(ir ~= N+1)
  rojp1=r(ir+1);
 else
  rojp1=r(ir);
 end
        
 rojp12=(rojp1+roj)/2;
 if(ir == N+1)
   rojp12=roj;
 end
 
 rojm12=(roj+rojm1)/2;
 if(ir == 1)
  rojp12=roj;
 end
	   
 hj=sqrt((rojp12*rojp12-rojm12*rojm12)/2);
	   
 vhj(ir)=hj;
 
 if(ir <= N-1)
  rojp2=r(ir+2);
 else
  rojp2=r(ir)+2*dr(3);
 end
	 
 rojp12p1=(rojp2+rojp1)/2;
 rojm12p1=(rojp1+roj)/2;

 hjp1=sqrt((rojp12p1*rojp12p1-rojm12p1*rojm12p1)/2);

 if(ir == 1)     
  muj=0;
 else
  muj=(roj/m(ir)+rojm1/m(ir-1))/(2*(roj-rojm1));
 end
 if(ir == N+1)
  muj=0;
 end  
 if(ir ~= N+1)
  mujp1=(rojp1/m(ir+1)+roj/m(ir))/(2*(rojp1-roj));
 end
 if(ir == N)
  mujp1=0;
 end
 if(ir == N+1)
  mujp1=0;
 end

 k = k+1;
 IH(k) = ir;
 JH(k) = ir;
 EH(k) = l2/(r(ir)^2) + V(ir) + (mujp1+muj)/(hj*hj);

 if(ir <= N)
  k = k+1;
  IH(k) = ir;
  JH(k) = ir+1;
  EH(k) = -mujp1/(hj*hjp1);	

  k = k+1;
  IH(k) = ir+1;
  JH(k) = ir;
  EH(k) = -mujp1/(hj*hjp1);
 end
	
end

H = sparse(IH,JH,EH);

[wf, E] = eigs(H,nstates,-5);

En = sort(diag(E));

En(1)*Ry
En(2)*Ry
En(3)*Ry
