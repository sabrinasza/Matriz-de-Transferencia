% Sabrina da Silva Sá
%Orientador: Andrey Chaves 
clc; clear all;

N = 1000; %numero de pontos no eixo da posiÃ§Ã£o
nstates = 5;

Ry = 13.605*10^3; %energia de Rydberg em meV
a0 = 0.5292; %raio de Bohr em Angstron
wf = zeros(N,N);
E = zeros(N,N);
hbar = 0.6582; %meV*ps

RR = 150/a0;  %raio do fio - retire o % se quiser um fio comum
%RRi = 150/a0;  %raio interno da casca - retire o % se quiser uma "casca"
RRe = 250/a0;  %raio externo da casca - retire o % se quiser uma "casca"
V0 = 100000/Ry; %altura da barreira de potencial

m1 = 1.0; %massa dentro do fio
m2  = 1.0; %massa fora do fio

dr = 0.25/a0; %tamanho do passo em x

rmin = 10^-8; %melhor nao usar o ponto r = 0, pode gerar casos do tipo "1/0"

for ir = 1:N+1
  r(ir) = ir*dr;	
end

for ir = 1:N+1
  if(r(ir) < RR) %use esse "if" se quiser um fio comum
  %if(r(ir) > RRi && r(ir) < RRe) %use esse "if" se quiser uma casca     
   V(ir) = 0;
   m(ir) = m1;
  else
   V(ir) = V0;
   m(ir) = m2;
  end	
end

file1 = fopen('potEmassa.dat','w');
file2 = fopen('ExB.dat','w');

for ir = 1:N+1
 fprintf(file1,'%8.6E \t %8.6E \t %8.6E\n',r(ir)*a0,V(ir)*Ry,m(ir));
end

l = -3; %momentum angular
l2 = l^2;

for iB = 0:20

B = iB*1.5

omegac = 0.176*B; %em 1/ps

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
  rojp2=r(ir)+2*dr;
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
 EH(k) = l2/(r(ir)^2) + V(ir) + (mujp1+muj)/(hj*hj) + 0.5*l*hbar*omegac/Ry + 0.7108*omegac^2*roj^2*a0^2*10^-5/Ry;

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

[wf, E] = eigs(H,nstates,0.001/Ry);

En = sort(diag(E));

fprintf(file2,'%8.6E \t %8.6E \t %8.6E\t %8.6E\n',B,En(1)*Ry,En(2)*Ry,En(3)*Ry);

end %end do loop no campo magnÃ©tico

En(1)*Ry
En(2)*Ry
En(3)*Ry 

disp('o que deveriam ser:')
Ry*(2.4048D0/RR)^2 %para l = 0
Ry*(5.5201D0/RR)^2 %para l = 0
Ry*(8.6537D0/RR)^2 %para l = 0
Ry*(3.8317D0/RR)^2 %para l = +-1
Ry*(7.0156D0/RR)^2 %para l = +-1
Ry*(10.1735D0/RR)^2 %para l = +-1
