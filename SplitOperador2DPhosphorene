program SplitOperator2DMultilayerPhosphorene  
!######################################################
! Programa "Dinâmica de pacotes de onda em fosforeno
! em um sistema bidimensional submetido a um potencial
! externo modelado com o Método de Split Operator"
!######################################################
!                  Versão 1.0
!Diego Rabelo da Costa & Gabriel Oliveira de Sousa
! Orientadores: Gil de Aquino Farias & Andrey Chaves
!######################################################

parameter (l=2, nx=512, ny=512) ! 'nx' e 'ny' sÃ£o o nÃºmero de passos em x e y, que deve ser definido como (2.d0)**N;
parameter (Nlayer = 1) ! number of BP layer
implicit real*8 (a,b,d-h,o-z) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o reais por definiÃ§Ã£o
implicit complex*16 (c) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o complexas por definiÃ§Ã£o
dimension X(nx), Y(ny), POT(nx,ny) !coordenadas x, y e potencial
dimension AKX(nx), AKY(ny), nn(l) 
real*8 S(nx,ny),Sx(nx,ny),Sy(ny),f(nx,ny)
double complex CS1(nx,ny),CS2(nx,ny)
double complex CPSI(2*nx), CEXPV(nx,ny) !funÃ§Ãµes de onda
double complex CPSIA(nx,ny), CPSIB(nx,ny) !Spinor |Psi> = [Phi_A,Phi_B]**t
double complex CMKA(nx,ny), CMKB(nx,ny), CPSIKA(nx,ny), CPSIKB(nx,ny)
integer i, nxy, nxy2, kx, ky
real*8 DKX, DKY, DKXY, AKX0, AKY0
 CHARACTER(len=40)          :: Filename, FilenameR, FilenameI
 CHARACTER(len=8)           :: middle

common / constants / dx, dy, cdt, ci, Ry, a0, hc, pi !constantes que serÃ£o utilizadas em comum nas subroutinas
common / constantsE / u0, eta_x, eta_y, delta, xi, gamma_x, gamma_y !constantes que serÃ£o utilizadas para o cÃ¡lculo da energia
common / constantsPOT / WIDTH, ipropagation !constantes para o potential

open(unit = 32, file = 'xyavg.dat')
open(unit = 33, file = 'prob.dat')

310 FORMAT(10F16.6)

imovie = 0
ipot = 1
iwave = 1
iwaveK = 0
ipropagation = 0 ! se 0 entÃ£o a propagaÃ§Ã£o Ã© ao longo da direÃ§Ã£o X, se 1 Ã© ao longo da direÃ§Ã£o Y

pi = 4.D0*datan(1.D0)
TWOPI = 2.0D0*pi !ParÃ¢metro que vai ser usado para calcular os ' k's '
a0 = 0.5292D0 ! 'a0' Ã© o raio de Bohr (dado em Angstron)
Ry = 13.6057D0 ! 'Ry' Ã© a energia de Rydberg (dada em eV)
hc = 1.054589D0/1.602189D0 ! 'hc' Ã© a constante h cortado (dada em eV*fs)
x0 = 0.D0/a0 ! PosiÃ§Ã£o x inicial do pacote de onda
y0 = 0.D0/a0 ! PosiÃ§Ã£o y inicial do pacote de onda

deltaAA = -0.338D0 ! eV
deltaAB = -2.912D0 ! eV
deltaAC = 3.831D0 ! eV
deltaAD = -0.076D0 ! eV
etaAA = 1.161D0  ! eV.Angs^2
etaAB = 2.05D0 ! eV.Angs^2
etaAC = 0.460D0 ! eV.Angs^2
etaAD = 0.104D0 ! eV.Angs^2
gamAA = -1.563D0 ! eV.Angs^2
gamAB = 3.607D0 ! eV.Angs^2
gamAC = -1.572D0 ! eV.Angs^2
gamAD = 0.179D0 ! eV.Angs^2
akhiAB = 3.688D0 ! eV.Angs
akhiAC = 2.208D0 ! eV.Angs
deltaACprime = 0.712D0 ! eV
deltaADprime = -0.132D0 ! eV
etaACprime = -0.9765D0 ! eV.Angs^2
etaADprime = 2.699D0 ! eV.Angs^2
gamACprime = 2.443D0 ! eV.Angs^2
gamADprime = 0.364D0 ! eV.Angs^2
akhiACprime = 2.071D0 ! eV.Angs

u0 = deltaAA + deltaAD
etaX = etaAA + etaAD
etaY = gamAA + gamAD
delta = deltaAB + deltaAC
gamX = etaAB + etaAC
gamY = gamAB + gamAC
akhi = akhiAB + akhiAC

alambN = dcos(pi/(Nlayer+1.D0)) ! alambN = cos(n*pi/(Nlayer+1)); verificar o valor de n para Nlayer > 1, onde n=1,Nlayer

u0n = u0 + alambN*deltaADprime
etaXn = etaX + alambN*etaADprime
etaYn = etaY + alambN*gamADprime
deltan = delta + alambN*deltaACprime
gamXn = gamX + alambN*etaACprime
gamYn = gamY + alambN*gamACprime
akhin = akhi + alambN*akhiACprime

u0 = u0n/Ry
eta_x = etaXn/(Ry*a0*a0)
eta_y = etaYn/(Ry*a0*a0)
delta = deltan/Ry
xi = akhin/(Ry*a0)
gamma_x = gamXn/(Ry*a0*a0)
gamma_y = gamYn/(Ry*a0*a0)

 ci = (0.D0,1.D0)  ! 'ci' Ã© o nÃºmero imaginÃ¡rio puro i
 cdt = 0.1D0 ! 'cdt' Ã© o intervalo de tempo (dado em fs - femtosegundos) ===PASSO===
nt = 20000.D0 ! 'nt' Ã© o nÃºmero de iteraÃ§Ãµes temporais

 CCOEF = -ci*cdt*Ry/hc !coeficiente no operador de evoluÃ§Ã£o temporal. NÃ³s mutiplicamos por Ry para fazer isso adimensional

WX = 2560.D0/a0 !tamanho do sistema
dx = WX/nx  !dx = WX/(nx-1) !tamanho do passo
WY = 2560.D0/a0 !tamanho do sistema
dy = WY/ny 

 nn(1) = nx
 nn(2) = ny

nxy = nx*ny 
nxy2 = 2*nxy

DKX = TWOPI/WX
DKY = TWOPI/WY
DKXY= DKX*DKY

AKBX = 0.D0
AKBY = 0.D0

AKX0 = AKBX * DKX / 2.D0
AKY0 = AKBY * DKY / 2.D0

nnx=nx/2
nny=ny/2

 DO 717 KX = 1, nx
  K = KX-1
  Q = DKX*K
  IF(K.GT.nnx) Q = DKX*(K-nx)
  AKX(KX) = Q + AKX0
717   CONTINUE
 DO 718 KY = 1, ny
  K = KY-1
  Q = DKY*K
  IF(K.GT.nny) Q = DKY*(K-ny)
  AKY(KY) = Q + AKY0
718   CONTINUE

do i = 1,nx
  X(i) = -WX/2.D0 + (i-1)*dx  !O sistema pertence ao intervalo [-WX/2,WX/2] 
end do
do j = 1,ny
  Y(j) = -WY/2.D0 + (j-1)*dy  !O sistema pertence ao intervalo [-WX/2,WX/2] 
end do

xsig = 300.D0/a0 !largura do pacote de onda inicial
ysig = 300.D0/a0 !largura do pacote de onda inicial

 call wavei(x0,xsig,X,y0,ysig,Y,CPSIA,CPSIB) !construindo o pacote de onda inicial

if(iwave.eq.1) then
 open(unit = 10, file = 'psi0.dat')
 do i = 1,nx
  do j = 1, ny
   write(10,310) X(i)*a0, Y(j)*a0, cdabs(CPSIA(i,j))**2/a0/a0+cdabs(CPSIB(i,j))**2/a0/a0 !escrevendo o mÃ³dulo quadrado do pacote de onda inicial
  enddo
   write(10,310) !Pular linha para plotar no gnuplot
 enddo
endif
 
V0 = 0.0D0/Ry  !Altura do potencial barreira
WIDTH = 0.0D0/a0

 call POTEXT(V0,X,Y,POT,it)

if(ipot.eq.1) then
 open(unit = 21, file = 'POT.dat')
  do i = 1,nx
   do j = 1,ny
      write(21,310)x(i)*a0,y(j)*a0,POT(i,j)*Ry
   enddo
  enddo
endif

!============================================================================================================
do it = 1, nt 
                                                                      !ComeÃ§o do cÃ¡lculo da evoluÃ§Ã£o temporal
!============================================================================================================

  ti = it*cdt

  write(*,*) it

  do i = 1,nx
   do j = 1,ny
    CEXPV(i,j) = cdexp(CCOEF*POT(i,j)/2.D0) !exponencial do termo externo
   enddo
  enddo

  do i = 1, nx 
   do j = 1, ny
    CPSIA(i,j) = CEXPV(i,j)*CPSIA(i,j) !Primeira multiplicaÃ§Ã£o final do termo da exponencial para Psi_A
    CPSIB(i,j) = CEXPV(i,j)*CPSIB(i,j) !Primeira multiplicaÃ§Ã£o final do termo da exponencial para Psi_B
   enddo
  enddo

!============================================================================================================
                                                         ! Realizando a Transformada Direta de Fourier RÃ¡pida 
!============================================================================================================

        isign = 1
        
!Para CPSIA

        call FFT(CPSIA,nx,ny,nn,2,isign,nxy2)

     do kx = 1, nx
       do ky = 1, ny    
       CPSIKA(kx,ky) = CPSIA(kx,ky) 
       end do
     end do   

!Para CPSIB

        call FFT(CPSIB,nx,ny,nn,2,isign,nxy2) 
 
     do kx = 1, nx
       do ky = 1, ny    
       CPSIKB(kx,ky) = CPSIB(kx,ky)  
       end do
     end do 
!============================================================================================================

     do kx = 1, nx
      do ky = 1, ny
 
         f(kx,ky) = u0 + eta_x*AKX(kx)*AKX(kx) + eta_y*AKY(ky)*AKY(ky) 
         S(kx,ky) = (dreal(cdt)/hc)*dsqrt((delta + gamma_x*AKX(kx)*AKX(kx) + gamma_y*AKY(ky)*AKY(ky))**2 + (xi*AKY(ky))**2) 
         Sx(kx,ky) = (dreal(cdt)/hc)*(delta + gamma_x*AKX(kx)*AKX(kx) + gamma_y*AKY(ky)*AKY(ky))
         Sy(ky) = -(dreal(cdt)/hc)*(xi*AKY(ky))
         CS1(kx,ky) = (Sx(kx,ky)-ci*Sy(ky))/S(kx,ky)
         CS2(kx,ky) = (Sx(kx,ky)+ci*Sy(ky))/S(kx,ky)
     
         CMKA(kx,ky)=cdexp(-ci*f(kx,ky)*dreal(cdt)/hc)*(dcos(S(kx,ky))*CPSIKA(kx,ky)-ci*(dsin(S(kx,ky))*CS1(kx,ky)*CPSIKB(kx,ky)))
         CMKB(kx,ky)=cdexp(-ci*f(kx,ky)*dreal(cdt)/hc)*(dcos(S(kx,ky))*CPSIKB(kx,ky)-ci*(dsin(S(kx,ky))*CS2(kx,ky)*CPSIKA(kx,ky)))
 
      enddo
     enddo

!============================================================================================================
                         ! Realizando a Transformada Direta de Fourier RÃ¡pida do Termo exp[CCOEF*POT(i)/2.D0]
!============================================================================================================ 
if(iwave.eq.1.and.iwaveK.eq.1) then
 open(unit = 11, file = 'psi0k.dat')   
 do kx = 1, nx
  do ky = 1, ny
    if(it.eq.1)write(11,310) AKX(kx)/a0, AKY(ky)/a0, cdabs(CMKA(kx,ky))**2 + cdabs(CMKB(kx,ky))**2
  end do
   write(11,*) !Pular linha para plotar no gnuplot 
 end do
end if

     isign = -1

!Para CPSIA

     call FFT(CMKA,nx,ny,nn,2,isign,nxy2)

     do i = 1, nx
       do j = 1, ny    
       CPSIA(i,j) = CMKA(i,j)  
       end do
     end do

!Para CPSIB

     call FFT(CMKB,nx,ny,nn,2,isign,nxy2)

     do i = 1, nx
       do j = 1, ny    
       CPSIB(i,j) = CMKB(i,j)  
       end do
     end do
!============================================================================================================ 

  do i = 1, nx 
   do j = 1, ny
    CPSIA(i,j) = CEXPV(i,j)*CPSIA(i,j) !MultiplicaÃ§Ã£o final do termo da exponencial para Psi_A
    CPSIB(i,j) = CEXPV(i,j)*CPSIB(i,j) !MultiplicaÃ§Ã£o final do termo da exponencial para Psi_B
   enddo
  enddo

      call centermass(CPSIA,CPSIB,X,Y,XCM,XCMr,XCMt,YCM,YCMr,YCMt) 
      write(32,310) cdabs(it*cdt), XCM, YCM, XCMr, YCMr, XCMt, YCMt
      call prob(X,Y,CPSIA,CPSIB,PBEF,PINT,PAFT)
      write(33,310) cdabs(it*cdt), PBEF, PINT, PAFT   

     if(imovie.eq.1) then

      if(mod(it,1000).eq.0) then

       WRITE(middle, '(I8.8)') it 
       Filename = 'PsiT'// middle //'.dat'
       FilenameR = 'PsiR'// middle //'.dat'
       FilenameI = 'PsiI'// middle //'.dat'
       OPEN(UNIT=13, FILE=Filename, STATUS='REPLACE')
       OPEN(UNIT=14, FILE=FilenameR, STATUS='REPLACE')
       OPEN(UNIT=15, FILE=FilenameI, STATUS='REPLACE')

       do i = 1,nx
        do j = 1,ny
           write(13,310) X(i)*a0, Y(j)*a0, cdabs(CPSIA(i,j))**2/a0/a0+cdabs(CPSIB(i,j))**2/a0/a0 ! Escrevendo o mÃ³dulo quadrado.
           write(14,310) X(i)*a0, Y(j)*a0, dreal(CPSIA(i,j)/a0/a0+CPSIB(i,j)/a0/a0) 
           write(15,310) X(i)*a0, Y(j)*a0, dimag(CPSIA(i,j)/a0/a0+CPSIB(i,j)/a0/a0) 
        enddo
        write(13,310) !Pular linha para plotar no gnuplot
        write(14,310) !Pular linha para plotar no gnuplot
        write(15,310) !Pular linha para plotar no gnuplot
       enddo

      CLOSE(13)
      CLOSE(14)
      CLOSE(15)

      endif

     endif
      
  end do !final do laÃ§o do tempo
 
end program

!============================================================================================================
subroutine wavei(x0,xsig,X,y0,ysig,Y,CPSIA,CPSIB)
                                                                        !Construindo o pacote de onda inicial
!============================================================================================================
parameter (nx=512, ny=512)
implicit real*8 (a,b,d-h,o-z) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o reais por definiÃ§Ã£o
implicit complex*16 (c) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o complexas por definiÃ§Ã£o
dimension X(nx),Y(ny)
double complex CPSIA(nx,ny),CPSIB(nx,ny)

common / constants / dx, dy, cdt, ci, Ry, a0, hc, pi !constantes que serÃ£o utilizadas em comum nas subroutinas
common / constantsE / u0, eta_x, eta_y, delta, xi, gamma_x, gamma_y !constantes que serÃ£o utilizadas para o cÃ¡lculo da energia

ak0 = 0.1!027002D0
theta = (90.D0*pi)/180.D0 
akx0 = dsin(theta)*ak0
aky0 = dcos(theta)*ak0

E = u0 + eta_x*(akx0**2) + eta_y*(aky0**2) + dsqrt((delta + gamma_x*(akx0**2) + gamma_y*(aky0**2))**2 + (xi**2)*(aky0**2)) 
write(*,*) 'Energy =', E*Ry

cac1 = 1.D0
cac2 = 0.D0 

do i = 1,nx
 do j = 1,ny

  CPSIA(i,j) = cac1*cdexp(ci*akx0*X(i))*cdexp(ci*aky0*Y(j))*dexp(-(X(i)-x0)**2/(2.D0*(xsig)**2)-(Y(j)-y0)**2/(2.D0*(ysig)**2)) !pacote de onda inicial Ã© uma Gaussiana
  CPSIB(i,j) = cac2*cdexp(ci*akx0*X(i))*cdexp(ci*aky0*Y(j))*dexp(-(X(i)-x0)**2/(2.D0*(xsig)**2)-(Y(j)-y0)**2/(2.D0*(ysig)**2)) !pacote de onda inicial Ã© uma Gaussiana

 enddo 
enddo

!normalizaÃ§Ã£o
 ASUM = 0.D0
 do i = 1,nx !integral
  do j = 1,ny
   ASUM = ASUM + (cdabs(CPSIA(i,j))**2 + cdabs(CPSIB(i,j))**2)
  enddo
 enddo

 ASUM = dsqrt(ASUM*dx*dy)
!====================NORMALIZAÃ‡ÃƒO===========================================
	
 if(ASUM.eq.0.D0)then !NORMALIZING ONLY NON-ZERO WAVEFUNCTIONS
  AINV = 0.D0
 else
  AINV = 1.D0 / ASUM
 endif

 do i = 1,nx
  do j = 1,ny
   CPSIA(i,j) = CPSIA(i,j)*AINV  !Dividindo pelo fator de normalizaÃ§Ã£o para Psi_A 
   CPSIB(i,j) = CPSIB(i,j)*AINV  !Dividindo pelo fator de normalizaÃ§Ã£o para Psi_B
  enddo
 enddo
!====================NORMALIZAÃ‡ÃƒO===========================================

return
end

!============================================================================================================
subroutine POTEXT(V0,X,Y,POT,it)
                                                                                     !Construindo o potencial
!============================================================================================================
parameter (nx=512, ny=512)
implicit real*8 (a,b,d-h,o-z) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o reais por definiÃ§Ã£o
implicit complex*16 (c) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o comlexas por definiÃ§Ã£o
dimension X(nx), POT(nx,ny), Y(ny)
common / constants / dx, dy, cdt, ci, Ry, a0, hc, pi !constantes que serÃ£o utilizadas em comum nas subroutinas
common / constantsPOT / WIDTH, ipropagation !constantes para o potential

       do i = 1, nx
         do j = 1, ny  
         POT(i,j) = 0.D0

         !!!!!!!!!!! barreira inclinada paralela ao eixo y !!!!!!!!!!!!!
!         if(X(i).ge.0.D0.and.X(i).le.WIDTH)then
!           POT(i,j) = (V0/WIDTH)*X(i)
!         endif
!         if(X(i).gt.WIDTH)then
!           POT(i,j) = V0
!         endif

         !!!!!!!!!!! barreira abrupta paralela ao eixo y !!!!!!!!!!!!!
!         if(X(i).ge.0.D0)then
!           POT(i,j) = V0
!         endif

         !!!!!!!!!!! barreira inclinada paralela ao eixo x !!!!!!!!!!!!!
!         if(Y(j).ge.0.D0.and.Y(j).le.WIDTH)then
!           POT(i,j) = (V0/WIDTH)*Y(j)
!         endif
!         if(Y(j).gt.WIDTH)then
!           POT(i,j) = V0
!         endif

         !!!!!!!!!!! barreira abrupta paralela ao eixo x !!!!!!!!!!!!!
!         if(Y(j).ge.0.D0)then
!           POT(i,j) = V0
!         endif

         end do
       end do

return
   
end
!==================================================================================================
subroutine centermass(CPSIA,CPSIB,X,Y,XCM,XCMr,XCMt,YCM,YCMr,YCMt)
                                                                                         !<x> e <y>
!==================================================================================================
parameter (nx=512, ny=512)
implicit real*8 (a,b,d-h,o-z) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o reais por definiÃ§Ã£o
implicit complex*16 (c) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o complexas por definiÃ§Ã£o
dimension X(nx),Y(ny)
double complex CPSIA(nx,ny),CPSIB(nx,ny)
common / constants / dx, dy, cdt, ci, Ry, a0, hc, pi !constantes que serÃ£o utilizadas em comum nas subroutinas
common / constantsPOT / WIDTH, ipropagation !constantes para o potential

      XCM = 0.D0 ! Calculando a mÃ©dia de x (<x>) no espaÃ§o real
      XCMr = 0.D0
      XCMt = 0.D0
      YCM = 0.D0 ! Calculando a mÃ©dia de y (<y>) no espaÃ§o real
      YCMr = 0.D0
      YCMt = 0.D0
      do i = 1, nx
       do j = 1, ny
        XCM = XCM + dconjg(CPSIA(i,j))*X(i)*CPSIA(i,j) + dconjg(CPSIB(i,j))*X(i)*CPSIB(i,j)
        YCM = YCM + dconjg(CPSIA(i,j))*Y(j)*CPSIA(i,j) + dconjg(CPSIB(i,j))*Y(j)*CPSIB(i,j) 
        if(ipropagation.eq.0)then
          if(X(i).lt.0.D0)then
            XCMr = XCMr + dconjg(CPSIA(i,j))*X(i)*CPSIA(i,j) + dconjg(CPSIB(i,j))*X(i)*CPSIB(i,j)
          endif
          if(X(i).ge.0.D0)then
            XCMt = XCMt + dconjg(CPSIA(i,j))*X(i)*CPSIA(i,j) + dconjg(CPSIB(i,j))*X(i)*CPSIB(i,j)
	      endif
        endif
        if(ipropagation.eq.1)then
          if(Y(j).lt.0.D0)then
            YCMr = YCMr + dconjg(CPSIA(i,j))*Y(j)*CPSIA(i,j) + dconjg(CPSIB(i,j))*Y(j)*CPSIB(i,j) 
	      endif
          if(Y(j).ge.0.D0)then
            YCMt = YCMt + dconjg(CPSIA(i,j))*Y(j)*CPSIA(i,j) + dconjg(CPSIB(i,j))*Y(j)*CPSIB(i,j) 
	      endif
        endif 
       enddo
      enddo

      XCM = XCM*dx*dy*a0
      XCMr = XCMr*dx*dy*a0
      XCMt = XCMt*dx*dy*a0
      YCM = YCM*dx*dy*a0
      YCMr = YCMr*dx*dy*a0
      YCMt = YCMt*dx*dy*a0

return

end
!============================================================================================================
subroutine prob(X,Y,CPSIA,CPSIB,PBEF,PINT,PAFT)
                                                                                   !Calcula as probabilidades
!============================================================================================================
parameter (nx=512, ny=512)
implicit real*8 (a,b,d-h,o-z) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o reais por definiÃ§Ã£o
implicit complex*16 (c) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o complexas por definiÃ§Ã£o
dimension X(nx),Y(ny),POT(nx,ny)
double complex CPSIA(nx,ny),CPSIB(nx,ny)
common / constants / dx, dy, cdt, ci, Ry, a0, hc, pi !constantes que serÃ£o utilizadas em comum nas subroutinas
common / constantsPOT / WIDTH, ipropagation !constantes para o potential

     SUMB = 0.D0
     SUMI = 0.D0
     SUMA = 0.D0
     do 1 j = 1, ny
       do 1 i = 1, nx
         if(ipropagation.eq.0)then
           if(X(i).lt.0.D0)then
   	     SUMB = SUMB + dconjg(CPSIA(i,j))*CPSIA(i,j) + dconjg(CPSIB(i,j))*CPSIB(i,j)
	       endif
!           if(X(i).ge.0.D0.and.X(i).le.WIDTH)then
! 	     SUMI = SUMI + dconjg(CPSIA(i,j))*CPSIA(i,j) + dconjg(CPSIB(i,j))*CPSIB(i,j)
!	       endif
           if(X(i).ge.0.D0)then !           if(X(i).gt.WIDTH)then
 	     SUMA = SUMA + dconjg(CPSIA(i,j))*CPSIA(i,j) + dconjg(CPSIB(i,j))*CPSIB(i,j)
	       endif
         endif
         if(ipropagation.eq.1)then
           if(Y(j).lt.0.D0)then
   	     SUMB = SUMB + dconjg(CPSIA(i,j))*CPSIA(i,j) + dconjg(CPSIB(i,j))*CPSIB(i,j)
	       endif
!           if(Y(j).ge.0.D0.and.Y(j).le.WIDTH)then
! 	     SUMI = SUMI + dconjg(CPSIA(i,j))*CPSIA(i,j) + dconjg(CPSIB(i,j))*CPSIB(i,j)
!	       endif
           if(Y(j).ge.0.D0)then!           if(Y(j).gt.WIDTH)then
 	     SUMA = SUMA + dconjg(CPSIA(i,j))*CPSIA(i,j) + dconjg(CPSIB(i,j))*CPSIB(i,j)
	       endif
         endif 
1     continue
      
      PBEF = SUMB*dx*dy
      PINT = SUMI*dx*dy
      PAFT = SUMA*dx*dy

return
end

!============================================================================================================
subroutine ENERGY(POT,CPSI,E)
                                                                                       !calculates the energy
!============================================================================================================
parameter (nx=2**12)
implicit real*8 (a,b,d-h,o-z) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o reais por definiÃ§Ã£o
implicit complex*16 (c) !as variÃ¡veis que comeÃ§am com essa letra sÃ£o complexas por definiÃ§Ã£o
double complex CPSI(NX), CKIN(NX), CA(NX), CB(NX), CC(NX)
real*8 E, POT(NX)

common / constants / dx, dy, cdt, ci, Ry, a0, hc, pi !constantes que serÃ£o utilizadas em comum nas subroutinas


 SUMV = 0.D0
 do ix = 1, nx
  SUMV = SUMV + DCONJG(CPSI(ix))*POT(ix)*CPSI(ix)*dx
 enddo

 SUMK = 0.D0
 do ix = 1, nx
  if(ix.eq.nx)then
   CA(ix) = 0.D0
  else
   CA(ix) = -CPSI(ix+1)/2.D0/dx**2 
  endif

  CB(ix)=CPSI(ix)/2.D0/dx**2 

  if(ix.eq.1)then
   CC(ix) = 0.D0
  else
   CC(ix) = -CPSI(ix-1)/2.D0/dx**2 
  endif

  CKIN(ix) = CA(ix) + CB(ix) + CC(ix)

  SUMK = SUMK + DCONJG(CPSI(ix))*CKIN(ix)*dx
 enddo

 E = (SUMV + SUMK)*Ry


return
end

!============================================================================================================
                                                                                     !Transformada de Fourier
!============================================================================================================

!==================================================================
SUBROUTINE FFT(CWAV,NX,NY,NN,NDIM,ISIGN,NXY2)
!                         Fast Fourier Transform in both directions
!==================================================================
      IMPLICIT REAL(8)(A,B,D-H,O-Z)
      IMPLICIT DOUBLE COMPLEX (C)
      DOUBLE COMPLEX CWAV(NX,NY)
      REAL(8) WR, WI, WPR, WPI, WTEMP, THETA
      DIMENSION NN(NDIM), DAT(NXY2)
      NXY = NX * NY
	CI = (0.D0,1.D0)
      DO 21 IY = 1, NY
       I = ( IY - 1 ) * 2
       DO 22 IX = 1, NX
        INX = I * NX
        J = INX + 2 * IX - 1
        K = INX + 2 * IX
        DAT(J) = DREAL(CWAV(IX,IY))
        DAT(K) = DIMAG(CWAV(IX,IY))
22     CONTINUE
21    CONTINUE
      NTOT = 1
      DO 11 IDIM = 1, NDIM
       NTOT = NTOT * NN(IDIM)
11    CONTINUE
      NPREV = 1
      DO 18 IDIM = 1, NDIM
       N = NN(IDIM)
       NREM = NTOT / ( N * NPREV )
       IP1 = 2 * NPREV
       IP2 = IP1 * N
       IP3 = IP2 * NREM
       I2REV = 1
       DO 14 I2 = 1, IP2, IP1
       IF(I2.LT.I2REV) THEN
        DO 13 I1 = I2, I2 + IP1 - 2, 2
         DO 12 I3 = I1, IP3, IP2
         I3REV = I2REV + I3 - I2
         TEMPR = DAT(I3)
         TEMPI = DAT(I3+1)
         DAT(I3) = DAT(I3REV)
         DAT(I3+1) = DAT(I3REV+1)
         DAT(I3REV) = TEMPR
         DAT(I3REV+1) = TEMPI
12       CONTINUE
13      CONTINUE
       END IF
      IBIT = IP2 / 2
1     IF((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
       I2REV = I2REV - IBIT
       IBIT = IBIT / 2
       GO TO 1
      END IF
      I2REV = I2REV + IBIT
14    CONTINUE
      IFP1 = IP1
2     IF(IFP1.LT.IP2) THEN
       IFP2 = 2 * IFP1
       THETA = ISIGN * 6.28318530717959D0 / (IFP2 / IP1 )
       WPR = -2.D0 * DSIN( 0.5D0 * THETA )**2
       WPI = DSIN(THETA)
       WR = 1.D0
       WI = 0.D0
       DO 17 I3 = 1, IFP1, IP1
        DO 16 I1 = I3, I3 + IP1 - 2, 2
         DO 15 I2 = I1, IP3, IFP2
         K1 = I2
         K2 = K1 + IFP1
         TEMPR = WR*DAT(K2)-WI*DAT(K2+1)
         TEMPI = WR*DAT(K2+1)+WI*DAT(K2)
         DAT(K2) = DAT(K1) - TEMPR
         DAT(K2+1) = DAT(K1+1) - TEMPI
         DAT(K1) = DAT(K1) + TEMPR
         DAT(K1+1) = DAT(K1+1) + TEMPI
15       CONTINUE
16      CONTINUE
        WTEMP = WR
        WR = WR * WPR - WI * WPI + WR
        WI = WI * WPR + WTEMP * WPI + WI
17     CONTINUE
       IFP1 = IFP2
       GO TO 2
      END IF
      NPREV = N * NPREV
18    CONTINUE

      ANINV = 1.D0
      IF( ISIGN.EQ.1 ) ANINV = 1.D0 / DFLOAT(NXY)

      DO 19 IY = 1, NY
       I = ( IY - 1 ) * 2
       DO 20 IX = 1, NX
        INX = I * NX
        J = INX + 2 * IX - 1
        K = INX + 2 * IX
        CWAV(IX,IY) = ANINV * ( DAT(J) + CI * DAT(K) )
20     CONTINUE
19    CONTINUE
      RETURN
      END

