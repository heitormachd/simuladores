%Discretização
Lx=1e-2;
Ly=1e-2;%comprimento de propagação em x e y em [m]
dx= 10e-5;%tamanho do passo em [m / pixel]
dy=dx; %superfície regular
nx=fix(Lx/dx);%n total de pontos em [pixel]
ny=nx;
x=linspace(0, Lx, nx);
y=linspace(0, Ly, ny); %cria um vetor de pontos discretos para x e y comecando em 0 indo ate L com um total de n pontos
T=10e-6; %tempo total em [s]
dt = 8e-9; % discretização temporal em [s / amostra] ou [s / timestep]
nt=(T/dt); % n total de amostras temporais ou timesteps


%variaveis de campo
pn=zeros(nx,ny);%gera um grid para armazenar a pressao no tempo
pnf=pn;% p(n+1)
pnb=pn;% p(n-1)

%parametros
c=5500;  % velocidade de propagação no meio em [m/s]
CFL = c * dt / dx; %CFL=0.5; %CFL = c.dt/dx
%dt=CFL*dx/c;

%fonte
Tt=linspace(0, T-dt, nt); % como começa em zero tem que descontar dt do total do vetor de tempo
freq = 5e6; % frequencia em [Hz]
delay = 5e-7; % atraso do sinal em [s], se 0 (zero) o pico no pulso é no t = 0
bandwidth = 0.82; % largura de banda em [%]
% sn= gauspuls(Tt,1000,0.5, -6);
sn= gauspuls(Tt-delay, freq, bandwidth);
figure();plot(Tt, sn);



%loop
tt=1;

while(tt < nt)
    % troca as variáveis
    pnb=pn;
    pn=pnf; %deixa salvo p(n-1), p(n) e p(n+1)

    for i=2:nx-1
        for j=2:ny-1
            pnf(i,j)= 2*pn(i,j) - pnb(i,j) ...
                + CFL^2 * (pn(i+1,j) + pn(i,j+1) - 4*pn(i,j) + pn(i-1,j) + pn(i,j-1));
        end
    end
    
    %condicoes de fronteira(reflexao)
    %pn(:,[1 end])=0;
    %pn([1 end],:)=0;

    %condicoes de frontira(absorcao)
    pnf(1,:)=pn(2,:) + ((CFL-1)/(CFL+1))*(pnf(2,:)-pn(1,:));
    pnf(end, :)=pn(end-1,:) + ((CFL-1)/(CFL+1))*(pnf(end-1,:)-pn(end,:));
    pnf(:,1)=pn(:,2) + ((CFL-1)/(CFL+1))*(pnf(:,2)-pn(:,1));
    pnf(:,end)=pn(:,end-1) + ((CFL-1)/(CFL+1))*(pnf(:,end-1)-pn(:,end));

    %source
    % pn(50,50)=dt^2*20*sin(30*pi*t/20);
    pnf(fix(ny/2), fix(nx/2)) = pnf(fix(ny/2), fix(nx/2)) + sn(tt);
    % pnf(50, 50) = pnf(100, 100) + sn(tt);

    %plot
    clf;
    subplot(2,1,1);
    imagesc(x, y, pn'); colorbar; clim([-0.35 0.35])
    title(sprintf('t = %d', tt));
    subplot(2,1,2);
    mesh(x,y,pn'); colorbar; clim([-0.35 0.35])
    axis([0 Lx 0 Ly -0.35 0.35]);
    shg; pause(0.001);
    
    % incrementa o t (passos devem ser inteiros)
    tt = tt + 1;
end