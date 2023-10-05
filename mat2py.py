##import

import numpy as np
from numpy import ndarray
from math import *
import matplotlib.pyplot as plt
from scipy import ndimage, datasets, signal
from tqdm import tqdm

##variaveis

Lx = 1.0 - 0.1
Ly = 0.8 - 0.1 # comprimento de propagação em x e y em [m]
nx = 91
ny = 71
dx = Lx/nx
dy = Ly/ny
# dx= 10e-5 #tamanho do passo em [m / pixel]
# dy=dx #superfície regular
# nx= int(Lx/dx) #n total de pontos em [pixel]
# ny=nx





x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)  #cria um vetor de pontos discretos para x e y comecando em 0 indo ate L com um total de n pontos
T = 3 # 10e-6  #tempo total em [s]
nt = 4301  # #n total de amostras temporais ou timesteps
dt = (T/nt)  # #discretização temporal em [s / amostra] ou [s / timestep]

 #variaveis de campo
p = (ny, nx)
pn = np.zeros(p)   #gera um grid para armazenar a pressao no tempo
pnf = pn   # p(n+1)
pnb = pn   # p(n-1)

 #parametros
c = 2  # 5500     # velocidade de propagação no meio em [m/s]
CFL = c * dt / dx    #CFL=0.5    #CFL = c.dt/dx
#dt=CFL*dx/c


  ##fonte
Tt= np.linspace(0, T-dt, nt)      # como começa em zero tem que descontar dt do total do vetor de tempo
# freq = 5e6      # frequencia em [Hz]
# delay = 5e-7      # atraso do sinal em [s], se 0 (zero) o pico no pulso é no t = 0
# bandwidth = 0.82      # largura de banda em [  #]
# sn= gauspuls(Tt,1000,0.5, -6)
#sn= signal.gausspulse(Tt-delay, freq, bandwidth)

t= Tt-T/20
a = 1e-3
def ricker(t, a):
    #x = np.diff(np.diff(np.exp(-t**2/a))/dt)/dt
    x = np.exp(-t**2/a)*(4*t**2/a**2-2/a)
    return x/np.max(np.abs(x))

sn = ricker(t,a)

# sn = signal.ricker(nt,100.0)
# sn[10:7100] = sn[11670:18670]
# sn[11670:18670] = 0

plt.figure()
plt.plot(sn)
plt.show()



 #shots
nr = np.arange(nx) ##n do receptor

pos_rec_x = nr.copy()
pos_rec_y = np.zeros_like(nr) ##matriz posição do receptor

#pos_rec = pos_rec.astype(int)

act = np.zeros((nt, nx))  # matriz das aquisicoes no tempo

for tt in tqdm(range(nt)):
    #   troca as variáveis
    pnb = pn
    pn = pnf #   deixa salvo  p(n - 1), p(n)  e p(n + 1)

    gaussian= ndimage.laplace(pnf)

    pnf = (CFL**2) * gaussian
    pnf += 2 * pn - pnb


    #condicoes de fronteira(reflexao)
    # pn[:,ny-1]=0;
    # pn[nx-1,:]=0;

    # #condicoes de frontira(absorcao)
    pnf[0,:]=pn[1,:] + ((CFL-1)/(CFL+1))*(pnf[1,:]-pn[0,:])
    pnf[ny-1, :]=pn[ny-2,:] + ((CFL-1)/(CFL+1))*(pnf[ny-2,:]-pn[ny-1,:])
    pnf[:,0]=pn[:,1] + ((CFL-1)/(CFL+1))*(pnf[:,1]-pn[:,0])
    pnf[:,nx-1]=pn[:,nx-2] + ((CFL-1)/(CFL+1))*(pnf[:,nx-2]-pn[:,nx-1])

    #source

    # pnf[int(ny / 2), int(nx / 2)] = pnf[int(ny / 2), int(nx / 2)] + sn[tt]
    pnf[10 , 45] +=  sn[tt]
    #pnf[ny//2, nx//2] += sn[tt]

    #aquisições
    # for n in range(nx):
    #     act[tt,n] = pnf[pos_rec_y[n], pos_rec_x[n]]
    act[tt, :] = pnf[pos_rec_y, pos_rec_x]


##plot
# ax = plt.axes(projection='3d')
# X, Y = np.meshgrid(x, y)
# ax.plot_surface(X, Y, pnf, cmap='plasma')
# plt.show()
# plt.figure(1)
# plt.imshow(pnf)
# plt.colorbar()
# plt.show()
# plt.pause(1e-12)
# plt.close(fig=1)



##plot shot
plt.figure()
plt.imshow(act, origin='upper', aspect='auto')
plt.colorbar()
plt.show()





print("FIM!")