import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

n = []
t_dense = []
t_csr = []
p = []

with open('output.txt', 'r') as f:
    for line in f.readlines():
        if (line != '\n'):
            [pp, nn, tt_dense, tt_csr] = [float(i) for i in line.split()]
            p.append(pp)
            n.append(nn)
            t_dense.append(tt_dense)
            t_csr.append(tt_csr)
            
fig, ax = plt.subplots(figsize=(8, 5), dpi=400)
ax.axis(xmin = 0, xmax = 1100, ymin = 0, ymax=1*10**9)
ax.set_xlabel('n ($n^2$ элементов)')
ax.set_ylabel('Время, нс')
ax.set_title('Сравнение скоростей матрично-векторного умножения')
dense, = ax.plot(n[0:10], t_dense[0:10], 'o', label='dense')
csr, = ax.plot(n[0:10], t_csr[0:10], 'o', label='csr')
ax.minorticks_on()
ax.grid(which='major')
ax.grid(which='minor', linestyle = ':')
text = ax.text(0.3, 0.7, f'density={p[0]}', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
ax.legend(loc="upper left")
def animate(i):
    #ax.clear()
    dense.set_ydata(t_dense[10*i:10*i+10])
    csr.set_ydata(t_csr[10*i:10*i+10])
    text.set_text(f'density={round(p[10*i]*100, 3)}%')
    ax.set_ylim(0, t_csr[10*i+9]*1.1)
    ax.autoscale()
    #fig.savefig(f'pics/graph_{i}.png')
    print(i)
    return (dense, csr)

ani = animation.FuncAnimation(fig, animate, repeat=True,
                                    frames=len(n)//10, interval=10)
writer = animation.PillowWriter(fps=5,
                                metadata=dict(artist='Me'),
                                bitrate=1800)
ani.save('fixed_scale.gif', writer=writer)