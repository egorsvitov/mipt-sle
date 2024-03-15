import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

ns = []
simple_norms = []
jacobi_norms = []
gauss_seidel_norms = []
cheb_norms = []

ns.append(0)
simple_norms.append(np.sqrt(3))
jacobi_norms.append(np.sqrt(3))
gauss_seidel_norms.append(np.sqrt(3))
cheb_norms.append(np.sqrt(3))

with open('by_iterations.txt', 'r') as f:
    for line in f.readlines():
        if (line != '\n'):
            [_ns, _simple_norms, _jacobi_norms, _gauss_seidel_norms, _cheb_norms] = [float(i) for i in line.split()]
            ns.append(_ns)
            simple_norms.append(_simple_norms)
            jacobi_norms.append(_jacobi_norms)
            gauss_seidel_norms.append(_gauss_seidel_norms)
            cheb_norms.append(_cheb_norms)

fig, ax = plt.subplots(figsize=(8, 5), dpi=400)
# ax.axis(xmin = 0, xmax = max(t), ymin = 0, ymax=max(u))
ax.set_xlabel('Количество итераций')
ax.set_ylabel('Норма вектора результата')
# ax.set_title('Зависимость напряжения от времени')
ax.plot(ns, simple_norms, label='simple')
ax.plot(ns, jacobi_norms, label='jacobi')
ax.plot(ns, gauss_seidel_norms, label='gauss_seidel')
ax.plot(ns, cheb_norms, label='chebyshev_simple')
#ax.minorticks_on()
ax.grid(which='major')
#ax.grid(which='minor', linestyle = ':')
ax.legend()
ax.set_xticks(ns)
# ax.text(80, 1.2, 'Время зарядки конденсатора 138 с')
# ax.text(80, 1, 'Время разрядки конденсатора 122 с')
fig.savefig('by_iterations.png')
# pyplot.show()

tolerances = []
simple_times = []
jacobi_times = []
gauss_seidel_times = []
cheb_times = []

with open('by_time.txt', 'r') as f:
    for line in f.readlines():
        if (line != '\n'):
            [_tolerances, _simple_times, _jacobi_times, _gauss_seidel_times, _cheb_times] = [float(i) for i in line.split()]
            tolerances.append(_tolerances)
            simple_times.append(_simple_times)
            jacobi_times.append(_jacobi_times)
            gauss_seidel_times.append(_gauss_seidel_times)
            cheb_times.append(_cheb_times)

fig, ax = plt.subplots(figsize=(8, 5), dpi=400)
# ax.axis(xmin = 0, xmax = max(t), ymin = 0, ymax=max(u))
ax.set_xlabel('tolerance')
ax.set_ylabel('time')
# ax.set_title('Зависимость напряжения от времени')
ax.plot(tolerances, simple_times, label='simple')
ax.plot(tolerances, jacobi_times, label='jacobi')
ax.plot(tolerances, gauss_seidel_times, label='gauss_seidel')
ax.plot(tolerances, cheb_times, label='chebyshev_simple')
ax.minorticks_on()
ax.grid(which='major')
ax.grid(which='minor', linestyle = ':')
ax.legend()
# ax.text(80, 1.2, 'Время зарядки конденсатора 138 с')
# ax.text(80, 1, 'Время разрядки конденсатора 122 с')
fig.savefig('by_time.png')
# pyplot.show()

ns, tolerances, simple_times, jacobi_times, gauss_seidel_times, cheb_times = np.array(ns), np.array(tolerances), np.array(simple_times), np.array(jacobi_times), np.array(gauss_seidel_times), np.array(cheb_times)
# ns = np.log(ns)
tolerances = np.log(tolerances)
simple_times, jacobi_times, gauss_seidel_times = (simple_times- simple_times.min() + 1), (jacobi_times - jacobi_times.min() + 1), (gauss_seidel_times - gauss_seidel_times.min() + 1)
simple_times = np.log(simple_times)
jacobi_times = np.log(jacobi_times)
gauss_seidel_times = np.log(gauss_seidel_times)
cheb_times = np.log(cheb_times)

fig, ax = plt.subplots(figsize=(8, 5), dpi=400)
# ax.axis(xmin = 0, xmax = max(t), ymin = 0, ymax=max(u))
ax.set_xlabel('$\ln(tolerance)$')
ax.set_ylabel('$\ln(time)$')
# ax.set_title('Зависимость напряжения от времени')
ax.plot(tolerances, simple_times, 'o-', label='simple')
ax.plot(tolerances, jacobi_times, 'o-', label='jacobi')
ax.plot(tolerances, gauss_seidel_times, 'o-', label='gauss_seidel')
ax.plot(tolerances, cheb_times, 'o-', label='chebyshev_simple')
ax.minorticks_on()
ax.grid(which='major')
ax.grid(which='minor', linestyle = ':')
ax.legend()
# ax.text(80, 1.2, 'Время зарядки конденсатора 138 с')
# ax.text(80, 1, 'Время разрядки конденсатора 122 с')
fig.savefig('ln_by_time.png')
# pyplot.show()