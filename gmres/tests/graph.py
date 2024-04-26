import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

tolerances = []
simple_times = []
jacobi_times = []
gauss_seidel_times = []
cheb_times = []
new_cheb_times = []
fastest_grad_times = []
sym_GS_times = []
conj_times = []
gmres_times = []

with open('by_time.txt', 'r') as f:
    for line in f.readlines():
        if (line != '\n'):
            [_tolerances, _simple_times, _jacobi_times, _gauss_seidel_times, _cheb_times, _new_cheb_times, _fastest_grad_times, _sym_GS_times, _conj_times, _gmres_times] = [float(i) for i in line.split()]
            tolerances.append(_tolerances)
            simple_times.append(_simple_times)
            jacobi_times.append(_jacobi_times)
            gauss_seidel_times.append(_gauss_seidel_times)
            cheb_times.append(_cheb_times)
            new_cheb_times.append(_new_cheb_times)
            fastest_grad_times.append(_fastest_grad_times)
            sym_GS_times.append(_sym_GS_times)
            conj_times.append(_conj_times)
            gmres_times.append(_gmres_times)

fig, ax = plt.subplots(figsize=(8, 5), dpi=400)
# ax.axis(xmin = 0, xmax = max(t), ymin = 0, ymax=max(u))
ax.set_xlabel('tolerance')
ax.set_ylabel('time')
# ax.set_title('Зависимость напряжения от времени')
ax.plot(tolerances, simple_times, 'o-', label='simple')
ax.plot(tolerances, jacobi_times, 'o-', label='jacobi')
ax.plot(tolerances, gauss_seidel_times, 'o-', label='gauss_seidel')
ax.plot(tolerances, cheb_times, 'o-', label='chebyshev_simple')
ax.plot(tolerances, new_cheb_times, 'o-', label='new_chebyshev')
ax.plot(tolerances, fastest_grad_times, 'o-', label='fastest_gradient_descent')
ax.plot(tolerances, sym_GS_times, 'o-', label='sym_GS')
ax.plot(tolerances, conj_times, 'o-', label='conjugate_gradient')
ax.plot(tolerances, gmres_times, 'o-', label='gmres')
ax.minorticks_on()
ax.grid(which='major')
ax.grid(which='minor', linestyle = ':')
ax.legend()
# ax.text(80, 1.2, 'Время зарядки конденсатора 138 с')
# ax.text(80, 1, 'Время разрядки конденсатора 122 с')
fig.savefig('by_time.png')
# pyplot.show()

ns, tolerances, simple_times, jacobi_times, gauss_seidel_times = np.array(ns), np.array(tolerances), np.array(simple_times), np.array(jacobi_times), np.array(gauss_seidel_times)
# ns = np.log(ns)
tolerances = np.log(tolerances)
simple_times, jacobi_times, gauss_seidel_times = (simple_times- simple_times.min() + 1), (jacobi_times - jacobi_times.min() + 1), (gauss_seidel_times - gauss_seidel_times.min() + 1)
simple_times = np.log(simple_times)
jacobi_times = np.log(jacobi_times)
gauss_seidel_times = np.log(gauss_seidel_times)

fig, ax = plt.subplots(figsize=(8, 5), dpi=400)
# ax.axis(xmin = 0, xmax = max(t), ymin = 0, ymax=max(u))
ax.set_xlabel('$\ln(tolerance)$')
ax.set_ylabel('$\ln(time)$')
# ax.set_title('Зависимость напряжения от времени')
ax.plot(tolerances, simple_times, 'o-', label='simple')
ax.plot(tolerances, jacobi_times, 'o-', label='jacobi')
ax.plot(tolerances, gauss_seidel_times, 'o-', label='gauss_seidel')
ax.plot(tolerances, cheb_times, 'o-', label='chebyshev_simple')
ax.plot(tolerances, new_cheb_times, 'o-', label='new_chebyshev')
ax.minorticks_on()
ax.grid(which='major')
ax.grid(which='minor', linestyle = ':')
ax.legend()
# ax.text(80, 1.2, 'Время зарядки конденсатора 138 с')
# ax.text(80, 1, 'Время разрядки конденсатора 122 с')
fig.savefig('ln_by_time.png')
# pyplot.show()