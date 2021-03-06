__author__ = 'fiodar'

from interpolation import *
import matplotlib.pyplot as plt


nodes = [-5., -4.69965, -3.64222, -1.33749, -1.0145, 2.46799, 3.]
y = lambda x: x * np.sin(x)
data = (nodes, y(nodes))
f_newton = newton_p(data)
f_spline = cubic_spline(data)
X = np.linspace(-6, 5, 300)
inner = np.random.uniform(nodes[0], nodes[-1])
outer = np.random.uniform(nodes[-1], nodes[-1] + 2)
fig1 = plt.figure('Lab2', figsize=plt.figaspect(0.5))
ax1, ax2 = fig1.add_subplot(121), fig1.add_subplot(122)
ax1.set_title('Newton polynom')
ax2.set_title('Cubic spline')
ax1.plot(X, y(X), 'r-', label='function')
ax2.plot(X, y(X), 'r-', label='function')
ax1.plot(X, f_newton(X), label='interpolant')
ax2.plot(X, f_spline(X), 'g-', label='interpolant')
ax1.plot(nodes, y(nodes), 'ro')
ax2.plot(nodes, y(nodes), 'ro')
ax1.text(-5.9, 1, r'$\varphi(%.1f)=%.1f$' % (inner, f_newton(inner)))
ax1.text(-5.9, 0.5, r'$\varphi(%.1f)=%.1f$' % (outer, f_newton(outer)))
ax2.text(-5.9, 1, r'$\varphi(%.1f)=%.1f$' % (inner, f_spline(inner)))
ax2.text(-5.9, 0.5, r'$\varphi(%.1f)=%.1f$' % (outer, f_spline(outer)))
ax1.grid()
ax2.grid()
ax1.legend(loc='best')
ax2.legend(loc='best')
ax1.set_ylim(-6, 3)
ax2.set_ylim(-6, 3)
fig1.tight_layout()
plt.show()
