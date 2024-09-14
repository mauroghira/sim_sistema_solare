import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import RK23

#non va bene perché posso fre soluzione soloo in 1d
def gravity(x, t, M, m, G):
	r1 = np.array(x[0:2])
	r2 = np.array(x[4:6])
	d = r1-r2
	app = -G/(np.linalg.norm(d))**3
	x1dot = x[2] #sostituzione nel sistema lineare
	y1dot = x[3]
	v1dot = app*M*d #accelerazione
	x2dot = x[6] #sostituzione nel sistema lineare
	y2dot = x[7]
	v2dot = app*m*d
	return x1dot, y1dot, v1dot[0], v1dot[1], x2dot, y2dot, v2dot[0], v2dot[1]
	
def traie(x, t):
	M = 1.98e30
	G = 6.67e-11
	m = 5.69e26
	r1 = np.array(x[0:2])
	r2 = np.array(x[4:6])
	d = r1-r2
	app = -G/(np.linalg.norm(d))**3
	x1dot = x[2] #sostituzione nel sistema lineare
	y1dot = x[3]
	v1dot = app*M*d #accelerazione
	x2dot = x[6] #sostituzione nel sistema lineare
	y2dot = x[7]
	v2dot = app*m*d
	return x1dot, y1dot, v1dot[0], v1dot[1], x2dot, y2dot, v2dot[0], v2dot[1]

halfpi, pi, twopi = [f*np.pi for f in (0.5, 1, 2)]

M = 1.98e30
G = 6.67e-11
m = 5.69e26

e = 0.0560
a = 1427.4e9
mu = G*(M+m)
r_peri, r_apo = a*(1.-e), a*(1.+e)
v_peri, v_apo = [np.sqrt(2./r - 1./a) for r in (r_peri, r_apo)]
T = twopi * np.sqrt(a**3/mu)

print(T/(3600*24), 10765.7)
print(r_peri, 1352.5E9)
print(r_apo, 1514.5E9)

X0 = np.array([r_peri, 0, 0, v_peri, 0,0,0,0])
#X0 = np.array([-r_apo, 0, 0, -v_apo])
times = np.linspace(0, T, 100000)

answer, info = odeint(gravity, X0, times, args=(M, m, G), full_output=True)

print(answer)

x1 = np.array([i[0] for i in answer])
x2 = np.array([i[2] for i in answer])
print(x1)
y1 = np.array([i[1] for i in answer])
y2 = np.array([i[3] for i in answer])

"""
r1, r2 = answer[:,0], answer[;,2]
print(r1,r2)
x1, x2 = np.array(r1[:,0]), np.array(r2[:,0])
y1, y2 = np.array(r1[:,1]). np.array(r1[:,1])

solve = RK23(traie, 0, X0, T)
answer = []
for t in range(0,T,100000):
	solve.steo()
	answer.append(solve.y)

r1, r2 = np.array([i[0] for i in answer]), np.array([i[2] for i in answer])
print(r1,r2)
x1, x2 = np.array([i[0] for i in r1]), np.array([i[0] for i in r2])
x1, x2 = np.array([i[1] for i in r1]), np.array([i[1] for i in r2])
"""

theta = np.arctan2(y1, x1)
E = 2. * np.arctan(np.sqrt((1.-e)/(1.+e)) * np.tan(theta/2))
t = a * np.sqrt(a/mu) * (E - e * np.sin(E))

if True:
    plt.figure()
    #plt.subplot(2, 1, 1)
    plt.plot(x1, y1, label="Pianeta")
    plt.plot(x2, y2, label="Sole")
    #plt.plot([0], [0], 'ok')
    #plt.gca().set_aspect('equal')
    plt.title('y vs. x numerical')
    plt.legend()
    plt.plot()
    """
    plt.subplot(2, 1, 2)
    plt.plot(times[1:-1], x)
    plt.plot(times[1:-1], y)
    plt.xlim(-pi, pi)
    plt.title('x(t) and y(t) numerical')
    plt.show()

    plt.subplot(2, 2, 1)
    plt.title('theta(t_numerical)')
    plt.plot(times[1:-1], theta)
    plt.xlim(-pi, pi)
    plt.ylim(-pi, pi)
    plt.gca().set_aspect('equal')
    plt.subplot(2, 2, 2)
    plt.title('E_analytic(theta_numerical)')
    plt.plot(E, theta)
    plt.xlim(-pi, pi)
    plt.ylim(-pi, pi)
    plt.gca().set_aspect('equal')
    plt.subplot(2, 2, 3)
    plt.title('theta(t_analytic)')
    plt.plot(t, theta)
    plt.xlim(-pi, pi)
    plt.ylim(-pi, pi)
    plt.gca().set_aspect('equal')
    """
    #plt.subplot(2, 2, 4)
    plt.title('t_analytic(t_numerical)')
    plt.plot(t, times)
    plt.xlim(-pi, pi)
    plt.ylim(-pi, pi)
    #plt.gca().set_aspect('equal')
    
    plt.show()
