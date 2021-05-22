import matplotlib.pyplot as plt
import numpy
from math import cos
from math import sin
from math import pi

# T1'(t) = f1(t,T1,T2,p1,p2)
# P1'(t) = g1(t,T1,T2,p1,p2)
# T2'(t) = f2(t,T1,T2,p1,p2)
# P2'(t) = g2(t,T1,T2,p1,p2)
# (x0,p0)=(x(t0),p(t0))

tin = 0 #a
dt = 0.01
tmax = 1
lar = int(tmax / dt)
cu1 = 2.0
cu2 = 2.0
m = 0.1
g = 9.81


def denom(t1, t2):
    return 1+(sin(t1-t2))**2


def f1(ti, t1, t2, p1, p2):
    return (cu2**2*p1-cu1*cu2*p2*cos(t1-t2))/(cu1**2 * cu2**2 * m * denom(t1, t2))


def g1(ti, t1, t2, p1, p2):
    return -(p1*p2*sin(t1-t2))/(cu1 * cu2 * m * denom(t1, t2)) \
           + ((cu2**2 * p1**2 + 2 * cu1**2 * p2**2 -
            2 * cu1*cu2*p1*p2*cos(t1-t2))/(2*cu1**2 * cu2**2 * denom(t1, t2)**2))*sin(2*t1-2*t2)\
           - 2 * m * g*cu1*sin(t1)


def f2(ti, t1, t2, p1, p2):
    return (2 * cu1**2 * p2 -
            cu1 * cu2 * p1 * cos(t1-t2))/(2 * cu1**2 * cu2**2 * m * (1+(sin(t1-t2))**2))


def g2(ti, t1, t2, p1, p2):
    return p1*p2*sin(t1-t2)/(cu1*cu2*m*(1+(sin(t1-t2))**2)) \
           - ((cu2**2 * p1**2 + 2 * cu1**2 * p2**2 -
            2 * cu1*cu2*p1*p2*cos(t1-t2))/(2 * cu1**2 * cu2**2
                                           * (1+(sin(t1-t2))**2)**2))*sin(2*(t1-t2)) \
           - m*cu2*g*sin(t2)


A = [[0, 0.5], [0, 1], [0, 1.5], [pi/4, 0], [pi/2, 0], [0, -0.5], [0, -1], [0, -1.5], [-pi/4, 0], [-pi/2, 0]]
# A = [[1.5, 0]]
print(len(A))
for j in range(0, len(A)):
    T1 = numpy.zeros((lar,), dtype=float)
    P1 = numpy.zeros((lar,), dtype=float)
    T2 = numpy.zeros((lar,), dtype=float)
    P2 = numpy.zeros((lar,), dtype=float)
    Td1 = numpy.zeros((lar,), dtype=float)
    Td2 = numpy.zeros((lar,), dtype=float)
    t = numpy.linspace(tin, tmax, lar)
    T1[0] = A[j][0]
    P1[0] = A[j][1]
    T2[0] = A[j][0]
    P2[0] = A[j][1]
    Td1[0] = f1(t[0], T1[0], T2[0], P1[0], P2[0])


    for i in range(0, lar - 1):
        k1 = f1(t[i], T1[i], T2[i], P1[i], P2[i])
        l1 = g1(t[i], T1[i], T2[i], P1[i], P2[i])
        o1 = f2(t[i], T1[i], T2[i], P1[i], P2[i])
        r1 = g2(t[i], T1[i], T2[i], P1[i], P2[i])

        k2 = f1(t[i] + dt / 2, T1[i] + k1 / 2, T2[i] + l1 / 2, P1[i] + o1 / 2, P2[i] + r1 / 2)
        l2 = g1(t[i] + dt / 2, T1[i] + k1 / 2, T2[i] + l1 / 2, P1[i] + o1 / 2, P2[i] + r1 / 2)
        o2 = f2(t[i] + dt / 2, T1[i] + k1 / 2, T2[i] + l1 / 2, P1[i] + o1 / 2, P2[i] + r1 / 2)
        r2 = g2(t[i] + dt / 2, T1[i] + k1 / 2, T2[i] + l1 / 2, P1[i] + o1 / 2, P2[i] + r1 / 2)

        k3 = f1(t[i] + dt / 2, T1[i] + k2 / 2, T2[i] + l2 / 2, P1[i] + o2 / 2, P2[i] + r2 / 2)
        l3 = g1(t[i] + dt / 2, T1[i] + k2 / 2, T2[i] + l2 / 2, P1[i] + o2 / 2, P2[i] + r2 / 2)
        o3 = f2(t[i] + dt / 2, T1[i] + k2 / 2, T2[i] + l2 / 2, P1[i] + o2 / 2, P2[i] + r2 / 2)
        r3 = g2(t[i] + dt / 2, T1[i] + k2 / 2, T2[i] + l2 / 2, P1[i] + o2 / 2, P2[i] + r2 / 2)

        k4 = f1(t[i] + dt, T1[i] + k3, T2[i] + l3, P1[i] + o3, P2[i] + r3)
        l4 = g1(t[i] + dt, T1[i] + k3, T2[i] + l3, P1[i] + o3, P2[i] + r3)
        o4 = f2(t[i] + dt, T1[i] + k3, T2[i] + l3, P1[i] + o3, P2[i] + r3)
        r4 = g2(t[i] + dt, T1[i] + k3, T2[i] + l3, P1[i] + o3, P2[i] + r3)

        T1[i + 1] = T1[i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        P1[i + 1] = P1[i] + dt * (l1 + 2 * l2 + 2 * l3 + l4) / 6
        T2[i + 1] = T2[i] + dt * (o1 + 2 * o2 + 2 * o3 + o4) / 6
        P2[i + 1] = P2[i] + dt * (r1 + 2 * r2 + 2 * r3 + r4) / 6
        Td1[i + 1] = f1(t[i + 1], T1[i + 1], T2[i + 1], P1[i + 1], P2[i + 1])
        Td2[i + 1] = f1(t[i + 1], T1[i + 1], T2[i + 1], P1[i + 1], P2[i + 1])
    print(j + 1)
    plt.figure(1)
    plt.plot(T1, Td1, 'r')
    plt.figure(2)
    plt.plot(T2, Td2, 'b')

plt.figure(1)
plt.xlabel('T_1(t)')
plt.ylabel('Td_1(t)')
plt.grid(b=None, which='major', axis='both')
plt.figure(2)
plt.xlabel('T_2(t)')
plt.ylabel('Td_2(t)')
plt.grid(b=None, which='major', axis='both')
plt.show()
