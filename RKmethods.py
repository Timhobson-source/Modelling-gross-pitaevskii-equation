from nodepy import runge_kutta_method as rk
import numpy as np
rk44=rk.loadRKM("RK44") #a famous RK method

#Make my own RK method
A = np.array([[0,0],[0.5,0]])
b=np.array([0,1])
rk22=rk.ExplicitRungeKuttaMethod(A,b) #used Explicit here for reasons...

# can get c-vector via rk22.c (or rkmethodname.c generically!)

# can get orders with the order attribute
print(rk22.order()) #has order 2 !
print(rk44.order()) #has order 4 !

