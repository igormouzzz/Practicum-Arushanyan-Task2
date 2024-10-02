import numpy as np
from math import exp, sin, cos
import matplotlib.pyplot as plt

with open('out.txt', 'r') as f:
    x,y = list(map(float, f.readline().strip().split())), list(map(float, f.readline().strip().split()))
    #x2,y2 = list(map(float, f.readline().strip().split())), list(map(float, f.readline().strip().split()))
x_real = np.linspace(0,1,100)
#y_real = [2*xx -1 for xx in x_real]        
#y_real = [2*xx for xx in x_real]
y_real = [xx + exp(-xx) for xx in x_real]   #11
#y_real = [1.0]*len(x_real)                 #12
#y_real = [1 + 4/9*xx for xx in x_real]     #17
#y_real = [xx for xx in x_real]             #18
plt.plot(x_real, y_real, label = 'Реальное решение',color = 'green')
plt.scatter(x,y, color = 'red')
#plt.scatter(x,y2, color = 'red')
plt.plot(x,y, label = 'Численное решение', color = 'blue')
#plt.plot(x2,y2, label = 'Численное решение 2', color = 'purple')

plt.title('Графики функций')  # Заголовок графика
plt.xlabel('X')  # Подпись оси X
plt.ylabel('Y')  # Подпись оси Y
plt.grid()
plt.legend()
plt.show()
