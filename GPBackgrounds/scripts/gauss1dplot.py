from matplotlib import pyplot as plt
import numpy as np

def gaussian(x, mu, sig, fact):
    return fact * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x = np.linspace(-3,3,200)
#for mu, sig, fact, color in [(0, 0.3, 1),(-1,0.5, 0.75)]:
#    plt.plot(x, gaussian(x, mu, sig, fact))

plt.plot(x,gaussian(x, 0, 0.3, 1), color='purple')
plt.plot(x,gaussian(x, -1, 0.5, 0.75), color='green')
plt.plot(x,gaussian(x, -2, 0.7, 0.5), color='yellow')

plt.savefig('Gauss1D.pdf')
#plt.show()
