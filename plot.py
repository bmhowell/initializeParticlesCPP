import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("particlePosition.dat")
data = np.asarray(df.values)
xpos = data[:, 0]
ypos = data[:, 1]
diam = data[:, 2]
nparticle = np.arange(0, len(data[:, 0]))

xmin = 3.98571e-5
xmax = 0.000797143
ymin = 6.2e-5
ymax = 0.001674

print(len(xpos))

# plot particle locations
plt.figure(figsize=(5, 5))
plt.scatter(xpos, ypos)
plt.plot([xmin, xmin], [ymin, ymax], 'k')
plt.plot([xmin, xmax], [ymin, ymin], 'k')
plt.plot([xmax, xmax], [ymin, ymax], 'k')
plt.plot([xmin, xmax], [ymax, ymax], 'k')
plt.show()

# plot profile particle diameters
plt.figure()
plt.plot(nparticle, diam)
plt.show()

