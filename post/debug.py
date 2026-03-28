import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("debug.csv")

plt.figure()
plt.plot(data['WallPointX'], data['WallPointY'], 'o')
plt.quiver(data['WallPointX'], data['WallPointY'], data['NormalX'], data['NormalY'], angles='xy', scale_units='xy', scale=1, color='r')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.show()