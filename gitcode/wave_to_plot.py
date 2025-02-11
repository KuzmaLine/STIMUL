import matplotlib.pyplot as plt
#import plotly.graph_objects as go
import pandas as pd
import seaborn as sns
import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)
w = np.linspace(0.001, 0.1, 32)
g = np.linspace(0.001, 0.1, 32)
data = pd.read_csv("wave_time.csv", header=None).to_numpy()
#fig = go.figure(data=[go.Surface(z=data.values, x=g, y=w)])
#fig.show()

ggrid, wgrid = np.meshgrid(g, w)
print(data)

fig = plt.figure(figsize=(19, int(10)))
#plt.axvline(x = id[-1], color='r', label='Конец превалирования перехода А')
axes = fig.add_subplot(projection='3d')
axes.plot_surface(ggrid, wgrid, data, cmap='Spectral', alpha=1)
axes.contour(ggrid, wgrid, data, zdir='x', offset=0.1, cmap='coolwarm')
axes.contour(ggrid, wgrid, data, zdir='y', offset=0.1, cmap='coolwarm')
axes.set_xlabel("Интенсивность утечки")
axes.set_ylabel("Интенсивность волновода")
axes.set_zlabel("Момент времени")
plt.title("Момент времени начала превалирования перехода А")
plt.grid()
plt.legend()
plt.savefig("wave_time.png")
plt.show()

data = pd.read_csv("wave_prob_diff.csv", header=None)

fig = plt.figure(figsize=(19, int(10)))
axes = fig.add_subplot(projection='3d')
axes.plot_surface(ggrid, wgrid, data, cmap='Spectral')
axes.contour(ggrid, wgrid, data, zdir='x', offset=0, cmap='coolwarm')
axes.contour(ggrid, wgrid, data, zdir='y', offset=0, cmap='coolwarm')
axes.set_xlabel("Интенсивность утечки")
axes.set_ylabel("Интенсивность волновода")
axes.set_zlabel("Разница вероятностей")
plt.title("Разница вероятностей между переходом А и B")
plt.grid()
plt.legend()
plt.savefig("wave_prob_diff.png")
plt.show()