import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

data = pd.read_csv("w_time.csv", header=None)
data.index = np.linspace(0.01, 30, 40)
id = data[data[0] != 0.0].index
data.columns = ["Момент времени"]

fig = plt.figure(figsize=(19, int(10)))
plt.axvline(x = id[-1], color='r', label='Конец превалирования перехода А')
sns.lineplot(data=data)
plt.title("Момент времени начала превалирования перехода А")
plt.xlabel("Разница между частотами")
plt.ylabel("Время")
plt.grid()
plt.legend()
plt.savefig("w_time.png")

data = pd.read_csv("w_prob_diff.csv", header=None)
data.index = np.linspace(0.01, 30, 40)
data.columns = ["Разница вероятностей"]

fig = plt.figure(figsize=(19, int(10)))
plt.axvline(x = id[-1], color='r', label='Конец превалирования перехода А')
sns.lineplot(data=data)
plt.title("Разница вероятностей между переходами А и В")
plt.xlabel("Разница между частотами")
plt.ylabel("Вероятность")
plt.grid()
plt.legend()
plt.savefig("w_prob_diff.png")