import numpy as np
import matplotlib.pyplot as plt


def density_dependent_growth_with_detritus(X0, mu, mr, K, t):
    """
    模拟密度依赖死亡率的细胞丰度变化，并考虑碎屑部分
    X0: 初始细胞丰度
    mu: 最大生长速率
    mr: 基准死亡率
    K: 环境承载力
    t: 时间数组
    """
    X = np.zeros_like(t)  # 细胞丰度数组
    D = np.zeros_like(t)  # 碎屑丰度数组
    X[0] = X0  # 设置初始细胞丰度

    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]

        # calculate the growth rate and mortality rate 计算生长速率和密度依赖的死亡速率
        growth_rate = mu * X[i - 1]
        # death_rate = mr * X[i - 1] * (1 + X[i - 1] / K)
        death_rate = mr * X[i - 1] ** 2

        # update the cell abundance 更新细胞丰度
        X[i] = X[i - 1] + (growth_rate - death_rate) * dt

        # 确保丰度不为负数
        if X[i] < 0:
            X[i] = 0

        # 碎屑部分由死亡部分积累
        D[i] = D[i - 1] + death_rate * dt  # 碎屑由死亡部分积累

    return X, D


# 参数设置
X0 = 10  # 初始细胞丰度
mu = 0.6  # 最大生长速率
mr = 0.02  # 基准死亡率
K = 1000  # 环境承载力
t = np.linspace(0, 100, 500)  # 时间数组

# 模拟细胞丰度和碎屑变化
X, D = density_dependent_growth_with_detritus(X0, mu, mr, K, t)

# 绘制细胞丰度和碎屑的变化
plt.plot(t, X, label="Cell Abundance")
plt.plot(t, D, label="Detritus (Dead Cells)")
# 设置 x 和 y 轴的范围
plt.xlim(0, 100)  # 设置 x 轴的范围为 2 到 8
plt.ylim(0, 100000)  # 设置 y 轴的范围为 -0.5 到 0.5
plt.xlabel("Time")
plt.ylabel("Abundance (cells/L)")
plt.title("Cell Growth and Detritus with Density-dependent Mortality")
plt.legend()
plt.show()
