
# 增加浮游植物死亡速率和 detritus 转化的参数
mortality_rate = 0.01  # 细胞死亡速率常数 (s⁻¹)，可以根据实际情况调整
degradation_rate = 0.001  # Detritus 中大分子物质降解速率 (s⁻¹)
detritus = {'carbohydrate': 0, 'protein': 0, 'lipid': 0}  # Detritus 中大分子物质初始组成

# 更新浮游植物生长模型，考虑细胞死亡并转化为 detritus 的过程
def update_phytoplankton_growth(phytoplankton, time_step):
    # 模拟浮游植物生长（原模型中的生长逻辑）
    growth = phytoplankton['growth_rate'] * phytoplankton['biomass'] * time_step
    
    # 引入细胞死亡
    mortality = mortality_rate * phytoplankton['biomass'] * time_step
    phytoplankton['biomass'] -= mortality  # 减去死亡细胞量
    
    # 死亡细胞形成 detritus
    detritus['carbohydrate'] += mortality * phytoplankton['carbohydrate_ratio']
    detritus['protein'] += mortality * phytoplankton['protein_ratio']
    detritus['lipid'] += mortality * phytoplankton['lipid_ratio']
    
    # 更新浮游植物生物量
    phytoplankton['biomass'] += growth - mortality
    return phytoplankton, detritus

# 定义 detritus 随时间降解的函数
def update_detritus_composition(detritus, time_step):
    # Detritus 中大分子物质的降解过程
    for key in detritus:
        detritus[key] *= (1 - degradation_rate * time_step)  # 一阶衰减模型
    return detritus

# 示例使用该更新函数
phytoplankton = {'biomass': 1e-6, 'growth_rate': 0.05, 'carbohydrate_ratio': 0.5, 'protein_ratio': 0.3, 'lipid_ratio': 0.2}
time_step = 3600  # 以秒为单位（例如1小时）

# 更新浮游植物生长并计算 detritus 形成
phytoplankton, detritus = update_phytoplankton_growth(phytoplankton, time_step)

# 更新 detritus 随时间的降解
detritus = update_detritus_composition(detritus, time_step)

print("Updated Phytoplankton Biomass:", phytoplankton['biomass'])
print("Detritus Composition after 1 hour:", detritus)
