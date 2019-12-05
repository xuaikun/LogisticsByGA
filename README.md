# LogisticsByGA
利用GA算法解决物流中的问题[考虑最大行驶距离限制、配送货物量]
针对1个仓库的优化
只用了一条DNA：
1 0 x x x x
第0位为货车类型
第1位位仓库节点
之后的xxxx表示接获点的标号，如：1 0 1 2 3 4 5

返回结果为：
best_Solve = [Truck1_Num, Truck2_Num, SolveSum, TimeSum, N_GENERATIONS]
包括：
Truck1_Num:第一种类型货车的数量
Truck2_Num:第二种类型货车的数量
SolveSum:返回当前迭代过程中的最佳解决方案{优化路径，最佳适度值，平均适度值，所花时间，所需货物量}
TimeSum：总的所需时间
