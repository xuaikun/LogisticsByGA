# encoding: utf-8
import os
import math
import random
import time
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# 文件保存路径(根据你电脑的不同可进行修改)
path = "I:\\graduatestudy\\coding\\paper2"
# 保存的文件名
popDNAData_txt = os.path.join(path, 'popDNA.txt')

GenerationAndFitness_txt = os.path.join(path, 'GenerationAndFitness.txt')

# 迭代次数和种群大小会影响最终所得结果
NodeNumber = 10  # 节点数目
DNA_SIZE = 11  # DNA 长度，例如：1 0 1 2 3 4 5 6 7 8 9
# 第一个bit表示货车类型：0 载重量为2t，1载重量为5t，如果类型增加，可适当增加位数
# 第二个bit表示仓库，用标号0表示，一直不变的，若仓库数目增加，可改变位数
# 案例给出9个接货点，所以用1 2 3 4 5 6 7 8 9表示，若接货点增加，可以按照接货点的数量增加位数
POP_SIZE = 100  # population size 种群大小
CROSS_RATE = 0.8  # mating probability (DNA crossover)
MUTATION_RATE = 0.003  # mutation probability

# 迭代的总次数
N_GENERATIONS_Sum = 500

# 起始迭代次数
Start_N_GENERATIONS = 1

# 迭代步数
Step_N_GENERATIONS = 50


# 输入：节点初始距离矩阵(少量已知距离)，节点标号
# 输出：节点完整距离矩阵
def mydijkstra(routingArray, nodeNum):
    # routingArray数组为方阵
    n = routingArray.shape[0]  # 网络中节点个数
    # print "n =", n
    visited = np.zeros([1, n], dtype=np.int)  # 未被访问过的节点
    distance = np.empty([1, n], dtype=np.float32)  # 用于保存起点到各个顶点之间的距离，一行n列
    for i in range(0, n):
        distance[0][i] = float('inf')
    distance[0][nodeNum] = 0  # 起点到起点的距离为0
    for i in range(0, n):
        temp = []
        for j in range(0, n):
            temp.append(distance[0][j])
        # enumerate为枚举，visited_num为索引，visited_value为索引对应的值
        # 对应所得id1对应的节点已被标号
        id1 = [visited_num for visited_num, visited_value in enumerate(visited[0]) if visited_value == 1]
        # print "id1 =", id1
        # print "len(id1) =", len(id1)
        for p in range(0, len(id1)):
            temp[id1[p]] = float('inf')  # 已标号节点距离替换为无穷
        u = temp.index(min(temp))  # 找到标号值最小的节点
        # print "u =", u
        # np,show()
        visited[0][u] = 1  # 标记已标号的节点
        # 查找未标号的顶点
        id2 = [visited_num for visited_num, visited_value in enumerate(visited[0]) if visited_value == 0]
        # print "id2 =", id2
        for v in id2:
            if routingArray[u][v] + distance[0][u] < distance[0][v]:
                distance[0][v] = distance[0][u] + routingArray[u][v]  # 修改标号值
    # 返回结果
    # print "distance =\n", distance
    return distance


# 输入：节点数，初始化相关数据
# 输出：节点之间的距离邻接矩阵，值为节点之间的直接距离
def Init(nodeNumber):
    # 货车类型及载重量
    Truck1 = 2.0  # 货车1，载重量为2t
    Truck2 = 5.0  # 货车2，载重量为5t

    # 车辆一次最大行驶距离
    MaxTravelDis = 35.0  # 最大行驶距离为35km

    # 卸货时间
    UnloadTime = 5.0  # 卸货时间固定为5min

    # 车辆行驶速度
    DrivingSpeed = 10.0  # 行驶速度为10km/h

    # 每个派送人员的工作时间
    WorkTime = 8.0  # 派送人员工作时间为8小时
    # 将初始化数据保存到DataList中，以便返回使用
    DataList = [Truck1, Truck2, MaxTravelDis, UnloadTime, DrivingSpeed, WorkTime]

    # 初始化邻接矩阵，值为节点之间的直接距离
    NodeDistanceMatrix = np.zeros([nodeNumber, nodeNumber])
    # 初始化接货点的货物量
    CargoWeight = np.zeros([1, nodeNumber])

    # 案例
    CargoWeight[0] = [0, 1.7, 0.8, 1.3, 2.8, 1.9, 3.5, 0.9, 0.3, 1.2]

    NodeDistanceMatrix[0][1] = 5  # PA
    NodeDistanceMatrix[0][2] = 8  # PB
    NodeDistanceMatrix[0][3] = 7  # PC
    NodeDistanceMatrix[0][4] = 7  # PD
    NodeDistanceMatrix[0][5] = 4  # PE
    NodeDistanceMatrix[0][6] = 12  # PF
    NodeDistanceMatrix[0][7] = 9  # PG
    NodeDistanceMatrix[0][8] = 12  # PH
    NodeDistanceMatrix[0][9] = 6  # PI

    NodeDistanceMatrix[1][2] = 4  # AB
    NodeDistanceMatrix[1][9] = 3  # AI
    NodeDistanceMatrix[2][3] = 3  # BC
    NodeDistanceMatrix[3][4] = 4  # CD
    NodeDistanceMatrix[3][5] = 7  # CE
    NodeDistanceMatrix[4][5] = 3  # DE
    NodeDistanceMatrix[5][6] = 10  # EF
    NodeDistanceMatrix[6][7] = 4  # FG
    NodeDistanceMatrix[6][8] = 7  # FH
    NodeDistanceMatrix[7][8] = 5  # GH
    NodeDistanceMatrix[8][9] = 9  # HI
    # 节点之间的距离初始化
    for i in range(nodeNumber):
        for j in range(i, nodeNumber):
            # 求最小距离前，对距离矩阵进行处理：保持除了0以外的值不变，对角线不变，其余值均为inf
            if i != j and NodeDistanceMatrix[i][j] == 0:
                NodeDistanceMatrix[i][j] = np.inf
            if NodeDistanceMatrix[i][j] != 0:
                NodeDistanceMatrix[j][i] = NodeDistanceMatrix[i][j]

    # 求出节点之间的最小距离
    # print "NodeDistanceMatrix =\n", NodeDistanceMatrix
    minRouting = np.empty([nodeNumber, nodeNumber], dtype=list)
    for i in range(nodeNumber):
        rou = mydijkstra(NodeDistanceMatrix, i)
        minRouting[i] = rou
    # print "minRouting =\n", minRouting
    return DataList, CargoWeight, minRouting


# 最重要
# 到给定点最小距离和公式
# 适度函数，适度评分，每一个节点都有适度评分
def get_fitness(predDNA, NodeDistanceMatrix):
    # 使用np.zeros()时， NumSum每次初始化为0了,其中的值不可能为随机值
    # 之前使用np.empty()初始化NumSum导致矩阵中产生随机值，使得计算结果不正确
    # 用于保存每个点的适度值
    NumSum = np.zeros([1, len(predDNA)], dtype=float)
    for i in range(len(predDNA)):
        for j in range(1, len(predDNA[0]) - 1):
            NumSum[0][i] = NumSum[0][i] + NodeDistanceMatrix[predDNA[i][j]][predDNA[i][j + 1]]
        # 从最后一个节点返回仓库
        NumSum[0][i] = NumSum[0][i] + NodeDistanceMatrix[predDNA[i][1]][predDNA[i][len(predDNA[0]) - 1]]
    return NumSum[0]


# nature selection wrt pop's fitness
# 使用的是随机抽样，不过不是等概率的抽样，而是根据适度值的大小与适度值之和进行比较
# 通过群体的适度函数进行自然选择操作
# 本程序当前适用于求最小值的类型，若日后遇到求最大值的，则对最大值的处理就进行取消
def select(popDNA, fitness, pOP_SIZE):
    #  np.random.choice(a, size=3, replace=False, p=None) 表示抽样选择
    #  从a（a = np.arange(a)-->a个随机数）中以p的概率选择size个不相同的数，replace=False 表示抽出后不放回，表示不会出现重复数据
    #  replace=True表示抽出后继续放回，会出现重复数据， p=None 表示 概率一致性， p =[0.1,0, 0.3, 0.6, 0]选中每一个数的概率不相同
    #  返回的结果为选中的数据在a中的位置【有size个id】

    # {
    # 	最求最小值的处理
    # 	idx = np.random.choice(np.arange(pOP_SIZE), size=pOP_SIZE, replace=True,
    #  	                     p = fitness/fitness.sum())
    # }

    # {
    # 求最大值的处理

    # 定义最大fitness值
    maxfitness = np.zeros([1, len(fitness)], dtype=float)
    # 只取第一个元素的值，并且需要＋1e-3加个小的数不至于新的fitness值出现0-》导致概率p等于0(错误)
    # 1e-3 = 1X10^-3 = 1/1000 = 0.001
    maxfitness[:] = fitness[np.argmax(fitness)] + 1e-3
    # 选择的概率，目前选择概率=是当前节点的适应度值/适应度总和，本文实验适应度值是选择低的
    # 导致适应度低的节点没有选择😂，与实际结果相反了，

    # 解决方法，用最大适应度 - 当前适应度/（最大适应度 - 当前适应度）总和
    # 产生的点附近还会有更多的点（昨晚情况相反，适应度低的点，周围没什么点）
    # 修改fitness，得到新的fitness值
    fitness = maxfitness[0] - fitness
    # p 为更新后的概率
    # p = fitness/fitness.sum()
    # 日后如果遇到求最大值的就不需要以上处理
    idx = np.random.choice(np.arange(pOP_SIZE), size=pOP_SIZE, replace=True,
                           p=fitness / fitness.sum())
    # }
    # print "idx =", idx
    # print "popDNA[idx] =\n", popDNA[idx]
    # 选出暂时存活的个体
    return popDNA[idx]


# 交叉(交配)过程
def crossover(popDNA_m, popDNA_copy, dNA_SIZE, pOP_SIZE, RemainNodeList):
    # print "parent =\n", parent
    # print "pop =\n", pop
    # 交配概率
    if np.random.rand() < CROSS_RATE:
        # select another individual from pop
        # 从群体中选择另一个个体
        # size 表示生成几个数
        i_ = np.random.randint(0, pOP_SIZE, size=1)
        # choose crossover points
        # 选择交叉的节点,以True or False 形式存在
        # size 表示生成几个数
        # print "dNA_SIZE =", dNA_SIZE
        cross_points = np.random.randint(0, 2, size=dNA_SIZE).astype(np.bool)
        cross_points[1] = True
        # mating and produce one child
        # 生成孩子，作为下一代的父母
        # 将pop[i_, cross_points]赋值给parent[cross_points]
        popDNA_m[cross_points] = popDNA_copy[i_, cross_points]
        delete = popDNA_m[2: dNA_SIZE]
        delete_temp = []
        for i in range(len(delete)):
            delete_temp.append(delete[i])
        # 能够保证每个接货点都还存在
        x = random.sample(RemainNodeList, dNA_SIZE - 2)
        for i in range(len(x)):
            if x[i] not in delete_temp:
                for j in range(len(delete_temp)):
                    if delete_temp.count(delete_temp[j]) >= 2:
                        delete_temp[delete_temp.index(delete_temp[j])] = x[i]
        popDNA_m[2: dNA_SIZE] = delete_temp
    # 生成孩子
    # child = parent
    return popDNA_m


# 变异过程
def mutate(childDNA, dNA_SIZE, RemainNodeList):
    # DNA中任意一个点
    # childDNA = [1, 0, 6, 9, 3, 8 ,4, 7,2,5, 1]
    for point in range(dNA_SIZE):
        # 从DNA中突变某一节点，MUTATION_RATE突变概率
        # 0 变 1 ， 1变 0
        # 确保所突变的点不能为仓库,仓库位于基因的第1个位置
        if (np.random.rand() < MUTATION_RATE) and (point != 1):
            # 所选用的货车，在第0位，用0表示载重量为2t的货车
            if point == 0:
                childDNA[point] = 1 if childDNA[point] == 0 else 0
            # 其余的point点为接货点
            else:
                # 随机抽样接货点,返回闭区间的值
                # childDNA[point] = random.randint(1, dNA_SIZE - 2)
                childDNA_temp = random.sample(RemainNodeList, 1)
                print "childDNA_temp =", childDNA_temp
                childDNA[point] = childDNA_temp[0]
    childDNAList = childDNA[2:]
    childDNAList_temp = []
    for c_i in range(len(childDNAList)):
        childDNAList_temp.append(childDNAList[c_i])
    # 能够保证每个接货点都还存在
    x = random.sample(RemainNodeList, dNA_SIZE - 2)
    # 保证接货点的完整
    for i in range(len(x)):
        if x[i] not in childDNAList_temp:
            for j in range(len(childDNAList_temp)):
                if childDNAList_temp.count(childDNAList_temp[j]) >= 2:
                    childDNAList_temp[childDNAList_temp.index(childDNAList_temp[j])] = x[i]
    childDNA[2: dNA_SIZE] = childDNAList_temp
    # 孩子成长
    return childDNA


# 创建种群popDNA
def CreatPopDNA(usebackupFlag, pOP_SIZE, dNA_SIZE, nodeNumber, RemainNodeList, stateFlag):
    # 备份剩余的节点列表
    if stateFlag == "init":
        # 初始化popDNA
        if usebackupFlag is False:
            # 定义种群
            popDNA = np.zeros([pOP_SIZE, dNA_SIZE], dtype=list)
            for i in range(pOP_SIZE):
                # 产生不重复的数据（接货点从1开始）
                x = random.sample(range(1, nodeNumber), dNA_SIZE - 2)
                # 第1位为0，表示仓库
                x.insert(0, 0)
                # 第0为随机生成货车类型，类型0为2t货车，类型1为5t货车
                x.insert(0, np.random.randint(2))
                # 存入第i个值中
                popDNA[i] = x
            # x轴的DNA
            # 保存整数
            np.savetxt(popDNAData_txt, popDNA, fmt='%d')

        else:  # 使用备份数据
            # 导出群--整型
            popDNA = np.loadtxt(popDNAData_txt, dtype=np.int)
        x = random.sample(range(1, nodeNumber), dNA_SIZE - 2)
    else:
        # 定义种群
        popDNA = np.zeros([pOP_SIZE, dNA_SIZE], dtype=list)
        for i in range(pOP_SIZE):
            # 产生不重复的数据（接货点从1开始）
            # x = random.sample(range(1, nodeNumber), dNA_SIZE - 2)
            # 从剩余的节点中抽取dNA_SIZE - 2个节点用于生成种群
            # print "RemainNodeList =", RemainNodeList
            # print "len(RemainNodeList) =", len(RemainNodeList)
            # print "dNA_SIZE - 2 =", dNA_SIZE - 2
            x = random.sample(RemainNodeList, dNA_SIZE - 2)
            # 第1位为0，表示仓库
            x.insert(0, 0)
            # 第0为随机生成货车类型，类型0为2t货车，类型1为5t货车
            x.insert(0, np.random.randint(2))
            # print "x =", x
            # 存入第i个值中
            popDNA[i] = x
        # x轴的DNA
        # 保存整数
        np.savetxt('test.txt', popDNA, fmt='%d')
        x = RemainNodeList
    NodeList = []
    for i in range(len(x)):
        NodeList.append(x[i])
    return NodeList, popDNA


# GA算法
def GA(N_GENERATIONS, dNA_SIZE, pOP_SIZE, nodeNumber):
    # initialize the pop DNA
    # 种群初始化
    # pop为pOP_SIZE X dNA_SIZE的矩阵
    # DNA的长度
    # dNA_SIZE = 11  # DNA length

    # 初始化相关数据
    DataList, CargoWeight, NodeDistanceMatrix = Init(nodeNumber)
    # 创建种群
    NodeList = []
    # 保存0值，用做判断
    NodeList.append(0)
    NodeListTemp, popDNA = CreatPopDNA(True, pOP_SIZE, dNA_SIZE, nodeNumber, NodeList, 'init')
    NodeList = NodeListTemp
    print "NodeList =", NodeList

    bestDNA = 0  # 当前最佳解决方案
    bestfitness = 0  # 当前的最佳适度值
    avefitness = 0  # 平均适度值
    SolveSum = []  # 总的解决方案
    solveFlag = False  # 获得解决方案的标志
    While_Count = 0  # While中迭代的次数
    while solveFlag is False:
        # 迭代次数为N_GENERATIONS
        While_Count += 1
        print "While_Count =", While_Count
        for Step in range(N_GENERATIONS):

            # GA part (evolution)
            # 获得每个节点的适度值
            # 对种群的DNA进行评估
            fitness = get_fitness(popDNA, NodeDistanceMatrix)
            print "(While_Count, Step) =", (While_Count, Step)
            # print "fitness =\n", fitness
            # np.argmin(a) 找出a的最小值索引
            # 获得解决方案
            SolveNumber = np.argmin(fitness)  # 最小适度函数对应的第几个解决方案
            FitnessValue = fitness[np.argmin(fitness)]  # 最小适度值的值
            Solve = popDNA[np.argmin(fitness), :]  # 最小适度函数对应的详细解决方案
            print "np.argmin(fitness) =", SolveNumber
            print "min(fitness) =", FitnessValue
            print "Most fitted xDNA: ", Solve

            # 找到当前代的最小距离和
            bestfitness = fitness[np.argmin(fitness)]
            avefitness = fitness.sum() / len(fitness)
            bestDNA = np.argmin(fitness)

            # 根据每个个体的适度值以及自然选择的概率，选择存活的个体，组成新的群体
            # 要知道种群的大小在自然选择时，总的数量并没有改变
            # 对x, y进行自然选择
            popDNA = select(popDNA, fitness, pOP_SIZE)

            # 复制群体
            popDNA_copy = popDNA.copy()

            # 从全体中选父母用于产生后代群体
            for m in range(0, len(popDNA)):
                # 交叉
                childx = crossover(popDNA[m], popDNA_copy, dNA_SIZE, pOP_SIZE, NodeList)
                # 突变
                childx = mutate(childx, dNA_SIZE, NodeList)
                # parent is replaced by its child
                # 孩子代替父母
                # DNA的每一位都在被替换了
                popDNA[m][:] = childx
        # 将初始化数据保存到DataList中，以便返回使用
        # DataList = [Truck1, Truck2, MaxTravelDis, UnloadTime, DrivingSpeed, WorkTime]
        # 每辆车的最大行驶距离不能超过35KM
        MaxTravelDis = DataList[2]
        # 如果求出的最小适度值(最小哈密顿回路距离)大于车辆行驶最大距离，则需要减少当前哈密顿回路中节点的数量，每次减1个
        if FitnessValue > MaxTravelDis:
            print "MaxTravelDis =", MaxTravelDis
            print "FitnessValue =", FitnessValue

            # 接货点过多，超过了车辆一次行驶的最大距离
            # 尝试减少节点的数量
            dNA_SIZE = dNA_SIZE - 1
            NodeList, popDNA = CreatPopDNA(False, pOP_SIZE, dNA_SIZE, nodeNumber, NodeList, 'update')
        # 当前哈密顿回路的距离小于最大行驶距离
        # 现在需要判断该哈密顿回路中的接货点的货物量总和是否大于当前解决方案的货车载重量
        else:
            # 遍历解决方案中的元素
            GoodDemand = 0
            # 定义判断关注需求量时，当前解决方案是否合理的标志
            jugeFlag = True  # 初始值为合理

            print "Solve =", Solve
            for S_i in range(2, len(Solve)):
                # 计算当前哈密顿回路的货物需求量
                GoodDemand = GoodDemand + CargoWeight[0][Solve[S_i]]
            print "GoodDemand =", GoodDemand
            # 判断当前货车的载重量
            if Solve[0] == 0:  # 第一种类型的货车
                if DataList[Solve[0]] < GoodDemand:  # 当前货车不满足条件，需要换车
                    Solve[0] = 1
            print "check Solve =", Solve
            if Solve[0] == 1:  # 第二种类型的货车
                if DataList[0] >= GoodDemand:  # 如果使用载重量小的车也能满足送货需求，可换小车
                    Solve[0] = 0
                if DataList[Solve[0]] < GoodDemand:
                    # 说明当前解决方案不合理，需要新的解决方案
                    jugeFlag = False
                    print "当前货车满足条件，需要减少该方案中的节点个数或更新解决方案"
                    print "GoodDemand =", GoodDemand
                    # 接货点的需求量过多，超过了载重量最大的车辆的服务能力
                    # 尝试减少节点的数量
                    dNA_SIZE = dNA_SIZE - 1
                    print "NodeList =", NodeList
                    # 主要改变种群和基因，NodeList不改变
                    NodeList, popDNA = CreatPopDNA(False, pOP_SIZE, dNA_SIZE, nodeNumber, NodeList, 'update')
                    print "NodeList =", NodeList

            # 当前解决方案合理，直接应用
            if jugeFlag is True:
                '''
                print "NodeList =", NodeList
                print "Solve =", Solve
                print "SolveSum =", SolveSum
                '''
                InserSolveSum = []
                InserSolveSum.append(Solve)
                InserSolveSum.append(FitnessValue)
                InserSolveSum.append(round(avefitness, 2))
                # 卸货时间=停留点个数*(5/60)h
                UnloadTime = (len(Solve) - 2) * (DataList[3] / 60.0)
                # 车辆行驶时间=距离/速度
                DrivingTime = FitnessValue / DataList[4]
                # 总耗时
                TotalTime = UnloadTime + DrivingTime
                InserSolveSum.append(round(TotalTime, 2))
                # 哈密顿回路货物需求量
                GoodDemand = 0
                for S_i in range(2, len(Solve)):
                    # 计算当前哈密顿回路的货物需求量
                    GoodDemand = GoodDemand + CargoWeight[0][Solve[S_i]]

                InserSolveSum.append(round(GoodDemand, 2))
                SolveSum.append(InserSolveSum)
                print "before Solve =", Solve
                print "before NodeList =", NodeList

                for i in range(2, len(Solve)):
                    # 已有所属哈密顿图的节点需要删除，剩余未操作的节点
                    NodeList.remove(Solve[i])
                Solve = []
                print "after Solve =", Solve
                print "after NodeList =", NodeList
                # m = input()
                nodeNumber = len(NodeList)
                dNA_SIZE = nodeNumber + 2
                # 还有接货点没有加入哈密顿回路
                if nodeNumber != 0:  # 仍然剩余接货点
                    NodeList, popDNA = CreatPopDNA(False, pOP_SIZE, dNA_SIZE, nodeNumber, NodeList, 'update')
                # 节点全部加入哈密顿回路，退出while
                else:  # 结束"
                    break
    # 打印解决方案
    print "SolveSum =", SolveSum
    print "len(SolveSum) =", len(SolveSum)
    print "SolveSum[0] =", SolveSum[0]
    print "SolveSum[1] =", SolveSum[1]
    print "SolveSum[0][0] =", SolveSum[0][0]
    print "SolveSum[0][1] =", SolveSum[0][1]
    Truck1_Num = 0
    Truck2_Num = 0
    FitnessSum = 0
    avefitnessSum = 0
    TimeSum = 0
    GoodDemandSum = 0
    People = 0
    SolveList = []
    for i in range(len(SolveSum)):
        # 第0个元素为解决方案[1, 0, xxxxx],统计车辆数
        Truck2_Num = Truck2_Num + SolveSum[i][0][0]
        # 剩余部分分别为最佳适度值，平均适度值，所花时间，需求货物量
        FitnessSum = FitnessSum + SolveSum[i][1]
        avefitnessSum = avefitnessSum + SolveSum[i][2]
        TimeSum = TimeSum + SolveSum[i][3]
        GoodDemandSum = GoodDemandSum + SolveSum[i][4]
    # Truck1_Num 的数量
    Truck1_Num = len(SolveSum) - Truck2_Num
    print "Truck2_Num =", Truck2_Num
    print "TimeSum =", TimeSum
    People = math.ceil(TimeSum / DataList[5])
    # 返回最小距离和最佳解决方案，平均距离
    resultList = [SolveSum, Truck1_Num, Truck2_Num, FitnessSum, round(avefitnessSum, 2), round(TimeSum, 2),
                  round(GoodDemandSum, 2), People, N_GENERATIONS]
    print "resultList =", resultList

    # result = [bestfitness, bestDNA, avefitness]
    return resultList


# 主函数
def main():
    plt.rcParams['figure.figsize'] = (7.0, 6.0)
    # 主要用于画图中进行操作，线条的颜色
    LineColor = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    # LineColor =['b']
    # 线条的风格
    LineStyle = ['-', '--', '-.', ':']
    # 线条的标志
    LineLogo = ['.', 'o', 'v', '^', '>', '<', '1', '2', '3', '4', 's', 'p', '*']
    CSL_string0 = LineColor[0] + LineStyle[0] + LineLogo[1]
    CSL_string1 = LineColor[1] + LineStyle[1] + LineLogo[2]
    CSL_string2 = LineColor[2] + LineStyle[3] + LineLogo[3]

    Generationlist = []
    AveDislist = []
    MinDislist = []
    TimeList = []
    # 找到最佳的解决方案
    best_MinDis = np.inf
    best_Solve = []
    for N_GENERATIONS in range(Start_N_GENERATIONS, N_GENERATIONS_Sum, Step_N_GENERATIONS):
        GA_Start = time.time()
        # MinDis, bestXY, AveDis = GA(N_GENERATIONS, DNA_SIZE, POP_SIZE, NodeNumber)
        SolveSum, Truck1_Num, Truck2_Num, MinDis, AveDis, TimeSum, GoodDemandSum, People, N_GENERATIONS = GA(
            N_GENERATIONS, DNA_SIZE, POP_SIZE, NodeNumber)
        GA_End = time.time()
        print "AveDis =\n", AveDis
        print "MinDis =\n", MinDis
        if best_MinDis > MinDis:
            best_MinDis = MinDis
            best_Solve = [Truck1_Num, Truck2_Num, SolveSum, TimeSum, N_GENERATIONS]

        Generationlist.append(N_GENERATIONS)
        AveDislist.append(AveDis)
        MinDislist.append(MinDis)
        TimeList.append(round((GA_End - GA_Start), 2))

    print "best_MinDis =", best_MinDis
    print "best_Solve =", best_Solve
    # plt.show()
    datalist = []
    datalist.append(Generationlist)
    datalist.append(AveDislist)
    datalist.append(MinDislist)
    datalist.append(TimeList)

    np.savetxt(GenerationAndFitness_txt, datalist, fmt='%0.2f')

    plt.plot(Generationlist, AveDislist, CSL_string0, label="AveDis")
    plt.plot(Generationlist, MinDislist, CSL_string1, label="MinDis")
    # plt.legend(loc='upper right', edgecolor='black')

    # 设置x轴、y轴名称
    ax = plt.gca()
    ax.set_xlabel('Generation')
    ax.set_ylabel('Fitness')
    ax.xaxis.set_major_locator(MultipleLocator(20))
    pl.xticks(rotation=90)
    # ax.yaxis.set_major_locator(MultipleLocator(kedu))
    graph_path = os.path.join(path, 'GenerationAndFitness.png')
    plt.savefig(graph_path)
    plt.show()

    plt.plot(Generationlist, TimeList, CSL_string0, label="Time")
    # plt.legend(loc='upper right', edgecolor='black')
    # 设置x轴、y轴名称
    ax = plt.gca()
    ax.set_xlabel('Generation')
    ax.set_ylabel('Time')
    ax.xaxis.set_major_locator(MultipleLocator(20))
    pl.xticks(rotation=90)
    # ax.yaxis.set_major_locator(MultipleLocator(kedu))
    graph_path = os.path.join(path, 'GenerationAndTime.png')
    plt.savefig(graph_path)
    plt.show()


# 程序运行入口
if __name__ == "__main__":
    main()
