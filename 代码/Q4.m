clear;clc;
%% 数据初始化
N1=4;
N2=8;
F_N1=dftmtx(N1);
F_N2=dftmtx(N2);
F_N=kron(F_N1,F_N2);%Kronecker积

k=3;%分解产生的矩阵数
q=3;

% 随机种子序列
seeds = zeros(1, 2*q+1);
for i = 1:q
    seeds(i) = 2^(i-1);
    seeds(q+i) = -2^(i-1);
end

% 分别对Fn1 Fn2进行近似稀疏分解
%% 第一步：分解Fn1

% 遗传算法参数设置
populationSize = 200;         % 种群数量
chromosomeLength = k*N1*N1+1;  % 染色体长度
max_iter = 1000;              % 最大迭代次数
mutationRate = 0.1;         % 变异率
crossoverRate = 0.85;         % 交叉率
tournamentSize =5;          % 锦标赛选择法中每轮锦标赛的参与个体数

% 存储每次迭代的最优适应度值
bestFitnessHistory1 = zeros(max_iter, 1);

% 初始化种群
population = zeros(populationSize, chromosomeLength);
A=cell(1,k+1);
for i = 1:populationSize
    for j=1:k
        A{j}=Matrix_generate(N1,seeds);
        population(i,(j-1)*N1*N1+1:j*N1*N1)=reshape(A{j},1,N1*N1);
    end
    population(i,end)=rand();
end

for iter = 1:max_iter
    % 计算适应度
    fitness = zeros(populationSize, 1);
    for i = 1:populationSize
        fitness(i) = norm_F(population(i, :), F_N1, k, N1);
    end

    % 记录当前迭代的最优适应度值以及最优解
    [bestFitness1,index] = min(fitness);
    bestSolution1 = population(index,:);
    bestFitnessHistory1(iter) = bestFitness1;

    % 选择
    % 锦标赛选择法
    tournamentPopulation = zeros(populationSize, chromosomeLength);
    for i = 1:populationSize
        tournamentIndices = randperm(populationSize, tournamentSize);%从所有个体中抽取出参与锦标赛的
        tournament = population(tournamentIndices, :);%出列
        tournamentFitness(i) = min(fitness(tournamentIndices));%让竞赛者竞争产生赢家
        %保留这个赢家的数据
        tournamentPopulation(i, :) = tournament(find(fitness(tournamentIndices) == tournamentFitness(i), 1, 'first'), :);
    end

    % 交叉
    for i = 1:2:populationSize-1
        r = rand();
        if r <= crossoverRate
            % 使用单点交叉法进行交叉，同时不破坏约束一
            if k<=2
                crossPoint = randi([1,N1-1])*N1;
            else
                crossPoint = min(randi([2,k])*N1*N1+1,randi([2,k-1])*N1*N1+1+randi([1,N1-1])*N1);
            end
            tournamentPopulation([i, i+1], crossPoint:end) = tournamentPopulation([i+1, i], crossPoint:end);
        end
    end

    % 变异
    for i = 1:populationSize
        for j = 1:chromosomeLength-1
            r = rand();%骰点
            if r < mutationRate
                r=rand();%骰点
                % 对染色体进行变异
                if r <=1 && r>0.1   %策略一
                    % 生成一个随机值替换当前基因
                    if population(i,j)~=0
                        population(i,j)=seeds(randperm(size(seeds,2),1))+1i*seeds(randperm(size(seeds,2),1));
                    end
                else                %策略二
                    % 敲除当前基因
                    population(i, j) = 0;
                end
            end
        end
    end

    % 更新种群
    population = tournamentPopulation;
end

%解码染色体
for i=1:k
    A{i}=reshape(bestSolution1(N1^2*(i-1)+1:N1^2*i),N1,N1);
end
A{k+1}=bestSolution1(end);
%% 第二步 对Fn2进行分解
% 遗传算法参数设置
populationSize = 200;         % 种群数量
chromosomeLength = k*N2*N2+1;  % 染色体长度
max_iter = 250;              % 最大迭代次数
mutationRate = 0.1;         % 变异率
crossoverRate = 0.85;         % 交叉率
tournamentSize =5;          % 锦标赛选择法中每轮锦标赛的参与个体数

% 存储每次迭代的最优适应度值
bestFitnessHistory2 = zeros(max_iter, 1);

% 初始化种群
population = zeros(populationSize, chromosomeLength);
B=cell(1,k+1);
for i = 1:populationSize
    for j=1:k
        B{j}=Matrix_generate(N2,seeds);
        population(i,(j-1)*N2*N2+1:j*N2*N2)=reshape(B{j},1,N2*N2);
    end
    population(i,end)=rand();
end

for iter = 1:max_iter
    % 计算适应度
    fitness = zeros(populationSize, 1);
    for i = 1:populationSize
        fitness(i) = norm_F(population(i, :), F_N2, k, N2);
    end

    % 记录当前迭代的最优适应度值以及最优解
    [bestFitness2,index] = min(fitness);
    bestSolution2 = population(index,:);
    bestFitnessHistory2(iter) = bestFitness2;

    % 选择
    % 锦标赛选择法
    tournamentPopulation = zeros(populationSize, chromosomeLength);
    for i = 1:populationSize
        tournamentIndices = randperm(populationSize, tournamentSize);%从所有个体中抽取出参与锦标赛的
        tournament = population(tournamentIndices, :);%出列
        tournamentFitness(i) = min(fitness(tournamentIndices));%让竞赛者竞争产生赢家
        %保留这个赢家的数据
        tournamentPopulation(i, :) = tournament(find(fitness(tournamentIndices) == tournamentFitness(i), 1, 'first'), :);
    end

    % 交叉
    for i = 1:2:populationSize-1
        r = rand();
        if r <= crossoverRate
            % 使用单点交叉法进行交叉，同时不破坏约束一
            if k<=2
                crossPoint = randi([1,N2-1])*N2;
            else
                crossPoint = min(randi([2,k])*N2*N2+1,randi([2,k-1])*N2*N2+1+randi([1,N2-1])*N2);
            end
            tournamentPopulation([i, i+1], crossPoint:end) = tournamentPopulation([i+1, i], crossPoint:end);
        end
    end

    % 变异
    for i = 1:populationSize
        for j = 1:chromosomeLength-1
            r = rand();%骰点
            if r < mutationRate
                r=rand();%骰点
                % 对染色体进行变异
                if r <=1 && r>0.1   %策略一
                    % 生成一个随机值替换当前基因
                    if population(i,j)~=0
                        population(i,j)=seeds(randperm(size(seeds,2),1))+1i*seeds(randperm(size(seeds,2),1));
                    end
                else                %策略二
                    % 敲除当前基因
                    population(i, j) = 0;
                end
            end
        end
    end

    % 更新种群
    population = tournamentPopulation;
end

%解码染色体
for i=1:k
    B{i}=reshape(bestSolution2(N2^2*(i-1)+1:N2^2*i),N2,N2);
end
B{k+1}=bestSolution2(end);

%% 将分解后的矩阵按照Kronecker性质进行相乘，并计算硬件复杂度
l=0;FN_approx=eye(N1*N2);
for i=1:k
    l=l+computeComplexity(FN_approx,kron(A{i},B{i}));
    FN_approx=FN_approx*kron(A{i},B{i});
end
C=q*l;

disp(['C=' num2str(C)]);
disp(['bestFitness1：' num2str(bestFitness1)]);
disp(['bestFitness2：' num2str(bestFitness2)]);
disp(['RMSE：' num2str(1/(N1*N2)*norm(F_N-FN_approx,'fro'))])

%% 函数定义
%生成随机矩阵
function M=Matrix_generate(N,x)
    M=zeros(N);
    for i=1:N
        x1=zeros(1,N);
        x1(1)=x(randperm(size(x,2),1))+1i*x(randperm(size(x,2),1));
        x1(2)=x(randperm(size(x,2),1))+1i*x(randperm(size(x,2),1));
        x1=x1(randperm(size(x1,2),size(x1,2)));
        M(i,:)=x1;
    end
end

%适应度-目标函数
function f = norm_F(x,F_N,k,N)
    A=cell(1,k);
    for i=1:k
        A{i}=reshape(x(N^2*(i-1)+1:N^2*i),N,N);
    end
    Ak=eye(N);
    for i=1:k
        Ak=Ak*A{i};
    end
    f=1/N*norm(x(end)*F_N-Ak,'fro');
end