% This function initialize the first population of search agents
function [ X ]=initialization_latin(N,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1                               % 二维覆盖中上下界不是均为1维向量，因而边界数为1
    X=latin(N,dim).*(ub-lb)+lb;               % 初始化种群
end

% If each variable has a different lb and ub    每维度边界值不同
if Boundary_no>1
    X_temp = latin(N,dim);
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=X_temp(:,i).*(ub_i-lb_i)+lb_i;
    end
end

end

function [X] = latin(n,d)

% n - 样本数量
% d - 维度数量

X = zeros(n, d);
for i = 1:d
    % 为每个维度生成 [1/n, 2/n, ..., 1] 的点
    X(:, i) = (randperm(n) - 0.5) / n;
end

% 随机排序
for i = 1:d
    % 为每个维度生成随机排列
    X(:, i) = X(randperm(n), i);
end
end
