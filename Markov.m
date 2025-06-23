% 三个节点的马尔科夫链

% L矩阵特征值 0 -1.5 -1.5   零特征值的特征向量(1,1,1)与概率平面正交，因此整体只有一个不动点，遍历
% 遍历就是说概率平面上的所有点都收敛到这个全局唯一的不动点上
% A1 = [0,0.5,0.5 ; 0.5,0,0.5; 0.5,0.5,0];
% test3(A1,Step=5,Visualize="Scatter");

% L矩阵特征值 0 -0.46 -1.44     和例一不同，收敛到非中心点，但仍然唯一，所以还是遍历的
% A2 = [1/2,1/2,0 ; 1/3,0,2/3; 0,2/5,3/5];
% test3(A2,Step=10,Visualize="Scatter");

% 遍历针对一个链的所有初始状态，平稳是针对特定初始状态而言的
% 前者看整个中心流形和概率平面的交集，后者看单条运动轨迹


% 猜测：常返态对应脑电的周期成分，滑过态对应非周期成分

% 常返与吸引子？

% 周期


% 吸收态
% L矩阵的特征值为 0 -1/3 -1
% A3 = [1/3,1/3,1/3 ; 0,1,0 ; 1/3,1/3,1/3];
% test3(A3,Step=10,Visualize="Scatter"); 

% A4 = [1/2,1/2,0 ; 1/2,0,1/2; 0,0,1];
% test3(A4,Step=20,Visualize="Scatter");

% L矩阵特征值 0 0 -1  ， 收缩的特征向量是[-1,1,0]  ,沿着这个方向空间变成一个竖着的切面，与概率平面交出一条直线
% 这条直线上的点都是平稳而非遍历的马尔科夫链
A5 = [1/2,1/2,0 ; 1/2,1/2,0 ; 0,0,1];
test3(A5,Step=3,Visualize="Scatter");

function [] = test3(A,Config)
    
    arguments 
        A (3,3) { mustBeNonnegative,mustBeNumeric,mustBeNonempty }
        Config.Step { mustBePositive,mustBeInteger } = 1
        Config.Visualize { mustBeMember(Config.Visualize,["Plot","Scatter"]) } = "Plot"
    end

    % 在棱锥面上生成概率空间，然后看迭代的结果
    delta = 0.1;

    X = zeros(0,3);
    N = 0;
    for x=0:delta:1
        for y=0:delta:1-x
            N = N+1;
            X(N,:) = [x,y,max(1-x-y,0)]; % 浮点数运算可能导致<0
        end
    end

    Matrix = zeros(N,Config.Step+1,3);
    for i=1:N
        if(Config.Visualize=="Plot")
            Markov3(A,X(i,:),Step=Config.Step,Plot=true);
            hold on;
        elseif(Config.Visualize=="Scatter")
            Matrix(i,:,:) = Markov3(A,X(i,:),Step=Config.Step,Plot=false);
        end
    end

    if(Config.Visualize=="Scatter")
        for t=1:Config.Step
            xs = squeeze(Matrix(:,t,1));
            ys = squeeze(Matrix(:,t,2));
            zs = squeeze(Matrix(:,t,3));
            scatter3(xs,ys,zs);
            xlim([0,1]);
            ylim([0,1]);
            zlim([0,1]);
            drawnow
            pause(1);
        end
    end

end

function [X] = Markov3(A,x,Config)
    arguments 
        A (3,3) { mustBeNonnegative,mustBeNumeric,mustBeNonempty }
        x { mustBeVector,mustBeNumeric,mustBeNonnegative, mustBeNonempty }
        Config.Step { mustBePositive,mustBeInteger } = 1
        Config.Plot { mustBeNumericOrLogical } = true
    end

    X = zeros(Config.Step+1,length(x));
    X(1,:) = x;

    for i=2:Config.Step+1
        X(i,:) = x * A; % 不是微分方程，是连续的线性变换
        x = X(i,:);
    end
    
    if(Config.Plot)
        plot3(X(:,1),X(:,2),X(:,3),'o-');
        xlim([0,1]);
        ylim([0,1]);
        zlim([0,1]);
        grid on;
    end

    % disp(X(1,:));
    % disp(X(end,:));
    % disp("=====");


end