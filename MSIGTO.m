% Initialization of Latin Hypercube && Levy Flight Factor && Cauchy Inverse Cumulative Distribution Factor
function [Silverback_Score,Silverback,conv_curve]=MSIGTO(pop_size,Max_iter,lb,ub,dim,fobj)
disp('IGTO is optimizing your problem');
%% LHS
% initialize Silverback 
Silverback = [];
Silverback_Score = inf;                   % inf 
X=initialization_latin(pop_size,dim,ub,lb); 
conv_curve=zeros(1,Max_iter);   
for i=1:pop_size
    Pop_Fit(i)=fobj(X(i,:));   
    if Pop_Fit(i) < Silverback_Score
        Silverback_Score=Pop_Fit(i);     
        Silverback=X(i,:);
    end
end

GX = X(:,:);                             
lb = ones(1,dim).*lb;                    
ub = ones(1,dim).*ub;

%%  Controlling parameter 

p = 0.03;
w = 0.8;

%% Main loop 
for It=1:Max_iter

    C = (cos(2*rand)+1)*(1-It/Max_iter);         % 'cos(2*rand)+1' serves as the F parameter 
    L = C*(2*rand-1);                            % L Simulation of Silverback Leadership

    %% Exploration: 
    for i=1:pop_size
        if rand<p
            GX(i,:) =(ub-lb)*rand+lb;                 % rand < p, migrate to an unknown area
        else
            if rand >= 0.5                            % rand >= 0.5, migrate to the position of another gorilla
                Z = unifrnd(-C,C,1,dim);     
                H = Z.*X(i,:);                        
                GX(i,:) = (rand-C)*X(randi([1,pop_size]),:)+L.*H;
            else                                      % rand < 0.5, migrate to the known area
                GX(i,:) = X(i,:)-L.*(L*(X(i,:)- GX(randi([1,pop_size]),:))+rand*(X(i,:)-GX(randi([1,pop_size]),:)));
            end
        end
    end
    % boundry check and population update
    GX = boundarycheck(GX, lb, ub);      
    for i=1:pop_size
        New_Fit = fobj(GX(i,:));        
        if New_Fit < Pop_Fit(i)
            Pop_Fit(i) = New_Fit;        
            X(i,:) = GX(i,:);             
        end
        if New_Fit < Silverback_Score   
            Silverback_Score = New_Fit;
            Silverback = GX(i,:);
        end
    end

    conv_curve(It)=Silverback_Score;

    %% Exploition 
    for i=1:pop_size
        if C >= w           % follow the silverback
            g = 2^L;
            M = (abs(mean(GX)).^g).^(1/g);             
            
            % introduce LF method
            v = levy(1,dim,1.5); 
            GX(i,:) = L*M.*(X(i,:)-Silverback)*0.01.*v + X(i,:);
        else                % computation for female
            r1=rand;
            % introduce CICDO
            p1 = rand(1,dim);           
            GX(i,:) = Silverback-(Silverback*(2*r1-1)-X(i,:)*(2*r1-1)).*tan(pi*(p1-0.5));
        end
    end

    % boundry check and population update
    GX = boundarycheck(GX, lb, ub);
    for i=1:pop_size
        New_Fit= fobj(GX(i,:));
        if New_Fit < Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
        end
        if New_Fit < Silverback_Score
            Silverback_Score=New_Fit;
            Silverback=GX(i,:);
        end
    end

    conv_curve(It)=Silverback_Score;
    %fprintf("In Iteration %d, best estimation of the global optimum is %4.4f \n ", It,Silverback_Score );
end