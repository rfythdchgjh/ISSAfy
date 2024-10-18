



%%% Designed and Developed by Pavel Trojovský and Mohammad Dehghani %%%


function[Best_score,Best_pos,Curve_OOA]=OOA(pop,M,lb,ub,dim,fobj)
lb=ones(1,dim).*(lb);                              % Lower limit for variables
ub=ones(1,dim).*(ub);                              % Upper limit for variables

%% INITIALIZATION
for i=1:dim
    X(:,i) = lb(i)+rand(pop,1).*(ub(i) - lb(i));                          % Initial population
end

for i =1:pop
    L=X(i,:);
    fit(i)=fobj(L);
end
%%

for t=1:M  % algorithm iteration
    
    %%  update: BEST proposed solution
    [Fbest , blocation]=min(fit);
    
    if t==1
        xbest=X(blocation,:);                                           % Optimal location
        fbest=Fbest;                                           % The optimization objective function
    elseif Fbest<fbest
        fbest=Fbest;
        xbest=X(blocation,:);
    end
    %%
    %%
    for i=1:pop
        %% Phase 1: : POSITION IDENTIFICATION AND HUNTING THE FISH (EXPLORATION)
        fish_position=find(fit<fit(i));% Eq(4)
        if size(fish_position,2)==0
            selected_fish=xbest;
        else
            if rand <0.5
                selected_fish=xbest;
            else
                k=randperm(size(fish_position,2),1);
                selected_fish=X(fish_position(k));
            end
        end
        %
        I=round(1+rand);
        X_new_P1=X(i,:)+rand(1,1).*(selected_fish-I.*X(i,:));%Eq(5)
        X_new_P1 = max(X_new_P1,lb);X_new_P1 = min(X_new_P1,ub);
        
        % update position based on Eq (6)
        L=X_new_P1;
        fit_new_P1=fobj(L);
        if fit_new_P1<fit(i)
            X(i,:) = X_new_P1;
            fit(i) = fit_new_P1;
        end
        %% END Phase 1
        
        %%
        %% PHASE 2: CARRYING THE FISH TO THE SUITABLE POSITION (EXPLOITATION)
        X_new_P1=X(i,:)+(lb+rand*(ub-lb))/t;%Eq(7)
        X_new_P1 = max(X_new_P1,lb);X_new_P1 = min(X_new_P1,ub);
        % update position based on Eq (8)
        L=X_new_P1;
        fit_new_P1=fobj(L);
        if fit_new_P1<fit(i)
            X(i,:) = X_new_P1;
            fit(i) = fit_new_P1;
        end
        %% END Phase 2
        %%
    end
    %%
    
    best_so_far(t)=fbest;
    average(t) = mean (fit);
    
end
Best_score=fbest;
Best_pos=xbest;
Curve_OOA=best_so_far;
end

