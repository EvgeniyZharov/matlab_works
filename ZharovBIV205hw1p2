function genetic_TSP1
    clc;close all;clear all;
    N=100;%Количество городов
    city_coordinate = rand(N, 2)*10;
    s=100;%Количество образцов
    c=100;%Номер замены
    pc=0.9;%Вероятность кроссовера
    pm=0.2;%Вероятность мутации
    times=5000;                      %Максимальное количество итераций
    time=0;                         %Фактическое количество итераций
    pop=zeros(s,N+1);%Начальная популяция+адаптируемость
    pop_fit_aver=[];%Тотальный фитнес
    min_dis=[];%Наименьшее расстояние
    pop_min=[];%Ген кратчайшего расстояния
    for i=1:s   %Initialize
        pop(i,1:N)=randperm(N);
    end
    clf
    plot(city_coordinate(:,1),city_coordinate(:,2),'ro');%Нарисуйте график средней пригодности
    for i=1:N
        test_t=num2str(i);
        text(city_coordinate(i,1),city_coordinate(i,2),test_t);%метка
    end
    grid on;
    title("City map");
    xlabel('x');
    ylabel('y');
    
    city_distance=CityDistance(city_coordinate,N);%Междугороднее расстояние
    [individual_fit,sum,min1,min_index]=GroupFit(city_distance,N,pop,s);%адаптируемость
    sumP=sum;
    pop_fit_aver=[pop_fit_aver;sum];
    min_dis=[min_dis;min1];
    pop(:,N+1)=individual_fit;
    pop_min=[pop_min;pop(min_index,:)];
    pop=ChooseParents(pop,N,s,c);%Выберите родителя
    for i=1:times
        time=time+1;
        E_new_new=[];  
        for j=1:s/2    
            a=rand(1);
            b=rand(1);
                if a>pc         
                    ;
                else
                    crosspoint=rand(1,2);
                    crosspoint=floor(crosspoint*N)+1;
                    [pop(j,:),pop(j+s/2,:)]=CrossVariation(pop(j,:),pop(j+s/2,:),crosspoint,N);       
                end
                if b>pm
                    ;
                else
                    pop(j,:)=Mutation(pop(j,:),N);
                    pop(j+s/2,:)=Mutation(pop(j+s/2,:),N);
                end
                E_new_new=[E_new_new;pop(j,:);pop(j+s/2,:)];
        end
        [individual_fit,sum,min1,min_index]=GroupFit(city_distance,N,E_new_new,s);
        sumS=sum;
        pop_fit_aver=[pop_fit_aver;sum];
        min_dis=[min_dis;min1];
        E_new_new(:,N+1)=individual_fit;
        pop_min=[pop_min;E_new_new(min_index,:)];
        if(abs(sumS-sumP)<0.001)%Условия выхода
            break;
        end
        pop=ChooseParents(E_new_new,N,s,c);
    end
    [a,min_index]=min(min_dis);
    a;
    time1=1:time+1;

    figure%Нарисуйте график средней пригодности
    plot(time1,min_dis,'k.');
    disp(min_dis(size(min_dis,1)));
    grid on;
    title("Разброс графика минимального значения на поколение");
    xlabel("Количество итераций");
    ylabel("Кратчайшее расстояние");

    figure%Нарисуй оптимальный путь

    DrawPath(city_coordinate,pop_min,min_index,N)

    grid on;
    title("Оптимальная схема пути");
    xlabel('x');
    ylabel('y');

end

function [city_distance] = CityDistance(city_coordinate,N)%Городская матрица расстояний
    city_distance=zeros(N,N);
    for i=1:N
        for j=1:N
            city_distance(i,j)=((city_coordinate(i,1)-city_coordinate(j,1))^2+...
                (city_coordinate(i,2)-city_coordinate(j,2))^2)^0.5;
        end
    end
end

function [individual_fit,num,min_distance,a] = GroupFit(city_distance,N,pop,s)%Популяционный фитнес
    individual_distance=zeros(s,1);
    for j=1:s
        sum_distance=0;
        for i=1:N-1
            sum_distance=sum_distance+city_distance(pop(j,i),pop(j,i+1));
        end
        sum_distance=sum_distance+city_distance(pop(j,N),pop(j,1));
        individual_distance(j,1)=sum_distance;
    end
    [min_distance,a]=min(individual_distance);
    individual_fit=1./individual_distance;
    num=0;
    for i=1:s
      num=num+individual_fit(i,1);
    end
end

function [pop_ok]=ChooseParents(pop,N,s,c)%Выберите родителя
    pop=sortrows(pop,N+1);
    for i=1:c
        pop(i,:)=pop(s+1-i,:);
    end
    randIndex=randperm(size(pop,1));
    pop=pop(randIndex,:);
    pop_ok=pop;
end

function [a,b]=SwapRepeat(tbl,pop1,pop2,c1,c2,N)%Генная дедупликация
    i=100/N;
    for k=1:(c1-1)
        if tbl(pop1(k),3)>i
            kk=find(pop1(c1:c2)==pop1(k))+c1-1;
            kkk=pop1(k);
            pop1(k)=pop2(kk);
            pop2(kk)=kkk;
        end
    end
    for k=c2+1:N
        if tbl(pop1(k),3)>i
            kk=find(pop1(c1:c2)==pop1(k))+c1-1;
            kkk=pop1(k);
            pop1(k)=pop2(kk);
            pop2(kk)=kkk;
        end
    end
    a=pop1;
    b=pop2;
end

function [a,b]=CrossVariation(pop1,pop2,crosspoint,N)%Кроссовер
    A=pop1;
    if(crosspoint(:,1)<crosspoint(:,2))
        pop1(crosspoint(:,1):crosspoint(:,2))=pop2(crosspoint(:,1):crosspoint(:,2));
        pop2(crosspoint(:,1):crosspoint(:,2))=A(1,crosspoint(:,1):crosspoint(:,2));
        while 1
            tbl = tabulate(pop1(1:N));
            if (tbl(:,3)<=(100/N))
                break;
            end
            [pop1,pop2]=SwapRepeat(tbl,pop1,pop2,crosspoint(:,1),crosspoint(:,2),N);
        end
    else
        pop1(crosspoint(:,2):crosspoint(:,1))=pop2(crosspoint(:,2):crosspoint(:,1));
        pop2(crosspoint(:,2):crosspoint(:,1))=A(1,crosspoint(:,2):crosspoint(:,1));
        while 1
            tbl = tabulate(pop1(1:N));
            if (tbl(:,3)<=(100/N))
                break;
            end
            [pop1,pop2]=SwapRepeat(tbl,pop1,pop2,crosspoint(:,2),crosspoint(:,1),N);
        end
    end
    a=pop1;b=pop2;
end

function [a]=SwapGene(sub,c1,c2)%обмен
    kk=ceil((c2-c1)/2);
    kkk=(c2-c1)+2;
    for k=1:kk
        kkkk=sub(k);
        sub(k)=sub(kkk-k);
        sub(kkk-k)=kkkk;
    end
    a=sub;
end

function [a]=Mutation(pop0,N)%Мутации
    crosspoint=rand(1,2);
    crosspoint=floor(crosspoint*N)+1;
    if(crosspoint(:,1)<crosspoint(:,2))
        sub=pop0(crosspoint(:,1):crosspoint(:,2));
        sub=SwapGene(sub,crosspoint(:,1),crosspoint(:,2));
        pop0(crosspoint(:,1):crosspoint(:,2))=sub;
    else
        sub=pop0(crosspoint(:,2):crosspoint(:,1));
        sub=SwapGene(sub,crosspoint(:,2),crosspoint(:,1));
        pop0(crosspoint(:,2):crosspoint(:,1))=sub;
    end
    a=pop0;
end

function DrawPath(city_coordinate,E_new_new,min_index,N)%Нарисовать карту пути
    k=E_new_new(min_index,1:N);
    %plot(kkk(:,1),kkk(:,2),'b');%Нарисуйте график средней пригодности
    plot(city_coordinate(:,1),city_coordinate(:,2),'bo');
    hold on;
    for i=1:N-1
        plot([city_coordinate(k(i),1),city_coordinate(k(i+1),1)],[city_coordinate(k(i),2),city_coordinate(k(i+1),2)],'r','LineWidth',2);
        test_t=num2str(i);
        % text(city_coordinate(k(i),1),city_coordinate(k(i),2),test_t);
        hold on;
    end
    test_t=[num2str(N)];
    text(city_coordinate(k(N),1),city_coordinate(k(N),2),test_t);
end
