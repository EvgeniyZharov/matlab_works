function SimulatedAnnealing()
    nCities = 100;
    initialTemperature = 100;
    endTemperature = 0;
    
    cities = rand(nCities, 2)*10; % setup initial cities positions as 2-dimentional array (x,y)

    figure  % Create new canvas for plot
    plot(cities(:,1),cities(:,2),"b--o" ); % draw initial route. The colon alone, without start or end values, specifies all of the elements in that dimension.
    title('Initial route')
    
    state = OptimiseRoute(cities, initialTemperature, endTemperature); % call our optimization function
    
    figure
    plot(cities(state,1),cities(state,2),"r--o"); % draw final route
    title('Optimized route')
end

function [ state ] = OptimiseRoute(cities, initialTemperature, endTemperature)
    nCities = size(cities,1);
    state = [1:nCities]'; % setup initial cities visit order as numbered column-vector (transposed row)

    currentEnergy = CalculateEnergy(state, cities); % calculate the energy for the initial condition
    disp('Initial route length: ');
    disp(currentEnergy);
    T = initialTemperature;
    
    for k = 1:100000 % main loop
        stateCandidate = GenerateStateCandidate(state); % create a new order for visiting cities
        candidateEnergy = CalculateEnergy(stateCandidate, cities); % calculate its energy
        
        if(candidateEnergy < currentEnergy) % if the new order has less energy
            state = stateCandidate; % it became the new order
            currentEnergy = candidateEnergy;
        else
            p = GetTransitionProbability(candidateEnergy-currentEnergy, T); % otherwise, calculate the probability
            if (IsTransition(p)) % if the transition occurs with a given probability
                state = stateCandidate; % accept the new order
                currentEnergy = candidateEnergy;
            end
        end

        T = DecreaseTemperature(initialTemperature, k);
        
        if(T <= endTemperature) % exit condition
            break;
        end
    end    
    disp('Final route length: ');
    disp(currentEnergy);
end

function [ E ] = CalculateEnergy(sequence, cities) % calculate route length
    n = size(sequence,1); % get size of first dimention (row count)
    E = 0;
    for i = 1:n-1
        E = E + Metric(cities(sequence(i),:), cities(sequence(i+1),:));
    end
    % add distance between finish and start to return to initial point
    E = E + Metric(cities(sequence(end),:), cities(sequence(1),:));
end

function [ distance ] = Metric( A, B ) % calculate distance between 2 points
    distance = (A - B).^2;
    distance = sqrt(distance);
    distance = sum(distance);
end

function [ T ] = DecreaseTemperature( initialTemperature, k)
    T = initialTemperature * 0.1 / k; 
end

function [ P ] = GetTransitionProbability( dE, T )
    P = exp(-dE/T);
end

function [ a ] = IsTransition( probability )
    if(rand(1) <= probability)
        a = 1;
    else
        a = 0; 
    end
end

function [ seq ] = GenerateStateCandidate(seq)
    n = size(seq, 1); % get size of cities indexes array
    i = randi(n); % get a pseudorandom index between 1 and n.
    j = randi(n);

    % 2 opt swap points

    if j < i
        k = j;
        j = i;
        i = k;
    end
    
    new_seq = seq(i:j);

    new_seq = new_seq.';
    new_seq = fliplr(new_seq);
    new_seq = new_seq.';
    seq(i:j) = new_seq;

%     t = seq(i);
%     seq(i) = seq(j);
%     seq(j) = t;
end
