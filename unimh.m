function [edgedist,time]= unimh(n,p,iterations)
    
    %Allocate space for vectors
    %G = spalloc(n,n,n^2*p); %Sparse matrix with approximate number of edges
    edgesum = zeros(p*n^2,1);
    time=zeros(iterations,1);
    windowWidth = floor(2);
    %Metropolis-Hastings algorithm
    G=sprandsym(n,p);
    G(:,:)=logical(G(:,:)); % Set nonzero elements to 1

    
    for l = 1:iterations
        Gprop = G; % Proposed graph
        
        %%%%%%%%%%%%%%%%%% 
        %Link Swapping - swapping uniform random number of links
        delta=randi(windowWidth);
        for k=1:floor(delta/2)
            i(k) = randi(n); j(k)=randi(n); % Uniformly choose node pair
            % Force pair to be distinct
            while i(k)==j(k)
                j(k) = randi(n);
            end
        
            % Create/Remove link between (i,j), force symmetry
            Gprop(i(k),j(k)) = 1-G(i(k),j(k));
            Gprop(j(k),i(k)) = Gprop(i(k),j(k));
        end
        %%%%%%%%%%%%%%%%%%
        
        %Metropolis Rejection however we must now consider Q...x   
        if log(rand) < H(Gprop,n,p)-H(G,n,p)+Q(Gprop,n,delta)-Q(G,n,delta)
            G = Gprop;
        end
        
        m = full(sum(sum(G)))/2; % Count edges
        edgesum(m+1) = edgesum(m+1) + 1;
        time(l+1) = m;

    end
    
    % Normalise to create mass function
    edgedist = edgesum/sum(edgesum);

    function hamiltonian=H(G,n,p)
        m=sum(sum(G))/2; %number of edges on graph 
        beta=log(p/(1-p)); %log odds p
        hamiltonian=(m*beta);
    end

    function propdist = Q(G,n,delta)
        m=sum(sum(G))/2;
        propdist=-gammaln(0.5*n*(n-1)-m+1)+gammaln(0.5*n*(n-1)-m-delta+1) + gammaln(delta+1);
        if m<delta
            propdist = 0;
        end 
    end

end
        
