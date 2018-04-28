%% M-H Algorithm for producing random graphs
%%     on n vertices with predefined link probability p
%%     by changing one link at a time.

% Seed a chain at convergence then run loopLength times to assess mixing,
% output a matrix of timeseries.
function time=atModeLOOPonelink(n,p,iterations,loopLength)
    tic 
    function hamiltonian=H(G,p)
        m=sum(sum(G))/2; %number of edges on graph
        beta=log(p/(1-p)); %log odds p
        hamiltonian=(m*beta);
    end
    time=zeros(iterations,loopLength);
    %%%% Seed Chain at convergence.
    G=spalloc(n,n,(n^2)*p); %preallocate sparse matrix 
    for k=1:10^5 %Approximate stopping time upper bound
        Gprop=G; %Proposed G
        i=randi(n);
        j=randi(n);
        while i==j  %make sure not changing diagonal entries
            i=randi(n);
            j=randi(n);
        end
        Gprop(i,j)=1-G(i,j); %Break a link if there is one, make one if there isn't
        Gprop(j,i)=Gprop(i,j); %Force symmetry in the matrix
        if log(rand)<=min(0,(H(Gprop,p)-H(G,p)))
            G=Gprop;
        end
        m=full(sum(sum(G))/2);%count number of edges
        if m>floor(n*(n-1)*0.5*p)
            Gseed=G;
            break
        end
    end
    
    %End Seeding Chain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start loop at convergence   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for l=1:loopLength
        G=Gseed; %preallocate sparse matrix
        time(1,l)=full(sum(sum(G))/2);
        for k=1:iterations
            Gprop=G; %Proposed G
            i=randi(n);
            j=randi(n);
            while i==j  %make sure not changing diagonal entries
                i=randi(n);
                j=randi(n);
            end
            Gprop(i,j)=1-G(i,j); %Break a link if there is one, make one if there isn't
            Gprop(j,i)=Gprop(i,j); %Force symmetry in the matrix
            if log(rand)<=min(0,(H(Gprop,p)-H(G,p)))
                G=Gprop;
            end
            m=full(sum(sum(G))/2);%count number of edges    
            time(k+1,l)=m;
        end
    end
    
    toc
end


