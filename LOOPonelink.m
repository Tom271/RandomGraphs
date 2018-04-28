%% M-H Algorithm for producing random graphs
%%     on n vertices with predefined link probability p
%%     by changing one link at a time.

function [stoppingTime,time]=LOOPonelink(n,p,iterations,loopLength)
    tic 
    function hamiltonian=H(G,p)
        m=sum(sum(G))/2; %number of edges on graph
        beta=log(p/(1-p)); %log odds p
        hamiltonian=(m*beta);
    end
    time=zeros(iterations,loopLength);
    time(1,loopLength)=0;
    stoppingTime=zeros(loopLength,1);

    for l=1:loopLength
        G=spalloc(n,n,(n^2)*p); %preallocate sparse matrix
        %G=sprandsym(n,p);
        %G(:,:)=logical(G(:,:)); % Set nonzero elements to 1
        for k=1:iterations
            Gprop=G; %Proposed G
            i=randi(n);
            j=randi(n);
            while i==j  %make sure not changing diagonal entries
                i=randi(n);
            end
            Gprop(i,j)=1-G(i,j); %Break a link if there is one, make one if there isn't
            Gprop(j,i)=Gprop(i,j); %Force symmetry in the matrix
            if log(rand)<=min(0,(H(Gprop,p)-H(G,p)))
                G=Gprop;
            end
            m=full(sum(sum(G))/2);%count number of edges
            
            if m>floor(n*(n-1)*0.5*p)
                stoppingTime(l)=k;
                time(k:length(time),l)=floor(0.5*n*(n-1)*p);
                break
            end     
            time(k+1,l)=m;
        end
        hold on;
        figure(2);
        plt=plot(time(1:stoppingTime(l),l));
        plt.Color(4) = 0.35;
    end
    maxTime=max(stoppingTime);
    avgTimeSeries=sum(time,2)./sum(time~=0,2);
    figure(2);plot(avgTimeSeries(1:maxTime),'LineWidth',3,'color','g')
    toc
end

