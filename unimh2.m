function [edgedist,time]= unimh2(n,p,iterations)
tic
    edgesum = zeros(n^2,1);
    time=zeros(iterations,1);
    windowWidth = 2;
    G=sprandsym(n,p);
    G=logical(G);
	G(logical(eye(n)))=0;
    [row,col]=find(tril(G));
    edgeList=[row,col];
    
    for l = 1:iterations
        propedgeList=edgeList;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Edge Swapping
        
        delta =randi([-windowWidth/2,windowWidth/2]);
		if delta==0
            time(l) = m;
			continue
		end 
        if delta > 0 || size(propedgeList,1)<=-delta
            i=zeros(delta,1);
            j=i;
            for k=1:abs(delta)
                i(k)=randi(n);
                j(k)=randi(n);
                %Remove self edges and do not pick edges that already exist
                while i(k)==j(k) || sum(ismember(propedgeList,[i(k),j(k)],'rows'))>0 ...
				||sum(ismember(propedgeList,[j(k),i(k)],'rows'))>0
                    j(k) = randi(n);
                    i(k)= randi(n);
                end
                %Update list
                propedgeList=[propedgeList;[i(k),j(k)]];
            end
        end
        
        if delta < 0
            for k=1:abs(delta)
                %Choose random element from edgeList for removal
                edge2remove=randperm(size(propedgeList,1));
                vertexPair=propedgeList(edge2remove(1),:);
                %Remove (i,j) from edgeList
                remove=intersect(find(propedgeList(:,1)==vertexPair(1)),find(propedgeList(:,2)==vertexPair(2)));
                propedgeList(remove,:)=[];
                %Remove (j,i) link just in case
                remove=intersect(find(propedgeList(:,1)==vertexPair(2)),find(propedgeList(:,2)==vertexPair(1)));
                propedgeList(remove,:)=[];
            end
        end
		
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Metropolis Hastings
        beta=log(p/(1-p));
        propadj = sparse(propedgeList(:, 1), propedgeList(:, 2), 1, n,n);
        adj = sparse(edgeList(:, 1),edgeList(:, 2), 1, n,n);
        %deltaH=beta*(sum(sum(propadj))/2-sum(sum(adj))/2);
        deltaH=beta*(size(propedgeList,1)-size(edgeList,1));
        
        if log(rand) <= deltaH %+deltaQ(size(edgeList,1),n,delta)
		edgeList=propedgeList;			
        end
        
        m = size(edgeList,1); % Count edges
        edgesum(m+1) = edgesum(m+1) + 1;
        time(l) = m;
    end
    
    % Normalise to create mass function
    edgedist = edgesum;%/sum(edgesum);
	plot(graph(table(edgeList,'VariableNames',{'EndNodes'})));
    %function propdist = deltaQ(m,n,delta)
        %propdist=gammaln(0.5*n*(n-1)-m+1)-2*gammaln(0.5*n*(n-1)-m-delta+1)+gammaln(0.5*n*(n-1)-m-2*delta+1);
        
        %If delta+-1 then propdist reduces to:  
        %propdist=log((0.5*n*(n-1)-m)/(0.5*n*(n-1)-m-1));
        %Is there a similar simple form for all delta? Saves invoking gammaln
%        if m<delta
%            propdist = 0;
%        endif 
    %end
toc
end