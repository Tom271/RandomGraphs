function [edgedist,time]= unimh2(n,p,iterations)

    edgesum = zeros(n^2,1);
    time=zeros(iterations,1);
    windowWidth = 2;
    G=sprandsym(n,p);
    G(:,:)=logical(G(:,:));
    [row,col]=find(G);
    edgeList=[row,col];
    
    for l = 1:iterations
        propedgeList=edgeList;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Edge Swapping
        
        delta = randi([-windowWidth/2,windowWidth/2]);
        if delta > 0 | size(propedgeList,1)<-delta
            for k=1:abs(delta)
                i(k)=randi(n);
                j(k)=randi(n);
                #Remove self edges and do no pick edges that already exist
                while i(k)==j(k) | sum(ismember(propedgeList,[i(k),j(k)],'rows'))>0
                    j(k) = randi(n);
                endwhile
                
                %Update list
                propedgeList=[propedgeList;[i(k),j(k)]];
            endfor
        endif
        
        if delta < 0
            for k=1:abs(delta)
                #Choose random element from edgeList for removal
                vertexPair=propedgeList(randperm(size(propedgeList,1))(1),:);
                #Remove (i,j) from edgeList
                remove=intersect(find(propedgeList(:,1)==vertexPair(1)),find(propedgeList(:,2)==vertexPair(2)));
                propedgeList(remove,:)=[];
                #Remove (j,i) link just in case
                remove=intersect(find(propedgeList(:,1)==vertexPair(2)),find(propedgeList(:,2)==vertexPair(1)));
                propedgeList(remove,:)=[];
            endfor
        endif    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Metropolis Hastings
    beta=log(p/(1-p));
    deltaH=beta*(size(propedgeList,1)-size(edgeList,1));
    
    if log(rand) <deltaH+deltaQ(size(edgeList,1),n,delta)
        edgeList=propedgeList;
    endif
    
    m = size(edgeList,1); % Count edges
    edgesum(m+1) = edgesum(m+1) + 1;
    time(l) = m;
    endfor
    
    % Normalise to create mass function
    edgedist = edgesum/sum(edgesum);

    function propdist = deltaQ(m,n,delta)
        %propdist=gammaln(0.5*n*(n-1)-m+1)-2*gammaln(0.5*n*(n-1)-m-delta+1)+gammaln(0.5*n*(n-1)-m-2*delta+1);
        propdist=log(0.5*n*(n-1)-m)-log(0.5*n*(n-1)-m-1);
%        if m<delta
%            propdist = 0;
%        endif 
    endfunction

endfunction
        
