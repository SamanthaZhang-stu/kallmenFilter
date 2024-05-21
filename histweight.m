function [plotx,ploty]=histweight(x,w,N,edges)


[x,I]=sort(x);
w=w(I);

if(nargin<4)
binwidth=(x(end)-x(1))/N;
centers=linspace(x(1),x(end),N);
lim=linspace(x(1)-binwidth/2,x(end)+binwidth/2,N+1);
else
    binwidth=(edges(2)-edges(1))/N;
centers=linspace(edges(1),edges(2),N);
lim=linspace(edges(1)-binwidth/2,edges(2)+binwidth/2,N+1);
lim(end)=x(end);%to prevent out of bounds on the next few rows
end
y=zeros(1,N);

k=2;%points at the upper bin limit
for i=1:length(x)
    if(x(i)>lim(k))
        %ok, we must jump to the next bin
        k=k+1;
    else
        %we are still in the right bin
    end
        y(k-1)=y(k-1)+w(i); %add the weight to the bin               
end

%scale everything so that it is a real pdf
area=sum(y)*binwidth;
plotx=x;
ploty=y;

