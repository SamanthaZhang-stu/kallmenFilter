function Inew=rsmp(pik,N)

akkupik=cumsum(pik);

Inew=zeros(N,1);
y=rand(1)/N;
j=1;
for i=1:N        
    while(akkupik(j)<y)%flush all values that are below the current level of integrated propability
        j=j+1;        
    end
    Inew(i)=j;
    y=y+1/N;
end
