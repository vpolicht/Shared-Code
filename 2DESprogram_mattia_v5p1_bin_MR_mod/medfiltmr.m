function [Y] = medfiltmr( X,n )
dim=size(X);
Y=zeros(dim(1),dim(2));
Y=X;

if mod(n,2)==0
    %smooth della matrice se l'indice è pari
    
      for i=1:dim(2)
      for k=(n/2)+1:1:dim(1)-(n/2)-1
      Y(k,i)=(sum(X(k-(n/2):k+(n/2)-1)))/n;
      end
      end
else
    %%smooth della matrice se l'indice è dispari
    for i=1:dim(2)
    for k=(n/2)+0.5:1:dim(1)-(n/2)-0.5
    Y(k,i)=(sum(X(k-(n-1)/2:k+(n-1)/2)))/n;
    end
     end
end

end

