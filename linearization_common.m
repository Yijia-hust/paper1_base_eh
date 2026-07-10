function [ A,B,C ] = linearization_common(c2,c1,c0,xmax,xmin,Q)
m=size(xmax,1);
interval = (xmax-xmin)./Q;
A = zeros(m,Q);
B = zeros(m,Q);
C = c2.*xmin.^2+c1.*xmin+c0;
X = zeros(m,Q+1);
Y = zeros(m,Q+1);
for i=1:m
    if (c2(i)==0)&&(c1(i)==0)&&(c0(i)==0)
        Y(i,:)=0;
    else
        for j=1:Q+1
            X(i,j)=xmin(i)+(j-1)*interval(i);
            Y(i,j)=c2(i).*X(i,j).^2+c1(i).*X(i,j)+c0(i);
        end
    end
end
for i=1:m
    if (c2(i)==0)&&(c1(i)==0)&&(c0(i)==0)
        A(i,:)=0;
        B(i,:)=0;
    else
        for j=1:Q
            A(i,j)=(Y(i,j+1)-Y(i,j))/interval(i);
            B(i,j)=Y(i,j)-A(i,j)*X(i,j);%-C(i)
        end
    end
end


end