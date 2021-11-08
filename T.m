function [ y ] = T( x )
%applies the teager operator, and returns the output
%y = T(x)

for i = 2:length(x)-1
    y(i) = x(i)^2 - x(i-1)*x(i+1);
    omega=acos((x(i-1)+x(i+1))/(2*x(i)));
    a=(sinc(omega)).^2;
    y(i)=y(i);%./a;
end
y=[0 y];
end

