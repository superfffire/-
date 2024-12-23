function Pasokh=bsm(y,a,b,eps)
syms x
X=(a+b)/2;
fx=subs(y,X);
E=abs(fx);
while E>eps
    fa=subs(y,a);
    fb=subs(y,b);
    fx=subs(y,X);
    if fa*fx<0
        b=X;
    else
        a=X;
    end
    X=(a+b)/2;
    fx=subs(y,X);
    E=abs(fx);
    %disp([' X = ' , num2str(X) , '       f(x) = ' , num2str(fx)])
end
Pasokh=X;