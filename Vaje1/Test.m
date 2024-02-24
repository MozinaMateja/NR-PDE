x = 0;
y = 0;
u = @(x,y) exp(-x.^2 - y.^2);
u1 = @(x,y) exp(x).*(x.*sin(y) - y.*cos(y));
u2 = @(x,y) exp(x).*sin(y);
h = 0.05;
LaplaceovOperator(u,x,y,h);
BiharmonicniOperator(u,x,y,h);
deltau = @(x,y) 4.*(x.^2 + y.^2 - 1).*exp(-x.^2 - y.^2);  % u_xx + u_yy
deltadeltau = @(x,y) 4.*(4*x.^4 - 12*x.^2 + 4*y.^4 - 12*y.^2 + 6 ...
    + 2*(2*y.^2 - 1).*(2*x.^2 - 1)).*exp(-x.^2 - y.^2); % u_xxxx + 2*u_xxyy + u_yyyy

abs(deltau(x,y) - LaplaceovOperator(u,x,y,h)); % razlika med delta u in priblizkom v (0,0)
abs(deltadeltau(x,y) - BiharmonicniOperator(u,x,y,h)); % razlika med dvojnim delta u in priblizkom v (0,0)


%funkcija napake v odvisnosti od h
napaka1 = @(h) abs(deltau(0,0) - LaplaceovOperator(u,x,y,h));  
napaka2 = @(h) abs(deltadeltau(0,0) - BiharmonicniOperator(u,x,y,h));


% grafa napake
h1 = 0.005:0.0001:1;
y1 = arrayfun(napaka1,h1);
y2 = arrayfun(napaka2,h1);
plot(h1,y1);
hold on
plot(h1,y2);
hold off


% grafa Laplaceovega operatorja in biharmonicnega operatorja
hn = 0.05;
t = -1:0.05:1;
[X,Y] = meshgrid(t,t);
Z1 = LaplaceovOperator(u,X,Y,hn);
Z2 = BiharmonicniOperator(u,X,Y,hn);
surf(X,Y,Z1);
surf(X,Y,Z2);


% grafa napake
Z3 = abs(deltau(X,Y) - LaplaceovOperator(u,X,Y,hn));
surf(X,Y,Z3);
Z4 = abs(deltadeltau(X,Y) - BiharmonicniOperator(u,X,Y,hn));
surf(X,Y,Z4);
