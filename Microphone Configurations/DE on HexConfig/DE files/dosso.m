function z = dosso(X,dummy1,dummy2)

x1 = X(1); x2 = X(2);
x3 = X(3); x4 = X(4);
x5 = X(5); x6 = X(6);

z = 4.8 + x1^2 + 5*x2^2 + 0.1*x3^2 + 0.05*x4^2 + x5^2 + x6^2 ...
    - 0.3*cos(1*pi*(x1-x2)) - 1.4*cos(1*pi*(x1+x2)) ...
    - 0.5*cos(10*pi*(0.05*x4 - 0.1*x3)) ...
    - 1.0*cos(10*pi*(0.05*x4+0.1*x3)) - 0.25*cos(2*pi*(x5-x6)) ...
    - 1.35*cos(2*pi*(x5+x6));
z = sqrt(z/20);

% end of program