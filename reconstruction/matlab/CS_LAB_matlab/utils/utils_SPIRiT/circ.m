function res = circ(scale,x,y,diameter,x_off,y_off)

%res = circ(scale,x,y,diameter,x_off,y_off)
diameter = diameter*scale;
x_off = x_off*scale;
y_off = y_off*scale;

res = 0.5*pi*diameter^2*exp(i*2*pi*(x_off*x+y_off*y)).*jinc(diameter*sqrt(x.^2+y.^2));
