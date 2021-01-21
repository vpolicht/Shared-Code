t = sym('t');
s = sym('s');
r = sym('r');
f = sym('f');
cexpcos = exp(-r*t^s)*cos(f*t);
cexpsin = exp(-r*t^s)*sin(f*t);
cexp = exp(-r*t^s);

%copy and paste the output of these

fprintf(1,'cexpsin = %s\n',func2str(matlabFunction(cexpsin)));
fprintf(1,'cexpcos = %s\n',func2str(matlabFunction(cexpcos)));
fprintf(1,'cexp = %s\n',func2str(matlabFunction(cexp)));

fprintf(1,'dcexpsinds = %s\n',func2str(matlabFunction(diff(cexpsin,s))));
fprintf(1,'dcexpsindr = %s\n',func2str(matlabFunction(diff(cexpsin,r))));
fprintf(1,'dcexpsindf = %s\n',func2str(matlabFunction(diff(cexpsin,f))));

fprintf(1,'dcexpcosds = %s\n',func2str(matlabFunction(diff(cexpcos,s))));
fprintf(1,'dcexpcosdr = %s\n',func2str(matlabFunction(diff(cexpcos,r))));
fprintf(1,'dcexpcosdf = %s\n',func2str(matlabFunction(diff(cexpcos,f))));

fprintf(1,'dcexpds = %s\n',func2str(matlabFunction(diff(cexp,s))));
fprintf(1,'dcexpdr = %s\n',func2str(matlabFunction(diff(cexp,r))));


