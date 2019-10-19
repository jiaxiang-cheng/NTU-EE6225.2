num=[-0.092];
den=[130,1];
G=tf(num,den)
G.iodelay=16;
margin(G);