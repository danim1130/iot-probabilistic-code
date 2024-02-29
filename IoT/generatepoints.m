%Simulation window parameters
xMin=0;xMax=10;
yMin=0;yMax=10;
xDelta=xMax-xMin;yDelta=yMax-yMin; %rectangle dimensions
areaTotal=xDelta*yDelta
 
%Point process parameters
lambda=2; %intensity (ie mean density) of the Poisson process
 
%Simulate Poisson point process
numbPoints=poissrnd(areaTotal*lambda);%Poisson number of points
numbPoints = 10;

xx=xDelta*(rand(numbPoints,1))+xMin;%x coordinates of Poisson points
yy=xDelta*(rand(numbPoints,1))+yMin;%y coordinates of Poisson points

%Plotting
scatter(xx,yy)
xlabel('x');ylabel('y');

A=[xx, yy]
squareform(pdist(A))