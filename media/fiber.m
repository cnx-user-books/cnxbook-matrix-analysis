%
% fiber.m
% Generate the A matrix for the 3 x 3 "tissue" of figure 2.3
% proceed to construct, solve, and plot A^TKA x = f
%

% The numbering of nodes and bars is the one shown in figure 2.3.
% number of nodes (nnodes) and bars/fibers (nbars) of the truss
nnodes = 9;
nbars  = 20;

data = [      % one row of data for each fiber, the
1  4  -pi/2   % first two columns are node numbers
1  5  -pi/4   % while the third is the fiber angle
1  2  0
2  4  -3*pi/4
2  5  -pi/2
2  6  -pi/4
2  3  0
3  5  -3*pi/4
3  6  -pi/2
4  7  -pi/2
4  8  -pi/4
4  5  0
5  7  -3*pi/4
5  8  -pi/2
5  9  -pi/4
5  6  0
6  8  -3*pi/4
6  9  -pi/2
7  8  0
8  9  0];

% xy-coordinates of the nodes
% one row for each node, the first column are the x-coordinates
% and the second column are the y coordinates.
xy = [  
10 10
20 10
30 10
10 20
20 20
30 20
10 30
20 30
30 30 ];



% stiffnesses  (all ones)
k = ones(nbars,1);

% external forces
f = [-1 1 0 1 1 1 -1 0 0 0 1 0 -1 -1 0 -1 1 -1]';


% initialize the A matrix and compute its entries based on
% the formula (2.9) in the class notes

A = zeros(nbars,2*nnodes);	

for j=1:nbars
    m     = data(j,1);  % left node of bar j
    n     = data(j,2);  % right node of bar j
    theta = data(j,3);  % angle of bar j

    A(j,2*n-1) = cos(theta);   % direct implementation of general formula
    A(j,2*m-1) = -cos(theta);
    A(j,2*n)   = sin(theta);
    A(j,2*m)   = -sin(theta);

end

% view the matrix A
% this A is too big to view conveniently. Therefore
% we use the symbolic toolbox of Matlab, round the
% zeros and ones and replace the bulky .7071 with `a'
syms a cA	
for i=1:nbars,
   for j=1:2*nnodes,  
     if abs(abs(A(i,j))-1/sqrt(2)) > .1,
        cA(i,j) = round(A(i,j));
     else
        if A(i,j) > .5,
           cA(i,j) = a;
        else
           cA(i,j) = -a;
        end
     end
   end
end
disp(cA)		
 
 
% form the stiffness matrix
S  = A'*diag(k)*A;

% compute the pseudo inverse of the stiffness matrix
pS = pinv(S);

% solve A'*diag(k)*A x = f
x  = pS*f;

% draw the nodes of the original and deformed specimens

% a reference circle
t  = 0:.01:1;		% a reference circle
mx = cos(2*pi*t);
my = sin(2*pi*t);

for i = 1:nnodes
    % node i in the undeformed truss
    plot(xy(i,1)+mx,xy(i,2)+my); hold on
    % node i in the deformed truss
    plot(x(2*i-1)+xy(i,1)+mx,-x(2*i)+xy(i,2)+my,'--')
    % node number 
    text(xy(i,1),xy(i,2),int2str(i));
end

axis equal
axis off
hold off

