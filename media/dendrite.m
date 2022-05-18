% dendrite.m
%
% M-file to hold and process matrices pertaining to the
% dendrite electrical modeling problem. It then uses
% the Backward-Euler method with average parameter values for
% the cell to obtain the cell potentials.  Finally, a plot
% of potentials is issued.
%
% Doug Daniels
% 06-07-01 

% average parameter values
dt = 1;      % time step for backward-euler method
ell = .1;        % fiber length (cm)
a = 0.025;      % fiber radius (cm)
rho_i = 100;     % cytoplasmic resistivity (Ohm cm)
rho_m = 1000;   % membrane resistivity (Ohm cm^2)
E_m = 0;      % resting potential mV
N = 8;         % number of compartments
c = .1;		% micro farad
R_i = rho_i*(ell/N)/(pi*a*a);   % axial resistance (Ohm)
R_m = rho_m/(2*pi*a*ell/N);     % lateral resistance (Ohm)

A_cb = 4*pi*(.05)^2;
R_cb = rho_m/A_cb;
C_cb = c*A_cb;
C_m = c*(2*pi*a*ell/N);

% node-edge adjacency matrix
A = [	-1 0 0 0 0 0 0 0 0 0;
     	-1 0 0 0 0 0 0 0 0 0;
 	-1 1 0 0 0 0 0 0 0 0;
        0 -1 0 0 0 0 0 0 0 0;
        0 -1 0 0 0 0 0 0 0 0; 
        0 -1 1 0 0 0 0 0 0 0;
        0 0 -1 0 0 0 0 0 0 0;
        0 0 -1 0 0 0 0 0 0 0;
      	0 0 -1 1 0 0 0 0 0 0;
  	0 0 0 1 -1 0 0 0 0 0;
	0 0 0 0 -1 0 0 0 0 0;
	0 0 0 0 -1 0 0 0 0 0;
	0 0 0 0 -1 1 0 0 0 0;
	0 0 0 0 0 -1 0 0 0 0;
	0 0 0 0 0 -1 0 0 0 0;
	0 0 0 0 0 -1 1 0 0 0;
	0 0 0 0 0 0 -1 0 0 0;
	0 0 0 0 0 0 -1 0 0 0;
	0 0 0 -1 0 0 0 1 0 0;
	0 0 0 0 0 0 0 -1 0 0;
	0 0 0 0 0 0 0 -1 0 0;
	0 0 0 0 0 0 0 -1 1 0;
	0 0 0 0 0 0 0 0 -1 0;
	0 0 0 0 0 0 0 0 -1 0;
	0 0 0 0 0 0 0 0 -1 1;
	0 0 0 0 0 0 0 0 0 -1;
	0 0 0 0 0 0 0 0 0 -1 ];

% vector of batteries
b = - E_m .* [0; 1; 0; 0; 1; 0; 0; 1; 0; 0; 0; 1; 0; 0; 1; 0; 0;
     1; 0; 0; 1; 0; 0; 1; 0; 0; 1 ];

% diagonal elements of C  
c_diags = [ C_cb; 0; 0; C_m; 0; 0; C_m; 0; 0; 0; C_m; 0; 0; 
	    C_m; 0; 0; C_m; 0; 0; C_m; 0; 0; C_m; 0; 0; C_m; 
	    0  ];

% C matrix
C = diag(c_diags);

% diagonal elements of G before subsequent transformation (see below)
% G_diags = [ 0; 1/R_cb; 1/R_i; 0; 1/R_m; 1/R_i; 0; 
%	    1/R_m; 1/R_i; 1/R_i; 0; 1/R_m; 1/R_i; 
%	    0; 1/R_m; 1/R_i; 0; 1/R_m; 1/R_i; 0; 
%	    1/R_m; 1/R_i; 0; 1/R_m; 1/R_i; 0; 1/R_m ];


% here we have replaced ( 1 / R_(foo) ) with ( G_(foo) ) 
% in the G_diags
G_cb = 1 / R_cb;
G_i = 1 / R_i;
G_m = 1 / R_m;

G_foo_diags = [ 0; G_cb; G_i; 0; G_m; G_i; 0;
	 	    G_m; G_i; G_i; 0; G_m; G_i;
		    0; G_m; G_i; 0; G_m; G_i; 0;
		    G_m; G_i; 0; G_m; G_i; 0; G_m ];


% G matrix, modified as described above
G = diag(G_foo_diags);

A_trans_C_A = A' * C * A;

A_trans_G_A = A' * G * A;

A_trans_G_b = A' * G * b;

% Backward-Euler calculations
D = A_trans_C_A;
E = - A_trans_G_A;

S = D - ( dt .* E );

x(:,1) = zeros(10,1);
t(1) = 0;

for j = 2:400/dt
   
   t(j) = (j-1) * dt;
   g = A_trans_G_b + [t(j)^3 * exp( -t(j) / 6 ) / 10000/C_cb 0 0 0 0 0 0 0 0 0]';
   x(:,j) = S \ (D * x(:,j-1) + ( g .* dt) );
   
end


hold off;
plot(t,x(1,:),'m-');
hold on;
plot(t,x(2,:),'g-');
plot(t,x(3,:),'c-');
plot(t,x(4,:),'r-');
plot(t,x(5,:),'b-');
plot(t,x(6,:),'k-');
plot(t,x(7,:),'m:');
plot(t,x(8,:),'g:');
plot(t,x(9,:),'r:');
plot(t,x(10,:),'b:');

legend('x_1', 'x_2', 'x_3', 'x_4', 'x_5', 'x_6', 'x_7', 'x_8', 'x_9', 'x_1_0');

xlabel('Time (msec)','fontsize',16);
ylabel('Potential (mV)','fontsize',16);
