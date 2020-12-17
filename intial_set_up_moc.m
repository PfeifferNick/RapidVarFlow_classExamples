function[H,V] = intial_set_up_moc(tsteps,H0,L,zeta,Av,V_initial)

%% set global variables to be used in other functions
global Ap a dt dx Tcl g f mode n D

%initialization 
if mode == "normalVersion"
H = zeros(n+1,tsteps);
V = zeros(n+1,tsteps);
 
elseif mode == "surgeTank"
H = zeros(n+1,tsteps);
V = zeros(n+2,tsteps);
  
elseif mode == "airChamber"
H = zeros(n+1,tsteps);
V = zeros(n+2,tsteps);
else 
    disp('no mode')
    return   
end



%% set initial conditions
% loop over lenght of pipe
if mode == "normalVersion"
for i = 1:n+1
    if i==1
        H(i,1) = H0;
    else
   H(i,1) =  H0-(i-1)*dx*f/D*V_initial^2/(2*g);    %H(i,1) =  % remember that the initial conditions come from steady
                    %state and from fluid mechanics you know how to calculate the
        %piezometric head 
    end
   V(i,1) = V_initial;
%   if i == left
%      H()
%   end
end
 
elseif mode == "surgeTank"
  for i = 1:n+1
    if i==1
        H(i,1) = H0;
    else
   H(i,1) =  H0-(i-1)*dx*f/D*V_initial^2/(2*g);    %H(i,1) =  % remember that the initial conditions come from steady
                    %state and from fluid mechanics you know how to calculate the
        %piezometric head 
    end
  end
  for i = 1:n+2
   V(i,1) = V_initial;
  end
  
   
elseif mode == "airChamber"
  for i = 1:n+1
    if i==1
        H(i,1) = H0;
    else
   H(i,1) =  H0-(i-1)*dx*f/D*V_initial^2/(2*g);    %H(i,1) =  % remember that the initial conditions come from steady
                    %state and from fluid mechanics you know how to calculate the
        %piezometric head 
    end
  end
  for i = 1:n+2
   V(i,1) = V_initial;
  end
end

end
