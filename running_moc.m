function[H,V,H_normal,H_transient,V_normal,V_transient] = running_moc(tsteps,theta,K_ut,K_ux,H0,V_initial,nsT,maxit,AsT,tol)

%% set global variables to be used in other functions
global Ap a dt dx Tcl g f mode n D V H

%% calculate H and V
% loop over time steps
%% normal version wihtout airchamber or surge tank
if mode == "normalVersion"
    for k = 1:2
    if k == 1
        friction = "normalFriction"
    else
        friction = "transientFriction"
    end
    
   
   
    
for j = 1:tsteps-1 
    
    % left boundary condition
    H(1,j+1) = H0; 
    V(1,j+1) = V(2,j)-g/a*(H(2,j)-H0)-f*dt/(2*D)*V(2,j)*abs(V(2,j));
    % calculation of pressure wave
    % loop over lenght of pipe
    if friction == "normalFriction"
        for i = 2:n
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)));
            V(i,j+1) = (V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)));
        end
    elseif friction == "transientFriction"
        for i = 2:n
            if j == 1
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-((f*dt*a)/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j)-V(i+1,j)+V(i+1,j)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
            
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j)+V(i+1,j)-V(i+1,j)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            
            else
            H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)))...
                        -(a/(2*g))*(K_ut*(1-theta)*(V(i-1,j)-V(i-1,j-1)-V(i+1,j)+V(i+1,j-1)))...
                        -((a*dt)/(2*g))*(K_ux*a*(sign(V(i-1,j))*abs((V(i,j)-V(i-1,j))/dx)-sign(V(i+1,j))*abs((V(i,j)-V(i+1,j))/dx)));
        
            V(i,j+1) = ((V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)))...
                        -dt/2*((K_ut/dt)*(-2*theta*V(i,j)+(1-theta)*(V(i-1,j)-V(i-1,j-1)+V(i+1,j)-V(i+1,j-1)))+((K_ux*a/dx))*...
                        (sign(V(i-1,j))*abs(V(i,j)-V(i-1,j))+sign(V(i+1,j))*abs(V(i,j)-V(i+1,j)))))*(1/(1+theta*K_ut));
            end
        end
    else
       disp("choose friction");
       return
    end
    
    % right boundary condition
    
    % the valve closing scheme has been already implemented in the function:
    % [V]=water_hammer_moc_boundary(type,timestep,initial fluid velocity)
    % This function only allows two types of closure: 10 for sinusoidal and 20 for linear. 
    % More on writing matlab functions in: https://de.mathworks.com/help/matlab/ref/function.html
    % More on calling functions in matlab: https://de.mathworks.com/help/matlab/learn_matlab/calling-functions.html
   help=(j-1)*dt/Tcl*pi();
    if help <= pi()
        V(n+1,j+1) = V_initial*0.5*(cos(help)+1);
    else
        V(n+1,j+1) = 0;
    end 
    H(n+1,j+1) = H(n,j)+a/g*(V(n,j)-V(n+1,j+1))-(a/g)*((f*dt)/(2*D))*V(n,j)*abs(V(n,j));
end
if k == 1
        H_normal = H;
        V_normal = V;
    else
        H_transient = H;
        V_transient = V;
    end
    end
    %% surge tank
elseif mode == "surgeTank"
    kcount=1;
    count = 1;
for j = 1:tsteps-1 
    
    % left boundary condition
    H(1,j+1) = H0; 
    V(1,j+1) = V(2,j)-g/a*(H(2,j)-H0)-f*dt/(2*D)*V(2,j)*abs(V(2,j));
   
    % calculation of pressure wave % loop over lenght of pipe
    for i = 2:nsT
        
        if i == nsT
        H(nsT,j+1) = H(nsT,j);%H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
        
        Hcheck = zeros (2,maxit);
        
            while true
        %Hcheck(1,count) = Hcheck(2,count-1);
        Hcheck(1,count) = H(nsT,j+1);
        V(nsT,j+1) = V(nsT-1,j)+g/a*(H(nsT-1,j)-H(nsT,j+1))-(f/(2*D))*(V(nsT-1,j)*abs(V(nsT-1,j))*dt);
        V(nsT+1,j+1)= V(nsT+2,j)-g/a*(H(nsT+1,j)-H(nsT,j+1))-(f/(2*D))*(V(nsT+2,j)*abs(V(nsT+2,j))*dt);
        H(nsT,j+1) = H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
                
        Hcheck(2,count) = H(nsT,j+1);
        
        diff = abs(Hcheck(1,count)-Hcheck(2,count));
        
          if (diff<=tol) || (count== maxit) 
            break;
          else
        count=count+1;
          end
       end   
              
        else
            
        H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)));
        V(i,j+1) = (V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)));
        
        end
    end
    for i = nsT+1:n
        %Note V is changed because it has n+2 values and H only n+1
        H(i,j+1) = (H (i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-(f*dt*a/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)));
        V(i+1,j+1) = (V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)));
        
    end
                
        % right boundary condition
    
    % the valve closing scheme has been already implemented in the function:
    % [V]=water_hammer_moc_boundary(type,timestep,initial fluid velocity)
    % This function only allows two types of closure: 10 for sinusoidal and
    % 20 for linear. 
    % More on writing matlab functions in: https://de.mathworks.com/help/matlab/ref/function.html
    % More on calling functions in matlab: https://de.mathworks.com/help/matlab/learn_matlab/calling-functions.html
   help=(j-1)*dt/Tcl*pi();
    if help <= pi()
        V(n+2,j+1) = V_initial*0.5*(cos(help)+1);
    else
        V(n+2,j+1) = 0;
    end 
     H(n+2,j+1) = H(n+1,j)+a/g*(V(n+1,j)-V(n+2,j+1))-(a/g)*((f*dt)/(2*D))*V(n+1,j)*abs(V(n+1,j));
end
 
%% air chamber  

elseif mode == "airChamber"
kcount=1;
count = 1;
for j = 1:tsteps-1 
    
    % left boundary condition
    H(1,j+1) = H0; 
    V(1,j+1) = V(2,j)-g/a*(H(2,j)-H0)-f*dt/(2*D)*V(2,j)*abs(V(2,j));
   
    % calculation of pressure wave % loop over lenght of pipe
    for i = 2:nsT
        
        if i == nsT
        H(nsT,j+1) = H(nsT,j);%H(nsT,j)+(Ap*dt)/(2*AsT)*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
        
        Hcheck = zeros (2,maxit);
        
            while true
        %Hcheck(1,count) = Hcheck(2,count-1);
        % Ask about the formular of H(nst,j+1), m etc. !
        %prepressurise the air cussions for intial setup set up intial
        %pressure atmospheric plus hydrostatic pressure times area of the
        %air chaber, politrophic coefficient is m here we use 1.2, turned around A is
        %volumen of air cushion.
        %
        Hcheck(1,count) = H(nsT,j+1);
        V(nsT,j+1) = V(nsT-1,j)+g/a*(H(nsT-1,j)-H(nsT,j+1))-(f/(2*D))*(V(nsT-1,j)*abs(V(nsT-1,j))*dt);
        V(nsT+1,j+1)= V(nsT+2,j)-g/a*(H(nsT+1,j)-H(nsT,j+1))-(f/(2*D))*(V(nsT+2,j)*abs(V(nsT+2,j))*dt);
        %H(nsT,j+1) = H(nsT,j)+((Ap*dt)/2)*((H(i,j)))*((V(nsT,j+1)-V(nsT+1,j+1))+(V(nsT,j)-V(nsT+1,j))); % not sure about which velocity to take
                
        Hcheck(2,count) = H(nsT,j+1);
        
        diff = abs(Hcheck(1,count)-Hcheck(2,count));
        
          if (diff<=tol) || (count== maxit) 
            break;
          else
        count=count+1;
             end
            end   
              
        else
            
        H(i,j+1) = (H(i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+1,j)-V(i-1,j))-(f*dt*a/(4*g*D))*(V(i-1,j)*abs(V(i-1,j))-V(i+1,j)*abs(V(i+1,j)));
        V(i,j+1) = (V(i-1,j)+V(i+1,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i-1,j)*abs(V(i-1,j))+V(i+1,j)*abs(V(i+1,j)));
        
        end
    end
    for i = nsT+1:n
        %Note V is changed because it has n+2 values and H only n+1
        H(i,j+1) = (H (i-1,j)+H(i+1,j))/2-(a/(2*g))*(V(i+2,j)-V(i,j))-(f*dt*a/(4*g*D))*(V(i,j)*abs(V(i,j))-V(i+2,j)*abs(V(i+2,j)));
        V(i+1,j+1) = (V(i,j)+V(i+2,j))/2-g/(2*a)*(H(i+1,j)-H(i-1,j))-(f*dt/(4*D))*(V(i,j)*abs(V(i,j))+V(i+2,j)*abs(V(i+2,j)));
        
    end
                
        % right boundary condition
    
    % the valve closing scheme has been already implemented in the function:
    % [V]=water_hammer_moc_boundary(type,timestep,initial fluid velocity)
    % This function only allows two types of closure: 10 for sinusoidal and
    % 20 for linear. 
    % More on writing matlab functions in: https://de.mathworks.com/help/matlab/ref/function.html
    % More on calling functions in matlab: https://de.mathworks.com/help/matlab/learn_matlab/calling-functions.html
   help=(j-1)*dt/Tcl*pi();
    if help <= pi()
        V(n+2,j+1) = V_initial*0.5*(cos(help)+1);
    else
        V(n+2,j+1) = 0;
    end 
     H(n+2,j+1) = H(n+1,j)+a/g*(V(n+1,j)-V(n+2,j+1))-(a/g)*((f*dt)/(2*D))*V(n+1,j)*abs(V(n+1,j));
 end
else 
    disp('no mode')
    return  
end
end