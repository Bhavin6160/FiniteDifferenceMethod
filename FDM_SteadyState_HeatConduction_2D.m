% This programme is to solve Laplace Equation by
% using various finite difference iteration methods.
% Written by Bhavin6160
clear;
repeat = 1;
while (repeat == 1)
clc;
clear;

fprintf("Please enter the choice to use the following methods of iteration to solve laplace equation.\n\n");
fprintf("Enter 1 for Jacobi Iteration Method\n");
fprintf("Enter 2 for Point Gauss Siedel Iteration Method\n");
fprintf("Enter 3 for Point Succesive Over & Under Relaxation Method\n\n");

m = input("Enter Your Choice:   ");

if m<=0 || m ~= fix(m) || ~ismember(m,[1 2 3])
    fprintf("\nPlease enter 1,2 or 3.\n\n")
    return;
else 
    l = 1.0;   % length of plate in x-direction 
    w = 1.0;   % length of plate in y-direction
    n = input("Enter the number of grid points required in x direction:  ");
    % Checking for the value of x-direction grid point
    if n~=fix(n) || n<=0
        fprintf("Please enter the positive integer value...!!\n\n");
        return;
    end

    o = input("Enter the number of grid points required in y direction:  ");
    % Checking for the value of y-direction grid point 
    if o~=fix(o) || o<=0
        fprintf("Please enter the positive integer value...!!\n\n");
        return;
    end



    x = linspace(0,l,n);   % Dividing Length in equal steps using user input
    y = linspace(0,w,o);   % Dividing width in equal steps using user input
    T = zeros(o,n);        % Creating n*o
    T(1,:)= input("Enter the temperature of the bottom of the plate:  ");           % Bottom side temperature
    T(o,:)= input("Enter the temperature of the top of the plate:  ");              % Top side temperature
    T(2:end-1,1) = input("Enter the temperature of the leftside of the plate:  ");  % Left side temperature
    T(2:end-1,n) = input("Enter the temperature of the rightside of the plate:  "); % Right side temperature
    dx=l/n;   % step size in x direction
    dy=w/o;   % step size in y direction
    beta = dx/dy;
    error_req = input("Enter the accuracy required in iteration:  ");  % Accuracy required for consecutive iteration
    fprintf("\n");
    error = 1;
    k=0;
end

switch m

    case 1
        % Jacobi Iteration Method
        fprintf("You have selected Jacobi Iteration Method.\n\n")    
        tic
        while error>error_req
            T0=T;
            for i=2:o-1
                for j=2:n-1
                    T(i,j)=(0.5/(1+beta^2)) * (T0(i+1,j) + T0(i-1,j) + (beta^2)*( T0(i,j+1) + T0(i,j-1)));            
                end
            end
            error = max(max(abs(T0-T)));
            k=k+1;
       end
       toc
       Toutput = T;
       fprintf("\nNumber of iterations to achieve given accuracy are: %f \n\n",k);
       fprintf("\n");
       subplot(2,1,1);
       contour(x,y,T);
       colormap hot;
       colorbar;
       title('Temrature Variation on Plate(Steady State) By Jacobi Iteration');
       xlabel('x');
       ylabel('y');
       subplot(2,1,2);
       pcolor(x,y,T);
       shading interp;
       colormap hot;
       colorbar;
       title('Temrature Variation on Plate(Steady State) By Jacobi Iteration');
       xlabel('x');
       ylabel('y');

    case 2
       % Point Gauss Siedel Iteration Method
       fprintf("You have selected Point Gauss Siedel Iteration Method.\n\n") 
       tic
       while error>error_req
           T0=T;
           for i=2:o-1
               for j=2:n-1
                   T(i,j)=(0.5/(1+beta^2)) * (T0(i+1,j) + T(i-1,j) + (beta^2)*( T0(i,j+1) + T(i,j-1)));
               end
           end
           error = max(max(abs(T0-T)));
           k=k+1;
      end
      toc
      fprintf("\nNumber of iterations to achieve given accuracy are: %f \n\n",k);
      subplot(2,1,1);
      contour(x,y,T);
      colormap hot;
      colorbar;
      title('Temrature Variation on Plate(Steady State) By Point Gauss Siedel Iteration');
      xlabel('x');
      ylabel('y');
      subplot(2,1,2);
      pcolor(x,y,T);
      shading interp;
      colormap hot;
      colorbar;
      title('Temrature Variation on Plate(Steady State) By Point Gauss Siedel Iteration');
      xlabel('x');
      ylabel('y');

    case 3
     % Point Succesive Over & Under Relaxation 
     fprintf("Would you like to enter the relaxation parameter...??\n");
     relax_para = input("Enter 1 for YES and 0 for No:  ");

     if relax_para~=fix(relax_para) || relax_para<0
        fprintf("Please enter the 0 or 1...!!\n\n");
        return;
     end
 
      
     if relax_para == 1
         omega = input("Please enter the value of relaxation parameter:  ");
         fprintf("\n");
         fprintf("Your entered relaxation parameter is %f",omega);
         fprintf("\n");
         if (1<omega && omega<2)
             fprintf("Your have selected Point Successive Over Relaxation method.\n");
         elseif (0<omega && omega<1)
             fprintf("Your have selected Point Successive Under Ralaxation method.\n");
         elseif (omega == 1)
             fprintf("Your have selected Point Gauss Siedel Iteration.\n");
         else
             fprintf("Please Enter the value of relaxation in between 0 and 2");
         return;
         end

     else  
        % Calculation for the optimum value of relaxation parameter 
        IM = n;
        JM = o;
        a = ( (cos(pi/(IM-1)) + (beta^2) * cos(pi/(JM-1))) / (1 + (beta^2))  )^2;
        omega =( 2 - (1-a)^0.5 )/a;
        fprintf("Your optimum calculated relaxtion parameter is %f .",omega);
        fprintf("\n");
     end

     tic
     while error>error_req
         T0=T;
         for i=2:o-1
             for j=2:n-1
                 T(i,j)= (1-omega)*T0(i,j) + ((0.5*omega)/(1+beta^2)) * (T0(i+1,j) + T(i-1,j) + (beta^2)*( T0(i,j+1) + T(i,j-1)));
             end
         end
         error = max(max(abs(T0-T)));
         k=k+1;
     end
     toc
     fprintf("\nNumber of iterations to achieve given accuracy are: %f \n\n",k);
     subplot(2,1,1);
     contour(x,y,T);
     colormap hot;
     colorbar;
     title("Temrature Variation on Plate (Steady State) Point Successive Relaxation");
     xlabel('x');
     ylabel('y');
     subplot(2,1,2);
     pcolor(x,y,T);
     shading interp;
     colormap hot;
     colorbar;
     title("Temrature Variation on Plate (Steady State)Point Successive Relaxation");
     xlabel('x');
     ylabel('y');

    otherwise
        fprintf("Please Enter a valid choice.\n");
end

repeat = input("Do you like to continue..?? Enter 1 for Yes and 0 for No :  ");  % If user want to run the program once again

end