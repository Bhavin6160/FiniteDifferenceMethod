# This Program is to solve laplace equation 2D steady state heat conduction using various FDM Methods.
# written by Bhavin6160, just for the fun.
import numpy as np
import sys
import math
from matplotlib import pyplot as plt

choice = float(input("\n\nThis program is written to solve Laplace Equation for 2D state heat conduction using various FDM Methods.\n\n"
      "Select any of the following method for the solution :\n"
      "1. Jacobi Iteration Method.\n"
      "2. Point Gauss Siedel Iteration Method\n"
      "3. Point Succesive Over & Under Relaxation Method\n"
      "Enter your choice  :  "
      ))

if (choice !=1) and (choice !=2) and (choice !=3):
    print("Please Enter the 1,2 or 3.")
    sys.exit("Invalid Choice")

length = float(input("Enter the length of the plate in X-direction : "))
b = float(input("Enter the length of the plate in Y-direction : "))

if (length <= 0) or (b <= 0):
    print("Please Enter the positive integer as length.\n")
    sys.exit("Invalid Input.")

x = int(input("Enter the number of grid points required in X-direction : "))
y = int(input("Enter the number of grid points required in y-direction : "))

if (x or y) <=0:
    print("Please Enter the positive integer as number of grid points.\n")
    sys.exit("Invalid Input.")

T = np.zeros([y, x])
bottom = float(input("Enter the temperature of the top side of the plate : "))
top = float(input("Enter the temperature of the bottom side of the plate : "))
right = float(input("Enter the temperature of the right side of the plate : "))
left = float(input("Enter the temperature of the left side of the plate : "))

T[[y-1], :] = bottom
T[[0], :] = top
T[:, 0] = left
T[:, [x-1]] = right

dx = length/x
dy = b/y
beta = (dx/dy)
error = 1
error_req = float(input("Enter the accuracy required in the solution : "))

xx = np.linspace(0, length, num=x, endpoint=True)
yy = np.linspace(0, b, num=y, endpoint=True)

# print(xx)
# print(yy)

if choice == 1:      # Jacobi Iteration Method
    while error > error_req:
        T1 = np.copy(T)
        for i in range(1, (y-1)):
            for j in range(1, (x-1)):
                T[i][j] = 0.5 / (1 + (beta ** 2)) * ((T1[i+1][j] + T1[i-1][j] + (beta ** 2) * (T1[i][j+1] + T1[i][j-1])))

        T_out = np.subtract(T1, T)
        error = np.amax(np.absolute(T_out))

elif choice == 2:    # Point Gauss Siedel Iteration Method
    while error > error_req:
        T1 = np.copy(T)
        for i in range(1, (y-1)):
            for j in range(1, (x-1)):
                T[i][j] = 0.5 / (1 + (beta ** 2)) * ((T1[i+1][j] + T[i-1][j] + (beta ** 2) * (T1[i][j+1] + T[i][j-1])))

        T_out = np.subtract(T1, T)
        error = np.amax(np.absolute(T_out))

elif choice == 3:   # Point Succesive Over & Under Relaxation Method

    print("Would you like to enter relaxation parameter..??")
    para = float(input("Enter 1 for YES and 0 for NO : "))

    if(para != int(para) or (para<0) or (para>1)):
        print("Please enter 1 or 0 only.")
        sys.exit("You have entered an invalid number")

    if (para == 1):
        omega = float(input("Please Enter the relaxation parameter : "))

        if(0 < omega < 1):
            print("You have chosen a Point Successive Under Relaxation method.")
        elif(1 < omega < 2):
            print("You have chosen a Point Successive Over Relaxation method.")
        elif( omega==1 ):
            print("You have chosen a Point Gauss Siedel Iteration method.")
        else:
            print("Please Enter the Relaxation Parameter in between 0 and 1.")
            sys.exit("Wrong Relaxation Parameter value.")
    else:
        a = (((math.cos(math.pi/(x-1))) + (beta ** 2) * (math.cos(math.pi/(y-1)))) / (1 + (beta ** 2))) ** 2
        omega = (2 - math.sqrt(1-a))/a
        print("Your calculated optimum relaxation parameter is ", omega)

    while error > error_req:
        T1 = np.copy(T)
        for i in range(1, (y - 1)):
            for j in range(1, (x - 1)):
                T[i][j] = ((1-omega) * T1[i][j]) + (((0.5*omega) / (1 + (beta ** 2))) * (T1[i+1][j] + T[i-1][j] + (beta ** 2) * (T1[i][j+1] + T[i][j-1])))

            T_out = np.subtract(T1, T)
            error = np.amax(np.absolute(T_out))

print(T)
fig1, ax1 = plt.subplots(1, 1, figsize=(6, 6))
# fig.tight_layout()
ax1.contour(xx, yy, T)
ax1.set_title('Temperature Distribution contour on the plate')
ax1.set_xlabel(r'$X\longrightarrow$')
ax1.set_ylabel(r'$Y\longrightarrow$')


fig2, ax2 = plt.subplots(1, 1, figsize=(6, 6))
# fig.tight_layout()
smooth = ax2.pcolormesh(xx, yy, T, shading='gouraud', vmin=T.min(), vmax=T.max())
ax2.set_title('Temperature Distribution on the plate')
ax2.set_xlabel(r'$X\longrightarrow$')
ax2.set_ylabel(r'$Y\longrightarrow$')
fig2.colorbar(smooth)
plt.show()
