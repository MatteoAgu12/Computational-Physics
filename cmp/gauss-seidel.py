import numpy as np

"""Calculates numerically the solution of a 2-variable system."""

# insert input parameters
a = float(input("a = "))
b = float(input("b = "))
c = float(input("c = "))
d = float(input("d = "))
e = float(input("e = "))
f = float(input("f = "))

matrix = np.array([a, b, c, d, e, f])

x = float(input("\nInitial guess for x = "))
y = float(input("Initial guess for y = "))

digit = int(input("\nInsert the number of significant digits you want: "))

# print system to solve
print("\n" + str(a) + "x + " + str(b) + "y = " + str(c) + "\n" +
      str(d) + "x + " + str(e) + "y = " + str(f) + "\n")

# define function for solving one iteration
def numerical_solver(x_0, y_0, matrix):
    try:
        x_1 = (matrix[2] - matrix[1] * y_0) / matrix[0]
        y_1 = (matrix[5] - matrix[3] * x_1) / matrix[4]
    except ZeroDivisionError:
        print("Encountered division by zero.")
        return np.array([x_0, y_0, float('inf'), float('inf')])
    
    return np.array([x_0, y_0, x_1, y_1])

# input for main loop
sol = numerical_solver(x, y, matrix)
print("Initial guess for solution:", sol)

# number of iterations
n = 0

# main loop
while (abs(sol[2] - sol[0]) > 10 ** -(digit+1)) and (abs(sol[3] - sol[1]) > 10 ** -(digit+1)):
    x_0 = sol[2]
    y_0 = sol[3]    
    n += 1
    sol = numerical_solver(x_0, y_0, matrix)
    print("At iteration " + str(n) + ", sol is: " + str(sol))

# print result
print("Final result: ", sol)
