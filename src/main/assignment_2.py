import numpy as np
from decimal import *

def nevilles_method(x_points, y_points, x):
    matrix = np.zeros((3, 3))
    res = 0

    for counter, row in enumerate(matrix):
        row[0] = y_points[counter]
   
    num_of_points = 3
    
    for i in range(1, num_of_points):
        for j in range(1, i+1):
            first_multiplication = float(x - x_points[i-j]) * matrix[i][j-1]
            second_multiplication = float(x - x_points[i]) * matrix[i-1][j-1]
            denominator = float(x_points[i] - x_points[i-j])
            # this is the value that we will find in the matrix
            coefficient = float(first_multiplication-second_multiplication)/denominator
            matrix[i][j] = coefficient
            res = matrix[i][j]
    
    return res


def divided_difference_table(x_points, y_points):
    # set up the matrix
    size: int = 4
    matrix: np.array = np.zeros((size, size))
        
    # fill the matrix
    for index, row in enumerate(matrix):
        row[0] = y_points[index]
    # populate the matrix (end points are based on matrix size and max operations we're using)
    for i in range(1, len(x_points)):
        for j in range(1, i+1):
            # the numerator are the immediate left and diagonal left indices...
            numerator = (matrix[i][j-1]) - (matrix[i-1][j-1])
            # the denominator is the X-SPAN...
            denominator = x_points[i] - x_points[i-j]
            operation = numerator / denominator
            # cut it off to view it more simpler
            matrix[i][j] = '{0:.7g}'.format(operation)
    
    
    print(f'[{Decimal(matrix[1][1])}, {Decimal(matrix[2][2])}, {Decimal(matrix[3][3])}]', end = '\n\n')
    return matrix


def get_approximate_result(matrix, x_points, value):
    # p0 is always y0 and we use a reoccuring x to avoid having to recalculate x 
    reoccuring_x_span = 1
    reoccuring_px_result = y_points[0]
    
    # we only need the diagonals...and that starts at the first row...
    for index in range(1, 4):
        polynomial_coefficient = matrix[index][index]
        # we use the previous index for x_points....
        reoccuring_x_span *= (value - x_points[index-1])
        
        # get a_of_x * the x_span
        mult_operation = polynomial_coefficient * reoccuring_x_span
        # add the reoccuring px result
        reoccuring_px_result += mult_operation
    
    # final result
    return reoccuring_px_result


np.set_printoptions(precision=7, suppress=True, linewidth=100)

def apply_div_dif(matrix: np.array):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, i+2):
            # skip if value is prefilled (we dont want to accidentally recalculate...)
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue
            
            # get left cell entry
            left: float = matrix[i][j-1]
            # get diagonal left entry
            diagonal_left: float = matrix[i-1][j-1]
            # order of numerator is SPECIFIC.
            numerator: float = left - diagonal_left
            # denominator is current i's x_val minus the starting i's x_val....
            denominator = matrix[i][0] - matrix[i-j+1][0]
            # something save into matrix
            operation = numerator / denominator
            matrix[i][j] = operation
    
    return matrix
def hermite_interpolation():
    x_points = [3.6, 3.8, 3.9]
    y_points = [1.675, 1.436, 1.318]
    slopes = [-1.195, -1.188, -1.182]
    # matrix size changes because of "doubling" up info for hermite 
    num_of_points = len(x_points)
    matrix = np.zeros((2*num_of_points, 2*num_of_points))
    # populate x values (make sure to fill every TWO rows)
    for x in range(num_of_points):
        matrix[2*x][0] = x_points[x]
        matrix[2*x+1][0] = x_points[x]
    
    # prepopulate y values (make sure to fill every TWO rows)
    for x in range(num_of_points):
        matrix[2*x][1] = y_points[x]
        matrix[2*x+1][1] = y_points[x]
        
    # prepopulate with derivates (make sure to fill every TWO rows. starting row CHANGES.)
    for x in range(num_of_points):
        matrix[2*x+1][2] = slopes[x]
        
    filled_matrix = apply_div_dif(matrix)
    print(filled_matrix, end='\n\n')


def cubic_spline():
    x_points = [2, 5, 8, 10]
    y_points = [3, 5, 7, 9]
    
    h = np.zeros(4)
    a_matrix: np.array = np.zeros((4, 4))
    b_matrix = np.zeros(4)
   
    # part a
    for i in range(3):
        h[i] = x_points[i+1] - x_points[i]
        
    a_matrix[0][0] = 1
    a_matrix[1][0] = h[0]
    a_matrix[1][1] = 2*(h[0]+ h[1])
    a_matrix[1][2] = h[1]
    a_matrix[2][1] = h[1]
    a_matrix[2][2] = 2*(h[1]+h[2])
    a_matrix[2][3] = h[2]
    a_matrix[3][3] = 1
    
    print(a_matrix, end = '\n\n')
    
    # part b
    b_matrix[1] = 3/h[1] * (y_points[2]-y_points[1]) -  3/h[0] * (y_points[1]-y_points[0])
    b_matrix[2] = 3/h[2] * (y_points[3]-y_points[2]) -  3/h[1] * (y_points[2]-y_points[1])
    print(b_matrix, end = '\n\n')
    
    # part c
    X = np.linalg.inv(a_matrix).dot(b_matrix)
    print(X)

if __name__ == "__main__":
    # problem 1
    x_points = [3.6, 3.8, 3.9]
    y_points = [1.675, 1.436, 1.318]
    approximating_value = 3.7
    print(nevilles_method(x_points, y_points, approximating_value), end = '\n\n')
    
    #problem 2
    x_points = [7.2, 7.4, 7.5, 7.6]
    y_points = [23.5492, 25.3913, 26.8224, 27.4589]
    divided_table = divided_difference_table(x_points, y_points)
    #print(divided_table[1][1], divided_table[2][2], divided_table[3][3], end = '\n\n')
    print(get_approximate_result(divided_table, x_points, 7.3), end = '\n\n')
    
    hermite_interpolation()
    
    cubic_spline()