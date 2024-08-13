import numpy as np
import matplotlib.pyplot as plt
from shapely import Point
from shapely.wkt import loads # for Point(loads('POINT (0.0 0.0)'))
from shapely import distance 
from shapely import get_x
from shapely import get_y
from shapely import transform

points_between = 10

def rotate(vector, angle):
    angle = np.radians(angle)
    c, s = np.cos(angle), np.sin(angle)
    Rot_matrix = np.array(((c, -s), (s, c)))
    return vector @ Rot_matrix


def main():
    input_filename = "examples/in_decagon.txt"
    output_filename = "examples/out_decagon.txt"

    file = open(input_filename, 'r')
    out_file = open(output_filename, 'w')
    in_amount_of_points = int(file.readline().strip("#"))
    out_amount_of_points = in_amount_of_points + (in_amount_of_points - 1) * points_between

    array_of_points = np.zeros((out_amount_of_points, ), dtype=Point)
    counter=0

    point_1 = Point(loads("POINT" + "(" + file.readline() + ")"))

    array_of_points[counter] = point_1
    counter += 1

    the_smallest_unit = 100 #bad, very bad
    the_biggest_unit = 0 

    for line in file:
        point_2 =  Point(loads("POINT" + "(" + line + ")"))

        unit_distance = distance(point_1, point_2)/(points_between + 1)

        if unit_distance > the_biggest_unit :
            the_biggest_unit = unit_distance
        if unit_distance < the_smallest_unit :
            the_smallest_unit = unit_distance

        unit_vector = Point(point_2.x - point_1.x, point_2.y - point_1.y)

        unit_vector = Point(unit_vector.x / distance(point_1, point_2), unit_vector.y / distance(point_1, point_2))

        for i in range(1, points_between + 1, 1):
            next_point = Point(point_1.x + i * unit_distance * unit_vector.x
                    ,point_1.y +  i * unit_distance * unit_vector.y);
            array_of_points[counter] = next_point
            counter += 1

        array_of_points[counter] = point_2
        counter += 1
        point_1 = point_2

    G1 = np.array([get_x(val) for val in array_of_points])
    G2 = np.array([get_y(val) for val in array_of_points])

    vector_to_shift = np.ndarray(shape=(2,), dtype=np.float64)
    vector_to_shift[0] = 1
    vector_to_shift[1] = 0
    
    to_cmp = array_of_points.copy()

    for index in range(len(array_of_points)):
        vector_to_shift = rotate(vector_to_shift, np.random.uniform(0, 360))
        vector_to_shift *= np.random.uniform(the_smallest_unit, the_biggest_unit)
        array_of_points[index] = Point(array_of_points[index].x + vector_to_shift[0]
                , array_of_points[index].y + vector_to_shift[1])
        vector_to_shift[0] = 1
        vector_to_shift[1] = 0
    
    plt.plot(G1, G2, color='blue', linewidth=1)
    plt.scatter(G1, G2, color='red', marker='o', s=10)

    G1 = np.array([get_x(val) for val in array_of_points])
    G2 = np.array([get_y(val) for val in array_of_points])

    plt.plot(G1, G2, color='green', linewidth=1)
    plt.scatter(G1, G2, color='orange', s=10)

    plt.title('Plot points')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    plt.show()  
    
    out_file.write("# " + str(len(array_of_points)) + "\n")
    for point in array_of_points:
        out_file.write(f"{point.x:.8f} " + f"{point.y:.8f}\n")

if __name__ == "__main__":
	main()


