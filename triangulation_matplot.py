# Vlachos Dimitris 1641

import sys
import matplotlib.pyplot as plt


# A simple Polygon Sample to kickstart the algorithm
sample_P = [
    (200.0, 180.0),
    (180.0, 190.0),
    (165.0, 210.0),
    (165.0, 235.0),
    (150.0, 250.0),
    (100.0, 170.0),
    (150.0, 100.0),
    (170.0, 140.0),
    (130.0, 150.0)
]


################################################################################
################################################################################
################# Triangulation of a simple polygon P ##########################
################################################################################
################################################################################


###########################################################################
######################## Check if  P is simple. ###########################
###########################################################################

############ Helper function for segments_intersection ##############
#####################################################################
# Check: if a point q lies on segment pr, then return True, else return False
def on_segment(p, q, r):
    # Coordinates for the 3 points:
    px, py = p
    qx, qy = q
    rx, ry = r
    
    # Point q belongs to straight line pr iff:
    # 
    #  Cond1: Check if the x‐coordinate of q is between the x-coordinates of p and r.
    x_within = False
    if (min(px, rx) <= qx) and (qx <= max(px, rx)):
        x_within = True
    #
    #  Cond2: Check if the y‐coordinate of q is between the y-coordinates of p and r.
    y_within = False
    if (min(py, ry) <= qy) and (qy <= max(py, ry)):
        y_within = True

    # Both conditions must be met.
    return x_within and y_within 


########### Helper function for segments_intersection ##############
####################################################################
# if p ==> q ==> r makes a left turn(Counter Clock Wise) then return positive ,else return negative.
def orientation(p, q, r):
    # Compute vector p->q
    dx1 = q[0] - p[0]
    dy1 = q[1] - p[1]
    # Compute vectot q->r
    dx2 = r[0] - p[0]
    dy2 = r[1] - p[1]

    # The cross product in 2D:
    cross_product = dx1 * dy2 - dy1 * dx2

    return cross_product # positive or negative number (zero if collinear)


############# Helper function for is_simple_polygon ##############
##################################################################
# Check if segments p1q1 and p2q2 intersect (we do not include shared endpoints)
def segments_intersection(p1, q1, p2, q2):
    
    ori1 = orientation(p1, q1, p2)
    ori2 = orientation(p1, q1, q2)
    ori3 = orientation(p2, q2, p1)
    ori4 = orientation(p2, q2, q1)
    # General case:
    # Step 1: Check if ori1 and ori2 have oposite sign(+/-)
    first_pair_opposite = False
    if (ori1 > 0 and ori2 < 0) or (ori1 < 0 and ori2 > 0):
        first_pair_opposite = True
    # Step 2: Check if ori3 and ori4 have oposite sign
    second_pair_opposite = False
    if (ori3 > 0 and ori4 < 0) or (ori3 < 0 and ori4 > 0):
        second_pair_opposite = True

    if first_pair_opposite and second_pair_opposite:
        return True

    # For simpicity we ignore collinear cases.  
    return False


# P: a list with points (x,y)
def is_simple_polygon(P):
    n = len(P) # number of vertices(or nodes)
    # Every outside loop works with an edge p1q1.
    for i, _ in enumerate(P):
        # p1: vertex with pointer i, q1: next vertex, adjacent to p1
        p1 = P[i]
        q1 = P[(i+1)%n] # if i == n-1, the next pointer is 0 (So we'll form a "cycle")
        # Compare every edge p1q1 with all the others p2q2
        for j in range(i+1, n):
            p2, q2 = P[j], P[(j+1)%n]
            # Case 1: Ignore edges that have in common a vertex. 
            if (i+1) % n == j:
                continue
            # Case 2: Ignore edges that are 1st and last (share end of "cycle").
            if i == (j+1) % n:
                continue
            # Case 3: Check if edges p1q1 and p2q2 are adjacent(common vertex). Ignore them if they are.
            if p1 == p2 or p1 == q2 or q1 == p2 or q1 == q2:
                continue
            # Case 4: # Do the orientation tests
            if segments_intersection(p1, q1, p2, q2):
                return False
    return True


###########################################################################
######################## Check if P is y-monotone #########################
###########################################################################
def check_if_y_monotone(P):
    n = len(P) # number of vertices

    # Init with the 1st vertex
    max_y = P[0][1]
    i_max = float('-inf')  # Smaller than any other number - the pointer with max y
    min_y = P[0][1] 
    i_min = float('inf')  # Larger than any other number - the pointer with min y

    # Compute i_max and i_min locally
    i_max = max(range(n), key=lambda i: P[i][1])
    i_min = min(range(n), key=lambda i: P[i][1])

    # Find index(or pointer) of vertex with max y and vertex with min y 
    for index, point in enumerate(P):
        y = point[1]
        # Update with the max y we've found
        if y > max_y:
            max_y = y
            i_max = index
        # Update with the min y we've found
        if y < min_y:
            min_y = y
            i_min = index

    # Create 2 chains: left_chain and right_chain
    left_chain = []
    i = i_max
    while True:
        left_chain.append(i)
        if i == i_min:
            break
        i = (i + 1) % n

    right_chain = []
    i = i_max
    while True:
        right_chain.append(i)
        if i == i_min:
            break
        i = (i - 1 + n) % n

    # Check if the left chain descends (y descending)
    for i in range(1, len(left_chain)):
        y_prev = P[left_chain[i - 1]][1]
        y_curr = P[left_chain[i]][1]
        if y_curr > y_prev: # if y ascends at any stepp, then it's not y-monotone 
            return False  

    # Check if the left chain descends (y descending)
    for i in range(1, len(right_chain)):
        y_prev = P[right_chain[i - 1]][1]
        y_curr = P[right_chain[i]][1]
        if y_curr > y_prev:
            return False

    # Iff both chains are descending  <====> P is y-monotone
    return True

##########################################################################    
####################### Triangulate the Polygon ##########################
##########################################################################
def triangulate_monotone_P(P):
    n = len(P) # number of vertices

    # Compute i_max and i_min locally
    i_max = max(range(n), key=lambda i: P[i][1])
    i_min = min(range(n), key=lambda i: P[i][1])

    # Create 2 chains: left_chain and right_chain
    left_chain = []
    i = i_max
    while True:
        left_chain.append(i)
        if i == i_min:
            break
        i = (i + 1) % n

    right_chain = []
    i = i_max
    while True:
        right_chain.append(i)
        if i == i_min:
            break
        i = (i - 1 + n) % n


    # Every vertex will belong to the left or right "chain".
    # So keep track for each one in a dictionary with a "label".
    chain_label = {}
    for index in left_chain:
        chain_label[index] = 'left'
    for index in right_chain:
        chain_label[index] = 'right'


    # Sort all vertices by decreasing y (if same y's, then increasing x).
    # Create a list of indices
    indices = list(range(n))
    # Sort the list
    indices.sort(key=lambda i: (-P[i][1], P[i][0]))


    # Kickstart the algorithm:
    # Initialize stack with first 2 vertices # 1st step always
    stack = [indices[0], indices[1]]
    triangles = []

    # Process remaining vertices in sorted order
    for k in range(2, n):
        u = indices[k]
        # If u belongs to a different chain than top of stack
        if chain_label[u] != chain_label[stack[-1]]:
            # Connect u to all vertices in the stack
            while len(stack) > 0: # while stack not empty
                v = stack.pop() # pop the "newer"(last insterted) element(vertex) & store it
                if stack: # if stack not empty after the pop
                    # Create a triange(u is the newest vertex ,v is the vertex just popped
                    # stack[-1] is the vertex which was below v before it was popped.
                    triangles.append((u, v, stack[-1]))
            # Push back the last processed and u
            stack = [indices[k - 1], u] # k: points at the current top of the list
        else:
            # u and the element at the top of stack are on the same chain(either left or right)
            v = stack.pop() # v: the most‐recently pushed vertex among those that remain on the chain
            while len(stack) > 0: # check at the new top of the stack (w) and check if (u, v, w) forms a “valid” triangle(T)
                w = stack[-1]
                orient = orientation(P[w], P[v], P[u])
                chain_label[u] == 'left'
                if (chain_label[u] == 'left' and orient > 0) or (chain_label[u] != 'left' and orient < 0):
                    triangles.append((w, v, u)) # append a valid triangle(T)
                    v = stack.pop()
                else:
                    break # exit the loop without "creating" more triangles
            stack.append(v) # push v back into the top of the stack
            stack.append(u) # also push u (it's the most‐recent vertex in our chain)
    
    # Convert triangle indices into actual coordinate triples
    triangle_coords = []
    for triangle in triangles:
        # trianle: in the form (p, q, r)
        p = P[triangle[0]]
        q = P[triangle[1]]
        r = P[triangle[2]]
        triangle_coords.append((p, q, r))

    return triangle_coords


###################################################################################
######################## Visualize the Triangulation ##############################
###################################################################################
def draw_polygon_and_triangles(P, triangles):
    plt.figure()
    ax = plt.gca() # get current axis (to "plot")
    ax.set_aspect(1.0) # the geometrical object(s) will be visualized correctly
    ax.set_title("Polygon Triangulation")

    # Plot the "border" of the polygon
    x_coordinates = []
    y_coordinates = []
    for p in P:
        x_coordinates.append(p[0])
        y_coordinates.append(p[1])

    # 'Close' the border of the P (from pn to p0)
    x_coordinates.append(P[0][0])
    y_coordinates.append(P[0][1])

    plt.plot(x_coordinates, y_coordinates, 'k-', label="Polygon") # k: for black

    # Design the triangles
    for tri in triangles:
        tri_x = [tri[0][0], tri[1][0], tri[2][0], tri[0][0]]
        tri_y = [tri[0][1], tri[1][1], tri[2][1], tri[0][1]]
        plt.plot(tri_x, tri_y, 'r')  # r: for red triangles

    # Plot points
    for (x, y) in P: # for every point in the list
        plt.plot(x, y, 'ro')  # red cycles
        plt.text(x + 2, y + 2, f"({x:.0f}, {y:.0f})", fontsize=10) # +2: so we can clearly see the coords

    plt.grid(True)
    plt.legend()
    plt.axis('equal')  # ensures both axes will have the same scale
    plt.show()


P = sample_P
###################################################################################
###################################################################################
#################################### main ########################################
###################################################################################


print("=========== Program to Triangulate simple polygons. ===========")
while True:

    # Check if a P is simple (iff: v1 ===> v2 ===> v3===>...===> vn-1 ===> vn===> v1)
    if not is_simple_polygon(P):
        print("P is not simple. ===> Termination...")
        sys.exit(1)

    # Check monotonicity for our P
    if not check_if_y_monotone(P):
        print("P is not y-monotone. ===> Skipping...")
        sys.exit(1)
    
    # Triangulate
    triangles = triangulate_monotone_P(P)
    print("Triangles each one contains 3 vertices):")
    for t in triangles:
        print(t)

    # Visualize the Triangulation
    draw_polygon_and_triangles(P, triangles)

    # Can continue giving Polygons...
    # Do you want to triangulate a new Polygon?
    answer = input("\n\nDo you want to insert a new Polygon? (Yes/No): ").strip()
    if answer not in ('Yes', 'yes', 'y'):
            print("=====> Termination....")
            break
    
    try:
        s = input("How many vertices do you want for P?").strip()
        n = int(s)
        if n < 3: # 3 vertices ===> Triangle
            print("P is invalid.")
            break
    except ValueError:
        print(f"'{s}' is not an integer.")        
        break


    # Read new coordinates
    valid_input = True
    P = [] # To store the coords.
    for i in range(n): 
        print("The coordinates should be given in Counter Clock wise or Clock wise turn!")    
        line = input(f"Please give coordinates {i+1} (x y): ").strip().split()
        if len(line) != 2:
            print("Error: Please give exactly 2 numbers (x y).")
            P = [] # Must init again.
            break
            
        try:
            x = float(line[0])
            y = float(line[1])
            P.append((x, y))
        except ValueError:
            print("Invalid input. Coordinates must be integers, floats or doubles.")
            valid_input = False
            break


    # If there was an error with loading the points (x,y), ask the user again
    if len(P) != n: # or not valid_input:
        print("There was a problem loading the points.Please try again.")
        P = []
        continue

    # Check if P is simple
    if not is_simple_polygon(P):
        print("The P is not simple. Please try again.")
        P = []
        continue

    # Check  monotonicity for user's P
    if not check_if_y_monotone(P):
        print("The P you gave is not y-monotone. Please try again.")
        print("Restarting polygon input...........................")
        P = []
        continue
