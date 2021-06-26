from functools import reduce
import random


def generate_dixon_cube_thm1(generating_set, k):
    """Generate cube with length k.

    Create cube following Theorem 1 of Di08:
    Dixon, "Generating random elements in finite groups", EJoC, 2008.
    See: https://people.math.carleton.ca/~jdixon/RandomGroupElts.pdf
    """
    s = type(generating_set[0])
    d = len(generating_set)
    # dixon_w = copy.deepcopy(generating_set)
    dixon_w = generating_set.copy()
    for i in range(d, k):
        x_rhs = reduce(s.operation, [x for x in dixon_w if bool(random.getrandbits(1))], s.identity)
        x_lhs = reduce(s.operation, [x.inverse() for x in dixon_w if bool(random.getrandbits(1))][::-1], s.identity)
        x_next = x_lhs @ x_rhs
        dixon_w.append(x_next)
        print(f"Generating Dixon cube: {round(100*i/(k+1))}%", end="\r")
    assert len(dixon_w) == k
    return dixon_w


def generate_dixon_cube_lemma13b(generating_set, l):
    """Generate cube following Lemma 13b from Di08.

    Create cube of length 2l+d, where d is len(generating_set).
    Idea of construction is that distributions of cubes X_k and Y_k
    "should be approximately independent so that arguments of Lemma 13b
    may possibly apply".
    See explanation on p. 11 of Di08.
    """
    s = type(generating_set[0])
    d = len(generating_set)
    # Set up cube Y_d
    cube_y_d = []
    new = s.identity
    for i, x_i in enumerate(generating_set):
        new = x_i @ new
        cube_y_d += [new]

    # Create cubes X_k and Y_k such that:
    # X_k = Cube(x1, .., xk) where xk is a random element from Y_k-1
    # Y_k = Cube(y1, .., yk) where yk is a random element from X_k-1
    cube_x_k = generating_set
    cube_y_k = cube_y_d
    for i in range(d, l+d):  # counting from index 0; TODO: verify if k goes to l or l+d (text seems to contradict)
        x_k = reduce(s.operation, [x for x in cube_y_k if bool(random.getrandbits(1))], s.identity)
        cube_x_k += [x_k]
        y_k = reduce(s.operation, [x for x in cube_x_k if bool(random.getrandbits(1))], s.identity)
        cube_y_k += [y_k]

    cube_z_m = cube_x_k[:l+d] + cube_y_k[d:l+d]
    assert len(cube_z_m) == 2*l+d
    return cube_z_m


def random_group_element(dixon_cube):
    """Generate next random element based on given Dixon cube."""
    # dix = generate_dixon_cube_thm1(dixon_cube, len(dixon_cube)+1)
    s = type(dixon_cube[0])
    x_k = reduce(s.operation, [x for x in dixon_cube if bool(random.getrandbits(1))], s.identity)
    return x_k