import copy



START_CUBE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 17, 18, 24, 25, 26, 32, 33, 34, 11, 12, 19, 20, 27, 28, 35, 36, 13, 14, 15, 21, 22, 23, 29, 30, 31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]

TEMPLATE = (""
    "                 +--------------+\n"
    "                 |              |\n"
    "                 | {:02d}   {:02d}   {:02d} |\n"
    "                 |              |\n"
    "                 | {:02d}  top   {:02d} |\n"
    "                 |              |\n"
    "                 | {:02d}   {:02d}   {:02d} |\n"
    "                 |              |\n"
    "  +--------------+--------------+--------------+--------------+\n"
    "  |              |              |              |              |\n"
    "  | {:02d}   {:02d}   {:02d} | {:02d}   {:02d}   {:02d} | {:02d}   {:02d}   {:02d} | {:02d}   {:02d}   {:02d} |\n"
    "  |              |              |              |              |\n"
    "  | {:02d}  left  {:02d} | {:02d} front  {:02d} | {:02d} right  {:02d} | {:02d}  rear  {:02d} |\n"
    "  |              |              |              |              |\n"
    "  | {:02d}   {:02d}   {:02d} | {:02d}   {:02d}   {:02d} | {:02d}   {:02d}   {:02d} | {:02d}   {:02d}   {:02d} |\n"
    "  |              |              |              |              |\n"
    "  +--------------+--------------+--------------+--------------+\n"
    "                 |              |\n"
    "                 | {:02d}   {:02d}   {:02d} |\n"
    "                 |              |\n"
    "                 | {:02d} bottom {:02d} |\n"
    "                 |              |\n"
    "                 | {:02d}   {:02d}   {:02d} |\n"
    "                 |              |\n"
    "                 +--------------+\n"
    )


def cycle_to_image_map(cycle, degree):
    map = list(range(degree))
    for i, c in enumerate(cycle):
        if i < len(cycle)-1:
            map[c] = cycle[i+1]
        else:
            map[c] = cycle[0]
    return map


def clockwise_rotations(s):
    # Counterclockwise cycles from https://www.gap-system.org/Doc/Examples/rubik.html
    top_ccw_cycles = [[ 0, 2, 7, 5],[ 1, 4, 6, 3],[ 8,32,24,16],[ 9,33,25,17],[10,34,26,18]]
    left_ccw_cycles = [[ 8,10,15,13],[ 9,12,14,11],[ 0,16,40,39],[ 3,19,43,36],[ 5,21,45,34]]
    front_ccw_cycles = [[16,18,23,21],[17,20,22,19],[ 5,24,42,15],[ 6,27,41,12],[ 7,29,40,10]]
    right_ccw_cycles = [[24,26,31,29],[25,28,30,27],[ 2,37,42,18],[ 4,35,44,20],[ 7,32,47,23]]
    rear_ccw_cycles = [[32,34,39,37],[33,36,38,35],[ 2, 8,45,31],[ 1,11,46,28],[ 0,13,47,26]]
    bottom_ccw_cycles = [[40,42,47,45],[41,44,46,43],[13,21,29,37],[14,22,30,38],[15,23,31,39]]    
    gset_cycles = [top_ccw_cycles, left_ccw_cycles, front_ccw_cycles, right_ccw_cycles, rear_ccw_cycles, bottom_ccw_cycles]

    generating_set = []
    for g in gset_cycles:
        new_g = copy.deepcopy(s.identity)
        for cycle in g:
            # Convert to image map representation.
            im = cycle_to_image_map([i for i in cycle], s.degree)
            new_g = new_g @ s(im)
        # Make clockwise.
        generating_set.append(new_g.inverse())

    return generating_set

