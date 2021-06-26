"""Shuffling of the Rubik's cube using Dixon cube.

See: https://people.math.carleton.ca/~jdixon/RandomGroupElts.pdf
Typical/offical number of moves in Rubik's cube tournament: 20-25-30
https://www.quora.com/How-are-Rubiks-Cubes-scrambled-at-a-competition
"""
from math import log
import os
import sys
project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


from sec_groups.fingroups import RubiksCube
from tools import sampling


def count_samples_in_buckets(samples):
    """Create 6 buckets, corresponding to 6 faces of the cube.

    Compare outcome of the histogram based on dice rolls. E.g.: 
    https://www.mathematik.tu-clausthal.de/en/mathematics-interactive/simulation/dice-roll-simulator/
    """
    # TODO: check for errors
    ctr = [0]*6
    for x in samples:
        for i in range(6):
            if i*8 <= x < (i+1)*8:
                ctr[i] += 1
    return ctr



if __name__ == '__main__':
    group = RubiksCube()
    
    print("Rubik's Cube in initial state; corresponds to identity permutation.")
    print(group.identity)
    print("Cube after 4 rotations of rear face (identical to initial state).")
    print(group.rear ^ 4)

    trials = {25: 1000, 1000: 1000}

    print("# Generate random elements of Rubik's cube group.")
    for k, n in trials.items():
        print("\nDixon cube length k=", k)
        print("Chi-squared test sample n=", n)

        # dix = sampling.generate_dixon_cube_thm1(group.generating_set, k)
        dix = sampling.generate_dixon_cube_lemma13b(group.generating_set, k)

        samples = []
        print(f"Start sampling.")
        for i in range(n):
            # Sample random element using Dixon cube
            r = sampling.random_group_element(dix)
            # Sample value of facet 0, the upper left facet on the top face.
            samples.append(r.indices_after_operation()[0])
            print(f"Generating samples: {round(100*i/(n+1))}%", end="\r")
        ctr = count_samples_in_buckets(samples)
        print("Frequencies=", ctr)
        print(f"Expected frequency per bucket, {n}/6=", round(n/6, 0))
        print("Mean squared error=", sum([(frq-n/6)**2 for frq in ctr])/6)
        print(f"chi^2= {round(sum([((c-n/6)**2)/(n/6) for c in ctr]), 2)}; critical value=11.07 for p-value of 0.05 and 5 degrees of freedom (less than critical value means statistically similar)")


    
# TODOs:
# TS-2: Expand RubiksCube api to n x n cubes, as per https://www.jaapsch.net/scramble_cube.htm?size=2&num=5&len=30&col=yobwrg&subbutton=Scramble%21
