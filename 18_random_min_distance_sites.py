#!/usr/bin/env python3

import bisect
import random
import sys
from collections import defaultdict

input_path, distance_text, seed_text, output_path = sys.argv[1:5]

distance = int(distance_text)
seed = int(seed_text)
rng = random.Random(seed)

chromosomes = []
positions = defaultdict(list)

with open(input_path) as handle:
    for line in handle:
        chrom, position = line.rstrip().split("\t")[:2]
        position = int(position)

        if chrom not in positions:
            chromosomes.append(chrom)

        positions[chrom].append(position)

selected = {}

for chrom in chromosomes:
    candidates = positions[chrom][:]
    rng.shuffle(candidates)

    retained = []

    for position in candidates:
        index = bisect.bisect_left(retained, position)

        too_close_left = (
            index > 0 and
            position - retained[index - 1] < distance
        )

        too_close_right = (
            index < len(retained) and
            retained[index] - position < distance
        )

        if not too_close_left and not too_close_right:
            bisect.insort(retained, position)

    selected[chrom] = retained

with open(output_path, "w") as handle:
    for chrom in chromosomes:
        for position in selected[chrom]:
            handle.write(f"{chrom}\t{position}\n")
