# Sketch distance calculation

Sketches represent biological samples for pairwise distance calculation and
nearest-neighbour selection.

## Language

**Sketch block**:
A contiguous group of sample sketches considered together for distance
calculation. Its sample count is bounded by the user-selected block size.
_Avoid_: Tile when referring only to a group of sketches.

**Distance tile**:
The comparisons between a pair of sketch blocks. A diagonal tile compares a
block with itself; an off-diagonal tile compares two different blocks.
_Avoid_: Batch when referring to a set of pairwise comparisons.

**Sample pair**:
Two distinct samples whose distance contributes a neighbour candidate for
each sample. The pair is unordered: comparing A with B represents the same
pair as comparing B with A.

**Diagonal tile**:
The unordered sample pairs within one sketch block. Self-pairs are excluded.

**Off-diagonal tile**:
All sample pairs formed by taking one sample from each of two different
sketch blocks.
