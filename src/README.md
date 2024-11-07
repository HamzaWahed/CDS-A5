## Description

`lf_movelet` takes an offset and run in the F column, finds their index in F, then returns the offset and run of that same index but in L.

`lf` takes an index in L and maps it to F (standard LF mapping).

`BWT_inverse_lf` computes the original string using the LF mapping technique used in `lf`.

`BWT_inverse_movelet` uses `lf` to map from L to F, then uses `lf_movelet` to continue the mapping while avoiding rank queries.

## Instructions

1. Run `make`.

2. Run `./lf <NAME OF .ri FILE>` and follow the prompt to map an index in L to an index in F using the approach described in part 1 of the assignment.

3. Run `./movelet <NAME OF .ri FILE>` and follow the prompt to find the run and offset of position LF(i) in B_L given its run and offset (part 2 of the assignment). 

4. Run `./BWT_inverse_lf <NAME of .ri FILE> <NAME of corresponding bwt file>` to do BWT inversion the "traditional way".

5. Run `./BWT_inverse_movelet <NAME of .ri FILE> <NAME of corresponding bwt file>` to do BWT inversion the cool way.