## Algorithm implementation

All algorithms presented below work with functions implemented as truth tables.
Truth table file consists of two lines. First line has one number denoting dimension, second line has 2^dimension numbers denoting function images.
Example truth table of Gold function (x^3) in dimension 6 can be found in file `example_function.tt`.

Contributors:
- Ivana Ivkovic
- Nikolay Stoyanov Kaleyski 

### Dependency
vbf.h - vectorial boolean function library header

vbf.c - vectorial boolean function library

latest implementation maintained at:
[VBF GitLab](https://git.app.uib.no/Nikolay.Kaleyski/vectorial-boolean-functions)

---

### Algorithm 1: Testing linear equivalence of two uniformly distributed functions
_Check whether two uniformly distributed 3-to-1 functions (C1 and C2) are linear-equivalents_

alg1.c

**to compile:**

`make alg1`

**or:**

`gcc alg1.c vbf.c -lm -o check_lin_eq_2x_uniform_3to1.o`

**to run:**

`./check_lin_eq_2x_uniform_3to1.o function1.tt function2.tt`

**result:**

Truth tables of the linear functions L1 and L2 such that L1 ◦ C1 ◦ L2 = C2

**or:**

"False"

---

### Algorithm 2: Testing linear equivalence of a 3-to-1 function to a uniformly distributed 3-to-1 function
_Check whether a given function (T) is linear-equivalent to a uniformly distributed 3-to-1 function (C)_

alg2.c

**to compile:**

`make alg2`

**or:**

`gcc alg2.c vbf.c -lm -o find_uniform_3to1_lin_eq_to_3to1.o`

**to run:**

`./find_uniform_3to1_lin_eq_to_3to1.o function.tt`

**result:**

Truth table of the uniformly distributed 3-to-1 function C such that T ◦ L = C where L is a linear function. 
_*It is possible to modify the code to return the truth table of L._


**or:**

"False"

---

## Extra Algorithms

### Algorithm E1: Finding linear self-equivalence classes of a uniformly distributed 3-to-1 function
_Find all linear functions L1 and L2 for which a given uniformly distributed 3-to-1 function is linearly equivalent to self_

alg_e1.c

**to compile:**

`make alg_e1`

**or:**

`gcc alg_e1.c vbf.c -lm -o lin_eq_to_self.o`

**to run:**

`./lin_eq_to_self.o function.tt`

**result:**

truth tables of all the linear functions L1 and L2 such that L1 ◦ C1 ◦ L2 = C1
