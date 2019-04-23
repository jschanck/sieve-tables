# Definitions
- C(d) is the normalized (d-1)-dimensional spherical measure of a spherical cap of angle pi/3.
- W(d,a,b) is the normalized (d-1)-dimensional spherical measure
of an intersection of spherical caps of angles a and b whose
centers are at an angular distance of pi/3.

# Table entries
- 2075.csv lists (d, 1/C(d))
- 2924.csv lists (d, cos(a), cos(b), 1/(C(d)W(d,a,b))) for a choice of a and b that is explained below.
- 3494.csv lists (d, cos(a), 1/(C(d)W(d,a,a))) for a choice of a that is explained below.

# How to interpret the tables

## 2075.csv
Sieving algorithms for the shortest vector problem attempt to find many
pairs of lattice vectors that are at an angular distance of pi/3 or
less.  The "kissing constant" in dimension d is the size of the
largest list of d-dimensional unit vectors that does not contain
a pair of vectors at angular distance < pi/3. The kissing
constant is known to be larger than 1/C(d) =
2^((0.2075....+o(1))d), but the best known lower bound only
exceeds 1/C(d) by a linear factor in d [JJP].

Some analyses of sieving algorithms, e.g. [BDGL], assume that many
pairs of lattice vectors with small angular distance can be
found in a random list of size 1/C(d). This assumption is consistent
with the largest sieving experiments performed to date [G6K].
However, to my knowledge, no one has performed experiments that
are specifically designed to estimate the minimum list size
necessary for sieving.

This table should be interpreted as a folklore estimate for
the number of lattice vectors that any sieving algorithm must
enumerate. It is neither an upper bound nor a lower bound on
this quantity.

## 2924.csv
Consider the following experiment with parameters d, n, t, a, b:
1. Sample a list L of n i.i.d. uniform points on the d-1 sphere,
2. Sample a list F of t i.i.d. uniform points on the d-1 sphere,
3. For each (v,w) in L x L, compare v and w if and only if there exists
f in F such that <v,f> > cos(a)  and  <w,f> > cos(b). (To compare v
and w is to test whether <v,w> > cos(pi/3).)

The fourth column of 2924.csv should be interpreted as the expected number of
comparisons in step 3 when n = 1/C(d), t = 1/W(d,a,b), and a=b=pi/3.

#### Relevance to sieving algorithms:

The experiment underlying the table also underlies the complexity analysis of
the BDGL sieve. The optimal choice of a and b for BDGL depends
on the cost of implementing near-neighbor queries, i.e. the cost
of implementing the map:

  `v -> { (v,w) : exists f in F with <v, f> > cos(a) and <w, f> > cos(b) }`

Using the data structures in [BDGL], the optimal, asymptotic,
choice of both a and b is pi/3. The expected number of
comparisons, in this case, tends to
1/(C(d)W(d,pi/3,pi/3)) ~ 2^((0.2924...+o(1))d) as d->oo.

## 3494.csv
Consider the following experiment with parameters d, n, a:
1. Sample a list L of n i.i.d. uniform points on the d-1 sphere,
2. Independently, draw a uniform point f on the d-1 sphere,
3. For each pair (v,w) in L^2, compare v and w if and only if
   <v,f> > cos(a)   and   <w,f> > cos(a).
4. repeat 2-3 until n pairs (v,w) with <v,w> < 1/2 are found.

The third column of 3494.csv should be interpreted as the expected number of
comparisons in this experiment when n = 1/C(d) and a is given by
the second column.

#### Relevance to sieving algorithms:

Step 3 of the experiment compares all pairs in

`L' = {v in L : <v,f> > cos(a)}`,

and the list L' can be constructed in time L.
With cos(a)=0 the algorithm makes all of its comparisons in one iteration and
performs ~n^2 comparisons in total. As a is decreased the number of iterations
increases but the number of comparisons per iteration decreases. The second
column of 3494.csv lists the value of cos(a) that minimizes the
total number of comparisons. The optimal, asymptotic, choice of cos(a) for
n=1/C(d) is ~ sqrt(1 - sqrt(3/4)). In this case, the expected number of
comparisons tends to 2^((0.3494...+o(1))d).

This sieving algorithm is described in Section 5.1.bgj1 of
[G6K]. Note that [G6K] describes an algorithm for Exact-SVP in
dimension d that calls this sieve in dimension < d. The data
reported in Figure 3 of [G6K] is for this Exact-SVP algorithm,
and not for a direct application of the sieve. Figure 3 reports 2^42
Xeon E5-2650v3 cycles for d=100 (on average), and calls the
sieve in dimension 82 (on average). The table estimates 2^35
comparisons in dimension 82. [G6k] reports that comparisons are
implemented using a procedure that "consists of about a dozen
x86 non-vectorised instructions for vectors of dimension
roughly one hundred." More systematic experiments are needed.


# Cap and Wedge volume calculations

- C(d) was computed using the function log2_cap_beta in spheres.gp.
- W(d,a,b) was computed using the function log2_wedge_intnum in spheres.gp.
The calculation is based on Case 8 of [LK].

The tables were computed using 57 digit (192 bit) precision with
GP/PARI 2.12.0 (alpha). The scripts will not work with PARI <
2.12.0 as they rely on new hypergeometric function support.


# References

[BDGL] Anja Becker, Léo Ducas, Nicolas Gama, Thijs Laarhoven.
"New directions in nearest neighbor searching with applications
to lattice sieving." SODA 2016.
[ePrint](https://eprint.iacr.org/2015/1128)

[G6K]
Martin R. Albrecht, Léo Ducas, Gottfried Herold, Elena
Kirshanova, Eamonn W. Postlethwaite, Marc Stevens.
"The General Sieve Kernel and New Records in Lattice Reduction."
2019.
 [ePrint](https://eprint.iacr.org/2019/089)

[JJP] Matthew Jenssen, Felix Joos, Will Perkins.
"On kissing numbers and spherical codes in high dimensions."
Advances in Mathematics, Volume 335, 2018.
[arXiv](https://arxiv.org/abs/1803.02702)

[LK] Yongjae Lee, Woo Chang Kim. "Concise formulas for the surface area of the
intersection of two hyperspherical caps." [pdf](https://ie.kaist.ac.kr/uploads/professor/tech_file/Concise+Formulas+for+the+Surface+Area+of+the+Intersection+of+Two+Hyperspherical+Caps.pdf)
