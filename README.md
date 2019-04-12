* b.2075.csv gives a lower bound on the kissing number in dimension b.
* b.2924.csv gives a lower bound on the cost of the BDGL sieve.
* b.3494.csv gives a lower bound on the cost of the bgj1 sieve.
* The \*.jjp.csv files use the JJP lower bound on the kissing number
instead of the cap volume lower bound.

Further information can be found in my Apr. 11, 2019 email to the pqc-forum,
reproduced below.

---------

Dear pqc-forum,

All of the second round lattice-based KEMs --- Frodo, Kyber, LAC, NewHope,
NTRU, NTRU Prime, Round5, Saber, ThreeBears --- make some security claim based
on the sieve algorithm of Becker--Ducas--Gama--Laarhoven [BDGL]. Moreover, they
use cost estimates like
  * sqrt(3/2)^b ~ 2^0.2925b operations, and
  * sqrt(3/2)^b ~ 2^0.2925b bits of memory, or
  * sqrt(4/3)^b ~ 2^0.2075b bits of memory.

One could view these numbers as approximations to the asymptotic complexity of
the BDGL sieve, i.e.
  * 2^((0.2925...+o(1))b) operations, and
  * 2^((0.2925...+o(1))b) bits of memory, or
  * 2^((0.2075...+o(1))b) bits of memory.

However, one can also view them as approximations to
  * A(b)/p(b) and
  * A(b),

where A(b) is the kissing number in dimension b and p(b) is the probability
that two random points on the b-1 sphere at angular distance pi/3 lie
in a random spherical cap of angle ~ pi/3. As b goes to infinity A(b) goes to
2^((0.2075...+o(1))b) and A(b)/p(b) goes to 2^((0.2925...+o(1))b).

This interpretation of the "2^0.2925b" and "2^0.2075b" values avoids the
discussion of subexponential terms that contribute to the complexity of
the sieve, e.g. those from
  * sampling lattice vectors,
  * computing inner products,
  * list decoding random spherical codes,
  * etc.

It also allows you to write down a lower bound on the cost of the
sieve in any fixed dimension using only
  1) The lower bound on the kissing number from [JJP], and
  2) The (exact) wedge volume formula from [BDGL].

Note that this is how the asymptotic cost of BDGL is derived, modulo arguments that
other parts of the algorithm have sub-exponential cost.

I'd like to suggest that we replace the inaccurate formulas "0.2925b"
and "0.2075b" with tables of log A(b)/p(b) and log A(b). I've posted some
preliminary tables at [1]. The repository includes scripts for generating
the tables. These scripts also produce tables for "bgj1" from [G6K].

This community might be surprised by the inaccuracy of the 0.2925b formula.
In dimension b the table value exceeds 0.2925b by
   * b=300 : 22.9
   * b=400 : 24.4
   * b=500 : 25.5
   * b=600 : 26.5
   * b=700 : 27.3
   * b=800 : 28.0

More pointedly: the attempts to estimate the sub-exponential terms
using formulas like "0.2925b + 16.4" have not accounted for the dominant
sub-exponential term in the relevant range.

Some technical notes follow.

Cheers,

John

- The Chabauty--Shannon--Wyner lower bound on the kissing number was
  recently improved by Jenssen, Joos, and Perkins [JJP]. The improvement
  is "only" linear in the dimension, but this is significant in our
  setting. I've used the JJP bound for the tables.

- In finite dimension, the angle of the cap in the definition of p(b)
  has to be optimized. With the notation from BDGL, the tables entries
  correspond to a local minimum along "alpha=beta" in the neighborhood
  of alpha=beta=1/2.

[1] https://github.com/jschanck/sieve-tables

[BDGL] https://eprint.iacr.org/2015/1128

[G6K] https://eprint.iacr.org/2019/089

[JJP] https://arxiv.org/abs/1803.02702
