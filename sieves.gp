\\ Spherical code size lower bounds, and related quantities, from
\\   [JJP18] Matthew Jenssen, Felix Joos, Will Perkins.
\\   "On kissing numbers and spherical codes in high dimensions."
\\   Advances in Mathematics, Volume 335, 2018.
\\   https://arxiv.org/abs/1803.02702
\\
\\ Wedge volume formula from
\\   [BDGL16] Anja Becker, Léo Ducas, Nicolas Gama, Thijs Laarhoven.
\\   "New directions in nearest neighbor searching with applications
\\    to lattice sieving." SODA 2016.
\\   https://eprint.iacr.org/2015/1128

log2(x) = {if(x == 0, 0, log(x)/log(2))};

log2_vol_cap(dim,angle) = {
  \\ Normalized spherical measure of a cap of the given angle.
  (dim-1)*log2(sin(angle)) - 1/2 * log2(2*Pi*dim) - log2(cos(angle));
}

CSW_log2_code_size_lower_bound(dim, angle) = {
  \\ Chabauty--Shannon--Wyner spherical code size lower bound.
  \\ A(dim,angle) >= 1/volume(spherical cap of given angle)
  \\ Explicitly [Eq. 2, JJP18]
  \\ A(dim,angle) >= (1+o(1))sqrt(2*Pi*dim)*cos(angle)/sin(angle)^(dim-1)
  \\ This function ignores the o(1).
  -log2_vol_cap(dim,angle);
}

JJP_sin_q(angle) = {
  my(r);
  \\ sin(q(angle)) where q is as defined in [Eq. 3, JJP18].
  r = (cos(angle) - 1)^2 * (1 + 2*cos(angle));
  sqrt(r) / sin(angle);
}

JJP_q(angle) = {
  \\ [Eq. 3, JJP18]
  \\ q(angle) = arcsin( sqrt(r) / sin(angle) )
  \\ where r = (cos(angle) - 1)^2 * (1 + 2*cos(angle))
  asin(JJP_sin_q(angle));
}

JJP_log2_code_size_lower_bound(dim, angle) = {
  my(log2_c);
  \\ Theorem 2 of [JJP18] (Linear improvement over CSW for any fixed angle).
  \\ A(dim,angle) >= (1+o(1)) * (CSW bound) * c(angle) * dim
  \\ where c(angle) = ln(sin(angle)/JJP_sin_q(angle))
  \\ This function ignores the o(1).
  log2_c = log2(log(sin(angle)) - log(JJP_sin_q(angle)));
  CSW_log2_code_size_lower_bound(dim, angle) + log2(dim) + log2_c;
}

BDGL_log2_wedge(n,a,b,t) = {
  my(g_sq, log2_A);
  \\ Lemma 2.2 of [BDGL].
  \\ Suppose <x,y> = cos(t), C1 is a spherical cap of
  \\ angle a with center x, and C2 is a spherical cap of
  \\ angle b with center y. The "wedge" C1 intersect C2
  \\ has volume given by this function.
  g_sq = (cos(a)^2 + cos(b)^2 - 2*cos(a)*cos(b)*cos(t))/sin(t)^2;
  log2_A = log2(g_sq) - 2*log2(1-g_sq);
  n * log2(sqrt(1-g_sq)) + log2_A - 2*log2(n);
}

NV_cost(dim, quantum=0, JJP=0) = {
  my(log2_N, log2_test_cost, log2_find_one_cost, log2_find_all_cost,\
     log2_ops, log2_space);
  \\ Nguyen-Vidick sieve
  \\ Start with a list of size N ~ kissing number for dimension
  if(JJP,
    log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3),
    log2_N = CSW_log2_code_size_lower_bound(dim, Pi/3));
  \\ Test all pairs
  log2_test_cost = 0;
  if(quantum,
    log2_find_one_cost = .5 * log2_N + log2_test_cost,
    log2_find_one_cost = log2_N + log2_test_cost);
  log2_find_all_cost = log2_N + log2_find_one_cost;

  log2_ops = log2_find_all_cost;
  log2_space = log2_N;
  [log2_ops, log2_space];
}

bgj1_cost(dim, t, quantum=0, JJP=0) = {
  my(log2_N, log2_test_cost, log2_fill_cost, log2_bucket_size,\
     log2_find_one_cost, log2_find_all_cost, log2_per_filter_cost,\
     log2_p, log2_ops, log2_space);
  \\ Start with a list of size N ~ kissing number for dimension
  if(JJP,
    log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3),
    log2_N = CSW_log2_code_size_lower_bound(dim, Pi/3));
  \\ Pick a random spherical cap of angle t.
  \\ Bucket the vectors that lie in this cap.
  \\ (Inclusion is an inner product test.)
  log2_test_cost = 0;
  log2_fill_cost = log2_N + log2_test_cost;
  log2_bucket_size = log2_N + log2_vol_cap(dim, t);
  \\ Fix some v in the bucket.
  \\ Find all w in the bucket that are close to v. (inner product test)
  \\ There are O(1) neighbors of v in bucket.
  if(quantum,
    log2_find_one_cost = .5 * log2_bucket_size + log2_test_cost, \\ Grover
    log2_find_one_cost = log2_bucket_size + log2_test_cost);     \\ Exhaustive
  \\ Repeat for each v.
  log2_find_all_cost = log2_bucket_size + log2_find_one_cost;
  \\ Cost per filter is maximum of filling cost and checking cost
  log2_per_filter_cost = max(log2_fill_cost, log2_find_all_cost);
  \\ The filter only detects a pair (v,w) at angle < Pi/3 if the filter lies
  \\ in the intersection of caps of angle t around v and w. This happens with
  \\ probability given by the BDGL wedge formula.
  log2_p = BDGL_log2_wedge(dim, t, t, Pi/3);
  \\ Repeat with random filters until all pairs are detected.
  log2_ops = log2_per_filter_cost - log2_p;
  log2_space = log2_N;
  [log2_ops, log2_space];
}

BDGL_cost(dim, a, b, quantum=0, JJP=0) = {
  my(log2_N, log2_t, log2_insert_cost, log2_bucket_size, \
     log2_relevant_filters, log2_find_one_cost, log2_find_all_cost,\
     log2_ops, log2_space);
  \\ Start with a list of size N ~ kissing number in dimension
  if(JJP,
    log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3),
    log2_N = CSW_log2_code_size_lower_bound(dim, Pi/3));
  \\ Initialize t := t(a,b) filters. Each filter defines a bucket.
  log2_t = -BDGL_log2_wedge(dim, a, b, Pi/3);
  \\ Each of the N vectors is inserted into a number of buckets
  \\ that is proportional to the measure of a cap of angle b.
  log2_insert_cost = log2_N + log2_vol_cap(dim, b) + log2_t;
  log2_bucket_size = log2_N + log2_vol_cap(dim, b);
  \\ Each vector is compared against a number of buckets
  \\ proportional to the measure of a cap of angle a.
  log2_relevant_filters = log2_vol_cap(dim, a) + log2_t;
  \\ The probability that a pair of vectors at angle Pi/3 is detected is
  \\ given by the measure of the intersection of caps of angle a and b at
  \\ angular distance Pi/3. Hence the definition of t.

  \\ Assume that the cost of identifying the relevant filters is roughly
  \\ equal to the number of relevant filters (!).
  \\ Cost to find a neighbor of v is proportional to aggregate bucket size 
  if(quantum,
    log2_find_one_cost = max(log2_relevant_filters, .5*(log2_relevant_filters + log2_bucket_size)),
    log2_find_one_cost = max(log2_relevant_filters,     log2_relevant_filters + log2_bucket_size));
  \\ This has to be done for each v.
  log2_find_all_cost = log2_N + log2_find_one_cost;

  log2_ops = max(log2_insert_cost, log2_find_all_cost);
  log2_space = max(log2_N, log2_t + log2_bucket_size);
  [log2_ops, log2_space];
}

local_min(f,x,A=1,B=5) = {
  \\ Search left f(x-ε) and right f(x+ε) for a local minimum near f(x).
  \\ ε will be determined to 10^-B, starting with steps of size 10^-A.
  my(y, k);
  y = f(x);
  for(k=A, B, e=0.1^k; while(f(x-e) < y, x=x-e; y = f(x)));
  for(k=A, B, e=0.1^k; while(f(x+e) < y, x=x+e; y = f(x)));
  x
}

bgj1_cost_min(dim,quantum=0,JJP=0) = {
  my(t);

  if(quantum,
    \\ Asymptotically optimal filter angle for quantum bgj1 is asin((3/4)^(1/6)) = 1.2635...
    t = local_min(x->(bgj1_cost(dim,x,1,JJP)[1]), asin((3/4)^(1/6))),
    \\ Asymptotically optimal filter angle for classical bgj1 is asin((3/4)^(1/4)) = 1.1960...
    t = local_min(x->(bgj1_cost(dim,x,0,JJP)[1]), asin((3/4)^(1/4))));

  bgj1_cost(dim,t,quantum)[1];
}

BDGL_cost_min(dim,quantum=0,JJP=0) = {
  my(t);

  if(quantum,
    \\ Asymptotically optimal filter angle for quantum BDGL is acos(sqrt(3/16)) = 1.2296...
    t = local_min(x->(BDGL_cost(dim,x,x,1,JJP)[1]), acos(sqrt(3/16))),
    \\ Asymptotically optimal filter angle for classical BDGL is Pi/3 = 1.0471...
    t = local_min(x->(BDGL_cost(dim,x,x,0,JJP)[1]), acos(1/2)));

  BDGL_cost(dim,t,t,quantum)[1];
}

gen_tables() = {
  my(bl,bh,s,b);

  bl=50;
  bh=1105;

  s = "";
  for(b=bl, bh,
    s = concat(s, Strprintf("%d, %.1f\n", b, BDGL_cost_min(b,quantum=0,JJP=0))));
  write("b.2924.csv", s);

  s = "";
  for(b=bl, bh,
    s = concat(s, Strprintf("%d, %.1f\n", b, BDGL_cost_min(b,quantum=0,JJP=1))));
  write("b.2924.jjp.csv", s);

  s = "";
  for(b=bl, bh,
    s = concat(s, Strprintf("%d, %.1f\n", b, bgj1_cost_min(b,quantum=0,JJP=0))));
  write("b.3494.csv", s);

  s = "";
  for(b=bl, bh,
    s = concat(s, Strprintf("%d, %.1f\n", b, bgj1_cost_min(b,quantum=0,JJP=1))));
  write("b.3494.jjp.csv", s);

  s = "";
  for(b=bl, bh,
    s = concat(s, Strprintf("%d, %.1f\n", b, CSW_log2_code_size_lower_bound(b, Pi/3))));
  write("b.2075.csv", s);

  s = "";
  for(b=bl, bh,
    s = concat(s, Strprintf("%d, %.1f\n", b, JJP_log2_code_size_lower_bound(b, Pi/3))));
  write("b.2075.jjp.csv", s);
}

;
