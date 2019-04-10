\\ Spherical code size lower bounds, and related quantities, from
\\   [JJP18] Matthew Jenssen, Felix Joos, Will Perkins.
\\   "On kissing numbers and spherical codes in high dimensions."
\\   Advances in Mathematics, Volume 335, 2018.
\\   https://arxiv.org/abs/1803.02702

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
  \\ Theorem 2 of [JJP18] (Linear improvement over CSW for any fixed angle).
  \\ A(dim,angle) >= (1+o(1)) * (CSW bound) * c(angle) * dim
  \\ where c(angle) = ln(sin(angle)/JJP_sin_q(angle))
  \\ This function ignores the o(1).
  log2_c = log2(log(sin(angle)) - log(JJP_sin_q(angle)));
  CSW_log2_code_size_lower_bound(dim, angle) + log2(dim) + log2_c;
}

BDGL_log2_wedge(n,a,b,t) = {
  \\ Suppose <x,y> = cos(t), C1 is a spherical cap of
  \\ angle a with center x, and C2 is a spherical cap of
  \\ angle b with center y. The "wedge" C1 intersect C2
  \\ has volume given by this function.
  g2 = (cos(a)^2 + cos(b)^2 - 2*cos(a)*cos(b)*cos(t))/sin(t)^2;
  n * log2(sqrt(1-g2)) - 2*log2(n);
}




NV_cost(dim, quantum=0) = {
  \\ Start with a list of size N ~ kissing number for dimension
  log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3);
  \\ Test all pairs
  log2_test_cost = 0;
  if(quantum,
    log2_find_one_cost = .5 * log2_N + log2_test_cost,
    log2_find_one_cost = log2_N);
  log2_find_all_cost = log2_N + log2_find_one_cost;

  log2_ops = log2_find_all_cost;
  log2_space = log2_N;
  [log2_ops, log2_space];
}

bgj1_cost(dim, t, quantum=0) = {
  \\ Start with a list of size N ~ kissing number for dimension
  log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3);
  \\ Pick a random spherical cap of angle t.
  \\ Make a bucket of vectors in N that lie in this cap.
  \\ (Inclusion in spherical cap is tested using inner product.)
  log2_test_cost = 0; \\ XXX: reasonable value here
  log2_fill_cost = log2_N + log2_test_cost;
  log2_bucket_size = log2_N + log2_vol_cap(dim, t);
  \\ For each v in bucket, find w in bucket close to v (another inner product test)
  \\ (There are O(1) neighbors of v in bucket)
  if(quantum,
    log2_find_one_cost = .5 * log2_bucket_size + log2_test_cost, \\ Grover
    log2_find_one_cost = log2_bucket_size + log2_test_cost);     \\ Exhaustive
  log2_find_all_cost = log2_bucket_size + log2_find_one_cost;
  \\ Cost per filter is maximum of filling cost and checking cost
  log2_per_filter_cost = max(log2_fill_cost, log2_find_all_cost);
  \\ The filter only detects a pair (v,w) with probability given by the
  \\ normalized spherical measure of an intersection of caps of angle t
  \\ at angular distance Pi/3.
  log2_p = BDGL_log2_wedge(dim, t, t, Pi/3);
  \\ Repeat with random filters until we detect all pairs.
  log2_ops = log2_per_filter_cost - log2_p;
  log2_space = log2_N;
  [log2_ops, log2_space];
}

bgj1_cost_min(dim) = {
  my(log2_N, a, b);
  log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3);
  \\ TODO: fine tuning
  a = asin((3/4)^(1/4));
  b = asin((3/4)^(1/6));
  printf("Classical filter θ=%f, Ops: %d, Mem: %d , Bucket size: %d\n", a, bgj1_cost(dim,a), log2_N, (log2_N + log2_vol_cap(dim, a)));
  printf("Quantum filter θ=%f, Ops: %d, Mem: %d, Bucket size: %d", b, bgj1_cost(dim,b,1), log2_N, (log2_N + log2_vol_cap(dim, b)));
  bgj1_cost(dim,a);
}

BDGL_cost(dim, a, b, quantum=0) = {
  \\ Start with a list of size N ~ kissing number in dimension
  log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3);
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
    log2_find_one_cost = log2_relevant_filters + max(0, log2_bucket_size));
  \\ This has to be done for each v.
  log2_find_all_cost = log2_N + log2_find_one_cost;
  log2_ops = max(log2_insert_cost, log2_find_all_cost);
  log2_space = max(log2_N, log2_t + log2_bucket_size);
  [log2_ops, log2_space];
}

BDGL_cost_min(dim) = {
  my(log2_N, a, b);
  log2_N = JJP_log2_code_size_lower_bound(dim, Pi/3);
  \\ TODO: fine tuning
  a1 = acos(1/2);
  a2 = acos(sqrt(3/16));
  printf("Classical Ops: %d, Mem: %d , Bucket size: %d\n", BDGL_cost(dim,a1,a1), log2_N, (log2_N + log2_vol_cap(dim, a1)));
  printf("Quantum Ops: %d, Mem: %d, Bucket size: %d", BDGL_cost(dim,a2,a2,1), log2_N, (log2_N + log2_vol_cap(dim, a2)));
  BDGL_cost(dim,a1,a1,0);
}

