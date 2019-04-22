read("spheres.gp")

\\ Digits of precision
\p 57

\\ Use exact cap and wedge volumes.
log2_cap    = log2_cap_beta;
log2_wedge  = log2_wedge_intnum;
\\ Use asymptotic volume estimates.
\\log2_cap    = log2_cap_asymp;
\\log2_wedge  = log2_wedge_bdgl;

\\ Use 1/(cap volume) for list size
log2_list_size = log2_CSW;
\\ Use a lower bound on the kissing constant for the list size
\\log2_list_size = log2_JJP;


NV_cost(dim, quantum=0) = {
  my(log2_N, log2_test_cost, log2_find_one_cost, log2_find_all_cost,\
     log2_ops, log2_space);
  \\ Nguyen-Vidick sieve
  \\ Start with a list of size N
  log2_N = log2_list_size(dim, Pi/3);
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

bgj1_cost(d, t, quantum=0) = {
  my(log2_N, log2_test_cost, log2_fill_cost, log2_bucket_size,\
     log2_find_one_cost, log2_find_all_cost, log2_per_filter_cost,\
     log2_p, log2_ops, log2_space);
  \\ Start with a list of size N
  log2_N = log2_list_size(d, Pi/3);
  \\ Pick a random spherical cap of angle t.
  \\ Bucket the vectors that lie in this cap.
  \\ (Inclusion is an inner product test.)
  log2_test_cost = 0;
  log2_fill_cost = log2_N + log2_test_cost;
  log2_bucket_size = log2_N + log2_cap(d, t);
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
  log2_p = log2_wedge(d, t, t);
  \\ Repeat with random filters until all pairs are detected.
  log2_ops = log2_per_filter_cost - log2_p;
  log2_space = log2_N;
  [log2_ops, log2_space];
}

BDGL_cost(d, a, b, quantum=0) = {
  my(log2_N, log2_t, log2_insert_cost, log2_bucket_size, \
     log2_relevant_filters, log2_find_one_cost, log2_find_all_cost,\
     log2_ops, log2_space);
  \\ Start with a list of size N
  log2_N = log2_list_size(d, Pi/3);
  \\ Initialize t := t(a,b) filters. Each filter defines a bucket.
  log2_t = -log2_wedge(d, a, b);
  \\ Each of the N vectors is inserted into a number of buckets
  \\ that is proportional to the measure of a cap of angle b.
  log2_insert_cost = log2_N + log2_cap(d, b) + log2_t;
  log2_bucket_size = log2_insert_cost - log2_t;
  \\ Each vector is compared against a number of buckets
  \\ proportional to the measure of a cap of angle a.
  log2_relevant_filters = log2_cap(d, a) + log2_t;
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

local_min(f,x,E=[1,5]) = {
  \\ Search left f(x-ε) and right f(x+ε) for a local minimum near f(x).
  \\ ε will be determined to 10^-E[2], starting with steps of size 10^-E[1].
  my(e,y, y2, k);
  y = f(x);
  for(k=E[1], E[2],
    e=0.1^k;
    while(y2 = f(x-e); y2 < y, x=x-e; y = y2);
    while(y2 = f(x+e); y2 < y, x=x+e; y = y2);
  );
  x
}

bgj1_cost_min(dim,guess=asin((3/4)^(1/4))) = {
  my(t, E);

  E = [3,3];

  t = local_min(x->(bgj1_cost(dim,x,0)[1]), guess, E);

  \\if(quantum,
  \\  \\ Asymptotically optimal filter angle for quantum bgj1 is asin((3/4)^(1/6)) = 1.2635...
  \\  t = local_min(x->(bgj1_cost(dim,x,1)[1]), asin((3/4)^(1/6)), E),
  \\  \\ Asymptotically optimal filter angle for classical bgj1 is asin((3/4)^(1/4)) = 1.1960...
  \\  t = local_min(x->(bgj1_cost(dim,x,0)[1]), asin((3/4)^(1/4)), E));

  [t, bgj1_cost(dim,t)[1]];
}

BDGL_cost_min(dim,guess=Pi/3) = {
  my(t, e1, e2);

  E = [3,3];
  t = Pi/3; \\local_min(x->(BDGL_cost(dim,x,x,0)[1]),guess, E);

  \\if(quantum,
  \\  \\ Asymptotically optimal filter angle for quantum BDGL is acos(sqrt(3/16)) = 1.2296...
  \\  t = local_min(x->(BDGL_cost(dim,x,x,1)[1]), acos(sqrt(3/16)), E),
  \\  \\ Asymptotically optimal filter angle for classical BDGL is Pi/3 = 1.0471...
  \\  t = local_min(x->(BDGL_cost(dim,x,x,0)[1]), acos(1/2), E));

  [t, BDGL_cost(dim,t,t)[1]];
}

gen_tables() = {
  my(bl,bh,s,b,t);

  bl = 50;
  bh = 850;
  s = "";
  for(b=bl,bh,
    [t, c] = BDGL_cost_min(b);
    s = concat(s, Strprintf("%d, %.4f, %.4f, %.1f\n", b, cos(t), cos(t), c));
    printf("%d, %.4f, %.4f, %.1f\n", b, cos(t), cos(t), c);
    );
  write("2924.csv", s);

  s = "";
  t = 1.238; \\ Initial guess for b=50
  for(b=bl,bh,
    [t, c] = bgj1_cost_min(b,t);
    s = concat(s, Strprintf("%d, %.4f, %.1f\n", b, cos(t), c));
    printf("%d, %.4f, %.1f\n", b, cos(t), c));
  write("3494.csv", s);

  s = "";
  for(b=bl,bh,
    s = concat(s, Strprintf("%d, %.1f\n", b, log2_list_size(b, Pi/3))));
  write("2075.csv", s);
}

;
