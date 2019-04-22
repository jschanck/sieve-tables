beta(a,b) = {
  gamma(a)*gamma(b)/gamma(a+b);
}

betainc(x,a,b) = {
  \\ https://dlmf.nist.gov/8.17#E7
  (x^a/a) * hypergeom([a,1-b],a+1,x) / beta(a,b);
  \\ https://dlmf.nist.gov/8.17#E8
  \\ (x^a*(1-x)^b/a) * hypergeom([a+b,1],a+1,x) / beta(a,b);
}

log2_gamma(d) = {
  lngamma(d)/log(2);
}

log2_sphere(d) = {
  (d/2*log2(Pi) + 1) - log2_gamma(d/2);
}

log2_cap_asymp(d,a) = {
  \\ Normalized spherical measure of a cap of angle a.
  (d-1)*log2(sin(a)) - 1/2 * log2(2*Pi*d) - log2(cos(a));
}

log2_cap_intnum(d,a) = {
  my(s);
  s = log2(intnum(x=0,a,sin(x)^(d-2)));
  s = s - log2(Pi)/2 + (log2_gamma(d/2) - log2_gamma((d-1)/2));
}

log2_cap_beta(d, a) = {
  log2(betainc(sin(a)^2, (d-1)/2, 1/2)) - 1;
}

log2_wedge_intnum(d,a,b,t=Pi/3) = {
  my(c,Sa,Sb);
  \\ The parameters a and b are the angles of spherical caps with
  \\ centers A and B that are at angular distance t (<A, B> = cos(t)).
  \\ Under suitable restrictions on a, b, and t, this function
  \\ computes the normalized (d-1)-dimensional volume of the
  \\ intersection of the caps. The conditions are: t < a+b;
  \\ 0 < a < Pi/2; 0 < b < Pi/2; and A and B must lie in opposite
  \\ half-spaces of the hyperplane orthogonal to (A/cos(a) - B/cos(b)).
  \\ The final condition can be restated as
  \\    (cos(b) - cos(a)cos(t))(cos(b)cos(t) - cos(a)) < 0.
  \\ This corresponds to Case 8 in
  \\ [LK14] Yongjae Lee, Woo Chang Kim.
  \\ "Concise formulas for the surface area of the intersection of two hyperspherical caps."
  \\ KAIST TECHNICAL REPORT.
  \\ https://ie.kaist.ac.kr/
  \\  uploads/professor/tech_file/
  \\  Concise+Formulas+for+the+Surface+Area+of+the+Intersection+of+Two+Hyperspherical+Caps.pdf

  \\ TODO: Verify input conditions.

  \\ The hyperplane orthogonal to (A/cos(a) - B/cos(b)) makes
  \\ an angle c with B where c is defined by
  \\ cos(b)*sec(c) = cos(a)*sec(t-c).
  c = atan(cos(a)/(cos(b)*sin(t)) - 1/tan(t));

  \\ For x between c and b, the portion of the wedge on the A side of the hyperplane
  \\ at height cos(x) along B is a (d-2)-dimensional spherical cap of radius sin(x)
  \\ and angle acos(tan(c) / tan(x)).
  if(d==2, Sb = intnum(x=c, b, 1),
           Sb = 1/2 * intnum(x=c, b, sin(x)^(d-2)*betainc(sin(acos(tan(c)/tan(x)))^2, (d-2)/2, 1/2)));
           \\Sb = intnum(x=c, b, sin(x)^(d-2)*2^log2_cap_asymp(d-1,acos(tan(c)/tan(x)))));

  \\ Likewise, for x between t-c and a, the portion of the wedge on the B side of the hyperplane
  \\ at height cos(x) along A is a (d-2)-dimensional spherical cap of radius sin(x)
  \\ and angle acos(tan(t-c) / tan(x)).
  if(d==2, Sa = intnum(x=t-c, a, 1),
           Sa = 1/2 * intnum(x=t-c, a, sin(x)^(d-2)*betainc(sin(acos(tan(t-c)/tan(x)))^2, (d-2)/2, 1/2)));
           \\Sa = intnum(x=t-c, a, sin(x)^(d-2)*2^log2_cap_asymp(d-1,acos(tan(t-c)/tan(x)))));

  \\ The wedge volume is Sa + Sb, up to normalization.
  log2(Sa + Sb) + log2_sphere(d-1) - log2_sphere(d);
}

log2_wedge_bdgl(d,a,b,t=Pi/3) = {
  my(g_sq, log2_A);
  \\ Wedge volume formula from Lemma 2.2 of
  \\   [BDGL16] Anja Becker, Léo Ducas, Nicolas Gama, Thijs Laarhoven.
  \\   "New directions in nearest neighbor searching with applications
  \\    to lattice sieving." SODA 2016.
  \\   https://eprint.iacr.org/2015/1128
  g_sq = (cos(a)^2 + cos(b)^2 - 2*cos(a)*cos(b)*cos(t))/sin(t)^2;
  log2_A = log2(g_sq) - 2*log2(1-g_sq);
  \\ XXX: Is (d-4) right?
  (d-4) * log2(sqrt(1-g_sq)) + log2_A - 2*log2(d-4) + log2_sphere(d-2) - log2_sphere(d);
}

log2_CSW(d, a) = {
  \\ Chabauty--Shannon--Wyner lower bound on the maximum number
  \\ of points that can be placed on the d-1 sphere such that
  \\ the angular distance between any two points is at least a.
  -log2_cap_beta(d,a);
}

log2_JJP(d,a) = {
  my(r, log2_c);
  \\ Jenssen--Joos--Perkins lower bound on the maximum number
  \\ of points that can be placed on the d-1 sphere such that
  \\ the angular distance between any two points is at least a.
  \\ Linear improvement over CSW lower bound for any fixed angle.
  \\   [JJP18] Matthew Jenssen, Felix Joos, Will Perkins.
  \\   "On kissing numbers and spherical codes in high dimensions."
  \\   Advances in Mathematics, Volume 335, 2018.
  \\   https://arxiv.org/abs/1803.02702
  \\ A(d,a) >= (1+o(1)) * ln(sin(a)/sin(q(a))) * d / cap_volume(d,a)
  \\ where q(a) = asin(sqrt(r) / sin(a)) for
  r = (cos(a) - 1)^2 * (1 + 2*cos(a));
  \\ is the angle of the smallest cap that contains the intersection
  \\ of caps of angle a with centers at angular distance a.
  \\ XXX: This function ignores the o(1), which comes from the asymptotics
  \\ of the lambert W function W(z) = log(z) - log(log(z)) + o(1).
  \\ TODO: compute the o(1) term.
  log2_c = log2(log(sin(a)) - log(sqrt(r) / sin(a)));
  -log2_cap_beta(d,a) + log2(d) + log2_c;
}

