#include <bits/stdc++.h>
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#define ll long long
#define int long long
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define ALL(a) (a).begin(), (a).end()
#define fore(i, a) for (auto &i : a)
#define MP make_pair
using namespace std;
using ld = long double;
using pll = pair<ll, ll>;
using pdd = pair<ld, ld>;
using Graph = vector<vector<ll>>;
using in = int;
using vin = vector<in>;
using vvin = vector<vector<int>>;
using PQS = priority_queue<int, vector<int>, greater<int>>;
using PQ = priority_queue<int>;
const ll MOD = 998244353;
const ll INF = 1LL << 60;
const string YYY = "YES";
const string yyy = "Yes";
const string NNN = "NO";
const string nnn = "No";
const int dx[4] = {1, 1, -1, -1};
const int dy[4] = {1, -1, 1, -1};
template <class T>
bool chmin(T &a, T b) {
  if (a > b) {
    a = b;
    return true;
  } else {
    return false;
  }
}
template <typename A, size_t N, typename T>
void Fill(A (&array)[N], const T &val) {
  std::fill((T *)array, (T *)(array + N), val);
}
template <class T>
bool chmax(T &a, T b) {
  if (a < b) {
    a = b;
    return true;
  } else {
    return false;
  }
}
bool cmp(const pll &a, const pll &b) { return a.second > b.second; }
ll GD(ll num) {  //各桁の和
  ll digit = 0;
  while (num != 0) {
    digit += num % 10;
    num /= 10;
  }
  return digit;
}
bool if_integer(ld x) {  //整数判定
  return std::floor(x) == x;
}
bool if_prime(ll x) {
  bool a = true;
  for (ll i = 2; i * i <= x; i++) {
    if (x % i == 0) {
      a = false;
      break;
    }
  }
  if (x == 1) a = false;
  return a;
}
ll gcd(ll x, ll y) {
  if (x == 0) {
    if (y == 0) {
      return 1;
    } else {
      return y;
    }
  }
  if (y == 0) {
    return x;
  }
  if (x % y == 0) {
    return (y);
  } else {
    return (gcd(y, x % y));
  }
}
ll lcm(int x, int y) { return x / gcd(x, y) * y; }
int GetDigit(int num) {
  int digit = 0;
  while (num != 0) {
    num /= 10;
    digit++;
  }
  return digit;
}
template <typename T>
vector<T> compress(vector<T> &X) {
  vector<T> vals = X;
  sort(vals.begin(), vals.end());
  vals.erase(unique(vals.begin(), vals.end()), vals.end());
  for (int i = 0; i < (int)X.size(); i++) {
    X[i] = lower_bound(vals.begin(), vals.end(), X[i]) - vals.begin();
  }
  return vals;
}

template <int MOD>
struct Fp {
  long long val;
  constexpr Fp(long long v = 0) noexcept : val(v % MOD) {
    if (val < 0) val += MOD;
  }
  constexpr int getmod() { return MOD; }
  constexpr Fp operator-() const noexcept { return val ? MOD - val : 0; }
  constexpr Fp operator+(const Fp &r) const noexcept { return Fp(*this) += r; }
  constexpr Fp operator-(const Fp &r) const noexcept { return Fp(*this) -= r; }
  constexpr Fp operator*(const Fp &r) const noexcept { return Fp(*this) *= r; }
  constexpr Fp operator/(const Fp &r) const noexcept { return Fp(*this) /= r; }
  constexpr Fp &operator+=(const Fp &r) noexcept {
    val += r.val;
    if (val >= MOD) val -= MOD;
    return *this;
  }
  constexpr Fp &operator-=(const Fp &r) noexcept {
    val -= r.val;
    if (val < 0) val += MOD;
    return *this;
  }
  constexpr Fp &operator*=(const Fp &r) noexcept {
    val = val * r.val % MOD;
    return *this;
  }
  constexpr Fp &operator/=(const Fp &r) noexcept {
    long long a = r.val, b = MOD, u = 1, v = 0;
    while (b) {
      long long t = a / b;
      a -= t * b;
      swap(a, b);
      u -= t * v;
      swap(u, v);
    }
    val = val * u % MOD;
    if (val < 0) val += MOD;
    return *this;
  }
  constexpr bool operator==(const Fp &r) const noexcept {
    return this->val == r.val;
  }
  constexpr bool operator!=(const Fp &r) const noexcept {
    return this->val != r.val;
  }
  friend constexpr ostream &operator<<(ostream &os, const Fp<MOD> &x) noexcept {
    return os << x.val;
  }
  friend constexpr Fp<MOD> modpow(const Fp<MOD> &a, long long n) noexcept {
    if (n == 0) return 1;
    auto t = modpow(a, n / 2);
    t = t * t;
    if (n & 1) t = t * a;
    return t;
  }
};
using mint = Fp<MOD>;

ll modpow(ll a, ll n) {
  ll res = 1;
  while (n) {
    if (n & 1) res = res * a % MOD;
    a = a * a % MOD;
    n >>= 1;
  }
  return res;
}
struct UnionFind {
  //自身が親であれば、その集合に属する頂点数に-1を掛けたもの
  //そうでなければ親のid
  vector<int> r;
  UnionFind(int N) { r = vector<int>(N, -1); }
  int root(int x) {
    if (r[x] < 0) return x;
    return r[x] = root(r[x]);
  }
  bool unite(int x, int y) {
    x = root(x);
    y = root(y);
    if (x == y) return false;
    if (r[x] > r[y]) swap(x, y);
    r[x] += r[y];
    r[y] = x;
    return true;
  }
  int size(int x) { return -r[root(x)]; }
  bool issame(int x, int y) { return root(x) == root(y); }
};
vector<pair<long long, long long>> prime_factorize(long long N) {
  vector<pair<long long, long long>> res;
  for (long long a = 2; a * a <= N; ++a) {
    if (N % a != 0) continue;
    long long ex = 0;
    while (N % a == 0) {
      ++ex;
      N /= a;
    }
    res.push_back({a, ex});
  }
  if (N != 1) res.push_back({N, 1});
  return res;
}
mint solve(in N, in K) {
  auto v = prime_factorize(K);
  in M = v.size();
  mint res = (1 + N) * N / 2;
  for (in i = 1; i < (1 << M); i++) {
    in num = 0;
    in l = 1;
    for (in j = i; j != 0; j >>= 1) num += j & 1;
    rep(j, M) {
      if (i >> j & 1) l *= v[j].first;
    }
    in x = N / l;
    mint p = (1 + x) * x / 2 * l;
    if (num % 2 == 1) {
      res -= p;
    } else {
      res += p;
    }
  }
  return res;
}
long long fac[300001], finv[300001], inv[300001];
void COMinit() {
  fac[0] = fac[1] = 1;
  finv[0] = finv[1] = 1;
  inv[1] = 1;
  for (int i = 2; i < 300001; i++) {
    fac[i] = fac[i - 1] * i % MOD;
    inv[i] = MOD - inv[MOD % i] * (MOD / i) % MOD;
    finv[i] = finv[i - 1] * inv[i] % MOD;
  }
}
ll comb(ll a, ll b) {
  if (a == 0 && b == 0) return 1;
  if (a < b || a < 0) return 0;
  return fac[a] * (finv[b] * finv[a - b] % MOD) % MOD;
}
ll evaluate(vin &c, vvin &d, vin &t, int K, int D) {
  in score = 0;
  vin last(26, -1);
  rep(i, t.size()) {
    score += d[i][t[i]];
    last[t[i]] = i;
    rep(j, 26) { score -= c[j] * (i - last[j]); }
  }
  for (in i = t.size(); i < min<long long>(D, t.size() + K); i++) {
    rep(j, 26) score -= (i - last[j]) * c[j];
  }
  return score;
}
ll get_score(vin &c, vvin &d, vin &t) {
  in score = 0;
  vin last(26, -1);
  rep(i, t.size()) {
    score += d[i][t[i]];
    last[t[i]] = i;
    rep(j, 26) score -= c[j] * (i - last[j]);
  }
  return score;
}
vin solve(int D, vin &c, vvin &d, int K) {
  vin t;
  rep(i, D) {
    int max_score = -INF;
    int p = -1;
    rep(j, 26) {
      t.push_back(j);
      in score = evaluate(c, d, t, K, D);
      if (chmax(max_score, score)) {
        p = j;
      }
      t.pop_back();
    }
    t.push_back(p);
  }
  return t;
}
signed main() {
  ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  cout << std::fixed << std::setprecision(50);
  ////////////////////////////
  ///////////////////////////
  in D;
  cin >> D;
  vin c(26);
  rep(i, 26) cin >> c[i];
  vvin d(D, vin(26));
  rep(i, D) rep(j, 26) { cin >> d[i][j]; }
  int ans_score = -INF;
  vin ans;
  rep(i, 26) {
    auto t = solve(D, c, d, i);
    in score = get_score(c, d, t);
    if (chmax(ans_score, score)) {
      ans = t;
    }
  }
  rep(i, D) { cout << ans[i] + 1 << endl; }
}
//
// rm -rf test/ && oj d
// g++ main.cpp
//  oj t --ignore-spaces-and-newline