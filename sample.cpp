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
constexpr int types = 26;
std::mt19937 mt{std::random_device{}()};

ll get_score(vin &c, vvin &S, vin &t) {
  int score = 0;
  vin prev(types, -1);
  for (int i = 0; i < t.size(); i++) {
    score += S[i][t[i]];
    prev[t[i]] = i;
    for (int j = 0; j < types; j++) score -= c[j] * (i - prev[j]);
  }
  return score;
}

ll evaluate(vin &c, vvin &S, vin &t, int k, int D) {
  ll score = 0;
  vin prev(types, -1);
  for (int i = 0; i < t.size(); i++) {
    score += S[i][t[i]];
    prev[t[i]] = i;
    for (int j = 0; j < types; j++) score -= c[j] * (i - prev[j]);
  }

  for (ll d = t.size(); d < min<int>(D, t.size() + k); d++) {
    for (int j = 0; j < types; j++) {
      score -= (d - prev[j]) * c[j];
    }
  }
  return score;
}

vin solve(int D, vin &c, vvin &S, int k) {
  vin out;

  for (int i = 0; i < D; i++) {
    int max_score = LLONG_MIN;
    int type = -1;

    for (int j = 0; j < types; j++) {
      out.push_back(j);
      int score = evaluate(c, S, out, k, D);
      if (score > max_score) {
        max_score = score;
        type = j;
      }
      out.pop_back();
    }
    out.push_back(type);
  }

  return out;
}

constexpr double TIME_LIMIT = 2000;

signed main() {
  auto start = chrono::system_clock::now();

  int D;
  cin >> D;

  vin c(types);
  rep(i, 26) cin >> c[i];

  vvin S(D, vin(types));
  rep(i, D) rep(j, 26) cin >> S[i][j];

  vin ans;
  int ans_score = LLONG_MIN;
  for (int k = 0; k < types; k++) {
    auto out = solve(D, c, S, k);
    int score = get_score(c, S, out);
    if (score > ans_score) {
      ans_score = score;
      ans = out;
    }
  }

  const double T0 = 2000;
  const double T1 = 600;
  double T = T0;

  int cnt = 0;

  // start = chrono::system_clock::now();
  //  search
  while (true) {
    auto now = chrono::system_clock::now();
    auto elapsed =
        chrono::duration_cast<chrono::milliseconds>(now - start).count();
    if (elapsed >= 1940) break;

    uniform_int_distribution<int> rand_day(0, D - 1);
    uniform_int_distribution<int> rand_type(0, types - 1);
    uniform_real_distribution<double> rand_01(0, 1);

    cnt++;

    if (cnt % 50 == 0) {
      auto t = elapsed / TIME_LIMIT;
      T = pow(T0, 1 - t) * pow(T1, t);
    }

    if (rand_01(mt) < 0.5) {
      auto day = rand_day(mt);
      auto type = rand_type(mt);
      auto prev_val = ans[day];
      ans[day] = type;
      int score = get_score(c, S, ans);
      auto p = exp((score - ans_score) / T);
      auto r = rand_01(mt);
      if (score > ans_score or r < p) {
        ans_score = score;
      } else {
        ans[day] = prev_val;
      }
    } else {
      auto day1 = max<ll>(0, rand_day(mt) - 1);
      uniform_int_distribution<int> rand_day2(day1 + 1, min(D - 1, day1 + 16));
      auto day2 = rand_day2(mt);
      swap(ans[day1], ans[day2]);
      int score = get_score(c, S, ans);
      auto p = exp((score - ans_score) / T);
      auto r = rand_01(mt);
      if (score > ans_score or r < p) {
        ans_score = score;
      } else {
        swap(ans[day1], ans[day2]);
      }
    }
  }

  for (auto x : ans) cout << x + 1 << endl;
}
//
// rm -rf test/ && oj d
// g++ main.cpp
//  oj t --ignore-spaces-and-newline