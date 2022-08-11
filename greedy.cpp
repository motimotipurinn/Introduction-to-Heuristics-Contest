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