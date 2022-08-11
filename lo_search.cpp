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
class Rand {
 private:
  // 32ビット版メルセンヌ・ツイスタ
  std::mt19937 mt;
  //非決定論的な乱数
  std::random_device rd;

 public:
  //コンストラクタ(初期化)
  Rand() { mt.seed(rd()); }

  //初期値
  void seed() { mt.seed(rd()); }
  void seed(const std::uint_fast32_t seed_) { mt.seed(seed_); }

  //通常の乱数
  std::uint_fast32_t operator()() { return mt(); }
  // 0～最大値-1 (余りの範囲の一様分布乱数)
  std::int_fast32_t operator()(const std::int_fast32_t max_) {
    std::uniform_int_distribution<> uid(
        0, ((max_ > 0) ? (std::int_fast32_t)max_ - 1 : 0));
    return uid(mt);
  }
  //最小値～最大値
  std::int_fast32_t operator()(const std::int_fast32_t min_,
                               const std::int_fast32_t max_) {
    std::uniform_int_distribution<> uid((min_ <= max_) ? min_ : max_,
                                        (min_ <= max_) ? max_ : min_);
    return uid(mt);
  }
  //確率
  bool randBool(const double probability_) {
    std::bernoulli_distribution uid(probability_);
    return uid(mt);
  }
  bool randBool() {
    std::uniform_int_distribution<> uid(0, 1);
    return ((uid(mt)) ? true : false);
  }
};
static thread_local Rand rnd;
signed main() {
  ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  cout << std::fixed << std::setprecision(50);
  ////////////////////////////
  ///////////////////////////
  std::chrono::system_clock::time_point start, end;
  start = std::chrono::system_clock::now();

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
  while (1) {
    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    if (elapsed < 1800) {
      int dice = rnd(2);
      if (dice == 0) {
        int p = rnd(D), q = rnd(26);
        int old = ans[p];
        ans[p] = q;
        int new_score = get_score(c, d, ans);
        if (!(chmax(ans_score, new_score))) {
          ans[p] = old;
        }
      } else {
        int p = rnd(D), q = rnd(1, 16);
        if (p + q < D) {
          swap(ans[p], ans[p + q]);
          int new_score = get_score(c, d, ans);
          if (!(chmax(ans_score, new_score))) {
            swap(ans[p], ans[p + q]);
          }
        }
      }
    } else {
      break;
    }
  }
  rep(i, D) { cout << ans[i] + 1 << endl; }
}