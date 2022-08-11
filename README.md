
# cppコードの説明

"greedy.cpp"は各日+k日後(0<=k<=26)におけるscoreの最大化を狙った貪欲。
https://atcoder.jp/contests/intro-heuristics/submissions/33926635
https://atcoder.jp/contests/intro-heuristics/submissions/33928087
"lo_search.cpp"は局所探索法で先ほどの貪欲解からそれぞれ1/2の確率でswapと各点changeを行います
https://atcoder.jp/contests/intro-heuristics/submissions/33929077
"main.cpp"も先ほどと同じ貪欲解をベースに焼きなましをしたものです
https://atcoder.jp/contests/intro-heuristics/submissions/33936800

# 各種コード使い方例

```
python3 generator.py 1234 > input.txt
g++ main.cpp
./a.out < input.txt > output.txt
python3 tester.py input.txt output.txt
firefox vis.html
```
