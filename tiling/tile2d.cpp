#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

//
// n0 = (1, 1, 0)
// n1 = (1, 0, 1)
// n2 = (1, -1, -1)
//

inline int roundUp(int val, int mod) {
	return (val + (mod - 1) * (val > 0)) / mod * mod;
}

inline int roundDown(int val, int mod) {
	return (val - (mod - 1) * (val < 0)) / mod * mod;
}

inline int mod(int val, int mod) {
	int res = val % mod;
	return res + (res < 0) * mod;
}

bool check(std::vector<int>& a, int i, int j, int t, int X, int Y) {
	bool b = (a[i + X * j] == t + 1 || a[i + X * j] == t);
	if (1 <= i) {
		b = b && (a[(i - 1) + X * j] == t + 1 || a[(i - 1) + X * j] == t);
	}
	if (i < X - 1) {
		b = b && (a[(i + 1) + X * j] == t + 1 || a[(i + 1) + X * j] == t);
	}
	if (1 <= j) {
		b = b && (a[i + X * (j - 1)] == t + 1 || a[i + X * (j - 1)] == t);
	}
	if (j < Y - 1) {
		b = b && (a[i + X * (j + 1)] == t + 1 || a[i + X * (j + 1)] == t);
	}
	a[i + X * j] += 1;
	return b;
}

int main() {
	int Lx = 10;
	int Ux = 28;
	int Ly = -8;
	int Uy = 11;
	int Lt = 5;
	int Ut = 15;
	int h = 7;

	std::vector<int> data((Ux - Lx) * (Uy - Ly), 0);

	for (int Ph = roundUp(3 * (Lt - h + 1), h); Ph <= roundDown(3 * Ut, h); Ph += h) {
		std::cerr << Ph << ":";
		int ah3_top = roundDown(3 * Ux + Ph + 2 * h - 2, 3 * h);
		int ah3_bot = roundUp(3 * Lx + Ph - 2 * h + 2, 3 * h);
		int bh3_top = roundDown(3 * Uy + Ph + 2 * h - 2, 3 * h);
		int bh3_bot = roundUp(3 * Ly + Ph - 2 * h + 2, 3 * h);
		for (int bh3 = bh3_bot; bh3 <= bh3_top; bh3 += 3 * h) {
			for (int ah3 = ah3_bot; ah3 <= ah3_top; ah3 += 3 * h) {
				assert(ah3 % (3 * h) == 0);
				int a = ah3 / (3 * h);
				assert(bh3 % (3 * h) == 0);
				int b = bh3 / (3 * h);
				assert(Ph % h == 0);
				int c = Ph / h - a - b;
				std::cerr << " (" << a << ',' << b << ',' << c << ')';
				//double t = (a + b + c) * h / 3.0;
				//double i = (2 * a - b - c) * h / 3.0;
				//double j = (2 * b - a - c) * h / 3.0;
				//std::cout << i << ' ' << j << ' ' << t << '\n';
				int t3_bot = std::max(roundUp((a + b + c) * h, 3), 3 * Lt);
				int t3_top = std::min(roundDown((a + b + c + 3) * h - 3, 3), 3 * (Ut - 1));
				//int t3_bot = roundUp((a + b + c) * h, 3);
				//int t3_top = roundDown((a + b + c + 3) * h - 3, 3);
				for (int t3 = t3_bot; t3 <= t3_top; t3 += 3) {
					assert(t3 % 3 == 0);
					int t = t3 / 3;
					int j_bot = std::max({-t + b * h, Ly});
					int j_top = std::min({-t + (b + 1) * h - 1, Uy - 1});
					for (int j = j_bot; j <= j_top; ++j) {
						//int i_bot = std::max(-t + a * h, t - j - (c + 1) * h + 1);
						//int i_top = std::min(-t + (a + 1) * h - 1, t - j  - c * h);
						int i_bot = std::max({-t + a * h, t - j - (c + 1) * h + 1, Lx});
						int i_top = std::min({-t + (a + 1) * h - 1, t - j  - c * h, Ux - 1});
						for (int i = i_bot; i <= i_top; ++i) {
							assert(a * h <= t + i);
							assert(t + i < (a + 1) * h);
							assert(b * h <= t + j);
							assert(t + j < (b + 1) * h);
							assert(c * h <= t - i - j);
							assert(t - i - j < (c + 1) * h);
							assert(check(data, i - Lx, j - Ly, t - Lt, Ux - Lx, Uy - Ly));
							std::cout << i << ' ' << j << ' ' << t << '\n';
						}
					}
				}
			}
		}
		std::cerr << '\n';
	}

	for (int i : data) {
		assert(i == Ut - Lt);
	}
}
