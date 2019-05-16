#ifndef GEN_H
#define GEN_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <array>
#include <chrono>

using Time = std::chrono::time_point<std::chrono::system_clock>;
using Clock = std::chrono::system_clock;

template <typename Duration>
auto as_ms(Duration const& duration) {
	return std::chrono::duration_cast<std::chrono::milliseconds>(duration);
}

using Int = int;

inline Int roundUp(Int value, Int mod) {
	return (value + (mod - 1) * (value > 0)) / mod;
}

inline Int roundDown(Int value, Int mod) {
	return (value - (mod - 1) * (value < 0)) / mod;
}

inline long long nextPower2(long long int value) {
	Int res = 1;
	while (res < value) {
		res *= 2;
	}
	return res;
}

inline bool check(std::vector<Int>& a, Int i, Int j, Int t, Int X, Int Y) {
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

struct RawGen {
	const Int T;
	const Int X;
	const Int Y;

	const Int N;
	const long long H;

	const long long total;

	RawGen(Int T, Int X, Int Y)
		: T(T)
		, X(X)
		, Y(Y)
		, N(roundUp(std::max(std::max(X - 1, Y - 1), T), 2))
		, H(nextPower2(6ll * N))
		, total(H * H * H / 3ll + (2ll - (H * H * H) % 3ll))
		, t(0)
		, i(0)
		, j(0)
	{
		assert(6ll * N <= H);
		assert(X <= 2 * N + 1);
		assert(Y <= 2 * N + 1);
		assert(T <= 2 * N);
		std::cerr << "can be used N = " << roundDown(H, 6) << std::endl;
	}

	bool enabled() {
		return
			-N <= i && i <= N &&
			-N <= j && j <= N &&
			-N <= t - H/2 && t - H/2 <= N;
	}

	bool inside() {
		return
			0 <= get_i() && get_i() < X &&
			0 <= get_j() && get_j() < Y &&
			0 <= get_t() && get_t() < T;
	}

	Int get_t() {
		return t - H/2 + N;
	}

	Int get_i() {
		return i + N;
	}

	Int get_j() {
		return j + N;
	}

	long long numTotal() const {
		return total;
	}

	void next() {
		Int h = 1;
		Int a = roundDown(t + i, h);
		Int b = roundDown(t + j, h);
		Int c = roundDown(t - i - j, h);
		while (1) {
			Int aa = a & 1;
			Int bb = b & 1;
			Int cc = c & 1;
			while (aa && bb && cc) {
				h *= 2;
				a = roundDown(t + i, h);
				b = roundDown(t + j, h);
				c = roundDown(t - i - j, h);
				aa = a & 1;
				bb = b & 1;
				cc = c & 1;
			}
			//std::cerr << a << ',' << b << ',' << c << ',' << h << " >>> ";

			if (!aa && !bb && !cc) {
				a += 1;
			} else if (aa && !bb && !cc) {
				a -= 1;
				b += 1;
			} else if (!aa && bb && !cc) {
				b -= 1;
				c += 1;
			} else if (!aa && !bb && cc) {
				b += 1;
			} else if (!aa && bb && cc) {
				b -= 1;
				a += 1;
			} else if (aa && !bb && cc) {
				c -= 1;
				b += 1;
			} else if (aa && bb && !cc) {
				c += 1;
			}

			while (h > 1) {
				h /= 2;
				a *= 2;
				b *= 2;
				c *= 2;
			}

			//Int nt = roundUp((a + b + c) * h, 3);
			//bool bt = nt <= roundDown((a + b + c + 3) * h, 3) - 1;
			//Int nj = -nt + b * h;
			//bool bj = nj <= -nt + (b + 1) * h - 1;
			//Int ni = std::max(-nt + a * h, nt - nj - (c + 1) * h + 1);
			//bool bi = ni <= std::min(-nt + (a + 1) * h - 1, nt - nj - c * h);

			if ((a + b + c) % 3) {
				continue;
			}
			Int nt = (a + b + c) / 3;
			Int nj = -nt + b;
			Int ni1 = -nt + a;
			Int ni2 = nt - nj - c;
			if (ni1 != ni2) {
				continue;
			}
			Int ni = ni1;
			//std::cerr << a << ',' << b << ',' << c << ',' << h << ' ' << ((bt && bj && bi) ? "ok" : "DROP") << std::endl;

			t = nt;
			i = ni;
			j = nj;
			break;
			//if (bt && bj && bi) {
			//	t = nt;
			//	i = ni;
			//	j = nj;
			//	break;
			//}
		}
	}

	Int t;
	Int i;
	Int j;
};

struct CachedGen {
	CachedGen(Int T, Int X, Int Y) : NUM(T * X * Y) {
		//fs.open("cache.bin", std::ios_base::binary | std::ios_base::out);
		//index = 0;
		//RawGen gen(T, X, Y);
		////std::vector<Int> test(X * Y, 0);
		////Int debugCounter = 0;
		//for (long long k = 0; k < gen.numTotal(); ++k) {
		//	if (gen.inside()) {
		//		write(gen.get_i());
		//		write(gen.get_j());
		//		write(gen.get_t());
		//		//assert(check(test, gen.get_i(), gen.get_j(), gen.get_t(), X, Y));
		//		//std::cout << gen.get_i() << ' ' << gen.get_j() << ' ' << gen.get_t() << std::endl;
		//		//debugCounter++;
		//	}
		//	gen.next();
		//}
		////std::cerr << debugCounter << " == " << NUM << std::endl;
		//flush();
		//fs.close();

		fs.open("cache.bin", std::ios_base::binary | std::ios_base::in);
		index = BUF_SIZE;
	}

	void next() {
		i = read();
		j = read();
		t = read();
	}
	Int read() {
		if (index == BUF_SIZE) {
			fs.read(reinterpret_cast<char*>(buf.data()), BUF_SIZE * sizeof(Int));
			index = 0;
		}
		return buf[index++];
	}
	void write(Int value) {
		buf[index++] = value;
		if (index == BUF_SIZE) {
			fs.write(reinterpret_cast<char*>(buf.data()), BUF_SIZE * sizeof(Int));
			index = 0;
		}
	};
	void flush() {
		if (index) {
			fs.write(reinterpret_cast<char*>(buf.data()), index * sizeof(Int));
		}
	}

	Int num() const {
		return NUM;
	}

	Int index = 0;

	std::fstream fs;

	static const Int BUF_SIZE = 1024;
	std::array<Int, BUF_SIZE> buf;

	const Int NUM;

	Int t = -1;
	Int i = -1;
	Int j = -1;
};

#endif
