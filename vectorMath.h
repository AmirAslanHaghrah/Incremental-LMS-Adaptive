#include <random>
#include <ctime>
#include <vector>
#include <cmath>
#include <ppl.hh>


double l2_norm(std::vector<std::vector<double>> const u);
std::vector<std::vector<double>> diag(std::vector<std::vector<double>> const u);
std::vector<std::vector<double>> eye(int m);
std::vector<std::vector<double>> ones(int m, int n);
std::vector<std::vector<double>> SQRT(std::vector<std::vector<double>> const u);
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> u);
std::vector<std::vector<double>> kron(const std::vector<std::vector<double>> A, const std::vector<std::vector<double>> B);

double l2_norm(std::vector<std::vector<double>> const u) {
	double accum = 0.0;
	for (int i = 0; i < (int)u.size(); i++) {
		for (int j = 0; j < (int)u[0].size(); j++) {
			accum += u[i][j] * u[i][j];
		}
	}
	return sqrt(accum);
}

std::vector<std::vector<double>> diag(std::vector<std::vector<double>> const u) {
	std::vector<std::vector<double>> result;
	if (u.size() == 1) {
		result.resize(u[0].size(), std::vector<double>(u[0].size(), 0.0));
		for (int i = 0; i < (int)u[0].size(); i++) {
			result[i][i] = u[0][i];
		}
	}
	else {
		result.resize(u.size(), std::vector<double>(u.size(), 0.0));
		for (int i = 0; i < (int)u.size(); i++) {
			result[i][i] = u[i][0];
		}
	}
	return result;
}

double trace(std::vector<std::vector<double>> const u) {
	int m = (int)u.size();
	double result = 0;
	for (int i = 0; i < m; i++) {
		result += u[i][i];
	}
	return result;
}

std::vector<std::vector<double>> eye(int m) {
	std::vector<std::vector<double>> result;
	result.resize(m, std::vector<double>(m, 0.0));


	for (int i = 0; i < m; i++) {
		result[i][i] = 1.0;
	}
	return result;
}

std::vector<std::vector<double>> ones(int m, int n) {
	std::vector<std::vector<double>> result;
	result.resize(m, std::vector<double>(n, 1.0));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result[i][j] = 1.0;
		}
	}
	return result;
}

std::vector<std::vector<double>> zeros(int m, int n) {
	std::vector<std::vector<double>> result;
	result.resize(m, std::vector<double>(n, 1.0));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result[i][j] = 0.0;
		}
	}
	return result;
}


std::vector<std::vector<double>> SQRT(std::vector<std::vector<double>> const u) {
	std::vector<std::vector<double>> result;
	result.resize(u.size(), std::vector<double>(u.size(), 0.0));
	
	#pragma omp parallel for
	for (register int i = 0; i < (int)u.size(); i++) {
		#pragma omp parallel for
		for (register int j = 0; j < (int)u[0].size(); j++) {
			result[i][j] = sqrt(u[i][j]);
		}
	}
	return result;
}

std::vector<std::vector<double>> operator + (const std::vector<std::vector<double>>& left, const std::vector<std::vector<double>>& right) {
	std::vector<std::vector<double>> result;
	result.resize(left.size(), std::vector<double>(left[0].size(), 0.0));
	
	#pragma omp parallel for
	for (register int i = 0; i < (int)left.size(); i++) {
		#pragma omp parallel for
		for (register int j = 0; j < (int)left[0].size(); j++) {
			result[i][j] = left[i][j] + right[i][j];
		}
	}
	return result;
}

std::vector<std::vector<double>> operator - (const std::vector<std::vector<double>>& left, const std::vector<std::vector<double>>& right) {
	std::vector<std::vector<double>> result;
	result.resize(left.size(), std::vector<double>(left[0].size(), 0.0));

	#pragma omp parallel for
	for (register int i = 0; i < (int)left.size(); i++) {
		#pragma omp parallel for
		for (register int j = 0; j < (int)left[0].size(); j++) {
			result[i][j] = left[i][j] - right[i][j];
		}
	}
	return result;
}

std::vector<std::vector<double>> operator * (const std::vector<std::vector<double>>& left, const std::vector<std::vector<double>>& right) {	
	std::vector<std::vector<double>> result;
	result.resize(left.size(), std::vector<double>(right[0].size(), 0.0));

	#pragma omp parallel for
	for (register int i = 0; i < (int)left.size(); i++) {
		for (register int j = 0; j < (int)right[0].size(); j++) {
			for (register int k = 0; k < (int)left[0].size(); k++) {
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return result;
}

std::vector<std::vector<double>> operator * (const double& left, const std::vector<std::vector<double>>& right) {
	std::vector<std::vector<double>> result;
	result.resize(right.size(), std::vector<double>(right[0].size(), 0.0));

	#pragma omp parallel for
	for (register int i = 0; i < (int)right.size(); i++) {
		#pragma omp parallel for
		for (register int j = 0; j < (int)right[0].size(); j++) {
			result[i][j] = left * right[i][j];
		}
	}
	return result;
}

std::vector<std::vector<double>> operator * (const std::vector<std::vector<double>>& left, const double& right) {
	std::vector<std::vector<double>> result;
	result.resize(left.size(), std::vector<double>(left[0].size(), 0.0));

	#pragma omp parallel for
	for (register int i = 0; i < (int)left.size(); i++) {
		#pragma omp parallel for
		for (register int j = 0; j < (int)left[0].size(); j++) {
			result[i][j] = left[i][j] * right;
		}
	}
	return result;
}

std::vector<std::vector<double>> operator / (const std::vector<std::vector<double>>& left, const double& right) {
	std::vector<std::vector<double>> result;
	result.resize(left.size(), std::vector<double>(left[0].size(), 0.0));

	#pragma omp parallel for
	for (register int i = 0; i < (int)left.size(); i++) {
		#pragma omp parallel for
		for (register int j = 0; j < (int)left[0].size(); j++) {
			result[i][j] = left[i][j] / right;
		}
	}
	return result;
}

std::vector<std::vector<double>> operator ^ (const std::vector<std::vector<double>>& left, const double& right) {
	std::vector<std::vector<double>> result;
	result.resize(left.size(), std::vector<double>(left[0].size(), 0.0));

	#pragma omp parallel for
	for (register int i = 0; i < (int)left.size(); i++) {
		#pragma omp parallel for
		for (register int j = 0; j < (int)left[0].size(); j++) {
			result[i][j] = pow(left[i][j], right);
		}
	}
	return result;
}

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> u) {
	std::vector<std::vector<double>> result;
	result.resize(u[0].size(), std::vector<double>(u.size(), 0.0));

	for (int i = 0; i < (int)u[0].size(); i++) {
		for (int j = 0; j < (int)u.size(); j++) {
			result[i][j] = u[j][i];
		}
	}
	return result;
}


std::vector<std::vector<double>> kron(const std::vector<std::vector<double>> A, const std::vector<std::vector<double>> B) {
	std::vector<std::vector<double>> result;
	result.resize(A.size() * B.size(), std::vector<double>(A[0].size() * B[0].size(), 0.0));
	int p = (int)B.size();
	int q = (int)B[0].size();
	int r = (int)result.size();
	int s = (int)result[0].size();
	
	#pragma omp parallel for
	for (register int i = 0; i < r; i++) {
		for (register int j = 0; j < s; j++) {
			result[i][j] = A[i / p][j / q] * B[i % p][j % q];
		}
	}

	return result;
}

//If A is a square matrix, then 'a' is its vectorized version.
//If A is a column vector, then 'a' is its corresponding matrix.
std::vector<std::vector<double>> vec(std::vector<std::vector<double>> u) {
	std::vector<std::vector<double>> result;

	if (((int)u[0].size()) == 1) {
		int p = (int)sqrt(u.size());
		result.resize(p, std::vector<double>(p, 0.0));
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < p; j++) {
				result[j][i] = u[i * p + j][0];
			}
		}
	}
	else {
		int p = (int)u.size();
		result.resize(p * p, std::vector<double>(1, 0.0));
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < p; j++) {
				result[i * p + j][0] = u[j][i];
			}
		}
	}
	return result;
}
