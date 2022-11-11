#include <iostream>
#include <string>
#include <fstream>
//#include <Windows.h>
#include <iomanip>
#include "vectorMath.h"


#define I 10000						// Number of iterations
#define K 20						// Node count
#define M 5							// Filter Taps

#define eyeMatrix eye(M)
#define eyeMatrixM2 eye(M * M)
#define onesMatrix ones(M, 1)
#define oneVec vec(eye(M))


std::default_random_engine generator((unsigned int)time(0));
std::normal_distribution<double> normal(0.0, 1.0);
std::uniform_real_distribution<double> uniformal(0.0, 1.0);



class Node {
public:
	int index = 0;
	double mu = 0.01;				// Step Size

	double var_n = 0;			// Guassian Noise
	double sqn = 0;				// sqrt(var_n);

	double rho = 5.0;				// Desired eigenvalue spread
	double lambda_min = 1.0;		// Generate positive eigenvalues in the interval
	double lambda_max = rho;		// [lambda_min; lambda_max];
	int delay = 0;
	int delayMean = 0;

	std::vector<std::vector<double>> eigenv;
	std::vector<std::vector<double>> R;						//Covariance matrix
	std::vector<std::vector<double>> sqR;
	std::vector<std::vector<double>> Lambda;				// Since R is diagonal; therefore U = I in this problem
	std::vector<std::vector<double>> lambda;				// a column of eigenvalues
	std::vector<std::vector<double>> msd_curve_theorical;
	std::vector<std::vector<double>> emse_curve_theorical;
	std::vector<std::vector<double>> F_bar;
	std::vector<std::vector<std::vector<double>>> Fa;
	std::vector<std::vector<double>> wo_bar;
	std::vector<std::vector<double>> u;
	std::vector<std::vector<double>> D;
	std::vector<std::vector<double>> Y;
	std::vector<std::vector<double>> E;
	std::vector<std::vector<std::vector<double>>> W;
	std::vector<std::vector<double>> nw;
	std::vector<std::vector<double>> msd_curve_simulation;
	std::vector<std::vector<double>> emse_curve_simulation;

public:
	Node() {
		wo_bar.resize(M, std::vector<double>(1, 0.0));
		eigenv.resize(M, std::vector<double>(1, 0.0));
		R.resize(M, std::vector<double>(M, 0.0));
		sqR.resize(M, std::vector<double>(M, 0.0));
		Lambda.resize(M, std::vector<double>(M, 0.0));
		lambda.resize(M, std::vector<double>(1, 0.0));
		F_bar.resize(M, std::vector<double>(M, 0.0));
		msd_curve_theorical.resize(1, std::vector<double>(I, 0.0));
		emse_curve_theorical.resize(1, std::vector<double>(I, 0.0));
		Fa.resize(I, std::vector<std::vector<double>>(M, std::vector<double>(M, 0.0)));
		u.resize(1, std::vector<double>(M, 0.0));
		D.resize(1, std::vector<double>(I, 0.0));
		Y.resize(1, std::vector<double>(I, 0.0));
		E.resize(1, std::vector<double>(I, 0.0));
		W.resize(I, std::vector<std::vector<double>>(1, std::vector<double>(M, 0.0)));
		nw.resize(I, std::vector<double>(1, 0.0));
		msd_curve_simulation.resize(1, std::vector<double>(I, 0.0));
		emse_curve_simulation.resize(1, std::vector<double>(I, 0.0));

		eigenv[0][0] = lambda_min;
		eigenv[M - 1][0] = lambda_max;
		for (int i = 1; i < M - 1; i++) {
			eigenv[i][0] = lambda_min + (lambda_max - lambda_min) * uniformal(generator);
		}
		R = diag(eigenv);		//Covariance matrix --- it is diagonal = Lambda
		sqR = SQRT(R);
		Lambda = R;
		lambda = eigenv;
		F_bar = (eyeMatrix - (2.0 * mu * Lambda) + ((mu * mu) * (Lambda * Lambda))) + ((mu * mu) * (lambda * transpose(lambda)));
		for (int i = 0; i < I; i++) {
			Fa[i] = eyeMatrix;
		}
		for (int i = 0; i < M; i++) {
			u[0][i] = normal(generator);
		}
	}

	void saveTheoricalMSDCurve() {
		std::ofstream temp;
		temp.open("results\\Node#" + std::to_string(index) + "_MSD_THE_Error_delay_" + std::to_string(delayMean) + ".txt");
		temp.precision(10);
		for (int i = 0; i < I; i++)
			temp << msd_curve_theorical[0][i] << "\n";

		temp.close();
	}
	void saveSimulationMSDCurve() {
		std::ofstream temp;
		temp.open("results\\Node#" + std::to_string(index) + "_MSD_EXP_Error_delay_" + std::to_string(delayMean) + ".txt");
		temp.precision(10);
		for (int i = 0; i < I; i++)
			temp << msd_curve_simulation[0][i] << "\n";

		temp.close();
	}
	void saveTheoricalEMSECurve() {
		std::ofstream temp;
		temp.open("results\\Node#" + std::to_string(index) + "_EMSE_THE_Error_delay_" + std::to_string(delayMean) +  ".txt");

		temp.precision(10);

		for (int i = 0; i < I; i++)
			temp << emse_curve_theorical[0][i] << "\n";

		temp.close();
	}
	void saveSimulationEMSECurve() {
		std::ofstream temp;
		temp.open("results\\Node#" + std::to_string(index) + "_EMSE_EXP_Error_delay_" + std::to_string(delayMean) + ".txt");
		temp.precision(10);
		for (int i = 0; i < I; i++)
			temp << emse_curve_simulation[0][i] << "\n";

		temp.close();
	}
	double calculateSteadyStateMSDValueTheoretical() {
		double avg = 0.0;
		for (int i = 0; i < 500; i++) {
			avg = avg + msd_curve_theorical[0][I - i - 1];
		}
		return avg / 500.0;
	}
	double calculateSteadyStateMSDValueSimulation() {
		double avg = 0.0;
		for (int i = 0; i < 500; i++) {
			avg = avg + msd_curve_simulation[0][I - i - 1];
		}
		return avg / 500.0;
	}
	double calculateSteadyStateEMSEValueTheoretical() {
		double avg = 0.0;
		for (int i = 0; i < 500; i++) {
			avg = avg + emse_curve_theorical[0][I - i - 1];
		}
		return avg / 500.0;
	}
	double calculateSteadyStateEMSEValueSimulation() {
		double avg = 0.0;
		for (int i = 0; i < 500; i++) {
			avg = avg + emse_curve_simulation[0][I - i - 1];
		}
		return avg / 500.0;
	}
};

void calculateTheoreticalMSDCurve(Node *node, std::vector<std::vector<double>> wo);
void calculateSimulationMSDCurve(Node *node, std::vector<std::vector<double>> wo, int L);
void calculateTheoreticalEMSECurve(Node *node, std::vector<std::vector<double>> wo);
void calculateSimulationEMSECurve(Node *node, std::vector<std::vector<double>> wo, int L);
void saveSteadyStateTheoreticalMSDValue(Node *node);
void saveSteadyStateSimulationMSDValue(Node *node);
void saveSteadyStateTheoreticalEMSEValue(Node *node);
void saveSteadyStateSimulationEMSEValue(Node *node);

class NonGaussianNode {
public:
	int index = 0;
	double mu = 0.01;				// Step Size

	double var_un = 0;
	double squn = 0;				// sqrt(var_un * 12);
	double alpha = 0;

	double rho = 5.0;				// Desired eigenvalue spread
	double lambda_min = 1.0;		// Generate positive eigenvalues in the interval
	double lambda_max = rho;		// [lambda_min; lambda_max];
	int delay = 0;
	int delayMean = 0;

	std::vector<std::vector<double>> eigenv;
	std::vector<std::vector<double>> R;						//Covariance matrix
	std::vector<std::vector<double>> sqR;
	std::vector<std::vector<double>> msd_curve_theorical;
	std::vector<std::vector<double>> emse_curve_theorical;

	std::vector<std::vector<double>> F_bar;
	std::vector<std::vector<std::vector<double>>> Fa;
	std::vector<std::vector<double>> wo_bar;

	double var_u = 0;
	std::vector<std::vector<double>> u;
	std::vector<std::vector<double>> D;
	std::vector<std::vector<double>> Y;
	std::vector<std::vector<double>> E;
	std::vector<std::vector<std::vector<double>>> W;
	std::vector<std::vector<double>> nw;
	std::vector<std::vector<double>> msd_curve_simulation;
	std::vector<std::vector<double>> emse_curve_simulation;

public:
	NonGaussianNode() {
		wo_bar.resize(M, std::vector<double>(1, 0.0));
		eigenv.resize(M, std::vector<double>(1, 0.0));
		R.resize(M, std::vector<double>(M, 0.0));
		sqR.resize(M, std::vector<double>(M, 0.0));
		F_bar.resize(M * M, std::vector<double>(M * M, 0.0));
		msd_curve_theorical.resize(1, std::vector<double>(I, 0.0));
		emse_curve_theorical.resize(1, std::vector<double>(I, 0.0));
		Fa.resize(I, std::vector<std::vector<double>>(M * M, std::vector<double>(M * M, 0.0)));
		u.resize(1, std::vector<double>(M, 0.0));
		D.resize(1, std::vector<double>(I, 0.0));
		Y.resize(1, std::vector<double>(I, 0.0));
		E.resize(1, std::vector<double>(I, 0.0));
		W.resize(I, std::vector<std::vector<double>>(1, std::vector<double>(M, 0.0)));
		nw.resize(I, std::vector<double>(1, 0.0));
		msd_curve_simulation.resize(1, std::vector<double>(I, 0.0));
		emse_curve_simulation.resize(1, std::vector<double>(I, 0.0));

		eigenv[0][0] = lambda_min;
		eigenv[M - 1][0] = lambda_max;
		for (int i = 1; i < M - 1; i++) {
			eigenv[i][0] = lambda_min + (lambda_max - lambda_min) * uniformal(generator);
		}
		R = diag(eigenv);		//Covariance matrix --- it is diagonal = Lambda
		sqR = SQRT(R);

		for (int i = 0; i < I; i++) {
			Fa[i] = eyeMatrixM2;
		}
	}

	void saveTheoricalMSDCurve() {
		std::ofstream temp;
		temp.precision(10);
		temp.open("results\\Node#" + std::to_string(index) + "_Non_Gaussian_MSD_THEORICAL_Error_" + std::to_string(delayMean) + ".txt");

		for (int i = 0; i < I; i++)
			temp << msd_curve_theorical[0][i] << "\n";

		temp.close();
	}
	void saveSimulationMSDCurve() {
		std::ofstream temp;
		temp.precision(10);
		temp.open("results\\Node#" + std::to_string(index) + "_Non_Gaussian_MSD_EXPERIMENTAL_Error_" + std::to_string(delayMean) + ".txt");

		for (int i = 0; i < I; i++)
			temp << msd_curve_simulation[0][i] << "\n";

		temp.close();
	}
	void saveTheoricalEMSECurve() {
		std::ofstream temp;
		temp.precision(10);
		temp.open("results\\Node#" + std::to_string(index) + "_Non_Gaussian_EMSE_THEORICAL_Error_" + std::to_string(delayMean) + ".txt");

		for (int i = 0; i < I; i++)
			temp << emse_curve_theorical[0][i] << "\n";

		temp.close();
	}
	void saveSimulationEMSECurve() {
		std::ofstream temp;
		temp.precision(10);
		temp.open("results\\Node#" + std::to_string(index) + "_Non_Gaussian_EMSE_EXPERIMENTAL_Error_" + std::to_string(delayMean) + ".txt");

		for (int i = 0; i < I; i++)
			temp << emse_curve_simulation[0][i] << "\n";

		temp.close();
	}
	double calculateSteadyStateMSDValueTheoretical() {
		double avg = 0;
		for (int i = 0; i < 500; i++) {
			avg = avg + msd_curve_theorical[0][I - i - 1];
		}
		return avg / 500.0;
	}
	double calculateSteadyStateMSDValueSimulation() {
		double avg = 0;
		for (int i = 0; i < 500; i++) {
			avg = avg + msd_curve_simulation[0][I - i - 1];
		}
		return avg / 500.0;
	}
	double calculateSteadyStateEMSEValueTheoretical() {
		double avg = 0;
		for (int i = 0; i < 500; i++) {
			avg = avg + emse_curve_theorical[0][I - i - 1];
		}
		return avg / 500.0;
	}
	double calculateSteadyStateEMSEValueSimulation() {
		double avg = 0;
		for (int i = 0; i < 500; i++) {
			avg = avg + emse_curve_simulation[0][I - i - 1];
		}
		return avg / 500.0;
	}
};

void calculateNonGaussianTheoreticalMSDCurve(NonGaussianNode *node, std::vector<std::vector<double>> wo);
void calculateNonGaussianSimulationMSDCurve(NonGaussianNode *node, std::vector<std::vector<double>> wo, int L);
void calculateNonGaussianTheoreticalEMSECurve(NonGaussianNode *node, std::vector<std::vector<double>> wo);
void calculateNonGaussianSimulationEMSECurve(NonGaussianNode *node, std::vector<std::vector<double>> wo, int L);
void saveSteadyStateNonGaussianTheoreticalMSDValue(NonGaussianNode *node);
void saveSteadyStateNonGaussianSimulationMSDValue(NonGaussianNode *node);
void saveSteadyStateNonGaussianTheoreticalEMSEValue(NonGaussianNode *node);
void saveSteadyStateNonGaussianSimulationEMSEValue(NonGaussianNode *node);

//void ShowConsoleCursor(bool showFlag);

int main() {
	//ShowConsoleCursor(false);
	clock_t start, end;
	start = clock();

	std::vector<std::vector<double>> wo;
	//Normally distributed pseudorandom numbers
	wo.resize(M, std::vector<double>(1, 0.0));
	for (int i = 0; i < (int)wo.size(); i++) {
		wo[i][0] =  uniformal(generator);
	}
	wo = wo / (l2_norm(wo));

	Node node[K];
	double varianceTemp[K] = {
		0.0592718,
		0.0453682,
		0.0305739,
		0.0807076,
		0.00811714,
		0.0162581,
		0.0617749,
		0.00357077,
		0.0863187,
		0.0626361,
		0.0973455,
		0.0231114,
		0.07729,
		0.0145409,
		0.0205027,
		0.0626162,
		0.0211048,
		0.0781336,
		0.0593297,
		0.0178738
	};

	for (int k = 0; k < K; k++) {
		node[k].index = k;
		node[k].var_n = varianceTemp[k];
		node[k].sqn = sqrt(node[k].var_n);
	}


	for (int k = 0; k < K; k++) {
		node[k].delayMean = 2;
		std::normal_distribution<double> normal_delay(node[k].delayMean, 1.0);
		node[k].delay = (int)abs(round(normal_delay(generator)));
	}
	calculateTheoreticalMSDCurve(node, wo);
	calculateSimulationMSDCurve(node, wo, 10);
	//calculateTheoreticalEMSECurve(node, wo);
	//calculateSimulationEMSECurve(node, wo, 10);
	for (int k = 0; k < K; k++) {
		node[k].saveTheoricalMSDCurve();
		node[k].saveSimulationMSDCurve();
		//node[k].saveTheoricalEMSECurve();
		//node[k].saveSimulationEMSECurve();
	}
	saveSteadyStateTheoreticalMSDValue(node);
	saveSteadyStateSimulationMSDValue(node);
	//saveSteadyStateTheoreticalEMSEValue(node);
	//saveSteadyStateSimulationEMSEValue(node);


	//---------------------------------------------------------------Non_Gaussian

	/*NonGaussianNode ngnode[K];
	double varianceTempNonGaussian[K] = {
		0.0364442,
		0.0190949,
		0.0256778,
		0.047748,
		0.0699879,
		0.0728704,
		0.0961621,
		0.0401297,
		0.0758908,
		0.0422516,
		0.0910473,
		0.00716829,
		0.075099,
		0.0880136,
		0.0661586,
		0.0886894,
		0.0372608,
		0.0468919,
		0.0939971,
		0.0408652
	};

	double alpha[K] = {
		0.344989,
		0.158934,
		0.290794,
		0.275534,
		0.264803,
		0.0159218,
		0.408175,
		0.159863,
		0.22325,
		0.277393,
		0.068108,
		0.0946405,
		0.0128532,
		0.311154,
		0.295435,
		0.147806,
		0.157287,
		0.185706,
		0.382797,
		0.106032
	};

	double var_un[K]{
		0.0437644,
		0.112294,
		0.0877024,
		0.374034,
		0.326542,
		0.467831,
		0.116653,
		0.405012,
		0.134277,
		0.324366,
		0.362124,
		0.107368,
		0.143197,
		0.398292,
		0.267924,
		0.281394,
		0.307436,
		0.0768452,
		0.371975,
		0.125418
	};

	for (int k = 0; k < K; k++) {
		ngnode[k].index = k;
		ngnode[k].var_un = varianceTempNonGaussian[k];
		ngnode[k].squn = sqrt(ngnode[k].var_un * 15);

		ngnode[k].alpha = alpha[k];
		ngnode[k].var_u = var_un[k];
	}

	for (int k = 0; k < K; k++) {
		ngnode[k].delayMean = 1;
		std::normal_distribution<double> normal_delay(ngnode[k].delayMean, 1.0);
		ngnode[k].delay = (int)abs(round(normal_delay(generator)));
	}
	calculateNonGaussianTheoreticalMSDCurve(ngnode, wo);
	calculateNonGaussianSimulationMSDCurve(ngnode, wo, 100);
	calculateNonGaussianTheoreticalEMSECurve(ngnode, wo);
	calculateNonGaussianSimulationEMSECurve(ngnode, wo, 100);
	for (int k = 0; k < K; k++) {
		ngnode[k].saveTheoricalMSDCurve();
		ngnode[k].saveSimulationMSDCurve();
		ngnode[k].saveTheoricalEMSECurve();
		ngnode[k].saveSimulationEMSECurve();
	}
	saveSteadyStateNonGaussianTheoreticalMSDValue(ngnode);
	saveSteadyStateNonGaussianSimulationMSDValue(ngnode);
	saveSteadyStateNonGaussianTheoreticalEMSEValue(ngnode);
	saveSteadyStateNonGaussianSimulationEMSEValue(ngnode);
*/

	end = clock();
	std::cout << "Time required for execution: " << (double)(end - start) / CLOCKS_PER_SEC << " seconds." << "\n\n";
	//system("pause");
	return 0;
}

void calculateTheoreticalMSDCurve(Node *node, std::vector<std::vector<double>> wo) {
	std::cout << "==================================\n";
	std::cout << "Theoretical MSD Curve Calculation\n";
	std::cout << "==================================\n";

	std::vector<std::vector<double>> temp = eye(M);
	double b = 0;
	double c = 0;

	for (int k = 0; k < K; k++) {
		for (int i = 0; i < I; i++) {
			node[k].Fa[i] = eyeMatrix;
		}
		node[k].index = k;
		node[k].msd_curve_theorical[0][0] = pow(l2_norm(wo), 2.0);		
		node[k].wo_bar = wo;
	}

	for (int nodeIndex = 0; nodeIndex < K; nodeIndex++) {		
		int localDelay = 0;
		for (int k = 0; k < nodeIndex; k++) {
			localDelay += node[k].delay;
		}
		for (int i = 0; i < localDelay ; i++) {
			node[nodeIndex].msd_curve_theorical[0][i] = node[nodeIndex].msd_curve_theorical[0][0];
		}

		for (int i = localDelay; i < I; i++) {
			std::cout << "progress:\t" << nodeIndex << "\t,\t" << i << "\r";
			int totalDelay = 0;
			for (int k = 0; k < K; k++) {
				totalDelay += node[k].delay;
			}

			b = 0;
			c = 0;

			if (i - 1 - totalDelay >= 0) {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrix;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b += (node[cl].mu * node[cl].mu * node[cl].var_n * transpose(node[cl].lambda) * temp * node[nodeIndex].Fa[i - 1 - totalDelay] * onesMatrix)[0][0];
				}

				temp = eyeMatrix;
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;
				}

				c = (transpose(node[nodeIndex].wo_bar) * (diag(node[nodeIndex].Fa[i - 1 - totalDelay] * (eyeMatrix - temp) * onesMatrix)) * node[nodeIndex].wo_bar)[0][0];

				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[i - 1 - totalDelay];
				node[nodeIndex].msd_curve_theorical[0][i] = node[nodeIndex].msd_curve_theorical[0][i - 1 - totalDelay] + b - c;
			}
			else {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrix;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b = b + (node[cl].mu * node[cl].mu * node[cl].var_n * transpose(node[cl].lambda) * temp * node[nodeIndex].Fa[0] * onesMatrix)[0][0];
				}

				temp = eyeMatrix;
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;
				}

				c = (transpose(node[nodeIndex].wo_bar) * (diag(node[nodeIndex].Fa[0] * (eyeMatrix - temp) * onesMatrix)) * node[nodeIndex].wo_bar)[0][0];

				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[0];
				node[nodeIndex].msd_curve_theorical[0][i] = node[nodeIndex].msd_curve_theorical[0][0] + b - c;
			}
		}
		std::cout << "\n";
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
void calculateSimulationMSDCurve(Node *node, std::vector<std::vector<double>> wo, int L) {
	std::cout << "==================================\n";
	std::cout << "Exprimental MSD Curve Calculation\n";
	std::cout << "==================================\n";
	for (int l = 0; l < L; l++) {
		std::cout << "progress:\t" << 1 + (int)(l  * 100.0f / L) << "%\r";
		//Beep(6000, 20);
		for (int k = 0; k < K; k++) {
			for (int i = 0; i < I; i++) {
				node[k].W[i] = zeros(1, M);
			}
			node[k].nw = zeros(I, 1);
		}

		for (int i = 0; i < I; i++) {
			for (int k = 0; k < K; k++) {				
				// Regressor
				for (int m = 0; m < M; m++) {
					node[k].u[0][m] = normal(generator);
				}
				node[k].u = node[k].u * node[k].sqR;

				// Reference
				node[k].D[0][i] = (node[k].u * wo)[0][0] + node[k].sqn * normal(generator);

				if (k > 0) {
					if (i - node[k].delay >= 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[k - 1].W[i - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k - 1].W[i - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}
				else {
					if (i - node[k].delay > 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[K - 1].W[i - 1 - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[K - 1].W[i - 1 - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}

				node[k].nw[i][0] = pow(l2_norm(wo - transpose(node[k].W[i])), 2.0);
			}
		}

		for (int k = 0; k < K; k++) {
			node[k].msd_curve_simulation = node[k].msd_curve_simulation + transpose(node[k].nw);
		}
	}

	for (int k = 0; k < K; k++) {
		node[k].msd_curve_simulation = node[k].msd_curve_simulation / L;
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
void calculateTheoreticalEMSECurve(Node *node, std::vector<std::vector<double>> wo) {
	std::cout << "==================================\n";
	std::cout << "Theoretical EMSE Curve Calculation\n";
	std::cout << "==================================\n";

	std::vector<std::vector<double>> temp = eye(M);
	std::vector<std::vector<double>> weightMatrix = eye(M);
	double b = 0;
	double c = 0;

	for (int k = 0; k < K; k++) {
		for (int i = 0; i < I; i++) {
			node[k].Fa[i] = eye(M);
		}
		node[k].index = k;		
		node[k].emse_curve_theorical[0][0] = (transpose(wo) * (node[k].Lambda / l2_norm(node[k].Lambda)) * wo)[0][0];		
		node[k].wo_bar = wo;
	}

	for (int nodeIndex = 0; nodeIndex < K; nodeIndex++) {
		int localDelay = 0;
		for (int k = 0; k < nodeIndex; k++) {
			localDelay += node[k].delay;
		}
		for (int i = 0; i < localDelay; i++) {
			node[nodeIndex].emse_curve_theorical[0][i] = node[nodeIndex].emse_curve_theorical[0][0];
		}

		weightMatrix = node[nodeIndex].lambda / l2_norm(node[nodeIndex].Lambda);
		for (int i = localDelay; i < I; i++) {
			std::cout << "progress:\t" << nodeIndex << "\t,\t" << i << "\r";

			int totalDelay = 0;
			for (int k = 0; k < K; k++) {
				totalDelay += node[k].delay;
			}
			
			b = 0;
			c = 0;						

			if (i - 1 - totalDelay >= 0) {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrix;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b += (node[cl].mu * node[cl].mu * node[cl].var_n * transpose(node[cl].lambda) * temp * node[nodeIndex].Fa[i - 1 - totalDelay] * weightMatrix)[0][0];
				}

				temp = eyeMatrix;				
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;					
				}


				c = (transpose(node[nodeIndex].wo_bar) * (diag(node[nodeIndex].Fa[i - 1 - totalDelay] * (eyeMatrix - temp) * weightMatrix)) * node[nodeIndex].wo_bar)[0][0];

				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[i - 1 - totalDelay];				
				node[nodeIndex].emse_curve_theorical[0][i] = node[nodeIndex].emse_curve_theorical[0][i - 1 - totalDelay] + b - c;
				
			}
			else {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrix;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b = b + (node[cl].mu * node[cl].mu * node[cl].var_n * transpose(node[cl].lambda) * temp * node[nodeIndex].Fa[0] * weightMatrix)[0][0];
				}

				temp = eyeMatrix;
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;
				}

				c = (transpose(node[nodeIndex].wo_bar) * (diag(node[nodeIndex].Fa[0] * (eyeMatrix - temp) * weightMatrix)) * node[nodeIndex].wo_bar)[0][0];
				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[0];
				node[nodeIndex].emse_curve_theorical[0][i] = node[nodeIndex].emse_curve_theorical[0][0] + b - c;
			}
		}
		std::cout << "\n";
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
void calculateSimulationEMSECurve(Node *node, std::vector<std::vector<double>> wo, int L) {
	std::cout << "==================================\n";
	std::cout << "Exprimental EMSE Curve Calculation\n";
	std::cout << "==================================\n";
	for (int l = 0; l < L; l++) {
		std::cout << "progress:\t" << 1 + (int)(l * 100.0f / L) << "%\r";
		//Beep(6000, 20);
		for (int k = 0; k < K; k++) {
			for (int i = 0; i < I; i++) {
				node[k].W[i] = zeros(1, M);
			}
			node[k].nw = zeros(I, 1);
		}

		for (int i = 0; i < I; i++) {
			for (int k = 0; k < K; k++) {
				// Regressor
				for (int m = 0; m < M; m++) {
					node[k].u[0][m] = normal(generator);
				}
				node[k].u = node[k].u * node[k].sqR;

				// Reference
				node[k].D[0][i] = (node[k].u * wo)[0][0] + node[k].sqn * normal(generator);

				if (k > 0) {
					if (i - node[k].delay >= 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[k - 1].W[i - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k - 1].W[i - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}
				else {
					if (i - node[k].delay > 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[K - 1].W[i - 1 - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[K - 1].W[i - 1 - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}

				node[k].nw[i][0] = (transpose(wo - transpose(node[k].W[i])) * (node[k].Lambda / l2_norm(node[k].Lambda)) * (wo - transpose(node[k].W[i])))[0][0];				
			}
		}

		for (int k = 0; k < K; k++) {
			node[k].emse_curve_simulation = node[k].emse_curve_simulation + transpose(node[k].nw);	
		}
	}

	for (int k = 0; k < K; k++) {
		node[k].emse_curve_simulation = node[k].emse_curve_simulation / L;
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
void saveSteadyStateTheoreticalMSDValue(Node *node) {
	std::ofstream temp;
	temp.open("results\\SS_MSD_THE_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp <<  node[k].calculateSteadyStateMSDValueTheoretical() << "\n";

	temp.close();
}
void saveSteadyStateSimulationMSDValue(Node *node) {
	std::ofstream temp;
	temp.open("results\\SS_MSD_EXP_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp << node[k].calculateSteadyStateMSDValueSimulation() << "\n";

	temp.close();
}
void saveSteadyStateTheoreticalEMSEValue(Node *node) {
	std::ofstream temp;
	temp.open("results\\SS_EMSE_THE_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp << node[k].calculateSteadyStateEMSEValueTheoretical() << "\n";

	temp.close();
}
void saveSteadyStateSimulationEMSEValue(Node *node) {
	std::ofstream temp;
	temp.open("results\\SS_EMSE_EXP_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp << node[k].calculateSteadyStateEMSEValueSimulation() << "\n";

	temp.close();
}

void calculateNonGaussianTheoreticalMSDCurve(NonGaussianNode *node, std::vector<std::vector<double>> wo) {
	std::cout << "==============================================\n";
	std::cout << "Theoretical MSD Curve Calculation Non-Gaussian\n";
	std::cout << "==============================================\n";
	
	std::vector<std::vector<double>> temp = eye(M * M);
	double b = 0;
	double c = 0;

	std::vector<std::vector<double>> temp0;
	temp0.resize(1, std::vector<double>(M, 0.0));
	std::vector<std::vector<double>> A;
	A.resize(M * M, std::vector<double>(M * M, 0.0));
	std::vector<std::vector<double>> B;
	B.resize(M * M, std::vector<double>(M * M, 0.0));

	std::vector<std::vector<double>> temp1;

	for (int k = 0; k < K; k++) {
		node[k].u = zeros(1, M);
		for (int i = 0; i < I; i++) {

			// Regressor
			for (register int j = 0; j < M; j++) {
				temp0[0][j] = (uniformal(generator) - 0.5);				
			}
			node[k].u = node[k].alpha * node[k].u + (1 - node[k].alpha) * sqrt(15) * temp0 * node[k].sqR;			

			temp1 = transpose(node[k].u) * node[k].u;
			B = B + kron(temp1, temp1);
			A = A + kron(eyeMatrix, temp1) + kron(temp1, eyeMatrix);			
		}
	}

	B = B / (I * K);
	A = A / (I * K);
 
	for (int k = 0; k < K; k++) {
		for (int i = 0; i < I; i++) {
			node[k].Fa[i] = eyeMatrixM2;
		}
		node[k].F_bar = eyeMatrixM2 - node[k].mu * A + (node[k].mu * node[k].mu) * B;
		node[k].index = k;
		node[k].msd_curve_theorical[0][0] = pow(l2_norm(wo), 2.0);
		node[k].wo_bar = wo;
	}


	for (int nodeIndex = 0; nodeIndex < K; nodeIndex++) {
		int localDelay = 0;
		for (int k = 0; k < nodeIndex; k++) {
			localDelay += node[k].delay;
		}
		for (int i = 0; i < localDelay; i++) {
			node[nodeIndex].msd_curve_theorical[0][i] = node[nodeIndex].msd_curve_theorical[0][0];
		}

		for (int i = localDelay; i < I; i++) {
			std::cout << "progress:\t" << nodeIndex << "\t,\t" << i << "\r";
			int totalDelay = 0;
			for (int k = 0; k < K; k++) {
				totalDelay += node[k].delay;
			}

			b = 0;
			c = 0;

			if (i - 1 - totalDelay >= 0) {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrixM2;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b += (node[cl].mu * node[cl].mu * node[cl].var_un * transpose(vec(node[cl].R)) * temp * node[nodeIndex].Fa[i - 1 - totalDelay] * oneVec)[0][0];
				}

				temp = eyeMatrixM2;
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;
				}

				c = (transpose(node[nodeIndex].wo_bar) * (vec(node[nodeIndex].Fa[i - 1 - totalDelay] * (eyeMatrixM2 - temp) * oneVec)) * node[nodeIndex].wo_bar)[0][0];

				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[i - 1 - totalDelay];
				node[nodeIndex].msd_curve_theorical[0][i] = node[nodeIndex].msd_curve_theorical[0][i - 1 - totalDelay] + b - c;
			}
			else {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrixM2;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b = b + (node[cl].mu * node[cl].mu * node[cl].var_un * transpose(vec(node[cl].R)) * temp * node[nodeIndex].Fa[0] * oneVec)[0][0];
				}

				temp = eyeMatrixM2;
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;
				}

				c = (transpose(node[nodeIndex].wo_bar) * (vec(node[nodeIndex].Fa[0] * (eyeMatrixM2 - temp) * oneVec)) * node[nodeIndex].wo_bar)[0][0];

				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[0];
				node[nodeIndex].msd_curve_theorical[0][i] = node[nodeIndex].msd_curve_theorical[0][0] + b - c;
			}
		}
		std::cout << "\n";
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
void calculateNonGaussianSimulationMSDCurve(NonGaussianNode *node, std::vector<std::vector<double>> wo, int L) {
	std::cout << "==============================================\n";
	std::cout << "Exprimental MSD Curve Calculation Non_Gaussian\n";
	std::cout << "==============================================\n";

	std::vector<std::vector<double>> temp0;
	temp0.resize(1, std::vector<double>(M, 0.0));


	for (int l = 0; l < L; l++) {
		std::cout << "progress:\t" << 1 + (int)(l  * 100.0f / L) << "%\r";
		//Beep(6000, 20);
		for (int k = 0; k < K; k++) {
			for (int i = 0; i < I; i++) {
				node[k].W[i] = zeros(1, M);
			}
			node[k].nw = zeros(I, 1);
			node[k].u = zeros(1, M);
		}

		for (int i = 0; i < I; i++) {
			for (int k = 0; k < K; k++) {

				// Regressor
				for (register int j = 0; j < M; j++) {
					temp0[0][j] = (uniformal(generator) - 0.5);					
				}
				node[k].u = node[k].alpha * node[k].u + (1 - node[k].alpha) * sqrt(15) * temp0 * node[k].sqR;				

				// Reference
				node[k].D[0][i] = (node[k].u * wo)[0][0] + node[k].squn * (uniformal(generator) - 0.5);

				if (k > 0) {
					if (i - node[k].delay >= 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[k - 1].W[i - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k - 1].W[i - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}
				else {
					if (i - node[k].delay > 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[K - 1].W[i - 1 - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[K - 1].W[i - 1 - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}

				node[k].nw[i][0] = pow(l2_norm(wo - transpose(node[k].W[i])), 2.0);
			}
		}

		for (int k = 0; k < K; k++) {
			node[k].msd_curve_simulation = node[k].msd_curve_simulation + transpose(node[k].nw);
		}
	}

	for (int k = 0; k < K; k++) {
		node[k].msd_curve_simulation = node[k].msd_curve_simulation / L;
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
void calculateNonGaussianTheoreticalEMSECurve(NonGaussianNode *node, std::vector<std::vector<double>> wo) {
	std::cout << "===============================================\n";
	std::cout << "Theoretical EMSE Curve Calculation Non-Gaussian\n";
	std::cout << "===============================================\n";

	std::vector<std::vector<double>> temp = eye(M * M);
	double b = 0;
	double c = 0;

	std::vector<std::vector<double>> temp0;
	temp0.resize(1, std::vector<double>(M, 0.0));
	std::vector<std::vector<double>> A;
	A.resize(M * M, std::vector<double>(M * M, 0.0));
	std::vector<std::vector<double>> B;
	B.resize(M * M, std::vector<double>(M * M, 0.0));

	std::vector<std::vector<double>> temp1;

	for (int k = 0; k < K; k++) {
		node[k].u = zeros(1, M);
		for (int i = 0; i < I; i++) {

			// Regressor
			for (register int j = 0; j < M; j++) {
				temp0[0][j] = (uniformal(generator) - 0.5);
			}
			node[k].u = node[k].alpha * node[k].u + (1 - node[k].alpha) * sqrt(15) * temp0 * node[k].sqR;

			temp1 = transpose(node[k].u) * node[k].u;
			B = B + kron(temp1, temp1);
			A = A + kron(eyeMatrix, temp1) + kron(temp1, eyeMatrix);
		}
	}

	B = B / (I * K);
	A = A / (I * K);

	for (int k = 0; k < K; k++) {
		for (int i = 0; i < I; i++) {
			node[k].Fa[i] = eyeMatrixM2;
		}
		node[k].F_bar = eyeMatrixM2 - node[k].mu * A + (node[k].mu * node[k].mu) * B;
		node[k].index = k;
		node[k].emse_curve_theorical[0][0] = (transpose(wo) * (node[k].R / l2_norm(node[k].R)) * wo)[0][0];
		node[k].wo_bar = wo;
	}


	for (int nodeIndex = 0; nodeIndex < K; nodeIndex++) {
		int localDelay = 0;
		for (int k = 0; k < nodeIndex; k++) {
			localDelay += node[k].delay;
		}
		for (int i = 0; i < localDelay; i++) {
			node[nodeIndex].emse_curve_theorical[0][i] = node[nodeIndex].emse_curve_theorical[0][0];
		}

		for (int i = localDelay; i < I; i++) {
			std::cout << "progress:\t" << nodeIndex << "\t,\t" << i << "\r";
			int totalDelay = 0;
			for (int k = 0; k < K; k++) {
				totalDelay += node[k].delay;
			}

			b = 0;
			c = 0;

			std::vector<std::vector<double>> weightMatrix = vec(node[nodeIndex].R / l2_norm(node[nodeIndex].R));

			if (i - 1 - totalDelay >= 0) {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrixM2;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b += (node[cl].mu * node[cl].mu * node[cl].var_un * transpose(vec(node[cl].R)) * temp * node[nodeIndex].Fa[i - 1 - totalDelay] * weightMatrix)[0][0];
				}

				temp = eyeMatrixM2;
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;
				}

				c = (transpose(node[nodeIndex].wo_bar) * (vec(node[nodeIndex].Fa[i - 1 - totalDelay] * (eyeMatrixM2 - temp) * weightMatrix)) * node[nodeIndex].wo_bar)[0][0];

				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[i - 1 - totalDelay];
				node[nodeIndex].emse_curve_theorical[0][i] = node[nodeIndex].emse_curve_theorical[0][i - 1 - totalDelay] + b - c;
			}
			else {
				for (int l = 0; l < K; l++) {
					int cl = (l + (nodeIndex + 1)) % K;
					temp = eyeMatrixM2;
					for (int k = K - 1; k > l; k--) {
						int ck = (k + (nodeIndex + 1)) % K;
						temp = node[ck].F_bar * temp;
					}
					b = b + (node[cl].mu * node[cl].mu * node[cl].var_un * transpose(vec(node[cl].R)) * temp * node[nodeIndex].Fa[0] * weightMatrix)[0][0];
				}

				temp = eyeMatrixM2;
				for (int k = K - 1; k >= 0; k--) {
					int ck = (k + (nodeIndex + 1)) % K;
					temp = node[ck].F_bar * temp;
				}

				c = (transpose(node[nodeIndex].wo_bar) * (vec(node[nodeIndex].Fa[0] * (eyeMatrixM2 - temp) * weightMatrix)) * node[nodeIndex].wo_bar)[0][0];

				node[nodeIndex].Fa[i] = temp * node[nodeIndex].Fa[0];
				node[nodeIndex].emse_curve_theorical[0][i] = node[nodeIndex].emse_curve_theorical[0][0] + b - c;
			}
		}
		std::cout << "\n";
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
void calculateNonGaussianSimulationEMSECurve(NonGaussianNode *node, std::vector<std::vector<double>> wo, int L) {
	std::cout << "===============================================\n";
	std::cout << "Exprimental EMSE Curve Calculation Non-Gaussain\n";
	std::cout << "===============================================\n";
	std::vector<std::vector<double>> temp;
	temp.resize(1, std::vector<double>(M, 0.0));
	
	for (int l = 0; l < L; l++) {
		std::cout << "progress:\t" << 1 + (int)(l  * 100.0f / L) << "%\r";
		//Beep(6000, 20);
		for (int k = 0; k < K; k++) {
			for (int i = 0; i < I; i++) {
				node[k].W[i] = zeros(1, M);
			}
			node[k].nw = zeros(I, 1);
			node[k].u = zeros(1, M);
		}

		for (int i = 0; i < I; i++) {
			for (int k = 0; k < K; k++) {
				// Regressor
				for (int j = 0; j < M; j++) {
					temp[0][j] = (uniformal(generator) - 0.5);
				}
				node[k].u = node[k].alpha * node[k].u + (1 - node[k].alpha) * sqrt(15) * temp * node[k].sqR;

				// Reference
				node[k].D[0][i] = (node[k].u * wo)[0][0] + node[k].squn * (uniformal(generator) - 0.5);

				if (k > 0) {
					if (i - node[k].delay >= 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[k - 1].W[i - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k - 1].W[i - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}
				else {
					if (i - node[k].delay > 0) {
						node[k].Y[0][i] = (node[k].u * transpose(node[K - 1].W[i - 1 - node[k].delay]))[0][0];
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[K - 1].W[i - 1 - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
					}
					else {
						node[k].Y[0][i] = 0.0;
						node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
						node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
					}
				}

				node[k].nw[i][0] = (transpose(wo - transpose(node[k].W[i])) * (node[k].R / l2_norm(node[k].R)) * (wo - transpose(node[k].W[i])))[0][0];
			}
		}

		for (int k = 0; k < K; k++) {
			node[k].emse_curve_simulation = node[k].emse_curve_simulation + transpose(node[k].nw);
		}
	}

	for (int k = 0; k < K; k++) {
		node[k].emse_curve_simulation = node[k].emse_curve_simulation / L;
	}
	std::cout << "\n";
	std::cout << "Done\n\n";
}
void saveSteadyStateNonGaussianTheoreticalMSDValue(NonGaussianNode *node) {
	std::ofstream temp;
	temp.open("results\\SS_NG_MSD_THE_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp << node[k].calculateSteadyStateMSDValueTheoretical() << "\n";

	temp.close();
}
void saveSteadyStateNonGaussianSimulationMSDValue(NonGaussianNode *node) {
	std::ofstream temp;
	temp.open("results\\SS_NG_MSD_EXP_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp << node[k].calculateSteadyStateMSDValueSimulation() << "\n";

	temp.close();
}
void saveSteadyStateNonGaussianTheoreticalEMSEValue(NonGaussianNode *node) {
	std::ofstream temp;
	temp.open("results\\SS_NG_EMSE_THE_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp << node[k].calculateSteadyStateEMSEValueTheoretical() << "\n";

	temp.close();
}
void saveSteadyStateNonGaussianSimulationEMSEValue(NonGaussianNode *node) {
	std::ofstream temp;
	temp.open("results\\SS_NG_EMSE_EXP_Error_delay_" + std::to_string(node[0].delayMean) + ".txt");
	temp.precision(10);

	for (int k = 0; k < K; k++)
		temp << node[k].calculateSteadyStateEMSEValueSimulation() << "\n";

	temp.close();
}


//void ShowConsoleCursor(bool showFlag) {
//	HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
//
//	CONSOLE_CURSOR_INFO     cursorInfo;
//
//	GetConsoleCursorInfo(out, &cursorInfo);
//	cursorInfo.bVisible = showFlag; // set the cursor visibility
//	SetConsoleCursorInfo(out, &cursorInfo);
//}



void calculateSimulationMSDCurve_NodeIdle(Node *node, std::vector<std::vector<double>> wo, int L) {
	std::cout << "==================================\n";
	std::cout << "Exprimental MSD Curve Calculation\n";
	std::cout << "==================================\n";
	for (int l = 0; l < L; l++) {
		std::cout << "progress:\t" << 1 + (int)(l  * 100.0f / L) << "%\r";
		//Beep(6000, 20);
		for (int k = 0; k < K; k++) {
			for (int i = 0; i < I; i++) {
				node[k].W[i] = zeros(1, M);
			}
			node[k].nw = zeros(I, 1);
		}

		for (int i = 0; i < I; i++) {
			for (int k = 0; k < K; k++) {
				// Regressor
				for (int m = 0; m < M; m++) {
					node[k].u[0][m] = normal(generator);
				}
				node[k].u = node[k].u * node[k].sqR;

				// Reference
				node[k].D[0][i] = (node[k].u * wo)[0][0] + node[k].sqn * normal(generator);
				if (i % 8 == 0) {
					if (k > 0) {
						if (i - node[k].delay >= 0) {
							node[k].Y[0][i] = (node[k].u * transpose(node[k - 1].W[i - node[k].delay]))[0][0];
							node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
							node[k].W[i] = node[k - 1].W[i - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
						}
						else {
							node[k].Y[0][i] = 0.0;
							node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
							node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
						}
					}
					else {
						if (i - node[k].delay > 0) {
							node[k].Y[0][i] = (node[k].u * transpose(node[K - 1].W[i - 1 - node[k].delay]))[0][0];
							node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
							node[k].W[i] = node[K - 1].W[i - 1 - node[k].delay] + node[k].mu * node[k].E[0][i] * node[k].u;
						}
						else {
							node[k].Y[0][i] = 0.0;
							node[k].E[0][i] = node[k].D[0][i] - node[k].Y[0][i];
							node[k].W[i] = node[k].mu * node[k].E[0][i] * node[k].u;
						}
					}
				}
				else {
					if (i > 0) {
						node[k].W[i] = node[k].W[i - 1];
					}
				}

				node[k].nw[i][0] = pow(l2_norm(wo - transpose(node[k].W[i])), 2.0);
			}
		}

		for (int k = 0; k < K; k++) {
			node[k].msd_curve_simulation = node[k].msd_curve_simulation + transpose(node[k].nw);
		}
	}

	for (int k = 0; k < K; k++) {
		node[k].msd_curve_simulation = node[k].msd_curve_simulation / L;
	}

	std::cout << "\n";
	std::cout << "Done\n\n";
}
