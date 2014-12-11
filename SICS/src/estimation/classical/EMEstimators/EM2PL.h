/*
 * EM2PL.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

#ifndef EM2PL_H_
#define EM2PL_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/TwoPLModel.h>
class EM2PL: public EMEstimator {
public:
	EM2PL() {
	}
	virtual ~EM2PL() {
	}

	virtual void transform(Model* m) {
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			double qa = (*m->getParameterModel()->getParameterSet()[a])(0, i);
			double qb = (*m->getParameterModel()->getParameterSet()[d])(0, i);
		}
	}

	virtual void untransform(Model* m) {
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			double qa = (*m->getParameterModel()->getParameterSet()[a])(0, i);
			double qb = (*m->getParameterModel()->getParameterSet()[d])(0, i);
		}
	}

	virtual void setInitialValues(map<Parameter, Matrix<double>*> parameterSet,
			Model* m) {
		m->getParameterModel()->setParameterSet(parameterSet);
	}

	virtual void setInitialValues(int method, Model* m) {
		//TODO MOVE ALGORITHMS TO ANOTHER FILE
		/*TODO
		 * Possible methods
		 * ANDRADE
		 * OSPINA
		 * RANDOM
		 *
		 * The default method is OSPINA
		 */
		if (!method == Constant::RANDOM) {
			std::srand(std::time(0)); // use current time as seed for random generator
			int items = m->getParameterModel()->getParameterSet()[a]->nC();
			for (int i = 0; i < items; i++) {
				(*m->getParameterModel()->getParameterSet()[a])(0, i) =
						randomd() * 2;
				//fill b
				(*m->getParameterModel()->getParameterSet()[d])(0, i) =
						randomd() * 4 - 2;
			}
		}

		if (!method == Constant::ANDRADE) {
			//Andrade method
			int items = m->getParameterModel()->getParameterSet()[a]->nC();
			//sums of the patterns
			int totalscores = 0;
			int *itemscores = new int[items];
			memset(itemscores, 0, sizeof(int) * items);
			double *covariances = new double[items];
			memset(covariances, 0, sizeof(double) * items);
			double variance = 0;
			PatternMatrix* data =
					dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
			double Ni = (double) data->countIndividuals();
			for (data->resetIterator(); !data->checkEnd(); data->iterate()) {
				double df = (double) data->getCurrentFrequency();
				double bs = (double) data->getCurrentBitSet().count();
				for (int i = 0; i < items; i++) {
					if (data->getCurrentBitSet()[i]) {
						itemscores[i] += df;
					}
				}
				totalscores += bs * df;
			}
			cout << "Score total : " << totalscores << endl;
			//calculate variances and covariances
			for (data->resetIterator(); !data->checkEnd(); data->iterate()) {
				double df = (double) data->getCurrentFrequency();
				double bs = (double) data->getCurrentBitSet().count();
				for (int i = 0; i < items; i++) {
					if (data->getCurrentBitSet()[i]) {
						covariances[i] += ((1 - itemscores[i] / Ni)
								* (1 - bs / items)) * df;
					}
				}
				variance += ((bs - ((double) totalscores / Ni))
						* (bs - ((double) totalscores / Ni))) * df;
			}
			variance /= Ni;
			for (int i = 0; i < items; i++) {
				covariances[i] /= (Ni - 1);
				cout << "cov : " << i << " " << covariances[i] << endl;
			}
			//Now calculate the standard deviations for the sums and the items
			long double*stddevs = new long double[items];
			memset(stddevs, 0, sizeof(long double) * items);
			long double*pearson = new long double[items];
			memset(pearson, 0, sizeof(long double) * items);
			long double*pis = new long double[items];
			memset(pis, 0, sizeof(long double) * items);
			for (int i = 0; i < items; i++) {
				double avg = totalscores / Ni;
				stddevs[i] = stdDev_bin(itemscores[i], Ni, avg);
				pis[i] = itemscores[i] / Ni;
				pearson[i] = (covariances[i]
						/ (stddevs[i] * std::sqrt(variance)));
				//fill a sqrt(pCoef * pCoef / (1.0 - pCoef * pCoef));
				(*m->getParameterModel()->getParameterSet()[a])(0, i) =
						std::sqrt(
								(pearson[i] * pearson[i])
										/ (1 / pearson[i] * pearson[i]));
				//fill b
				(*m->getParameterModel()->getParameterSet()[d])(0, i) =
						normalInverse(pis[i]);
			}
		}

		if (!method == Constant::CRISTIAN) {
			//fill a
			int items = m->getParameterModel()->getParameterSet()[a]->nC();
			for (int i = 0; i < items; i++) {
				(*m->getParameterModel()->getParameterSet()[a])(0, i) = 0.851;
			}
			//fill b
			(*m->getParameterModel()->getParameterSet()[d])(0, 0) =
					-0.566568124471702 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 1) =
					-1.14511732607529 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 2) =
					0.628823949705996 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 3) =
					1.98939265101279 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 4) =
					1.11965696346624 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 5) =
					-0.478666609627076 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 6) =
					-1.178623542548 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 7) =
					-1.06978229722843 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 8) =
					-1.4047686979683 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 9) =
					0.683068898124222 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 10) =
					-0.0532617736425396 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 11) =
					-1.10535471401831 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 12) =
					-0.193572680371671 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 13) =
					-0.0445766051573807 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 14) =
					-0.661327353196143 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 15) =
					-1.26249866491924 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 16) =
					0.555240135200819 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 17) =
					-1.36389252732882 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 18) =
					-1.46232353425884 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 19) =
					-0.319984341017887 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 20) =
					-1.97016609909551 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 21) =
					-2.99684557896768 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 22) =
					0.146555202506193 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 23) =
					0.23426909224306 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 24) =
					0.970850927690999 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 25) =
					1.90188814491673 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 26) =
					1.21769820870232 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 27) =
					-0.436750297269261 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 28) =
					0.743178112615275 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 29) =
					-0.445003669917367 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 30) =
					-1.85115355373292 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 31) =
					-1.5272787583475 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 32) =
					-0.588667339897699 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 33) =
					0.987951431011379 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 34) =
					-2.91798964751742 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 35) =
					-0.559412036351752 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 36) =
					-0.0746886675230089 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 37) =
					-0.860005617142715 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 38) =
					-0.751107888939827 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 39) =
					-1.96933326133621 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 40) =
					1.39850515263464 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 41) =
					0.934922614426893 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 42) =
					2.25800750013661 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 43) =
					0.675815082921858 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 44) =
					-0.0451555954167578 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 45) =
					-1.01338578098843 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 46) =
					-0.408494649436581 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 47) =
					1.04856058360661 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 48) =
					-0.216238170167062 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 49) =
					-0.311803829141431 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 50) =
					0.685488331328067 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 51) =
					-2.38365746189177 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 52) =
					-0.490501944430744 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 53) =
					2.12895158575715 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 54) =
					0.856284068445749 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 55) =
					-0.545115885651256 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 56) =
					0.710331472101276 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 57) =
					1.620011667106 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 58) =
					-1.619271389274 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 59) =
					-0.858144564233935 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 60) =
					1.54968254624868 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 61) =
					-1.79267801974796 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 62) =
					1.14642671851194 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 63) =
					-0.0231556281292887 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 64) =
					-2.39750585134654 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 65) =
					-1.83443982471499 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 66) =
					1.50001049162294 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 67) =
					-1.10146147008743 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 68) =
					0.8074863386241 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 69) =
					0.393215270239003 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 70) =
					0.397327228486701 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 71) =
					-0.442055389800962 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 72) =
					0.454443066816656 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 73) =
					1.12682418983781 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 74) =
					-0.837091198602913 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 75) =
					1.36458176012761 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 76) =
					-0.116983973327465 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 77) =
					1.56348620206992 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 78) =
					-0.419671630600554 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 79) =
					0.450902411677853 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 80) =
					-2.69627444177271 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 81) =
					-0.142495092212555 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 82) =
					0.153516452580127 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 83) =
					-0.371501010161694 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 84) =
					0.269209250206473 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 85) =
					-1.33505315586072 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 86) =
					-0.792733824801248 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 87) =
					-0.88487364983854 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 88) =
					-1.19116957722689 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 89) =
					-0.646865603719435 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 90) =
					-0.448542564109881 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 91) =
					0.595846146713524 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 92) =
					3.0878294021287 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 93) =
					-0.22554273716775 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 94) =
					0.151195884024596 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 95) =
					0.212749870413368 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 96) =
					0.222053103595936 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 97) =
					-0.300124140888755 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 98) =
					-0.107131397076053 * 0.851;
			(*m->getParameterModel()->getParameterSet()[d])(0, 99) =
					-0.700020844834876 * 0.851;
		}

		if (!method == Constant::CRISTIAN2) {
			//fill a
			int items = m->getParameterModel()->getParameterSet()[a]->nC();
			for (int i = 0; i < items; i++) {
				(*m->getParameterModel()->getParameterSet()[a])(0, i) = 0.851;
			}
			//fill b
			(*m->getParameterModel()->getParameterSet()[d])(0, 0) =
					-0.566568124471702;
			(*m->getParameterModel()->getParameterSet()[d])(0, 1) =
					-1.14511732607529;
			(*m->getParameterModel()->getParameterSet()[d])(0, 2) =
					0.628823949705996;
			(*m->getParameterModel()->getParameterSet()[d])(0, 3) =
					1.98939265101279;
			(*m->getParameterModel()->getParameterSet()[d])(0, 4) =
					1.11965696346624;
			(*m->getParameterModel()->getParameterSet()[d])(0, 5) =
					-0.478666609627076;
			(*m->getParameterModel()->getParameterSet()[d])(0, 6) =
					-1.178623542548;
			(*m->getParameterModel()->getParameterSet()[d])(0, 7) =
					-1.06978229722843;
			(*m->getParameterModel()->getParameterSet()[d])(0, 8) =
					-1.4047686979683;
			(*m->getParameterModel()->getParameterSet()[d])(0, 9) =
					0.683068898124222;
			(*m->getParameterModel()->getParameterSet()[d])(0, 10) =
					-0.0532617736425396;
			(*m->getParameterModel()->getParameterSet()[d])(0, 11) =
					-1.10535471401831;
			(*m->getParameterModel()->getParameterSet()[d])(0, 12) =
					-0.193572680371671;
			(*m->getParameterModel()->getParameterSet()[d])(0, 13) =
					-0.0445766051573807;
			(*m->getParameterModel()->getParameterSet()[d])(0, 14) =
					-0.661327353196143;
			(*m->getParameterModel()->getParameterSet()[d])(0, 15) =
					-1.26249866491924;
			(*m->getParameterModel()->getParameterSet()[d])(0, 16) =
					0.555240135200819;
			(*m->getParameterModel()->getParameterSet()[d])(0, 17) =
					-1.36389252732882;
			(*m->getParameterModel()->getParameterSet()[d])(0, 18) =
					-1.46232353425884;
			(*m->getParameterModel()->getParameterSet()[d])(0, 19) =
					-0.319984341017887;
			(*m->getParameterModel()->getParameterSet()[d])(0, 20) =
					-1.97016609909551;
			(*m->getParameterModel()->getParameterSet()[d])(0, 21) =
					-2.99684557896768;
			(*m->getParameterModel()->getParameterSet()[d])(0, 22) =
					0.146555202506193;
			(*m->getParameterModel()->getParameterSet()[d])(0, 23) =
					0.23426909224306;
			(*m->getParameterModel()->getParameterSet()[d])(0, 24) =
					0.970850927690999;
			(*m->getParameterModel()->getParameterSet()[d])(0, 25) =
					1.90188814491673;
			(*m->getParameterModel()->getParameterSet()[d])(0, 26) =
					1.21769820870232;
			(*m->getParameterModel()->getParameterSet()[d])(0, 27) =
					-0.436750297269261;
			(*m->getParameterModel()->getParameterSet()[d])(0, 28) =
					0.743178112615275;
			(*m->getParameterModel()->getParameterSet()[d])(0, 29) =
					-0.445003669917367;
			(*m->getParameterModel()->getParameterSet()[d])(0, 30) =
					-1.85115355373292;
			(*m->getParameterModel()->getParameterSet()[d])(0, 31) =
					-1.5272787583475;
			(*m->getParameterModel()->getParameterSet()[d])(0, 32) =
					-0.588667339897699;
			(*m->getParameterModel()->getParameterSet()[d])(0, 33) =
					0.987951431011379;
			(*m->getParameterModel()->getParameterSet()[d])(0, 34) =
					-2.91798964751742;
			(*m->getParameterModel()->getParameterSet()[d])(0, 35) =
					-0.559412036351752;
			(*m->getParameterModel()->getParameterSet()[d])(0, 36) =
					-0.0746886675230089;
			(*m->getParameterModel()->getParameterSet()[d])(0, 37) =
					-0.860005617142715;
			(*m->getParameterModel()->getParameterSet()[d])(0, 38) =
					-0.751107888939827;
			(*m->getParameterModel()->getParameterSet()[d])(0, 39) =
					-1.96933326133621;
			(*m->getParameterModel()->getParameterSet()[d])(0, 40) =
					1.39850515263464;
			(*m->getParameterModel()->getParameterSet()[d])(0, 41) =
					0.934922614426893;
			(*m->getParameterModel()->getParameterSet()[d])(0, 42) =
					2.25800750013661;
			(*m->getParameterModel()->getParameterSet()[d])(0, 43) =
					0.675815082921858;
			(*m->getParameterModel()->getParameterSet()[d])(0, 44) =
					-0.0451555954167578;
			(*m->getParameterModel()->getParameterSet()[d])(0, 45) =
					-1.01338578098843;
			(*m->getParameterModel()->getParameterSet()[d])(0, 46) =
					-0.408494649436581;
			(*m->getParameterModel()->getParameterSet()[d])(0, 47) =
					1.04856058360661;
			(*m->getParameterModel()->getParameterSet()[d])(0, 48) =
					-0.216238170167062;
			(*m->getParameterModel()->getParameterSet()[d])(0, 49) =
					-0.311803829141431;
			(*m->getParameterModel()->getParameterSet()[d])(0, 50) =
					0.685488331328067;
			(*m->getParameterModel()->getParameterSet()[d])(0, 51) =
					-2.38365746189177;
			(*m->getParameterModel()->getParameterSet()[d])(0, 52) =
					-0.490501944430744;
			(*m->getParameterModel()->getParameterSet()[d])(0, 53) =
					2.12895158575715;
			(*m->getParameterModel()->getParameterSet()[d])(0, 54) =
					0.856284068445749;
			(*m->getParameterModel()->getParameterSet()[d])(0, 55) =
					-0.545115885651256;
			(*m->getParameterModel()->getParameterSet()[d])(0, 56) =
					0.710331472101276;
			(*m->getParameterModel()->getParameterSet()[d])(0, 57) =
					1.620011667106;
			(*m->getParameterModel()->getParameterSet()[d])(0, 58) =
					-1.619271389274;
			(*m->getParameterModel()->getParameterSet()[d])(0, 59) =
					-0.858144564233935;
			(*m->getParameterModel()->getParameterSet()[d])(0, 60) =
					1.54968254624868;
			(*m->getParameterModel()->getParameterSet()[d])(0, 61) =
					-1.79267801974796;
			(*m->getParameterModel()->getParameterSet()[d])(0, 62) =
					1.14642671851194;
			(*m->getParameterModel()->getParameterSet()[d])(0, 63) =
					-0.0231556281292887;
			(*m->getParameterModel()->getParameterSet()[d])(0, 64) =
					-2.39750585134654;
			(*m->getParameterModel()->getParameterSet()[d])(0, 65) =
					-1.83443982471499;
			(*m->getParameterModel()->getParameterSet()[d])(0, 66) =
					1.50001049162294;
			(*m->getParameterModel()->getParameterSet()[d])(0, 67) =
					-1.10146147008743;
			(*m->getParameterModel()->getParameterSet()[d])(0, 68) =
					0.8074863386241;
			(*m->getParameterModel()->getParameterSet()[d])(0, 69) =
					0.393215270239003;
			(*m->getParameterModel()->getParameterSet()[d])(0, 70) =
					0.397327228486701;
			(*m->getParameterModel()->getParameterSet()[d])(0, 71) =
					-0.442055389800962;
			(*m->getParameterModel()->getParameterSet()[d])(0, 72) =
					0.454443066816656;
			(*m->getParameterModel()->getParameterSet()[d])(0, 73) =
					1.12682418983781;
			(*m->getParameterModel()->getParameterSet()[d])(0, 74) =
					-0.837091198602913;
			(*m->getParameterModel()->getParameterSet()[d])(0, 75) =
					1.36458176012761;
			(*m->getParameterModel()->getParameterSet()[d])(0, 76) =
					-0.116983973327465;
			(*m->getParameterModel()->getParameterSet()[d])(0, 77) =
					1.56348620206992;
			(*m->getParameterModel()->getParameterSet()[d])(0, 78) =
					-0.419671630600554;
			(*m->getParameterModel()->getParameterSet()[d])(0, 79) =
					0.450902411677853;
			(*m->getParameterModel()->getParameterSet()[d])(0, 80) =
					-2.69627444177271;
			(*m->getParameterModel()->getParameterSet()[d])(0, 81) =
					-0.142495092212555;
			(*m->getParameterModel()->getParameterSet()[d])(0, 82) =
					0.153516452580127;
			(*m->getParameterModel()->getParameterSet()[d])(0, 83) =
					-0.371501010161694;
			(*m->getParameterModel()->getParameterSet()[d])(0, 84) =
					0.269209250206473;
			(*m->getParameterModel()->getParameterSet()[d])(0, 85) =
					-1.33505315586072;
			(*m->getParameterModel()->getParameterSet()[d])(0, 86) =
					-0.792733824801248;
			(*m->getParameterModel()->getParameterSet()[d])(0, 87) =
					-0.88487364983854;
			(*m->getParameterModel()->getParameterSet()[d])(0, 88) =
					-1.19116957722689;
			(*m->getParameterModel()->getParameterSet()[d])(0, 89) =
					-0.646865603719435;
			(*m->getParameterModel()->getParameterSet()[d])(0, 90) =
					-0.448542564109881;
			(*m->getParameterModel()->getParameterSet()[d])(0, 91) =
					0.595846146713524;
			(*m->getParameterModel()->getParameterSet()[d])(0, 92) =
					3.0878294021287;
			(*m->getParameterModel()->getParameterSet()[d])(0, 93) =
					-0.22554273716775;
			(*m->getParameterModel()->getParameterSet()[d])(0, 94) =
					0.151195884024596;
			(*m->getParameterModel()->getParameterSet()[d])(0, 95) =
					0.212749870413368;
			(*m->getParameterModel()->getParameterSet()[d])(0, 96) =
					0.222053103595936;
			(*m->getParameterModel()->getParameterSet()[d])(0, 97) =
					-0.300124140888755;
			(*m->getParameterModel()->getParameterSet()[d])(0, 98) =
					-0.107131397076053;
			(*m->getParameterModel()->getParameterSet()[d])(0, 99) =
					-0.700020844834876;
		}
	}

	virtual void stepE(Model* m, Matrix<double>* f, Matrix<double>* r,
			QuadratureNodes* nodes) {

		/*
		 * What we need
		 * q
		 * a pattern iterator
		 * item number
		 * success probability matrix
		 * thetas
		 * weights
		 * parameter set
		 */
		//Dataset by patterns
		PatternMatrix* data =
				dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
		//Pattern iterator is data->iterator
		//Item number
		const double items = data->countItems();
		//Success probability matrix is obtained via pm->getProbability(int,int)
		ParameterModel* pm = m->getParameterModel();
		//Thetas
		Matrix<double>* thetas = nodes->getTheta();
		//Amount of nodes
		const int q = nodes->size();
		//Weights
		Matrix<double>* weights = nodes->getWeight();
		//A Matrix
		Matrix<double>* A = m->getParameterModel()->getParameterSet()[a];
		//B Matrix
		Matrix<double>* B = m->getParameterModel()->getParameterSet()[d];
		//C Matrix
		//Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];
		//Auxiliar array for the nodes
		long double faux[q];
		long double sum = 0.0;
		//Restart f and r to zero
		f->reset();
		r->reset();
		//Calculates the success probability for all the nodes.
		m->successProbability(nodes);

		int k, i;
		double prob;
		boost::dynamic_bitset<> current_bitset;

		//TODO CAREFULLY PARALLELIZE FOR
		for (data->resetIterator(); !data->checkEnd(); data->iterate()) {

			current_bitset = data->getCurrentBitSet();

			sum = 0.0;
			//Calculate g*(k) for all the k's
			//first calculate the P for each k and store it in the array f aux
			for (k = 0; k < q; k++) {
				faux[k] = (*weights)(0, k);
				//Calculate the p (iterate over the items in the productory)
				for (i = 0; i < items; i++) {
					prob = pm->getProbability(k, i);
					if (!current_bitset[items - i - 1]) {
						prob = 1 - prob;
					}
					faux[k] = faux[k] * prob;
				}
				//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
				//Now multiply by the weight
				sum += faux[k];
			}

			long double current_frequency =
					(long double) data->getCurrentFrequency();

			for (k = 0; k < q; k++) {
				faux[k] = faux[k] / sum; //This is g*_j_k
				//Multiply the f to the frequency of the pattern
				faux[k] = (current_frequency) * faux[k];
				(*f)(0, k) += faux[k];
				//Now selectively add the faux to the r
				for (i = 0; i < items; i++) {
					if (current_bitset[items - i - 1]) {
						(*r)(k, i) += faux[k];
					} // if
				} // for
			} // for
		}

	}

	virtual void stepM(Model* m, Matrix<double>* f, Matrix<double>* r,
			QuadratureNodes* nodes) {
		/*
		 */
		//Step M implementation using the BFGS Algorithm
		/*
		 * What we need
		 * fptr the pointer to loglik
		 * gprt the pointer to gradient
		 * hessptr the pointer to hessian matrix calculatrix
		 * args the a,b, and c
		 * pars, the other parameters q and stuff
		 * nargs, npars, sizes.
		 */
		//fptr
		double (*fptr)(double*, double*, int, int);
		void (*gptr)(double*, double*, int, int, double*);
		void (*hptr)(double*, double*, int, int, double*);
		fptr = &TwoPLModel::logLikelihood;
		gptr = &TwoPLModel::gradient;
		hptr = NULL;
		//cout<<"Address : "<<&gptr<<" "<<&hptr<<endl;
		int It = m->getItemModel()->getDataset()->countItems();
		int q = nodes->size();
		double args[3 * It];
		double pars[2 + 2 * q + q * It];
		int nargs = 3 * It;
		int npars = 2 + 2 * q + q * It;
		//filling args
		int nA = 0;
		// Obtain a
		//A Matrix
		Matrix<double>* A = m->getParameterModel()->getParameterSet()[a];
		//B Matrix
		Matrix<double>* B = m->getParameterModel()->getParameterSet()[d];
		//C Matrix
		//Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];

		Matrix<double> DA(*A);
		Matrix<double> DB(*B);
		//Matrix<double> DC(*C);

		for (int i = 0; i < It; i++) {
			args[nA] = (*A)(0, i);
			nA++;
		}

		// Obtain b
		for (int i = 0; i < It; i++) {
			args[nA] = (*B)(0, i);
			nA++;
		}

		//Filling pars
		int nP = 0;
		// Obtain q
		pars[nP] = q;
		nP++;
		// Obtain I
		pars[nP] = It;
		nP++;
		// Obtain theta
		//Thetas

		Matrix<double>* thetas = nodes->getTheta();
		for (int k = 0; k < q; k++) {
			pars[nP] = (*thetas)(0, k);	//TODO correct indexing on this and nearby matrices
			nP++;
		}
		// Obtain f
		for (int k = 0; k < q; k++) {
			pars[nP] = (*f)(0, k);
			nP++;
		}
		// Obtain r
		for (int k = 0; k < q; k++) {
			for (int i = 0; i < It; i++) {
				pars[nP] = (*r)(k, i);
				nP++;
			}
		}
		nargs = nA;
		npars = nP;
		/*
		 * Chooses the method
		 * method 1 is NR
		 * method 2 is BFGS
		 */
		Optimizer* optim;
		optim = new Optimizer();
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars);

		// Now pass the optimals to the Arrays.

		nA = 0;
		// Obtain a
		for (int i = 0; i < It; i++) {
			(*A)(0, i) = args[nA++];
			if (fabs((*A)(0, i)) > abs(10)) {
				(*A)(0, i) = 0.851;
				//			cout << "A reset." << endl;
			}

		}
		// Obtain b
		for (int i = 0; i < It; i++) {
			(*B)(0, i) = args[nA++];
			if (fabs((*B)(0, i)) > abs(-50)) {
				(*B)(0, i) = 0.5;
				//			cout << "B reset." << endl;
			}
		}

		//Boundary regularize the arguments
		//B = 0.5;
		//A = 0.851

		//Obtain the deltas
		//Perform substracts
		double maxDelta = 0;
		double meanDelta = 0;
		int DeltaC = 0;
		for (int v1 = 0; v1 < It; ++v1) {
			DA(0, v1) = DA(0, v1) - (*A)(0, v1);
			DB(0, v1) = DB(0, v1) - (*B)(0, v1);
			meanDelta = +fabs(DA(0, v1));
			meanDelta = +fabs(DB(0, v1));
			DeltaC += 3;
			if (fabs(DA(0, v1)) > maxDelta) {
				maxDelta = fabs(DA(0, v1));
			}
			if (fabs(DB(0, v1)) > maxDelta) {
				maxDelta = fabs(DB(0, v1));
			}
		}
		meanDelta = meanDelta / DeltaC;
		if (meanDelta < 0.0001 and maxDelta < 0.001) {
			m->itemParametersEstimated = true;
		}
		//And set the parameter sets
		map<Parameter, Matrix<double> *> parSet;
		parSet[a] = A;
		parSet[d] = B;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}
	;

};
#endif /* EM2PL_H_ */
