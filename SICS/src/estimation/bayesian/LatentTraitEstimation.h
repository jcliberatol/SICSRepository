/*
 * LatentTraitEstimation.h
 *
 *  Created on: Oct 16, 2014
 *      Author: jliberato
 */

#ifndef LATENTTRAITESTIMATION_H_
#define LATENTTRAITESTIMATION_H_
#include <model/Model.h>
#include <type/LatentTraits.h>
#include <sstream>
#include <type/Constant.h>
#include <optimizer/Brent_fmin.h>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#define t_zita t_model->parameterModel->parameterSet
#define gg t_model->parameterModel->successProbability

class LatentTraitEstimation {
public:

	struct parameters_logL {
		bool * pattern;
		int size;
		int node;
		Model *model;
	};

	LatentTraitEstimation() {
	}

	virtual ~LatentTraitEstimation() {
	}

	Model* model;
	QuadratureNodes * quadNodes;
	LatentTraits * lt;

	//Deprecated
	inline double patternProbabilities(vector<char> pattern, int node) {
		double p = 1;

		for (int i = 0; i < pattern.size(); i++) {
			if (pattern.at(i) > 0) {
				p *= (*model->parameterModel->probabilityMatrix)(node, i);
			} else {
				p *= 1 - (*model->parameterModel->probabilityMatrix)(node, i);
			}
		}
		return (p);
	}

	inline double patternProbabilities(bool * pattern, int size, int node) {
		double p = 1;

		for (int i = 0; i < size; i++) {
			if (pattern[i]) {
				p *= (*model->parameterModel->probabilityMatrix)(node, i);
			} else {
				p *= 1 - (*model->parameterModel->probabilityMatrix)(node, i);
			}
		}
		return (p);
	}

	//Deprecated
	static double patternProbabilities(double theta, vector<char> pattern,
			int node, Model * t_model) {
		double p = 1;

		for (int i = 0; i < pattern.size(); i++) {
			if (pattern.at(i) > 0) {
				p *= gg(theta, new double[3] { t_zita[0][0][i], t_zita[1][0][i],
				t_zita[2][0][i] });

			} else {
				p *= 1 - gg(theta, new double[3] { t_zita[0][0][i],
				t_zita[1][0][i], t_zita[2][0][i] });
			}
		}
		return (p);
	}

	static double patternProbabilities(double theta, bool * pattern, int size,
			int node, Model * t_model) {
		double p = 1;

		for (int i = 0; i < size; i++) {
			if (pattern[i] > 0) {
				p *= gg(theta, new double[3] { t_zita[0][0][i], t_zita[1][0][i],
				t_zita[2][0][i] });

			} else {
				p *= 1 - gg(theta, new double[3] { t_zita[0][0][i],
				t_zita[1][0][i], t_zita[2][0][i] });
			}
		}
		return (p);
	}

	static double patternProbabilities(double theta, bool * pattern, int size,
			int node, Model * t_model, double *** parSet) {
		double p = 1;

		for (int i = 0; i < size; i++) {
			if (pattern[i] > 0) {
				p *= gg(theta, new double[3] { parSet[0][0][i], parSet[1][0][i],
						parSet[2][0][i] });

			} else {
				p *= 1 - gg(theta, new double[3] { parSet[0][0][i],
						parSet[1][0][i], parSet[2][0][i] });
			}
		}
		return (p);
	}

	LatentTraits * getLatentTraits() {
		return (lt);
	}

	void setLatentTraits(LatentTraits * ltt) {
		lt = ltt;
	}

	void setModel(Model* m) {
		model = m;
	}

	//Deprecated
	void estimateLatentTraitsEAP_() {

		map<vector<char>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator begin = lt->pm->matrix.begin();
		map<vector<char>, int>::const_iterator end = lt->pm->matrix.end();

		int counter = 0;

		for (it = begin; it != end; ++it, ++counter) {
			double sum_num = 0;
			double sum_den = 0;

			for (int i = 0; i < quadNodes->size(); ++i) {
				double pp = patternProbabilities(it->first, i);
				sum_num += (*quadNodes->getTheta())(0, i)
						* ((*quadNodes->getWeight())(0, i)) * pp;
				sum_den += (*quadNodes->getWeight())(0, i) * pp;
			}

			(*lt->traits)(counter, lt->dim - 1) = sum_num / sum_den;
		}
	}

	void estimateLatentTraitsEAP() {

		bool ** pattern_list = lt->pm->getBitsetList();
		int size = lt->pm->matrix.size();

		int counter = 0;

		for (int index = 0; index < size; index++, ++counter) {
			double sum_num = 0;
			double sum_den = 0;

			for (int i = 0; i < quadNodes->size(); ++i) {
				double pp = patternProbabilities(pattern_list[index],
						lt->pm->size, i);
				sum_num += (*quadNodes->getTheta())(0, i)
						* ((*quadNodes->getWeight())(0, i)) * pp;
				sum_den += (*quadNodes->getWeight())(0, i) * pp;
			}

			(*lt->traits)(counter, lt->dim - 1) = sum_num / sum_den;
		}
	}

	//Deprecated
	static double logL_(double theta, vector<char> pattern, int node,
			Model *model) {
		return (-(log(patternProbabilities(theta, pattern, node, model))
				- ((theta * theta) / 2)));

	}

	static double logL(double theta, bool * pattern, int size, int node,
			Model *model) {
		return (-(log(patternProbabilities(theta, pattern, size, node, model))
				- ((theta * theta) / 2)));
	}

	static double logLP(double theta, bool * pattern, int size, int node,
			Model *model, double *** parSet) {
		return (-(log(
				patternProbabilities(theta, pattern, size, node, model, parSet))
				- ((theta * theta) / 2)));
	}

	static double logLR(double theta, void * params) {
		parameters_logL * temp = (struct parameters_logL*) params;

		return (-(log(
				patternProbabilities(theta, temp->pattern, temp->size,
						temp->node, temp->model)) - ((theta * theta) / 2)));
	}

	void printVectors() {
		bool ** pattern_list = lt->pm->getBitsetList();
		int * frequency_list = lt->pm->getFrequencyList();
		int size = lt->pm->matrix.size();

		int counter = 0;

		for (double i = -3.3; i <= 2.8; i += .01) {
			cout << i << ","
					<< logL(i, pattern_list[0], lt->pm->size, counter,
							this->model) << endl;
		}
	}

	void evaluate_theta(double theta, int pattern) {
		bool ** pattern_list = lt->pm->getBitsetList();
		int * frequency_list = lt->pm->getFrequencyList();
		int size = lt->pm->matrix.size();

		int counter = 0;

		cout << theta << ","
				<< logL(theta, pattern_list[pattern], lt->pm->size, counter,
						this->model) << endl;

	}

//Deprecated
	void estimateLatentTraitsMAP_() {
		map<vector<char>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator begin = lt->pm->matrix.begin();
		map<vector<char>, int>::const_iterator end = lt->pm->matrix.end();

		int counter = 0;

		for (it = begin; it != end; ++it, ++counter) {
			double (*function)(double, vector<char>, int, Model *) = &logL_;
			(*lt->traits)(counter, lt->dim - 1) = Brent_fmin(new double[2] { -5,
					5 }, 0.0001220703, function, it->first, counter,
					this->model, 1);
		}
	}

	void estimateLatentTraitsMAP() {
		bool ** pattern_list = lt->pm->getBitsetList();
		int * frequency_list = lt->pm->getFrequencyList();
		int size = lt->pm->matrix.size();

		int counter = 0;

		for (int index = 0; index < size; index++, ++counter) {
			double (*function)(double, bool *, int, int, Model *) = &logL;
			(*lt->traits)(counter, lt->dim - 1) = Brent_fmin(new double[2] { -5,
					5 }, 0.0001220703, function, pattern_list[index],
					lt->pm->size, counter, this->model, 1);
		}

		delete frequency_list;
	}

	void estimateLatentTraitsMAP(double *** parSet) {
		bool ** pattern_list = lt->pm->getBitsetList();
		int * frequency_list = lt->pm->getFrequencyList();
		int size = lt->pm->matrix.size();

		int counter = 0;

		for (int index = 0; index < size; index++, ++counter) {
			double (*function)(double, bool *, int, int, Model *, double ***) = &logLP;
			(*lt->traits)(counter, lt->dim - 1) = Brent_fmin(new double[2] { -5,
					5 }, 0.0001220703, function, pattern_list[index],
					lt->pm->size, counter, this->model, parSet, 1);
		}

		delete frequency_list;
	}

	void estimateLatentTraitsMAP_R() {
		bool ** pattern_list = lt->pm->getBitsetList();
		int * frequency_list = lt->pm->getFrequencyList();
		int size = lt->pm->matrix.size();
		parameters_logL temp;

		int counter = 0;

		for (int index = 0; index < size; index++, ++counter) {
			temp.model = model;
			temp.node = counter;
			temp.pattern = pattern_list[index];
			temp.size = lt->pm->size;
			double (*function)(double, void*) = &logLR;
			(*lt->traits)(counter, lt->dim - 1) = Brent_fmin(-5, 5, function,
					(void*) &temp, 0.0001220703);
		}

		for (int j = 0; j < size; j++) {
			delete pattern_list[j];
		}

		delete pattern_list;
		delete frequency_list;
	}

	void estimateLatentTraitsMAP_GSL() {
		int status;
		int iter = 0, max_iter = 100;
		const gsl_min_fminimizer_type *T;
		parameters_logL temp;
		gsl_min_fminimizer *s;
		gsl_function F;
		F.function = &logLR;
		F.params = &temp;

		T = gsl_min_fminimizer_quad_golden;
		s = gsl_min_fminimizer_alloc(T);

		double m = 2.0, m_expected = M_PI;
		double a = -5.0, b = 5.0;

		gsl_min_fminimizer_set(s, &F, m, a, b);

		bool ** pattern_list = lt->pm->getBitsetList();
		int * frequency_list = lt->pm->getFrequencyList();
		int size = lt->pm->matrix.size();

		int counter = 0;

		for (int index = 0; index < size; index++, ++counter) {
			cout << "." << endl;
			temp.model = model;
			temp.node = counter;
			temp.pattern = pattern_list[index];
			temp.size = lt->pm->size;

			F.params = (void*) &temp;

			do {
				status = gsl_min_fminimizer_iterate(s);

				m = gsl_min_fminimizer_x_minimum(s);
				a = gsl_min_fminimizer_x_lower(s);
				b = gsl_min_fminimizer_x_upper(s);

				status = gsl_min_test_interval(a, b, 0.001, 0.0);

			} while (status == GSL_CONTINUE && iter < max_iter);

			(*lt->traits)(counter, lt->dim - 1) = m;
		}

		for (int j = 0; j < size; j++) {
			delete pattern_list[j];
		}

		delete pattern_list;
		delete frequency_list;
	}

	void setQuadratureNodes(QuadratureNodes *nodes) {
		quadNodes = nodes;
		model->successProbability(quadNodes);
	}
};

#endif /* LATENTTRAITESTIMATION_H_ */
