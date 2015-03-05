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
#define t_zita t_model->parameterModel->parameterSet
#define gg t_model->parameterModel->successProbability

class LatentTraitEstimation {
public:
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

		for (int j = 0; j < size; j++) {
			delete pattern_list[j];
		}

		delete pattern_list;
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
