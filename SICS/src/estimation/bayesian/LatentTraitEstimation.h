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

class LatentTraitEstimation {
public:
	LatentTraitEstimation() {

	}
	virtual ~LatentTraitEstimation() {
	}

	Model* model;
	QuadratureNodes * quadNodes;
	LatentTraits * lt;

	inline double patternProbabilities(vector<char> pattern, int node) {
		double p = 1;

		for (int i = 0; i < pattern.size(); i++) {
			if (pattern.at(i) > 0) {
				p *= (*model->parameterModel->probabilityMatrix)(node, i);
			} else {
				p *= 1 - (*model->parameterModel->probabilityMatrix)(node, i);
			}
		}
		return p;
	}

	static double patternProbabilities(double theta, vector<char> pattern, int node,
			Model * t_model) {
		double p = 1;

		for (int i = 0; i < pattern.size(); i++) {
			if (pattern.at(i) > 0) {
				p *= t_model->parameterModel->successProbability(theta,
						new double[3] { t_zita[0][0][i], t_zita[1][0][i] });

			} else {
				p *= 1
						- t_model->parameterModel->successProbability(theta,
								new double[3] { t_zita[0][0][i],
								t_zita[1][0][i] });
			}
		}
		return p;
	}

	LatentTraits * getLatentTraits() {
		return lt;
	}

	void setLatentTraits(LatentTraits * ltt) {
		lt = ltt;
	}

	void setModel(Model* m) {
		model = m;
	}

	void estimateLatentTraitsEAP() {

		map<vector<char>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator begin = lt->pm->matrix.begin();
		map<vector<char>, int>::const_iterator end = lt->pm->matrix.end();

		int counter = 0;

		for (it = begin; it != end; ++it, ++counter) {
			double sum_num = 0;
			double sum_den = 0;

//			TODO export output to test with Liberato profiler.

//			for (int i = 0; i < 10; i++) {
//				cout << (int) it->first.at(i) << " ";
//			}
//			cout << endl;

			for (int i = 0; i < quadNodes->size(); ++i) {
				double pp = patternProbabilities(it->first, i);
				sum_num += (*quadNodes->getTheta())(0, i)
						* ((*quadNodes->getWeight())(0, i)) * pp;
				sum_den += (*quadNodes->getWeight())(0, i) * pp;
			}

			(*lt->traits)(counter, lt->dim - 1) = sum_num / sum_den;
		}
	}

	static double logL(double theta, vector<char> pattern, int node,
			Model *model) {
		return -(log(patternProbabilities(theta, pattern, node, model))
				- ((theta * theta) / 2));

	}

	void estimateLatentTraitsMAP() {
		map<vector<char>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator begin = lt->pm->matrix.begin();
		map<vector<char>, int>::const_iterator end = lt->pm->matrix.end();

		int counter = 0;

		for (it = begin; it != end; ++it, ++counter) {

			for (int i = 0; i < 10; i++) {
				cout << (int) it->first.at(i) << " ";
			}
			cout << endl;

			double (*function)(double, vector<char>, int, Model *) = &logL;
			(*lt->traits)(counter, lt->dim - 1) = Brent_fmin(new double[2] { -5,
					5 }, 0.0001220703, function, it->first, counter,
					this->model, 100000);
		}
	}

	void setQuadratureNodes(QuadratureNodes *nodes) {
		quadNodes = nodes;
		model->successProbability(quadNodes);
	}
};

#endif /* LATENTTRAITESTIMATION_H_ */
