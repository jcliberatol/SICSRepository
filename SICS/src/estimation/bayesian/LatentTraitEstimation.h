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

class LatentTraitEstimation {
public:
	LatentTraitEstimation();
	virtual ~LatentTraitEstimation();

	Model* model;
	QuadratureNodes * quadNodes;
	LatentTraits * lt;


	double patternProbabilities(vector<char> pattern, int node) {
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

	LatentTraits * getLatentTraits() {
		return lt;
	}

	void setLatentTraits(LatentTraits * ltt) {
		lt = ltt;
	}
	void setModel(Model* m){
			model = m;
		}
	void estimateLatentTraits() {

		map<vector<char>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator begin = lt->pm->matrix.begin();
		map<vector<char>, int>::const_iterator end = lt->pm->matrix.end();

		int counter = 0;

		for (it = begin; it != end; ++it, ++counter) {
			double sum_num = 0;
			double sum_den = 0;
			for(int i = 0; i < quadNodes->size() ; ++i){
				double pp = patternProbabilities(it->first, i);
				sum_num += (*quadNodes->getTheta())(0,i) * ((* quadNodes->getWeight())(0,i)) * pp;
				sum_den += (*quadNodes->getWeight())(0,i) * pp;
			}

			(*lt->traits)(counter, lt->dim -1) = sum_num/sum_den;
		}
	}

	void setQuadratureNodes(QuadratureNodes* nodes) {
		quadNodes = nodes;
		model->successProbability(quadNodes);
	}
};

#endif /* LATENTTRAITESTIMATION_H_ */
