/*
 * EMTest.cpp
 *
 *  Created on: 25 Jul 2014
 *      Author: jlgpisa
 */

#include "EMTest.h"

EMTest::EMTest() {

	pob = NULL;
	convergence = NULL;
	initialValues = NULL;
	pM = NULL;

}

EMTest::EMTest(string configFile) {
	this->configFile = configFile;
	loadConfiguration();
}

const string& EMTest::getConfigFile() const {
	return (configFile);
}

void EMTest::setConfigFile(const string& configFile) {
	this->configFile = configFile;
	loadConfiguration();
}

map<InTestFiletype, TestConfigFile*> EMTest::getInputConfig() {
	return (inputConfig);
}

void EMTest::setInputConfig(map<InTestFiletype, TestConfigFile*> inputConfig) {
	this->inputConfig = inputConfig;
}

void EMTest::loadConfiguration() {
	TestConfigFile * testConfigFile;
	string currentRead;
	ifstream pFile;

	pFile.open(configFile.c_str(), std::ifstream::in);

	while (pFile >> currentRead) {

		if (currentRead == TestConstant::kNewConf) {
			testConfigFile = new TestConfigFile();
		}

		if (currentRead == TestConstant::kFileType) {
			pFile >> currentRead;
			inputConfig[TestConfigFile::getFiletype(currentRead)] =
					testConfigFile;
		}

		if (currentRead == TestConstant::kFullPath) {
			pFile >> currentRead;
			testConfigFile->setFullPath(currentRead);
		}

		if (currentRead == TestConstant::kHeader) {
			int header;
			pFile >> header;
			testConfigFile->setHeader(header);
		}

		if (currentRead == TestConstant::kSep) {
			char sep;
			pFile >> sep;
			testConfigFile->setSeparator(sep);
		}

		if (currentRead == TestConstant::kNodeCount) {
			int nodes;
			pFile >> nodes;
			TestConstant::q = nodes;
		}
	}

	pFile.close();
}

void EMTest::loadInput(map<InTestFiletype, string> paths) {
	TestInput * testInput = new TestInput();
	int items = 0;

	if (inputConfig[DATASET] != NULL) {
		testInput->setConfig(inputConfig[DATASET]);
		testInput->setPath(paths[DATASET]);

		pM = testInput->loadPattern();
		items = pM->countItems();
	}

	if (inputConfig[INITIAL_VALUE] != NULL) {
		testInput->setConfig(inputConfig[DATASET]);
		testInput->setPath(paths[INITIAL_VALUE]);

		initialValues = testInput->loadMatrix(3, items);
	}

	if (inputConfig[CONVERGENCE] != NULL) {
		testInput->setConfig(inputConfig[CONVERGENCE]);
		testInput->setPath(paths[CONVERGENCE]);

		convergence = testInput->loadMatrix(3, items);
	}

	if (inputConfig[POB] != NULL) {
		testInput->setConfig(inputConfig[POB]);
		testInput->setPath(paths[POB]);

		pob = testInput->loadMatrixTransp(items, 3);
	}

}

void EMTest::runTest() {

	string reportFilename = "TestReport.json";
	TestReport report(reportFilename);

	// Path file streams
	ifstream dataSetF;
	ifstream convergenceF;
	ifstream initialValuesF;
	ifstream pobF;

	map<InTestFiletype, string> paths;

	loadConfiguration();

	dataSetF.open(inputConfig[DATASET]->getFullPath().c_str(), ifstream::in);
	convergenceF.open(inputConfig[CONVERGENCE]->getFullPath().c_str(),
			ifstream::in);
	initialValuesF.open(inputConfig[INITIAL_VALUE]->getFullPath().c_str(),
			ifstream::in);
	pobF.open(inputConfig[POB]->getFullPath().c_str(), ifstream::in);

	// load nodes
	Input input;
	Matrix<double> cuad(41, 2);
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);

	Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);
	for (int k = 0; k < cuad.nR(); k++) {
		(*theta)(0, k) = cuad(k, 0);
		(*weight)(0, k) = cuad(k, 1);
	}

	while (getline(dataSetF, paths[DATASET])) {

		// Obtain remaining file locations
		getline(initialValuesF, paths[INITIAL_VALUE]);
		getline(convergenceF, paths[CONVERGENCE]);
		getline(pobF, paths[POB]);

		loadInput(paths);

		// Create model
		Model *model = new Model();
		ModelFactory *modelFactory = new SICSGeneralModel();
		model->setModel(modelFactory, Constant::THREE_PL);

		// Set dataset to model
		model->getItemModel()->setDataset(pM);

		// Set nodes and weights to model
		model->getDimensionModel()->getLatentTraitSet()->setTheta(theta);
		model->getDimensionModel()->getLatentTraitSet()->setWeight(weight);

		// Build parameter set TODO: a model function that loads item and dimension
		model->getParameterModel()->buildParameterSet(model->getItemModel(),
				model->getDimensionModel());

		// Initial parameters
		int items = model->getItemModel()->countItems();
		for (int i = 0; i < items; i++) {
			(*model->getParameterModel()->getParameterSet()[a])(0, i) =
					(*initialValues)(0, i);
			(*model->getParameterModel()->getParameterSet()[d])(0, i) =
					(*initialValues)(1, i);
			(*model->getParameterModel()->getParameterSet()[c])(0, i) =
					(*initialValues)(2, i);
		}

		cout << *initialValues << endl;

		// Create estimation
		EMEstimation *em = new EMEstimation();
		em->setModel(model);

		report.startTime();
		em->estimate();
		report.endTime();

		report.setModel(model);
		report.setEmEstimation(em);
		report.setReportName(paths[DATASET]);
		report.addJsonResult(convergence, pob);

		delete modelFactory;
		//delete em;
		//delete model;

		delete pM;
		delete initialValues;
		delete convergence;
		delete pob;

	}

	report.reportJSON();

	dataSetF.close();
	convergenceF.close();
	initialValuesF.close();
	pobF.close();

	//delete theta;
	//delete weight;
}

void EMTest::runProcessor() {

	string reportFilename = "processReport.json";
	TestReport report(reportFilename);

	// Path file streams
	ifstream dataSetF;
	ifstream initialValuesF;

	map<InTestFiletype, string> paths;

	loadConfiguration();

	dataSetF.open(inputConfig[DATASET]->getFullPath().c_str(), ifstream::in);
	initialValuesF.open(inputConfig[INITIAL_VALUE]->getFullPath().c_str(),
			ifstream::in);

	// load nodes
	Input input;
	Matrix<double> cuad(41, 2);
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);

	Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);
	for (int k = 0; k < cuad.nR(); k++) {
		(*theta)(0, k) = cuad(k, 0);
		(*weight)(0, k) = cuad(k, 1);
	}

	while (getline(dataSetF, paths[DATASET])) {
		getline(initialValuesF, paths[INITIAL_VALUE]);

		loadInput(paths);

		// Create model
		Model *model = new Model();
		ModelFactory *modelFactory = new SICSGeneralModel();
		model->setModel(modelFactory, Constant::THREE_PL);

		// Set dataset to model
		model->getItemModel()->setDataset(pM);

		// Set nodes and weights to model
		model->getDimensionModel()->getLatentTraitSet()->setTheta(theta);
		model->getDimensionModel()->getLatentTraitSet()->setWeight(weight);

		// Build parameter set TODO: a model function that loads item and dimension
		model->getParameterModel()->buildParameterSet(model->getItemModel(),
				model->getDimensionModel());

		// Initial parameters
		int items = model->getItemModel()->countItems();
		for (int i = 0; i < items; i++) {
			(*model->getParameterModel()->getParameterSet()[a])(0, i) =
					(*initialValues)(0, i);
			(*model->getParameterModel()->getParameterSet()[d])(0, i) =
					(*initialValues)(1, i);
			(*model->getParameterModel()->getParameterSet()[c])(0, i) =
					(*initialValues)(2, i);
		}

		cout << *initialValues << endl;

		// Create estimation
		EMEstimation *em = new EMEstimation();
		em->setModel(model);

		report.startTime();
		em->estimate();
		report.endTime();

		report.setModel(model);
		report.setEmEstimation(em);
		report.setReportName(paths[DATASET]);
		report.reportProcess();

		delete modelFactory;
		//delete em;
		//delete model;

		delete pM;
		delete initialValues;

	}

	report.reportJSON();

	dataSetF.close();
	initialValuesF.close();
}

EMTest::~EMTest() {

}

