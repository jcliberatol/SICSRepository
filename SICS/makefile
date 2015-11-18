SRC =   src/type/Constant.cpp \
	src/type/PatternMatrix.cpp \
	src/type/QuadratureNodes.cpp \
	src/input/Input.cpp \
	src/estimation/classical/EMEstimation.cpp \
	src/model/Model.cpp \
	src/model/SICSGeneralModel.cpp \
	src/model/dimension/DimensionModel.cpp \
	src/model/dimension/MultidimensionalModel.cpp \
	src/model/dimension/MultiUniDimModel.cpp \
	src/model/dimension/UnidimensionalModel.cpp \
	src/model/item/DichotomousModel.cpp \
	src/model/item/ItemModel.cpp \
	src/model/item/PolytomousModel.cpp \
	src/model/parameter/OnePLACModel.cpp \
	src/model/parameter/OnePLModel.cpp \
	src/model/parameter/ThreePLModel.cpp \
	src/model/parameter/TwoPLModel.cpp \
	src/optimizer/Optimizer.cpp \
	src/output/Output.cpp \
	src/util/asa111.cpp \
	src/util/util.cpp
OBJ = $(SRC:.cpp=.o)
INCL = -I/home/cristian/git/CppAD/ -I./src/
CFLAGS = -std=c++11 -Wunused-function -g3 -fdiagnostics-color

all: SICS

$(OBJ): %.o : %.h

#.cpp.o:
src/%.o: src/%.cpp
	g++ $(CFLAGS) -Wall -c -I./src/ -o $@ $<

SICS: $(OBJ)
	g++ -Wall $(CFLAGS) -I./src/ $^ src/main/Main2.cpp -o $@2
desde_cero:
	g++ $(CFLAGS) -I./src/ -o SICS $(SRC) src/main/Main.cpp

debug:
	g++ $(CFLAGS) -ggdb -I./src/ -o SICS_dbg $(SRC) src/main/Main2.cpp
clean:
	rm $(OBJ) SICS2
