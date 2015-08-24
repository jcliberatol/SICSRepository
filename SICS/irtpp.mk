SRCS =   src/type/Constant.cpp \
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
	src/util/asa111.cpp

LIBRARY = irtpp
SRC_DIR = src
CPPFLAGS = -std=c++11 -march=native -O3 -Wall -fPIC
INCLUDES = -I./$(SRC_DIR)
OBJS = $(SRCS:.cpp=.o)

all : lib$(LIBRARY).a

$(OBJS): %.o : %.h

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ $(CPPFLAGS) -c $(INCLUDES) -o $@ $<

lib$(LIBRARY).a : $(OBJS)
	ar rsv $@ $^

clean:
	$(RM) $(OBJS)

