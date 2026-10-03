FC = gfortran
# Добавлен флаг -J obj/ для сохранения модулей .mod в папке obj
FFLAGS = -Wall -O3 -march=native -funroll-loops -ftree-vectorize -fopenmp -I. -J obj/
CUBA_DIR = /home/dmitry/Cuba-4.2.2
LIBS = -L$(CUBA_DIR) -lcubaq -lm

# Папки для исходников и объектных файлов
SRC_DIR = src
OBJ_DIR = obj

# 1. СПИСОК ОБЪЕКТОВ (С указанием пути к папке obj/)
OBJECTS = $(addprefix $(OBJ_DIR)/, \
          parameters.o ang.o angINT.o coeff_2004.o dirac_pnt.o \
          GAMMA_MOD.o Integration_Cuba.o others.o pot0_2004.o pot2_via_pot0.o \
          pot1_2004.o prep_to_1pot.o quadpack_generic.o \
          quantum_state.o spec_func.o main.o)

TARGET = main

all: $(TARGET)

# 2. СБОРКА (Перед компиляцией автоматически создается папка obj/)
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)
	@echo "✅ Сборка завершена!"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# 3. ЯВНЫЕ ЗАВИСИМОСТИ (С учетом путей к OBJ_DIR. Исправлена опечатка в pot2_via_pot0)
$(OBJ_DIR)/parameters.o:
$(OBJ_DIR)/GAMMA_MOD.o:
$(OBJ_DIR)/spec_func.o:

$(OBJ_DIR)/ang.o: $(OBJ_DIR)/parameters.o
$(OBJ_DIR)/coeff_2004.o: $(OBJ_DIR)/parameters.o
$(OBJ_DIR)/quadpack_generic.o: $(OBJ_DIR)/parameters.o
$(OBJ_DIR)/others.o: $(OBJ_DIR)/quantum_state.o
$(OBJ_DIR)/prep_to_1pot.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/others.o
$(OBJ_DIR)/angINT.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/ang.o
$(OBJ_DIR)/dirac_pnt.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/ang.o
$(OBJ_DIR)/quantum_state.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/dirac_pnt.o
$(OBJ_DIR)/Integration_Cuba.o: $(OBJ_DIR)/quantum_state.o $(OBJ_DIR)/parameters.o

# Модули расчетов
# Модули расчетов
$(OBJ_DIR)/pot0_2004.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/quantum_state.o $(OBJ_DIR)/coeff_2004.o $(OBJ_DIR)/angINT.o $(OBJ_DIR)/Integration_Cuba.o $(OBJ_DIR)/others.o
$(OBJ_DIR)/pot1_2004.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/prep_to_1pot.o $(OBJ_DIR)/quantum_state.o $(OBJ_DIR)/others.o $(OBJ_DIR)/Integration_Cuba.o $(OBJ_DIR)/angINT.o
$(OBJ_DIR)/pot2_via_pot0.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/quantum_state.o $(OBJ_DIR)/coeff_2004.o $(OBJ_DIR)/angINT.o $(OBJ_DIR)/Integration_Cuba.o $(OBJ_DIR)/others.o $(OBJ_DIR)/prep_to_1pot.o


# Главная программа
$(OBJ_DIR)/main.o: $(OBJ_DIR)/parameters.o $(OBJ_DIR)/quantum_state.o $(OBJ_DIR)/pot0_2004.o $(OBJ_DIR)/pot1_2004.o $(OBJ_DIR)/others.o $(OBJ_DIR)/pot2_via_pot0.o

# 4. КОМАНДЫ
.PHONY: all run clean rebuild

run: $(TARGET)
	CUBACORES=1 ./$(TARGET)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
	@echo "🧹 Очистка завершена!"

rebuild: clean all
