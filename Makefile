TARGET = main
FC = gfortran

PREF_SRC = ./src/
PREF_OBJ = ./obj/

# 1. АВТОМАТИЧЕСКИЙ ПОИСК ВСЕХ ФАЙЛОВ НА ДИСКЕ
SRC_F90 = $(wildcard $(PREF_SRC)*.f90)
SRC_F   = $(wildcard $(PREF_SRC)*.f)

OBJ = $(patsubst $(PREF_SRC)%.f90, $(PREF_OBJ)%.o, $(SRC_F90)) \
      $(patsubst $(PREF_SRC)%.f, $(PREF_OBJ)%.o, $(SRC_F))

# 2. ФЛАГИ КОМПИЛЯЦИИ И ЛИНКОВКИ
FFLAGS = -Wall -O3 -march=native -funroll-loops -ftree-vectorize -fopenmp -I. -J $(PREF_OBJ)
CUBA_DIR = /home/dmitry/Cuba-4.2.2
LIBS = -L$(CUBA_DIR) -lcubaq -lm

.PHONY: all run clean rebuild

all: $(TARGET)

# Линковка исполняемого файла
$(TARGET) : $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $(TARGET) $(LIBS)
	@echo "✅ Сборка завершена!"

# Правило для компиляции .f90 файлов (с автоматическим созданием папки obj/)
$(PREF_OBJ)%.o : $(PREF_SRC)%.f90
	@mkdir -p $(PREF_OBJ)
	$(FC) -c $(FFLAGS) $< -o $@

# Правило для компиляции .f файлов
$(PREF_OBJ)%.o : $(PREF_SRC)%.f
	@mkdir -p $(PREF_OBJ)
	$(FC) -c $(FFLAGS) $< -o $@

# Удобные команды для запуска и очистки
run: $(TARGET)
	CUBACORES=1 ./$(TARGET)

clean:
	rm -rf $(PREF_OBJ) $(TARGET)
	@echo "🧹 Очистка завершена!"

rebuild: clean all


# ==============================================================================
# 3. ЯВНЫЕ ЗАВИСИМОСТИ МЕЖДУ МОДУЛЯМИ
# ==============================================================================
# Здесь мы гарантируем правильный порядок компиляции Fortran-модулей.
# Если вы добавите новый файл, который использует, например, parameters,
# просто допишите в этот блок: $(PREF_OBJ)new_file.o: $(PREF_OBJ)parameters.o

$(PREF_OBJ)ang.o: $(PREF_OBJ)parameters.o
$(PREF_OBJ)coeff_2004.o: $(PREF_OBJ)parameters.o
$(PREF_OBJ)quadpack_generic.o: $(PREF_OBJ)parameters.o
$(PREF_OBJ)others.o: $(PREF_OBJ)quantum_state.o
$(PREF_OBJ)prep_to_1pot.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)others.o
$(PREF_OBJ)angINT.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)ang.o
$(PREF_OBJ)dirac_pnt.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)ang.o
$(PREF_OBJ)quantum_state.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)dirac_pnt.o
$(PREF_OBJ)Integration_Cuba.o: $(PREF_OBJ)quantum_state.o $(PREF_OBJ)parameters.o

# Модули расчетов
$(PREF_OBJ)pot0_2004.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)quantum_state.o $(PREF_OBJ)coeff_2004.o $(PREF_OBJ)angINT.o $(PREF_OBJ)Integration_Cuba.o $(PREF_OBJ)others.o
$(PREF_OBJ)pot1_2004.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)prep_to_1pot.o $(PREF_OBJ)quantum_state.o $(PREF_OBJ)others.o $(PREF_OBJ)Integration_Cuba.o $(PREF_OBJ)angINT.o
$(PREF_OBJ)pot2_via_pot0.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)quantum_state.o $(PREF_OBJ)coeff_2004.o $(PREF_OBJ)angINT.o $(PREF_OBJ)Integration_Cuba.o $(PREF_OBJ)others.o $(PREF_OBJ)prep_to_1pot.o

# Главная программа
$(PREF_OBJ)main.o: $(PREF_OBJ)parameters.o $(PREF_OBJ)quantum_state.o $(PREF_OBJ)pot0_2004.o $(PREF_OBJ)pot1_2004.o $(PREF_OBJ)others.o $(PREF_OBJ)pot2_via_pot0.o
