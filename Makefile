FC = gfortran
FFLAGS = -Wall -O3 -march=native -funroll-loops -ftree-vectorize -fopenmp -I.
CUBA_DIR = /home/dmitry/Cuba-4.2.2
LIBS = -L$(CUBA_DIR) -lcubaq -lm

# 1. СПИСОК ОБЪЕКТОВ (Только один раз, без wildcard)
OBJECTS = parameters.o ang.o angINT.o coeff_2004.o dirac_pnt.o \
          GAMMA_MOD.o Integration_Cuba.o others.o pot0_2004.o \
          pot1_2004.o prep_to_1pot.o quadpack_generic.o \
          quantum_state.o spec_func.o main.o

TARGET = main

all: $(TARGET)

# 2. СБОРКА (Добавлен FFLAGS для OpenMP)
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)
	@echo "✅ Сборка завершена!"

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# 3. ЯВНЫЕ ЗАВИСИМОСТИ (Каждая цель встречается СТРОГО один раз)
parameters.o:
GAMMA_MOD.o:
spec_func.o:

ang.o: parameters.o
coeff_2004.o: parameters.o
quadpack_generic.o: parameters.o
others.o: quantum_state.o
prep_to_1pot.o: parameters.o others.o
angINT.o: parameters.o ang.o
dirac_pnt.o: parameters.o ang.o
quantum_state.o: parameters.o dirac_pnt.o
Integration_Cuba.o: quantum_state.o parameters.o

# Модули расчетов
pot0_2004.o: parameters.o quantum_state.o coeff_2004.o angINT.o Integration_Cuba.o others.o
pot1_2004.o: parameters.o prep_to_1pot.o quantum_state.o others.o Integration_Cuba.o angINT.o

# Главная программа
main.o: parameters.o quantum_state.o pot0_2004.o pot1_2004.o others.o

# 4. КОМАНДЫ
run: $(TARGET)
	CUBACORES=1 ./$(TARGET)

clean:
	rm -f $(TARGET) *.o *.mod
	@echo "🧹 Очистка завершена!"

rebuild: clean all